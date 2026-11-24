"""Load HMT metabolome workbooks into tidy structures.

Group labels are taken from the 'Table4 Absolute Concentration' sheet, with one
correction: female rats FR1 and FR2 are swapped in all three tissues (see
README / report Methods). Heart & liver: FR2 -> Ctr, FR1 -> 30min.
Brain (sampled in reverse order): FR1 -> Ctr, FR2 -> 30min.
"""
from pathlib import Path
import re
import numpy as np
import pandas as pd
import openpyxl

ROOT = Path(__file__).resolve().parent.parent
FILES = {
    "Heart": ("OSLUH-HMT-011_Heart_ MetabolomeAnalysisData.xlsx", "OSLUH-HMT-011"),
    "Brain": ("OSLUH-HMT-012_Brain_MetabolomeAnalysisData.xlsx", "OSLUH-HMT-012"),
    "Liver": ("OSLUH-HMT-013_Liver_MetabolomeAnalysisData.xlsx", "OSLUH-HMT-013"),
}
TIMES = ["Ctr", "30min", "60min", "180min", "360min", "18hours"]
TIME_LABEL = {"Ctr": "Ctr", "30min": "30 min", "60min": "60 min", "180min": "180 min",
              "360min": "360 min", "18hours": "18 h"}
SHARED_TIMES = TIMES[:5]          # time points sampled in both species
FR_CORRECTION = {                 # sample -> corrected group
    "Heart": {"FR1H": "30min", "FR2H": "Ctr"},
    "Liver": {"FR1L": "30min", "FR2L": "Ctr"},
    "Brain": {"FR1B": "Ctr", "FR2B": "30min"},
}


def species_of(sample):
    if sample.startswith("FR") or re.fullmatch(r"\d+[hbL]", sample):
        return "Rat"
    return "Fish"


def sex_of(sample):
    if sample.startswith("FR"):
        return "F"
    if species_of(sample) == "Rat":
        return "M"
    return sample[-1]          # fish: trailing F/M


class Tissue:
    def __init__(self, name, correct_fr=True, male_rats_only=False):
        self.name = name
        fname, self.hmt_id = FILES[name]
        wb = openpyxl.load_workbook(ROOT / fname, read_only=True, data_only=True)
        rows = list(wb["Table4 Absolute Concentration"].iter_rows(values_only=True))
        hdr, names = rows[2][5:41], rows[3][5:41]
        grp, g = [], None
        for h in hdr:
            g = h or g
            grp.append("Ctr" if g.startswith("Ctr") else g)
        self.original_group = dict(zip(names, grp))
        self.group = dict(self.original_group)
        if correct_fr:
            self.group.update(FR_CORRECTION[name])
        recs, meta = {}, []
        for r in rows[4:]:
            if not (r[0] and r[1]):
                continue
            vals = r[5:41]
            recs[r[1]] = [v if isinstance(v, (int, float)) else np.nan for v in vals]
            nd = [isinstance(v, str) and v.strip().upper().startswith("N.D") for v in vals]
            meta.append(dict(name=r[1], id=r[0], label=r[2] or r[1], cid=r[3], hmdb=r[4],
                             is_param=(str(r[0]).strip() in ("-", "") or r[0] is None or not re.match(r"^[AC]_\d+", str(r[0]))),
                             n_nd=sum(nd)))
        self.data = pd.DataFrame(recs, index=names).T            # metabolite x sample
        self.meta = pd.DataFrame(meta).set_index("name")
        self.meta.loc[self.meta.is_param, "label"] = self.meta.index[self.meta.is_param]
        self.samples = pd.DataFrame({
            "sample": names,
            "species": [species_of(s) for s in names],
            "sex": [sex_of(s) for s in names],
            "group": [self.group[s] for s in names],
            "hmt_group": [self.original_group[s] for s in names],
        }).set_index("sample")
        # tissue weights from sample sheet
        ws = list(wb["Table 1 Sample Information"].iter_rows(values_only=True))
        w = {r[0]: r[1] for r in ws[2:] if r[0] in self.samples.index}
        self.samples["weight_mg"] = pd.Series(w)
        if male_rats_only:   # drop female rats (FR*) – used by the fish-vs-rat reports
            keep = [s for s in self.samples.index if not s.startswith("FR")]
            self.samples = self.samples.loc[keep]
            self.data = self.data[keep]
            self.group = {k: v for k, v in self.group.items() if k in keep}
            self.original_group = {k: v for k, v in self.original_group.items() if k in keep}

    @property
    def compounds(self):
        return self.meta.index[~self.meta.is_param]

    @property
    def params(self):
        return self.meta.index[self.meta.is_param]

    def cols(self, species=None, group=None, sex=None):
        s = self.samples
        m = pd.Series(True, index=s.index)
        if species: m &= s.species == species
        if group: m &= s.group == group
        if sex: m &= s.sex == sex
        return list(s.index[m])

    def log2_imputed(self, rows=None, cols=None):
        """log2 values with N.D. replaced by half the lowest detected value of that metabolite
        (across all samples of the tissue)."""
        X = self.data.loc[rows if rows is not None else self.compounds]
        mins = X.min(axis=1)
        X = X.apply(lambda r: r.fillna(mins[r.name] / 2), axis=1)
        X = np.log2(X.where(X > 0))
        return X if cols is None else X[cols]


def load_all(correct_fr=True):
    return {t: Tissue(t, correct_fr) for t in FILES}
