"""HTML for the 'Sampling times' part of the Methods (built from Rat_sampling_times.xlsx)."""
import numpy as np
import pandas as pd
import openpyxl
from sampling import load, ROOT
from report_html import esc

GL = {0: "Control", 30: "30 min", 60: "60 min", 180: "180 min", 360: "360 min"}
NOTES = {
    ("FR1", "Brain"): "brain taken first, within 1 min (exact time not recorded)",
    ("FR1", "Heart"): "~90 s (approximate); sampled in reverse of planned order",
    ("FR1", "Liver"): "~90 s (approximate); sampled in reverse of planned order",
    ("FR2", "Heart"): "sampled in reverse of planned order",
    ("FR2", "Liver"): "sampled in reverse of planned order",
    ("FR2", "Brain"): "sampled in reverse of planned order",
    ("Rat 9", "Brain"): "arrest time noted to the minute",
    ("Rat 10", "Brain"): "arrest time noted to the minute",
    ("Rat 10", "Liver"): "arrest time noted to the minute",
    ("Rat 13", "Heart"): "cut-to-freeze 7 s",
    ("Rat 14", "Brain"): "47 s direct procurement",
    ("Rat 14", "Liver"): "arrest derived from brain time",
    ("Rat 14", "Heart"): "arrest derived from brain time",
    ("Rat 27", "Heart"): "arrest written 1:30:02 (13:30:02)",
    ("Rat 27", "Liver"): "arrest written 1:30:02 (13:30:02)",
}
SELECTION = ("Six animals were sampled per group and organ, but for cost reasons only 3 (2 at 360 min) were sent "
             "for metabolomics; the samples with the shortest sampling times and no procedural errors were chosen. "
             "Circulatory arrest (time of death) was noted when the heart was cut out, and the clock time at which "
             "each organ was removed was recorded. Data: <i>Rat_sampling_times.xlsx</i>.")


def _fmt_s(v, rat, organ):
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return "&lt;60" if (rat, organ) == ("FR1", "Brain") else "–"
    return f"{int(round(v))}"


def rat_table(organ, groups=None):
    df = load()
    d = df[(df.organ == organ) & df.hmt.notna()].copy()
    if groups is not None:
        d = d[d.group.isin(groups)]
    d["gord"] = d.group.map(lambda g: list(GL).index(g) if g in GL else 99)
    d = d.sort_values(["gord", "sex", "num"])
    rows = []
    for _, r in d.iterrows():
        day = r.day if isinstance(r.day, str) else None
        date = r.date.strftime("%d %b %Y") if pd.notna(r.date) else (day or "not recorded")
        if pd.notna(r.date) and day and "June" not in day:
            date += f" ({day})"
        rows.append(f"<tr><td class='l'>{esc(r.hmt)}</td><td class='l' style='white-space:nowrap'>{esc(r.rat)}</td><td>{r.sex}</td>"
                    f"<td>{'–' if pd.isna(r.weight) else int(r.weight)}</td><td class='l' style='white-space:nowrap'>{esc(date)}</td>"
                    f"<td class='l'>{GL.get(r.group, r.group)}</td><td>{_fmt_s(r.seconds, r.rat, organ)}</td>"
                    f"<td class='l small'>{esc(NOTES.get((r.rat, organ), ''))}</td></tr>")
    return ("<table style='font-size:7pt'><thead><tr><th class='l'>HMT sample</th><th class='l'>Rat</th><th>Sex</th>"
            "<th>Weight (g)</th><th class='l'>Sampling date</th><th class='l'>Group</th>"
            f"<th>{organ} removed<br>(s after arrest)</th><th class='l'>Note</th></tr></thead><tbody>"
            + "".join(rows) + "</tbody></table>")


def control_summary(organ):
    df = load()
    d = df[(df.organ == organ) & df.hmt.notna() & (df.group == 0)]
    parts = []
    for sx, name in [("M", "males"), ("F", "females")]:
        s = d[d.sex == sx].seconds.dropna()
        extra = "; FR1 &lt;60 s, exact time not recorded" if organ == "Brain" and sx == "F" else ""
        if len(s):
            parts.append(f"{int(s.min())}–{int(s.max())} s ({name}{extra})")
    return f"{organ} control samples were removed {' and '.join(parts)} after circulatory arrest."


def age_weight_text():
    df = load()
    m = df[(df.organ == "Heart") & (df.sex == "M")].dropna(subset=["weight"])
    from scipy import stats
    rho, p = stats.spearmanr(m.num, m.weight)
    f = df[(df.organ == "Heart") & (df.sex == "F")]
    return (f"Rats were numbered in sampling order. Males were sampled between 7 April and May 2026 (days 1–6) and "
            f"weighed {int(m.weight.min())}–{int(m.weight.max())} g; weight increased with rat number (Spearman ρ = "
            f"{rho:.2f}, p = {p:.4f}), consistent with animals of the same age at arrival growing during the sampling "
            f"period. Females were sampled on 18–22 June 2026 and weighed {int(f.weight.min())}–{int(f.weight.max())} g. "
            "If all rats arrived at the same age (to be confirmed), the females were roughly 6–10 weeks older than the "
            "males at sampling.")


def fish_text(tissue):
    wb = openpyxl.load_workbook(ROOT / "Copy of Owerview Rat_Fish sampled_weight_sex_ischemic time_22062026.xlsx",
                                data_only=True)
    w = {}
    for r in list(wb["Crucian carp"].iter_rows(values_only=True))[2:]:
        if isinstance(r[0], int) and isinstance(r[1], (int, float)):
            w[r[0]] = r[1]
    txt = ("Crucian carp were sampled from 6 June 2025 over three consecutive days (room temperature 22.3–22.4 °C). "
           "Organ removal times were not recorded individually; for control fish 26 and 30 the heart was sampled "
           "first. Body weights of the fish analysed: "
           + ", ".join(f"{n}{s} {w[n]:.0f} g" for n, s in
                       [(26, "F"), (30, "F"), (38, "M"), (24, "F"), (25, "M"), (36, "F"), (27, "M"), (28, "F"), (29, "F"),
                        (3, "F"), (4, "F"), (21, "M"), (15, "F"), (16, "F"), (11, "F"), (13, "M")] if n in w) + ".")
    if tissue == "Brain":
        txt += (" Lab note: the brain labelled fish 21 (180 min) may be the brain of fish 22 (240 min, not otherwise "
                "analysed); the metabolite profile cannot distinguish 180 from 240 min.")
    return txt


DATE_NUMS = {  # hydroxyproline, nmol/g: early males (rats 9-14) / late males (rats 17-29) / females; late/early ratio
    "Heart": ("52–102", "13–53", "16–23", "0.37"),
    "Liver": ("32–39", "20–36", "12–17", "0.79"),
    "Brain": ("13–27", "8–14", "6–9", "0.72"),
}
ALLOC = {
    "Heart": "all heart 180 and 360 min samples were taken on the late days (4–6)",
    "Liver": "all liver 180 and 360 min samples were taken on the late days (4–6)",
    "Brain": "all brain control and 60 min samples were taken on the early days (2–3)",
}


def date_effect_note(tissue=None):
    ts = [tissue] if tissue else ["Heart", "Liver", "Brain"]
    nums = "; ".join(f"{t.lower()} {DATE_NUMS[t][0]} (early males), {DATE_NUMS[t][1]} (late males), "
                     f"{DATE_NUMS[t][2]} (females), late/early {DATE_NUMS[t][3]}×" for t in ts)
    alloc = "; ".join(ALLOC[t] for t in ts)
    return ("<div class='note warn'><b>Sampling-date effect (hydroxyproline).</b> Within each ischemia-time group, male "
            "rats sampled later had less hydroxyproline (p = 0.01–0.04), and the females, sampled last, continued the "
            f"same gradient: {nums} nmol/g. This fits a decline with age (slowing growth and collagen turnover) over "
            "the sampling period rather than a sex or ischemia effect. Time groups were not spread evenly over the "
            f"sampling days ({alloc}), so part of the change in rat hydroxyproline with ischemia time is probably a "
            "sampling-date effect. Other metabolites showed no consistent period effect after FDR correction "
            "(a tendency to ~35% lower levels in later male hearts rests mainly on two animals).</div>")


def species_section(tissue):
    return (f"<h3>Sampling times</h3><p>{SELECTION} {control_summary(tissue)}</p>"
            f"{rat_table(tissue)}<p class='small'>{age_weight_text()}</p>{date_effect_note(tissue)}"
            f"<p class='small'>{fish_text(tissue)}</p>")


def sex_section():
    parts = [f"<h3>Sampling times</h3><p>{SELECTION}</p><p>{age_weight_text()}</p>",
             "<p>" + " ".join(control_summary(o) for o in ["Heart", "Liver", "Brain"]) +
             " Female control hearts and livers were therefore removed later, not earlier, than male controls.</p>",
             date_effect_note()]
    for o in ["Heart", "Liver", "Brain"]:
        parts.append(f"<h3>{o}</h3>{rat_table(o, groups=[0, 30])}")
    return "".join(parts)
