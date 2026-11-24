"""HTML for the 'Sampling times' part of the Methods (built from Rat_sampling_times.xlsx)."""
import numpy as np
import pandas as pd
import openpyxl
from sampling import load, ROOT
from report_html import esc

BIRTH = pd.Timestamp("2026-01-15")   # birth date of all rats (males and females)


def age_label(date, day):
    """Age at sampling in weeks; estimated range when the sampling date was not recorded."""
    if pd.notna(date):
        return f"{(date - BIRTH).days / 7:.1f}"
    if isinstance(day, str) and day.startswith("Day 4"):
        return "12–16*"
    if isinstance(day, str) and day.startswith("Day 6"):
        return "≥16*"
    return "~22*"          # FR3: sampled 18 or 22 June


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


def rat_table(organ, groups=None, males_only=False):
    df = load()
    d = df[(df.organ == organ) & df.hmt.notna()].copy()
    if groups is not None:
        d = d[d.group.isin(groups)]
    if males_only:
        d = d[d.sex == "M"]
    d["gord"] = d.group.map(lambda g: list(GL).index(g) if g in GL else 99)
    d = d.sort_values(["gord", "sex", "num"])
    rows = []
    for _, r in d.iterrows():
        day = r.day if isinstance(r.day, str) else None
        date = r.date.strftime("%d %b %Y") if pd.notna(r.date) else (day or "not recorded")
        if pd.notna(r.date) and day and "June" not in day:
            date += f" ({day})"
        rows.append(f"<tr><td class='l'>{esc(r.hmt)}</td><td class='l' style='white-space:nowrap'>{esc(r.rat)}</td><td>{r.sex}</td>"
                    f"<td>{'–' if pd.isna(r.weight) else int(r.weight)}</td><td class='l' style='white-space:nowrap'>{esc(date)}</td><td>{age_label(r.date, r.day)}</td>"
                    f"<td class='l'>{GL.get(r.group, r.group)}</td><td>{_fmt_s(r.seconds, r.rat, organ)}</td>"
                    f"<td class='l small'>{esc(NOTES.get((r.rat, organ), ''))}</td></tr>")
    return ("<table style='font-size:7pt'><thead><tr><th class='l'>HMT sample</th><th class='l'>Rat</th><th>Sex</th>"
            "<th>Weight (g)</th><th class='l'>Sampling date</th><th>Age<br>(weeks)</th><th class='l'>Group</th>"
            f"<th>{organ} removed<br>(s after arrest)</th><th class='l'>Note</th></tr></thead><tbody>"
            + "".join(rows) + "</tbody></table><p class='small'>Age from the common birth date, 15 January 2026. "
            "</p>")


def control_summary(organ, males_only=False):
    df = load()
    d = df[(df.organ == organ) & df.hmt.notna() & (df.group == 0)]
    parts = []
    for sx, name in ([("M", "males")] if males_only else [("M", "males"), ("F", "females")]):
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
    return (f"All rats, males and females, were born on 15 January 2026. Rats were numbered in sampling order. Males "
            f"were sampled on 7–10 April (days 1–4, rats 1–20), 7 May (day 5) and 20 May 2026 (day 6), at 11.7–17.9 "
            f"weeks of age, and weighed "
            f"{int(m.weight.min())}–{int(m.weight.max())} g; weight increased with rat number (Spearman ρ = "
            f"{rho:.2f}, p = {p:.4f}) as the animals grew during the sampling period. Females were sampled on "
            f"18–22 June 2026 at 22.0–22.6 weeks of age and weighed {int(f.weight.min())}–{int(f.weight.max())} g, "
            "i.e. they were 4–11 weeks older than the males at sampling. Sex and age are therefore fully "
            "confounded in this data set.")


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


DATE_NUMS = {  # hydroxyproline median nmol/g: April males (~12 wk) / May males (16-18 wk) / females (22 wk);
    # change per week of age within ischemia groups (males), p
    "Heart": ("63", "31", "17", "−17%", "0.015"),
    "Liver": ("33", "24", "15", "−4%", "0.08"),
    "Brain": ("14", "12", "7", "−8%", "0.046"),
}
ALLOC = {
    "Heart": "heart 30 and 180 min samples were all taken in May and 360 min samples all in April",
    "Liver": "liver 180 min samples were all taken in May and 360 min samples all in April",
    "Brain": "brain control and 60 min samples were all taken in April",
}


def date_effect_note(tissue=None, males_only=False):
    ts = [tissue] if tissue else ["Heart", "Liver", "Brain"]
    if males_only:
        nums = "; ".join(f"{t.lower()} {DATE_NUMS[t][0]} → {DATE_NUMS[t][1]} nmol/g "
                         f"({DATE_NUMS[t][3]} per week of age, p = {DATE_NUMS[t][4]})" for t in ts)
        alloc = "; ".join(ALLOC[t] for t in ts)
        return ("<div class='note warn'><b>Sampling-date effect (free hydroxyproline).</b> Within each ischemia-time "
                "group, male rats sampled in May (16–18 weeks old) had less free hydroxyproline than those sampled in "
                f"April (~12 weeks): median {nums}. This fits a decline with age (slowing growth and collagen turnover). "
                "Free hydroxyproline reflects collagen turnover, not collagen content (total hydroxyproline after acid "
                f"hydrolysis), which in the heart increases with age. Time groups were not spread evenly over April and "
                f"May ({alloc}), so part of the differences in rat free hydroxyproline between ischemia-time groups "
                "probably reflects age at sampling. No other metabolite showed a significant sampling-date effect after "
                "FDR correction.</div>")
    nums = "; ".join(f"{t.lower()} {DATE_NUMS[t][0]} → {DATE_NUMS[t][1]} → {DATE_NUMS[t][2]} nmol/g "
                     f"({DATE_NUMS[t][3]} per week of age in males, p = {DATE_NUMS[t][4]})" for t in ts)
    alloc = "; ".join(ALLOC[t] for t in ts)
    return ("<div class='note warn'><b>Sampling-date effect (free hydroxyproline).</b> Within each ischemia-time group, older "
            "male rats had less free hydroxyproline, and the females, sampled last, continued the same gradient (median, "
            f"April males ~12 weeks → May males 16–18 weeks → females 22 weeks): {nums}. This fits a decline with age (slowing growth and collagen turnover) between "
            "~12 and ~22 weeks of age rather than a sex or ischemia effect. Free hydroxyproline, as measured here, reflects collagen turnover; "
            "it is not a measure of collagen content (total hydroxyproline after acid hydrolysis), which in the heart "
            "increases with age. Time groups were not spread evenly over April "
            f"and May ({alloc}), so part of the differences in rat free hydroxyproline between ischemia-time groups "
            "probably reflects age at sampling. No other metabolite showed a significant sampling-date effect after "
            "FDR correction (a tendency to ~35% lower levels in May male hearts rests mainly on two animals).</div>")


def species_section(tissue):
    df = load()
    m = df[(df.organ == "Heart") & (df.sex == "M")].dropna(subset=["weight"])
    age = (f"All rats were born on 15 January 2026. The male rats were numbered in sampling order and sampled on 7–10 "
           f"April (rats 1–20, ~12 weeks old), 7 May and 20 May 2026 (rats 21–30, 16–18 weeks old); they weighed "
           f"{int(m.weight.min())}–{int(m.weight.max())} g and weight increased with rat number.")
    return (f"<h3>Sampling times</h3><p>{SELECTION} {control_summary(tissue, males_only=True)}</p>"
            f"{rat_table(tissue, males_only=True)}<p class='small'>{age}</p>{date_effect_note(tissue, males_only=True)}"
            f"<p class='small'>{fish_text(tissue)}</p>{swelling_note(tissue)}")


def sex_section():
    parts = [f"<h3>Sampling times</h3><p>{SELECTION}</p><p>{age_weight_text()}</p>",
             "<p>" + " ".join(control_summary(o) for o in ["Heart", "Liver", "Brain"]) +
             " Female control hearts and livers were therefore removed later, not earlier, than male controls.</p>",
             date_effect_note()]
    for o in ["Heart", "Liver", "Brain"]:
        parts.append(f"<h3>{o}</h3>{rat_table(o, groups=[0, 30])}")
    return "".join(parts)


def swelling_note(tissue=None):
    """Methods paragraph on the tissue-weight (swelling) correction."""
    from data import SWELLING
    ts = [tissue] if tissue else ["Heart", "Brain", "Liver"]
    cols = ["60min", "180min", "360min", "18hours"]
    head = "".join(f"<th>{c.replace('min', ' min').replace('18hours', '18 h')}</th>" for c in cols)
    rows = []
    for t in ts:
        for sp, name in [("Rat", "Rat"), ("Fish", "Crucian carp")]:
            f = SWELLING.get((sp, t))
            cells = "".join(f"<td>{f[c]:.3f}</td>" if f and c in f else
                            ("<td>–</td>" if (sp == "Rat" and c == "18hours") else
                             "<td class='small'>1 (in situ)</td>" if not f else "<td>–</td>") for c in cols)
            rows.append(f"<tr><td class='l'>{t}</td><td class='l'>{name}</td>{cells}</tr>")
    table = (f"<table style='font-size:7pt'><thead><tr><th class='l'>Organ</th><th class='l'>Species</th>{head}</tr>"
             f"</thead><tbody>{''.join(rows)}</tbody></table>")
    return ("<h3>Correction for tissue swelling</h3><p>All rat organs and the crucian carp brain were removed and "
            "kept in buffer during ischemia, and change weight over time (rat organs gain weight; the crucian carp "
            "brain first loses and then gains weight). The crucian carp heart and liver stayed inside the fish "
            "during ischemia (in situ, room temperature), were never in buffer and need no correction. The tissue "
            "weights given to HMT were weighed after the ischemia period, and according to HMT the swelling rates "
            "listed in their Table 1 were not used when calculating nmol/g (to be confirmed in writing). HMT's values "
            "are therefore per gram of post-ischemia tissue. Organ weights were measured in a separate experiment "
            "before and after 1, 3 and 6 h (3 rats; 3 crucian carp) and after 18 h (3 other crucian carp), and all concentrations and "
            "concentration totals in this report were multiplied by the measured fold change in weight to express "
            "them per gram of original tissue (ratios such as energy charge are unaffected). Control and 30 min "
            "samples, and the in-situ crucian carp heart and liver, were not corrected. Factors applied:</p>"
            f"{table}<p class='small'>The factors submitted to HMT for crucian carp were entered for the heart "
            "instead of the brain, so HMT's Table 1 lists them for the carp heart; this report uses them for the "
            "brain. The corrections are at most ~25% (carp brain, 18 h) and ~15% (rat brain, 360 min), small "
            "compared with the ischemic changes.</p>")
