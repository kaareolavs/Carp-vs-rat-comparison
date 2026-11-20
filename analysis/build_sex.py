"""Male vs female rats at control and 30 min ischemia – one report covering heart, brain and liver."""
import datetime as dt
import warnings
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
import matplotlib.pyplot as plt
import figures as F
from data import Tissue, FR_CORRECTION
from stats import sex_comparison
from pathways import MAPS
from report_html import doc, esc, fmt_q, fc_cell, num, figure
import findings

XL = ["Ctr", "30 min"]
HALLMARK = ["ATP", "ADP", "AMP", "Phosphocreatine", "Lactic acid", "Succinic acid", "Ala", "γ-Aminobutyric acid",
            "Glycerol 3-phosphate", "Asp", "Glu", "Citric acid", "Inosine", "Hypoxanthine", "Glucose 6-phosphate",
            "Adenylate Energy Charge"]
SEXN = {"M": "Male", "F": "Female"}


def series_for(x):
    return [dict(name="Male rats", color=F.C1, offset=-0.06,
                 groups=[(0, x.cols("Rat", "Ctr", "M")), (1, x.cols("Rat", "30min", "M"))]),
            dict(name="Female rats", color=F.C2, offset=0.06,
                 groups=[(0, x.cols("Rat", "Ctr", "F")), (1, x.cols("Rat", "30min", "F"))])]


def tissue_section(tname, sec):
    x = Tissue(tname)
    R = sex_comparison(x)
    keep = R["keep"]
    label = lambda k: x.meta.loc[k, "label"] if k in x.meta.index else k
    comp = [c for c in x.compounds if keep[c]]
    par = [p for p in x.params if keep[p]]
    A = R["anova"]
    marks = {}
    for k in x.data.index:
        m = ""
        if A["sex"]["q"].get(k, np.nan) < 0.05: m += "*"
        if A["interaction"]["q"].get(k, np.nan) < 0.05: m += "†"
        if m: marks[k] = m
    series = series_for(x)
    S = R["summary"]
    mean = lambda sx, t, k: S.loc[k, (f"{sx}|{t}", "mean")]

    # counts
    cnt = {
        "Female vs male, control": int((R["sex_at"]["Ctr"].loc[comp].q < 0.05).sum()),
        "Female vs male, 30 min": int((R["sex_at"]["30min"].loc[comp].q < 0.05).sum()),
        "Male: 30 min vs control": int((R["response"]["M"].loc[comp].q < 0.05).sum()),
        "Female: 30 min vs control": int((R["response"]["F"].loc[comp].q < 0.05).sum()),
        "Sex main effect": int((A["sex"].loc[comp].q < 0.05).sum()),
        "Sex × ischemia interaction": int((A["interaction"].loc[comp].q < 0.05).sum()),
    }

    # PCA (rats Ctr + 30)
    cols = [c for v in R["groups"].values() for c in v]
    sc, ev = F.pca(R["Y"].loc[comp, cols].T)
    pser = [dict(name=SEXN[sx], color=c, groups=[("Ctr", R["groups"][(sx, "Ctr")]), ("30", R["groups"][(sx, "30min")])])
            for sx, c in [("M", F.C1), ("F", F.C2)]]
    fig, ax = plt.subplots(figsize=(4.2, 3.2))
    F.pca_panel(ax, sc, ev, x.samples, pser, "")
    for c in cols:
        ax.annotate(c, (sc.loc[c, "PC1"], sc.loc[c, "PC2"]), xytext=(2, -6), textcoords="offset points",
                    fontsize=4.8, color=F.MUTED)
    ax.legend(handles=F.legend_handles(pser), frameon=False, fontsize=6.5, loc="upper left", bbox_to_anchor=(0, 1.13), ncol=2)
    pca_svg = F.svg(fig)

    hm = {"F/M|Ctr": R["sex_at"]["Ctr"].loc[comp, "effect"], "F/M|30": R["sex_at"]["30min"].loc[comp, "effect"],
          "M|30 vs Ctr": R["response"]["M"].loc[comp, "effect"], "F|30 vs Ctr": R["response"]["F"].loc[comp, "effect"]}
    sg = {"F/M|Ctr": R["sex_at"]["Ctr"].loc[comp, "q"] < 0.05, "F/M|30": R["sex_at"]["30min"].loc[comp, "q"] < 0.05,
          "M|30 vs Ctr": R["response"]["M"].loc[comp, "q"] < 0.05, "F|30 vs Ctr": R["response"]["F"].loc[comp, "q"] < 0.05}
    heat_svg = F.heatmap(pd.DataFrame(hm), pd.DataFrame(sg),
                         [("Female / male", ["F/M|Ctr", "F/M|30"]), ("Ischemia response", ["M|30 vs Ctr", "F|30 vs Ctr"])],
                         label, cbar_label="log2 ratio")
    hall = [h for h in HALLMARK if h in x.data.index]
    hall_svg = F.panel_grid(x, hall, series, XL, label, marks, ncols=4, height=1.5,
                            ylabel_of=lambda k: "no unit" if "Charge" in k else "nmol/g")
    par_svg = F.panel_grid(x, par, series, XL, label, marks, ncols=5, height=1.25)
    maps = [(m["title"], F.pathway_map(m, x, series, XL, marks, label)) for m in MAPS[:-1]]

    def table(names):
        rows = []
        for k in names:
            lin = not R["islog"].get(k, True)
            cells = [f'<td class="l">{esc(label(k))}{" " + marks[k] if k in marks else ""}</td>']
            for sx in ["M", "F"]:
                for t in ["Ctr", "30min"]:
                    cells.append(f"<td>{num(mean(sx, t, k))}</td>")
            for d in [R["sex_at"]["Ctr"], R["sex_at"]["30min"], R["response"]["M"], R["response"]["F"]]:
                cells.append(fc_cell(d.loc[k, "effect"], d.loc[k, "q"], linear=lin))
            cells.append(f"<td>{fmt_q(A['sex'].loc[k, 'q'])}</td><td>{fmt_q(A['interaction'].loc[k, 'q'])}</td>")
            rows.append("<tr>" + "".join(cells) + "</tr>")
        head = ('<thead><tr><th class="l" rowspan="2">Metabolite</th><th colspan="2" style="text-align:center">Male '
                '(mean)</th><th colspan="2" style="text-align:center">Female (mean)</th><th colspan="2" '
                'style="text-align:center">log2 F/M</th><th colspan="2" style="text-align:center">log2 30 min/Ctr</th>'
                '<th rowspan="2">Sex<br>q</th><th rowspan="2">Sex×time<br>q</th></tr><tr><th>Ctr</th><th>30</th>'
                '<th>Ctr</th><th>30</th><th>Ctr</th><th>30</th><th>M</th><th>F</th></tr></thead>')
        return f'<table style="font-size:6.9pt">{head}<tbody>{"".join(rows)}</tbody></table>'

    sig = sorted({k for k in comp + par if (A["sex"].loc[k, "q"] < 0.05) or (A["interaction"].loc[k, "q"] < 0.05)
                  or R["sex_at"]["Ctr"].loc[k, "q"] < 0.05 or R["sex_at"]["30min"].loc[k, "q"] < 0.05},
                 key=lambda k: A["sex"].loc[k, "p"] if not np.isnan(A["sex"].loc[k, "p"]) else 1)

    ids = lambda sx, t: ", ".join(f"<b>{c}</b>†" if c in FR_CORRECTION[tname] else c for c in x.cols("Rat", t, sx))
    out = []
    out.append(f"""<section class="page"><h2>{sec}. {tname}</h2>
<table><thead><tr><th class="l"></th><th class="l">Control</th><th class="l">30 min ischemia</th></tr></thead><tbody>
<tr><td class="l">Male</td><td class="l">{ids("M", "Ctr")}</td><td class="l">{ids("M", "30min")}</td></tr>
<tr><td class="l">Female</td><td class="l">{ids("F", "Ctr")}</td><td class="l">{ids("F", "30min")}</td></tr></tbody></table>
<div class="kpis" style="grid-template-columns:repeat(3,1fr)">{''.join(f'<div class="kpi"><div class="v">{v}</div><div class="l">{esc(k)}</div></div>' for k, v in cnt.items())}</div>
<p class="small">Counts: metabolites with q &lt; 0.05 out of {len(comp)} tested. † corrected label (FR1/FR2, see Methods).</p>
<div style="display:grid;grid-template-columns:1.15fr 1fr;gap:10pt;align-items:start">
{figure(pca_svg, f"<b>Fig. {sec}.1. PCA</b> of rat {tname.lower()} samples (autoscaled log2, {len(comp)} metabolites). "
        "Large dots: group means, Ctr → 30 min.")}
<div>{''.join(f'<p>{p}</p>' for p in findings.SEX_TISSUE.get(tname, []))}</div></div></section>""")
    out.append(f'<section class="page">{figure(hall_svg, f"<b>Fig. {sec}.2. Hallmark ischemia metabolites</b>, male vs female rats (mean ± SD, individual animals). * sex main effect q &lt; 0.05; † sex × ischemia interaction q &lt; 0.05.")}</section>')
    out.append(f'<section class="page">{figure(heat_svg, f"<b>Fig. {sec}.3. Sex differences and ischemia responses.</b> Left: log2(female/male) at control and 30 min. Right: log2(30 min / control) in each sex. Dots: q &lt; 0.05.")}</section>')
    for i, (title, s) in enumerate(maps):
        out.append(f'<section class="page"><h3>{sec}.{i+1} {esc(title)} – {tname.lower()}</h3>{figure(s, "")}</section>')
    out.append(f'<section class="page"><h3>Metabolic parameters – {tname.lower()}</h3>{figure(par_svg, f"<b>Fig. {sec}.4. Derived metabolic parameters</b>, male vs female rats.")}</section>')
    out.append(f"""<section class="page"><h3>Metabolites with a sex difference – {tname.lower()}</h3>
<p class="small">Listed if q &lt; 0.05 for sex at control, sex at 30 min, the sex main effect or the sex × ischemia
interaction ({len(sig)} metabolites). Cells: log2 ratios, bold * = q &lt; 0.05.</p>
{table(sig) if sig else '<p>No metabolite reached q &lt; 0.05 for any sex comparison.</p>'}
<h3>Hallmark metabolites</h3>{table([h for h in hall if h in comp + par])}</section>""")

    xl = {}
    m = S.copy(); m.columns = [f"{a} {b}" for a, b in m.columns]
    rows = {}
    for t, d in R["sex_at"].items(): rows[f"log2_F_over_M_{t}"] = d.effect; rows[f"q_F_vs_M_{t}"] = d.q
    for sx, d in R["response"].items(): rows[f"log2_30min_over_Ctr_{sx}"] = d.effect; rows[f"q_response_{sx}"] = d.q
    for k, d in A.items(): rows[f"model_{k}_effect"] = d.effect; rows[f"model_{k}_q"] = d.q
    xl[f"{tname}_means"] = m
    xl[f"{tname}_stats"] = pd.DataFrame(rows)
    return "".join(out), xl, cnt


def build():
    today = dt.date.today().strftime("%d %B %Y")
    body = [f"""
<section class="cover">
  <div class="kicker">Comparative metabolomics · rat · heart, brain, liver</div>
  <h1>Male vs female rats:<br>metabolome at control and after 30 min ischemia</h1>
  <div class="sub">Re-analysis of HMT data sets OSLUH-HMT-011 (heart), -012 (brain) and -013 (liver)</div>
  <p><span class="swatch" style="background:{F.C1}"></span>Male rats &nbsp;
     <span class="swatch" style="background:{F.C2}"></span>Female rats</p>
  <ol class="toc"><li>1. Summary</li><li>2. Design and methods</li><li>3. Heart</li><li>4. Brain</li><li>5. Liver</li></ol>
  <div class="meta">Source data: HMT MetabolomeAnalysisData workbooks (CE-TOFMS absolute quantification). Report
  generated {today}. Female rats FR1 and FR2 are analysed with corrected group labels (section 2).</div>
</section>"""]
    sections, xl, counts = [], {}, {}
    for i, t in enumerate(["Heart", "Brain", "Liver"]):
        s, x_, c = tissue_section(t, i + 3)
        sections.append(s); xl.update(x_); counts[t] = c
    ct = pd.DataFrame(counts)
    body.append(f"""<section class="page"><h2>1. Summary</h2>
<h3>Key findings</h3><ul class="findings">{''.join(f'<li>{c}</li>' for c in findings.SEX)}</ul>
<h3>Number of metabolites with q &lt; 0.05</h3>
<table><thead><tr><th class="l">Comparison</th>{''.join(f'<th>{t}</th>' for t in ct.columns)}</tr></thead><tbody>
{''.join(f"<tr><td class='l'>{esc(i)}</td>" + ''.join(f'<td>{v}</td>' for v in r.values) + '</tr>' for i, r in ct.iterrows())}
</tbody></table>
<p class="small">Tested metabolites per tissue: those detected in at least three of the 12 rat samples at control and
30 min plus the derived parameters are reported in the tissue sections; counts here are for metabolites only.</p>
</section>""")
    body.append(f"""<section class="page"><h2>2. Design and methods</h2>
<p>In each tissue, 3 male and 3 female rats were sampled at control and 3 + 3 after 30 min ischemia (2 × 2 design,
n = 3 per cell). Later time points (60–360 min) contain only male rats and are therefore not part of this report.</p>
<p>Female rats gave heart and liver at one time point and brain at the other, because brain removal takes longer
(~45 s) than heart or liver (~15 s): for FR2, FR4 and FR5 the brain was sampled first (control) and heart/liver represent
30 min; for FR1, FR3 and FR6 the order was reversed. Male rat IDs differ between tissues (e.g. 13h, 11b, 12L), i.e.
group membership is defined per tissue as in the HMT files.</p>
<div class="note warn"><b>† Corrected labels for FR1 and FR2.</b> In all three tissues the metabolite profiles of FR1
and FR2 were the reverse of their labels (the "control" sample was fully ischemic and the "30 min" sample had a normal
control profile), consistent with the two animals' identities being exchanged. Used here: heart and liver – FR2 =
control, FR1 = 30 min; brain – FR1 = control, FR2 = 30 min.</div>
<h3>Statistics</h3>
<p>Concentrations were log2-transformed (energy charges kept linear); N.D. values were replaced by half the lowest
detected value of the metabolite in that tissue for testing, and metabolites detected in fewer than three samples were not
tested. For every metabolite:</p>
<ul><li><b>Female vs male</b> at control and at 30 min (Welch's t-test; effect = log2 female/male).</li>
<li><b>Ischemia response</b> (30 min vs control) within each sex (Welch's t-test; log2 fold change).</li>
<li><b>Two-way linear model</b> log2(concentration) ~ sex × time: the sex main effect (average female–male difference
over both time points, marked <b>*</b>) and the sex × ischemia interaction (does the 30 min response differ between
sexes, marked <b>†</b>).</li></ul>
<p>Benjamini–Hochberg FDR (q) was applied within each comparison and tissue, separately for metabolites and derived
parameters.</p>
<div class="note"><b>Limitations.</b> With n = 3 per cell only large, consistent differences can reach significance, and
a non-significant result does not show that the sexes are equal. Sex is confounded with anything else that differed
between the male and female cohorts (e.g. age, body weight, handling or the paired tissue-sampling scheme used for
females). Some male controls (e.g. heart 13h, 10h) show more AMP and inosine than female controls, which may reflect
sampling speed rather than sex.</div>
<h3>Figure conventions</h3><p><span class="swatch" style="background:{F.C1}"></span>male rats &nbsp;
<span class="swatch" style="background:{F.C2}"></span>female rats; mean ± SD; hollow marker = at least one animal
N.D. In pathway maps and panels: <b>*</b> sex main effect q &lt; 0.05, <b>†</b> sex × ischemia interaction q &lt; 0.05.</p>
</section>""")
    body += sections
    return doc("Male vs female rats", "".join(body)), xl
