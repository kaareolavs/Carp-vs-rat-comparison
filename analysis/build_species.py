"""Build the crucian carp vs rat report for one tissue (HTML; printed to PDF by print_pdf.js)."""
import datetime as dt
import warnings
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
import figures as F
from data import Tissue, TIMES, TIME_LABEL, SHARED_TIMES, FR_CORRECTION
from stats import fish_vs_rat
from pathways import MAPS
from report_html import doc, esc, fmt_q, fc_cell, num, figure
from sampling_section import species_section
import findings

XL = ["C", "30", "60", "180", "360", "18h"]
SP_NAME = {"Fish": "Crucian carp", "Rat": "Rat"}
HALLMARK = ["ATP", "ADP", "AMP", "Phosphocreatine", "Lactic acid", "Succinic acid", "Ala", "γ-Aminobutyric acid",
            "Glycerol 3-phosphate", "Asp", "Glu", "Citric acid", "Malic acid", "Inosine", "Hypoxanthine",
            "Glucose 6-phosphate"]
PARAM_UNIT_FREE = ("Charge", "/", "Ratio")


def series_for(x):
    return [dict(name="Crucian carp", color=F.C1, offset=-0.07,
                 groups=[(i, x.cols("Fish", t)) for i, t in enumerate(TIMES) if x.cols("Fish", t)]),
            dict(name="Rat", color=F.C2, offset=0.07,
                 groups=[(i, x.cols("Rat", t)) for i, t in enumerate(TIMES) if x.cols("Rat", t)])]


def build(tissue_name):
    x = Tissue(tissue_name, male_rats_only=True)
    R = fish_vs_rat(x)
    keep = R["keep"]
    label = lambda k: x.meta.loc[k, "label"] if k in x.meta.index else k
    comp = [c for c in x.compounds if keep[c]]
    par = [p for p in x.params if keep[p]]
    Fq = R["interaction_F"]["q"]
    marks = {k: "*" for k in x.data.index if Fq.get(k, np.nan) < 0.05}
    series = series_for(x)
    S = R["summary"]

    def mean(sp, t, k):
        return S.loc[k, (f"{sp}|{t}", "mean")] if (f"{sp}|{t}", "mean") in S.columns else np.nan

    # ------------------------------------------------------------------ counts
    counts = []
    for t in TIMES[1:]:
        row = {"Time": TIME_LABEL[t]}
        for sp in ["Fish", "Rat"]:
            d = R["response"].get((sp, t))
            if d is None:
                row[f"{sp}_up"] = row[f"{sp}_dn"] = None
                continue
            dc = d.loc[comp]
            row[f"{sp}_up"] = int(((dc.q < 0.05) & (dc.effect > 0)).sum())
            row[f"{sp}_dn"] = int(((dc.q < 0.05) & (dc.effect < 0)).sum())
        if t in SHARED_TIMES:
            row["species"] = int((R["species"][t].loc[comp].q < 0.05).sum())
            row["inter"] = int((R["interaction"][t].loc[comp].q < 0.05).sum())
        counts.append(row)
    counts = pd.DataFrame(counts).set_index("Time")
    base_diff = int((R["species"]["Ctr"].loc[comp].q < 0.05).sum())
    n_inter = int((Fq.loc[comp] < 0.05).sum())

    # ------------------------------------------------------------------ figures
    Y = R["Y"].loc[comp]
    cols_all = [c for (sp, t), cs in R["groups"].items() for c in cs]
    sc, ev = F.pca(Y[cols_all].T)
    # response space: centre each species on own control mean
    Yr = Y[cols_all].copy()
    for sp in ["Fish", "Rat"]:
        cs = [c for c in cols_all if x.samples.loc[c, "species"] == sp]
        Yr[cs] = Yr[cs].sub(Y[R["groups"][(sp, "Ctr")]].mean(axis=1), axis=0)
    Zr = Yr.T
    Zr = Zr.loc[:, Zr.std() > 0]
    U, Sv, Vt = np.linalg.svd((Zr / Zr.std()).values, full_matrices=False)
    evr = Sv ** 2 / (Sv ** 2).sum()
    scr = pd.DataFrame(U[:, :2] * Sv[:2], index=Zr.index, columns=["PC1", "PC2"])
    pser = [dict(name=SP_NAME[sp], color=c, groups=[(XL[i], R["groups"][(sp, t)]) for i, t in enumerate(TIMES)
                                                     if (sp, t) in R["groups"]]) for sp, c in [("Fish", F.C1), ("Rat", F.C2)]]
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(7.3, 3.3))
    F.pca_panel(axes[0], sc, ev, x.samples, pser, "A  Absolute levels (both species)")
    F.pca_panel(axes[1], scr, evr, x.samples, pser, "B  Response relative to own control")
    fig.legend(handles=F.legend_handles(pser), loc="upper right", ncol=2, frameon=False, fontsize=7,
               bbox_to_anchor=(1, 1.04))
    fig.tight_layout()
    pca_svg = F.svg(fig)

    ct = counts.rename(columns={})
    cb = pd.DataFrame({
        "Crucian carp vs own control": counts["Fish_up"].fillna(0) + counts["Fish_dn"].fillna(0),
        "Rat vs own control": counts["Rat_up"].fillna(0) + counts["Rat_dn"].fillna(0),
        "Carp vs rat (absolute)": counts["species"].fillna(0),
        "Response differs (interaction)": counts["inter"].fillna(0),
    })
    counts_svg = F.count_bars(cb, ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"])

    hm_cols, hm_sig, groups = {}, {}, []
    for sp in ["Fish", "Rat"]:
        cs = []
        for t in TIMES[1:]:
            d = R["response"].get((sp, t))
            if d is None:
                continue
            key = f"{sp}|{XL[TIMES.index(t)]}"
            hm_cols[key] = d.loc[comp, "effect"]; hm_sig[key] = d.loc[comp, "q"] < 0.05
            cs.append(key)
        groups.append((SP_NAME[sp], cs))
    heat_svg = F.heatmap(pd.DataFrame(hm_cols), pd.DataFrame(hm_sig), groups, label)

    # species ratio heat map (absolute differences)
    sr_cols, sr_sig = {}, {}
    for t in SHARED_TIMES:
        key = f"Carp/Rat|{XL[TIMES.index(t)]}"
        sr_cols[key] = R["species"][t].loc[comp, "effect"]; sr_sig[key] = R["species"][t].loc[comp, "q"] < 0.05
    ratio_svg = F.heatmap(pd.DataFrame(sr_cols), pd.DataFrame(sr_sig), [("log2 (crucian carp / rat)", list(sr_cols))],
                          label, cbar_label="log2 carp / rat", clip=5)

    hall = [h for h in HALLMARK if h in x.data.index]
    hall_svg = F.panel_grid(x, hall, series, XL, label, marks, ncols=4,
                            ylabel_of=lambda k: "nmol/g")
    param_pages = []
    for i in range(0, len(par), 12):
        chunk = par[i:i + 12]
        param_pages.append(F.panel_grid(x, chunk, series, XL, label, marks, ncols=3, height=1.75,
                                        ylabel_of=lambda k: ("ratio / no unit" if any(s in k for s in PARAM_UNIT_FREE)
                                                             else "nmol/g")))
    map_svgs = [(m["title"], F.pathway_map(m, x, series, XL, marks, label)) for m in MAPS]

    # ------------------------------------------------------------------ tables
    def resp_table(names, title_cols=True):
        rows = []
        for k in names:
            cells = [f'<td class="l">{esc(label(k))}{" *" if marks.get(k) else ""}</td>']
            lin = not R["islog"].get(k, True)
            for sp in ["Fish", "Rat"]:
                for t in SHARED_TIMES[1:]:
                    d = R["response"][(sp, t)]
                    cells.append(fc_cell(d.loc[k, "effect"], d.loc[k, "q"], linear=lin))
            cells.append(f"<td>{fmt_q(Fq.get(k))}</td>")
            rows.append("<tr>" + "".join(cells) + "</tr>")
        head = ('<thead><tr><th class="l" rowspan="2">Metabolite</th><th colspan="4" style="text-align:center">'
                'Crucian carp vs own control</th><th colspan="4" style="text-align:center">Rat vs own control</th>'
                '<th rowspan="2">Interaction<br>q</th></tr><tr>'
                + "".join(f"<th>{XL[TIMES.index(t)]}</th>" for t in SHARED_TIMES[1:]) * 2 + "</tr></thead>")
        return f"<table>{head}<tbody>{''.join(rows)}</tbody></table>"

    top_inter = Fq.loc[comp].dropna().sort_values().index[:40]

    def species_table(t, names):
        rows = []
        for k in names:
            d = R["species"][t]
            rows.append(f'<tr><td class="l">{esc(label(k))}</td><td>{num(mean("Fish", t, k))}</td>'
                        f'<td>{num(mean("Rat", t, k))}</td>{fc_cell(d.loc[k, "effect"], d.loc[k, "q"], clip=5)}'
                        f'<td>{fmt_q(d.loc[k, "q"])}</td></tr>')
        return ('<table><thead><tr><th class="l">Metabolite</th><th>Carp<br>nmol/g</th><th>Rat<br>nmol/g</th>'
                '<th>log2<br>carp/rat</th><th>q</th></tr></thead><tbody>' + "".join(rows) + "</tbody></table>")

    base = R["species"]["Ctr"].loc[comp]
    base_sig = base[base.q < 0.05].sort_values("effect", ascending=False).index

    def means_table(names, digits=None):
        rows = []
        for k in names:
            cells = [f'<td class="l">{esc(label(k))}</td>']
            for sp in ["Fish", "Rat"]:
                for t in TIMES:
                    if (f"{sp}|{t}", "mean") in S.columns:
                        cells.append(f"<td>{num(mean(sp, t, k), digits)}</td>")
            rows.append("<tr>" + "".join(cells) + "</tr>")
        ts_f = [XL[i] for i, t in enumerate(TIMES) if ("Fish", t) in R["groups"]]
        ts_r = [XL[i] for i, t in enumerate(TIMES) if ("Rat", t) in R["groups"]]
        head = ('<thead><tr><th class="l" rowspan="2">Metabolite</th>'
                f'<th colspan="{len(ts_f)}" style="text-align:center">Crucian carp</th>'
                f'<th colspan="{len(ts_r)}" style="text-align:center">Rat</th></tr><tr>'
                + "".join(f"<th>{t}</th>" for t in ts_f) + "".join(f"<th>{t}</th>" for t in ts_r) + "</tr></thead>")
        return f'<table style="font-size:6.6pt">{head}<tbody>{"".join(rows)}</tbody></table>'

    # sample table
    srows = []
    for t in TIMES:
        cells = [f"<td class='l'>{'Control' if t == 'Ctr' else TIME_LABEL[t]}</td>"]
        for sp in ["Fish", "Rat"]:
            cs = x.cols(sp, t)
            ids = []
            for c in cs:
                tag = c
                if c in FR_CORRECTION[tissue_name]:
                    tag = f"<b>{c}</b>†"
                ids.append(tag)
            cells.append(f"<td>{len(cs) or '–'}</td><td class='l'>{', '.join(ids) or 'not sampled'}</td>")
        srows.append("<tr>" + "".join(cells) + "</tr>")
    sample_table = ('<table><thead><tr><th class="l">Ischemia time</th><th>n</th><th class="l">Crucian carp '
                    '(F/M = sex)</th><th>n</th><th class="l">Rat (male)</th></tr></thead><tbody>'
                    + "".join(srows) + "</tbody></table>")

    # ------------------------------------------------------------------ auto findings
    def aec(sp, t):
        return mean(sp, t, "Adenylate Energy Charge")

    auto = []
    auto.append(f"Adenylate energy charge at control: crucian carp {aec('Fish', 'Ctr'):.2f}, rat {aec('Rat', 'Ctr'):.2f}. "
                f"After 30 min: carp {aec('Fish', '30min'):.2f}, rat {aec('Rat', '30min'):.2f}; after 360 min: carp "
                f"{aec('Fish', '360min'):.2f}, rat {aec('Rat', '360min'):.2f}.")
    for k in ["Lactic acid", "Succinic acid"]:
        f30 = R["response"][("Fish", "30min")].loc[k, "effect"]; r30 = R["response"][("Rat", "30min")].loc[k, "effect"]
        auto.append(f"{label(k)} at 30 min changed {2**f30:.1f}-fold in crucian carp and {2**r30:.1f}-fold in rat "
                    f"(relative to each species' own control).")
    auto.append(f"{base_diff} of {len(comp)} metabolites differ between species already at control (q < 0.05) – "
                f"absolute levels are therefore compared alongside the response relative to each species' own control.")
    auto.append(f"{n_inter} metabolites show a significantly different time-course between species "
                f"(species × time interaction, q < 0.05); these are marked with * in figures and tables.")
    curated = findings.SPECIES.get(tissue_name, [])

    today = dt.date.today().strftime("%d %B %Y")
    body = []
    body.append(f"""
<section class="cover">
  <div class="kicker">Comparative metabolomics · {esc(tissue_name)}</div>
  <h1>Crucian carp vs rat:<br>metabolic response of the {tissue_name.lower()} to ischemia</h1>
  <div class="sub">Re-analysis of HMT data set {x.hmt_id} with species-resolved statistics</div>
  <p><span class="swatch" style="background:{F.C1}"></span>Crucian carp (<i>Carassius carassius</i>) &nbsp;
     <span class="swatch" style="background:{F.C2}"></span>Rat</p>
  <ol class="toc">
    <li>1. Summary</li><li>2. Study design and samples</li><li>3. Methods</li>
    <li>4. Overall metabolome response</li><li>5. Hallmark ischemia metabolites</li>
    <li>6. Pathway maps</li><li>7. Metabolic parameters</li><li>8. Statistical tables</li>
    <li>Appendix. Group means for all metabolites</li>
  </ol>
  <div class="meta">Source data: <i>{esc(x.hmt_id)}_{tissue_name}_MetabolomeAnalysisData.xlsx</i> (Human Metabolome
  Technologies, CE-TOFMS absolute quantification). Report generated {today}. Rats: males only; female rats (control and
  30 min) are analysed in the separate male-vs-female report.</div>
</section>""")

    kpi = (f'<div class="kpis"><div class="kpi"><div class="v">{len(comp)}</div><div class="l">metabolites tested '
           f'(detected in ≥3 samples)</div></div><div class="kpi"><div class="v">{base_diff}</div><div class="l">differ '
           f'between species at control</div></div><div class="kpi"><div class="v">{n_inter}</div><div class="l">'
           f'species-specific time-course (interaction)</div></div><div class="kpi"><div class="v">'
           f'{int(counts.loc["30 min", "Rat_up"] + counts.loc["30 min", "Rat_dn"])} / '
           f'{int(counts.loc["30 min", "Fish_up"] + counts.loc["30 min", "Fish_dn"])}</div><div class="l">changed '
           f'after 30 min: rat / carp</div></div></div>')
    body.append(f"""<section class="page summary"><h2>1. Summary</h2>{kpi}
<h3>Key findings</h3><ul class="findings">{''.join(f'<li>{c}</li>' for c in curated)}</ul>
<h3>Numbers behind the findings</h3><ul class="findings">{''.join(f'<li>{a}</li>' for a in auto)}</ul>
<h3>Metabolites changed relative to each species' own control (q &lt; 0.05)</h3>
<table><thead><tr><th class="l">Time</th><th>Carp ↑</th><th>Carp ↓</th><th>Rat ↑</th><th>Rat ↓</th>
<th>Carp vs rat<br>(absolute)</th><th>Response differs<br>(interaction)</th></tr></thead><tbody>
{''.join(f"<tr><td class='l'>{i}</td>" + ''.join(f"<td>{'–' if pd.isna(v) else int(v)}</td>" for v in r[['Fish_up','Fish_dn','Rat_up','Rat_dn','species','inter']].values) + '</tr>' for i, r in counts.iterrows())}
</tbody></table>
<p class="small">Rats were not sampled at 18 h; the 18 h crucian carp group is shown in all figures and tested against
the carp control, but is not part of the between-species statistics.</p></section>""")

    body.append(f"""<section class="page"><h2>2. Study design and samples</h2>
<p>Crucian carp and rats were exposed to increasing ischemia times before the {tissue_name.lower()} was sampled and
frozen. The same crucian carp individuals provided heart, brain and liver at one time point. Each rat contributed
different organs to different time groups (e.g. rat 13: heart control, liver and brain 30 min).</p>
<p><b>Ischemia protocol.</b> All rat organs and the crucian carp brain were removed from the animal and kept in buffer
at room temperature (~20–23 °C) for the ischemia period. Ischemia was deliberately run at ~20 °C so that both species
were compared at the same tissue temperature (and thus a comparable temperature effect on metabolic rate), rather than
at their very different body temperatures. To bring large organs such as the rat liver quickly to this temperature,
organs were cut into 3–8 slices, in approximately the same way for every animal and in the separate swelling
experiment. The crucian carp heart and liver stayed inside the dead fish at
room temperature for the whole period (in-situ ischemia).</p>
{sample_table}
<div class="note"><b>Male rats only.</b> Female rats (FR1–FR6) were sampled only at control and 30 min and are
excluded from this report so that all rat time points come from males of the same cohort (3 per group, 2 at 360 min).
They are analysed, together with the males, in the separate male-vs-female report, which found no sex difference in the
ischemic response. Note that HMT's own group means and figures include the females.</div>
<p class="small">Tissue weights used for extraction ranged {x.samples.weight_mg.min():.0f}–{x.samples.weight_mg.max():.0f} mg.
Sample IDs are those used by HMT.</p></section>""")

    body.append(f"""<section class="page"><h2>3. Methods</h2>
<h3>Data</h3><p>Absolute concentrations (nmol/g tissue) of {len(x.compounds)} targeted metabolites and
{len(x.params)} derived metabolic parameters (energy charge, redox ratios, amino-acid sums, etc.) were taken from the
HMT workbook ({esc(x.hmt_id)}, sheet <i>Table4 Absolute Concentration</i>). Values reported as N.D. (not detected)
were treated as missing for means and figures; for statistics they were replaced by half of the lowest detected value
of that metabolite in the {tissue_name.lower()} data set. Metabolites detected in fewer than three samples were not
tested ({len(x.compounds) - len(comp)} metabolites).</p>
{species_section(tissue_name)}
<h3>Statistics</h3><p>All tests used log2-transformed concentrations (energy charges were analysed on their linear 0–1
scale). Three complementary questions were asked for every metabolite:</p>
<ol>
<li><b>Does each species change relative to its own control?</b> Welch's t-test of each ischemia time vs control,
separately for crucian carp and rat. The effect is the log2 fold change (log2FC).</li>
<li><b>Are absolute levels different between species?</b> Welch's t-test crucian carp vs rat at each shared time point
(control, 30, 60, 180, 360 min). Effect: log2(carp / rat).</li>
<li><b>Does the response differ between species?</b> A linear model log2(concentration) ~ species × time on the shared
time points. The interaction contrast at each time point equals log2FC<sub>carp</sub> − log2FC<sub>rat</sub>; the
overall species × time F-test (4 df) asks whether the time-courses differ anywhere. This is the main
fish-vs-rat comparison and is not affected by constant differences in absolute level between the species.</li>
</ol>
<p>P-values were adjusted with the Benjamini–Hochberg false discovery rate (q) within each comparison, separately for
metabolites and derived parameters; q &lt; 0.05 is called significant. PCA used autoscaled log2 values. In panel B of
the PCA each species was first centred on its own control mean, so the plot shows response trajectories rather than
species differences.</p>
<div class="note"><b>Limitations.</b> Group sizes are small (2–3 crucian carp and 2–3 male rats per time point), so only large effects reach significance and q-values should be read as screening
evidence. Where a metabolite is N.D. in one group, significance can partly reflect detection rather than a measured
concentration difference (hollow markers in figures). Absolute levels can differ between species for reasons unrelated
to ischemia (body temperature, tissue composition, diet), which is why the response-based comparison (3) is
emphasised. For heart and liver, crucian carp tissue underwent in-situ ischemia (inside the fish) whereas rat tissue was
kept in buffer, so these comparisons combine a species difference with a difference in ischemia mode; the brain was
kept in buffer in both species. Cutting creates damaged surfaces that can take up buffer, so part of the measured
swelling may come from the cutting itself; because cutting was done the same way throughout, this mainly affects
comparisons between organs rather than between time points.</div>
<h3>Figure conventions</h3>
<p><span class="swatch" style="background:{F.C1}"></span>crucian carp &nbsp;
<span class="swatch" style="background:{F.C2}"></span>rat. Time axis: C = control, 30–360 = minutes of ischemia,
18h = 18 hours (carp only). Points and lines show mean ± SD (large panels also show individual animals). Hollow
markers: at least one animal in the group was N.D.; the mean is of detected values only. <b>Bold name *</b> =
significant species × time interaction (q &lt; 0.05).</p></section>""")

    body.append(f"""<section class="page"><h2>4. Overall metabolome response</h2>
{figure(pca_svg, "<b>Fig. 1. Principal component analysis.</b> A: all samples on autoscaled log2 concentrations; "
        "small dots are animals, large dots group means connected in time order. B: the same data after centring "
        "each species on its own control mean, which removes baseline species differences and shows how far and in "
        f"which direction each species moves with ischemia time. {len(comp)} metabolites.")}
{figure(counts_svg, "<b>Fig. 2. Number of significantly changed metabolites</b> (q &lt; 0.05) at each ischemia time: "
        "change vs own control in each species, absolute species difference, and species × time interaction contrast.")}
</section>""")
    body.append(f"""<section class="page">{figure(heat_svg, "<b>Fig. 3. Response heat map.</b> log2 fold change vs "
        "own control for each species and time point (red = increase, blue = decrease, clipped at ±4). Rows clustered by "
        "response pattern. Dots: q &lt; 0.05 vs own control. White = not tested.")}</section>""")
    body.append(f"""<section class="page">{figure(ratio_svg, "<b>Fig. 4. Absolute species difference.</b> "
        "log2(crucian carp / rat) at each shared time point (red = higher in carp, blue = higher in rat, clipped at ±5). "
        "Dots: q &lt; 0.05.")}</section>""")

    body.append(f"""<section class="page"><h2>5. Hallmark ischemia metabolites</h2>
<p>Metabolites that report on energy status (adenylates, phosphocreatine), anaerobic end-products (lactate, succinate,
alanine, GABA, glycerol 3-phosphate) and purine breakdown (inosine, hypoxanthine).</p>
{figure(hall_svg, "<b>Fig. 5. Hallmark metabolites</b> (nmol/g, mean ± SD with individual animals). "
        "* species × time interaction q &lt; 0.05.")}</section>""")
    body.append(f"""<section class="page"><h3>Hallmark metabolites – response relative to own control (log2FC)</h3>
{resp_table(hall + [p for p in ['Adenylate Energy Charge', 'Lactate/Pyruvate', 'NADH/NAD+', 'Total Adenylate'] if p in par])}
<p class="small">Cells: log2 fold change vs the species' own control (energy charge: absolute difference). Bold with * :
q &lt; 0.05. Interaction q: overall species × time F-test.</p></section>""")

    for i, (title, s) in enumerate(map_svgs):
        intro = ("<h2>6. Pathway maps</h2><p>Each box shows the time-course of one metabolite (nmol/g, mean ± SD) in "
                 "crucian carp and rat on a simplified pathway. Grey italic names are intermediates not measured. "
                 "Dashed arrows: several reactions. Bold name with * : species × time interaction q &lt; 0.05.</p>"
                 if i == 0 else "")
        body.append(f'<section class="page">{intro}<h3>6.{i+1} {esc(title)}</h3>{figure(s, "")}</section>')

    for i, s in enumerate(param_pages):
        intro = ("<h2>7. Metabolic parameters</h2><p>Derived parameters calculated by HMT from the measured "
                 "concentrations (e.g. adenylate energy charge = (ATP + ½ADP)/(ATP + ADP + AMP)).</p>" if i == 0 else "")
        body.append(f'<section class="page">{intro}{figure(s, f"<b>Fig. {6+i}. Metabolic parameters</b> ({i+1}/{len(param_pages)}).")}</section>')
    body.append(f"""<section class="page"><h3>Metabolic parameters – response relative to own control</h3>
{resp_table(par)}</section>""")

    body.append(f"""<section class="page"><h2>8. Statistical tables</h2>
<h3>8.1 Metabolites whose ischemia response differs most between species</h3>
<p class="small">Top {len(top_inter)} by species × time interaction (F-test). Cells: log2FC vs own control.</p>
{resp_table(top_inter)}</section>""")
    half = (len(base_sig) + 1) // 2
    body.append(f"""<section class="page"><h3>8.2 Absolute species differences at control (q &lt; 0.05, {len(base_sig)} metabolites)</h3>
<div class="two">{species_table("Ctr", base_sig)}</div></section>""")
    body.append(f"""<section class="page"><h2>Appendix. Group means (nmol/g)</h2>
<p class="small">Mean of detected values; N.D. = not detected in any animal of the group. Time: C = control, minutes,
18h = 18 hours.</p>{means_table(list(x.compounds))}
<h3>Derived parameters</h3>{means_table(list(x.params))}</section>""")

    # Excel export of all statistics
    xl = {}
    m = S.copy(); m.columns = [f"{a} {b}" for a, b in m.columns]
    xl["group_means"] = m
    rows = {}
    for t in SHARED_TIMES:
        d = R["species"][t]; rows[f"carp_vs_rat_{t}_log2"] = d.effect; rows[f"carp_vs_rat_{t}_q"] = d.q
    xl["species_difference"] = pd.DataFrame(rows)
    rows = {}
    for (sp, t), d in R["response"].items():
        rows[f"{sp}_{t}_log2FC"] = d.effect; rows[f"{sp}_{t}_q"] = d.q
    xl["response_vs_control"] = pd.DataFrame(rows)
    rows = {"interaction_F_p": R["interaction_F"].p, "interaction_F_q": R["interaction_F"].q}
    for t, d in R["interaction"].items():
        rows[f"dlog2FC_carp_minus_rat_{t}"] = d.effect; rows[f"q_{t}"] = d.q
    xl["species_x_time"] = pd.DataFrame(rows)
    return doc(f"Crucian carp vs rat – {tissue_name}", "".join(body)), xl
