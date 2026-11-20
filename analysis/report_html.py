"""HTML shell and small helpers shared by the report builders."""
import html as _html
import numpy as np
import matplotlib.colors as mcolors
from figures import DIVERGING

CSS = r"""
@page { size: A4; margin: 16mm 15mm 17mm 15mm; }
:root { --ink:#0b0b0b; --ink2:#52514e; --muted:#8a8984; --rule:#d8d7d2; --soft:#f4f4f1;
        --fish:#2a78d6; --rat:#eb6834; }
html { -webkit-print-color-adjust: exact; print-color-adjust: exact; }
body { font-family: Inter, "DejaVu Sans", sans-serif; color: var(--ink); background:#fff;
       font-size: 9.2pt; line-height: 1.45; margin: 0; }
h1 { font-size: 22pt; line-height: 1.15; margin: 0 0 6pt; font-weight: 700; letter-spacing: -0.01em; }
h2 { font-size: 14pt; margin: 0 0 8pt; padding-bottom: 3pt; border-bottom: 1.2pt solid var(--ink);
     break-after: avoid; }
h3 { font-size: 10.5pt; margin: 12pt 0 4pt; break-after: avoid; }
p { margin: 0 0 6pt; }
.page { break-before: page; }
.cover { height: 245mm; display: flex; flex-direction: column; }
.cover .kicker { text-transform: uppercase; letter-spacing: .12em; font-size: 8pt; color: var(--ink2);
                 margin-top: 40mm; margin-bottom: 8pt; }
.cover .sub { font-size: 12pt; color: var(--ink2); margin-bottom: 18pt; }
.cover .meta { margin-top: auto; font-size: 8.5pt; color: var(--ink2); border-top: .6pt solid var(--rule);
               padding-top: 8pt; }
.swatch { display:inline-block; width: 9pt; height: 9pt; border-radius: 2pt; vertical-align: -1pt; margin-right: 3pt; }
.toc { columns: 1; margin: 0; padding: 0; list-style: none; }
.toc li { padding: 2pt 0; border-bottom: .4pt dotted var(--rule); }
figure { margin: 6pt 0 10pt; break-inside: avoid; }
figure svg { display: block; width: 100%; height: auto; }
figcaption { font-size: 7.8pt; color: var(--ink2); margin-top: 4pt; line-height: 1.4; }
figcaption b { color: var(--ink); }
table { border-collapse: collapse; width: 100%; font-size: 7.4pt; margin: 4pt 0 10pt;
        font-variant-numeric: tabular-nums; }
th, td { padding: 1.6pt 4pt; text-align: right; border-bottom: .4pt solid #ecebe7; }
th { font-weight: 600; color: var(--ink2); border-bottom: .8pt solid var(--ink2); vertical-align: bottom; }
td:first-child, th:first-child { text-align: left; }
td.l, th.l { text-align: left; }
tr { break-inside: avoid; }
thead { display: table-header-group; }
.sig { font-weight: 700; }
.note { background: var(--soft); border-left: 2.5pt solid var(--ink2); padding: 6pt 9pt; margin: 6pt 0 10pt;
        font-size: 8.4pt; break-inside: avoid; }
.note.warn { border-left-color: #c98500; }
.kpis { display: grid; grid-template-columns: repeat(4, 1fr); gap: 6pt; margin: 6pt 0 12pt; }
.kpi { border: .6pt solid var(--rule); border-radius: 4pt; padding: 6pt 8pt; }
.kpi .v { font-size: 15pt; font-weight: 700; line-height: 1.1; }
.kpi .l { font-size: 7.4pt; color: var(--ink2); }
ul.findings { padding-left: 13pt; margin: 4pt 0 8pt; }
ul.findings li { margin-bottom: 4pt; }
.small { font-size: 7.6pt; color: var(--ink2); }
.two { columns: 2; column-gap: 14pt; }
.two table { break-inside: auto; }
.tag { display:inline-block; font-size: 6.8pt; padding: 0 4pt; border-radius: 6pt; background: var(--soft);
       color: var(--ink2); border: .4pt solid var(--rule); }
"""


def esc(s):
    return _html.escape(str(s))


def doc(title, body):
    return f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1"><title>{esc(title)}</title>
<style>{CSS}</style></head><body>{body}</body></html>"""


def fmt_q(q):
    if q is None or (isinstance(q, float) and np.isnan(q)):
        return "–"
    if q < 0.001:
        return "&lt;0.001"
    return f"{q:.3f}"


def fc_cell(v, q=None, clip=4, digits=2, linear=False):
    """table cell with a faint diverging background; bold when q < 0.05."""
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return '<td class="small">–</td>'
    scale = 0.5 if linear else clip
    rgba = DIVERGING((np.clip(v, -scale, scale) + scale) / (2 * scale))
    bg = mcolors.to_hex((rgba[0] * 0.55 + 0.45, rgba[1] * 0.55 + 0.45, rgba[2] * 0.55 + 0.45))
    sig = q is not None and not np.isnan(q) and q < 0.05
    txt = f"{v:+.{digits}f}"
    return f'<td style="background:{bg}"{" class=sig" if sig else ""}>{txt}{"*" if sig else ""}</td>'


def num(v, digits=None):
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return "N.D."
    a = abs(v)
    if digits is not None:
        return f"{v:.{digits}f}"
    if a >= 1000:
        return f"{v:,.0f}"
    if a >= 100:
        return f"{v:.0f}"
    if a >= 10:
        return f"{v:.1f}"
    if a >= 1:
        return f"{v:.2f}"
    return f"{v:.3g}"


def figure(svg, caption):
    return f"<figure>{svg}<figcaption>{caption}</figcaption></figure>"
