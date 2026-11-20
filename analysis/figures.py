"""Matplotlib figures rendered to inline SVG."""
import io, re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyBboxPatch
from scipy.cluster.hierarchy import linkage, leaves_list

C1, C2 = "#2a78d6", "#eb6834"          # categorical slots 1-2 (validated palette)
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df"
DIVERGING = LinearSegmentedColormap.from_list("bgr", ["#1c5cab", "#86b6ef", "#f0efec", "#f0a3a2", "#b8302f"])
plt.rcParams.update({
    "font.family": ["Inter", "DejaVu Sans"], "font.size": 7, "axes.linewidth": 0.6, "axes.edgecolor": MUTED,
    "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
    "axes.spines.top": False, "axes.spines.right": False, "svg.fonttype": "none",
    "xtick.major.width": 0.5, "ytick.major.width": 0.5, "xtick.major.size": 2, "ytick.major.size": 2,
})


def svg(fig):
    buf = io.StringIO()
    fig.savefig(buf, format="svg", bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    s = buf.getvalue()
    s = s[s.index("<svg"):]
    return re.sub(r'(<svg[^>]*?) width="[^"]*" height="[^"]*"', r'\1 width="100%"', s, count=1)


def compact(v, _=None):
    a = abs(v)
    if a >= 1e4: return f"{v/1e3:.0f}k"
    if a >= 1e3: return f"{v/1e3:.1f}k".replace(".0k", "k")
    if a >= 10: return f"{v:.0f}"
    if a >= 1: return f"{v:.1f}".rstrip("0").rstrip(".")
    if a == 0: return "0"
    return f"{v:.2g}"


# ------------------------------------------------------------------ time-course panels
def draw_series(ax, data, series, xlabels, small=True, points=False):
    """series: list of dict(name,color,groups=[(x, cols)]). data: Series of one metabolite (raw)."""
    ymax = 0
    for s in series:
        xs, ms, sds, hollow = [], [], [], []
        for x, cols in s["groups"]:
            v = data[cols]
            d = v.dropna()
            if len(d) == 0:
                continue
            xs.append(x); ms.append(d.mean()); sds.append(d.std() if len(d) > 1 else 0)
            hollow.append(len(d) < len(v))
            ymax = max(ymax, (d.max() if points else d.mean() + (d.std() if len(d) > 1 else 0)))
            if points:
                jit = (np.arange(len(d)) - (len(d) - 1) / 2) * 0.06 + s.get("offset", 0)
                ax.scatter(x + jit, d.values, s=7, color=s["color"], alpha=0.55, lw=0, zorder=2)
        if not xs:
            continue
        off = s.get("offset", 0)
        xs = np.array(xs) + off
        lw = 1.0 if small else 1.6
        ax.errorbar(xs, ms, yerr=sds, color=s["color"], lw=lw, elinewidth=0.5 if small else 0.8,
                    capsize=0, zorder=3)
        for x, m, h in zip(xs, ms, hollow):
            ax.plot(x, m, "o", ms=2.6 if small else 4.2, mfc="white" if h else s["color"],
                    mec=s["color"], mew=0.7 if small else 1, zorder=4)
    ax.set_xticks(range(len(xlabels)))
    ax.set_xticklabels(xlabels)
    ax.set_xlim(-0.5, len(xlabels) - 0.5)
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(compact))
    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(3 if small else 4))
    ax.set_ylim(0, ymax * 1.12 if ymax > 0 else 1)
    ax.grid(axis="y", color=GRID, lw=0.4)
    ax.set_axisbelow(True)
    if ymax == 0:
        ax.text(0.5, 0.5, "N.D.", transform=ax.transAxes, ha="center", va="center", color=MUTED, fontsize=6)
        ax.set_yticks([])


def legend_handles(series):
    from matplotlib.lines import Line2D
    return [Line2D([], [], color=s["color"], marker="o", ms=4, lw=1.4, label=s["name"]) for s in series]


def pathway_map(pmap, tissue, series, xlabels, marks, label_of):
    """marks: dict compound -> annotation string (e.g. '*')."""
    cols, rows = pmap["cols"], pmap["rows"]
    W = 7.3
    cw = W / cols
    ch = 0.98
    H = rows * ch + 0.45
    fig = plt.figure(figsize=(W, H))
    bw, bh = cw * 0.86 / W, ch * 0.80 / H        # box size in fig fraction

    def centre(c, r):
        return ((c + 0.5) * cw / W, 1 - (0.45 + (r + 0.5) * ch) / H)

    boxes = {}
    for k, (c, r) in pmap["nodes"].items():
        x, y = centre(c, r)
        if k.startswith("~"):
            boxes[k] = (x, y, bw * 0.75, bh * 0.28)
        else:
            boxes[k] = (x, y, bw, bh)
    # edges under everything
    ov = fig.add_axes([0, 0, 1, 1], zorder=0); ov.set_axis_off(); ov.set_xlim(0, 1); ov.set_ylim(0, 1)
    for a, b, style in pmap["edges"]:
        (x1, y1, w1, h1), (x2, y2, w2, h2) = boxes[a], boxes[b]
        def clip(x, y, w, h, dx, dy):
            tx = (w / 2 + 0.006) / abs(dx) if dx else np.inf
            ty = (h / 2 + 0.006) / abs(dy) if dy else np.inf
            return min(tx, ty)
        dx, dy = x2 - x1, y2 - y1
        t1 = clip(x1, y1, w1, h1, dx, dy); t2 = clip(x2, y2, w2, h2, -dx, -dy)
        if t1 + t2 >= 1:
            continue
        ov.annotate("", xy=(x2 - dx * t2, y2 - dy * t2), xytext=(x1 + dx * t1, y1 + dy * t1),
                    arrowprops=dict(arrowstyle="-|>,head_length=0.35,head_width=0.18", color=MUTED, lw=0.7,
                                    ls=(0, (3, 2)) if style == "--" else "-", shrinkA=0, shrinkB=0))
    for k, (x, y, w, h) in boxes.items():
        if k.startswith("~"):
            ov.text(x, y, k[1:], ha="center", va="center", fontsize=5.6, style="italic", color=INK2,
                    bbox=dict(boxstyle="round,pad=0.25", fc="white", ec=GRID, lw=0.5))
            continue
        ov.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0.004,rounding_size=0.006",
                                    fc="#f7f7f5", ec="#d8d7d2", lw=0.5, zorder=1))
        ax = fig.add_axes([x - w / 2 + w * 0.2, y - h / 2 + h * 0.17, w * 0.76, h * 0.62], zorder=2)
        ax.patch.set_alpha(0)
        if k in tissue.data.index:
            draw_series(ax, tissue.data.loc[k], series, xlabels)
        ax.tick_params(labelsize=4.3, pad=1, length=1.5)
        title = label_of(k) + (" " + marks[k] if marks.get(k) else "")
        ov.text(x - w / 2 + w * 0.05, y + h / 2 - h * 0.05, title, ha="left", va="top", fontsize=5.8,
                fontweight="bold" if marks.get(k) else "normal", color=INK, zorder=5,
                clip_on=False)
    fig.legend(handles=legend_handles(series), loc="upper left", bbox_to_anchor=(0.0, 1.0), ncol=len(series),
               frameon=False, fontsize=7)
    return svg(fig)


def panel_grid(tissue, names, series, xlabels, label_of, marks=None, ncols=4, ylabel_of=None, height=1.55):
    marks = marks or {}
    n = len(names); nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(7.3, height * nrows + 0.3), squeeze=False)
    for i, ax in enumerate(axes.flat):
        if i >= n:
            ax.set_axis_off(); continue
        k = names[i]
        draw_series(ax, tissue.data.loc[k], series, xlabels, small=False, points=True)
        ax.set_title(label_of(k) + (" " + marks[k] if marks.get(k) else ""), fontsize=7.5, loc="left",
                     fontweight="bold" if marks.get(k) else "normal")
        ax.tick_params(labelsize=6)
        if ylabel_of:
            ax.set_ylabel(ylabel_of(k), fontsize=6)
    fig.legend(handles=legend_handles(series), loc="upper left", bbox_to_anchor=(0, 1.0 + 0.25 / (height * nrows)),
               ncol=len(series), frameon=False, fontsize=7.5)
    fig.tight_layout(h_pad=1.0, w_pad=0.8)
    return svg(fig)


# ------------------------------------------------------------------ PCA
def pca(X):
    """X: samples x features (already transformed). Autoscale and PCA."""
    Z = (X - X.mean()) / X.std().replace(0, 1)
    Z = Z.loc[:, Z.notna().all()]
    U, S, Vt = np.linalg.svd(Z.values, full_matrices=False)
    ev = S ** 2 / (S ** 2).sum()
    return pd.DataFrame(U[:, :3] * S[:3], index=X.index, columns=["PC1", "PC2", "PC3"]), ev[:3]


def pca_panel(ax, scores, ev, sample_info, series, title):
    """series: list of dict(name,color,groups=[(label, cols)]) - draws points + centroid trajectory."""
    for s in series:
        cx, cy, labs = [], [], []
        for lab, cols in s["groups"]:
            cols = [c for c in cols if c in scores.index]
            if not cols: continue
            ax.scatter(scores.loc[cols, "PC1"], scores.loc[cols, "PC2"], s=12, color=s["color"], alpha=0.45,
                       lw=0, zorder=2)
            cx.append(scores.loc[cols, "PC1"].mean()); cy.append(scores.loc[cols, "PC2"].mean()); labs.append(lab)
        ax.plot(cx, cy, "-", color=s["color"], lw=1.4, zorder=3)
        ax.scatter(cx, cy, s=34, color=s["color"], edgecolor="white", lw=1.2, zorder=4)
        for x, y, l in zip(cx, cy, labs):
            ax.annotate(l, (x, y), xytext=(4, 3), textcoords="offset points", fontsize=6, color=INK2, zorder=5)
    ax.axhline(0, color=GRID, lw=0.5, zorder=0); ax.axvline(0, color=GRID, lw=0.5, zorder=0)
    ax.set_xlabel(f"PC1 ({ev[0]*100:.0f}%)"); ax.set_ylabel(f"PC2 ({ev[1]*100:.0f}%)")
    ax.set_title(title, fontsize=8, loc="left")


# ------------------------------------------------------------------ heat map
def heatmap(M, sig, col_groups, label_of, title_note="", clip=4, row_h=None, cbar_label="log2 fold change"):
    """M: metabolites x columns (log2 values). sig: same shape bool. col_groups: list of (name, [cols])."""
    M = M.copy()
    filled = M.fillna(0).clip(-clip, clip)
    order = leaves_list(linkage(filled.values, "average", metric="euclidean")) if len(M) > 2 else np.arange(len(M))
    M, S = M.iloc[order], sig.iloc[order]
    allcols = [c for _, cs in col_groups for c in cs]
    gap = 0.6
    xs, x = [], 0
    for gi, (_, cs) in enumerate(col_groups):
        for c in cs:
            xs.append(x); x += 1
        x += gap
    nrow = len(M)
    if row_h is None:
        row_h = min(0.075, 7.4 / max(nrow, 1))
    fig = plt.figure(figsize=(7.0, row_h * nrow + 1.2))
    ax = fig.add_axes([0.26, 0.6 / (row_h * nrow + 1.2), 0.5, (row_h * nrow) / (row_h * nrow + 1.2)])
    for xi, c in zip(xs, allcols):
        vals = M[c].values.astype(float)
        col = DIVERGING((np.clip(np.nan_to_num(vals), -clip, clip) + clip) / (2 * clip))
        for r in range(nrow):
            fc = col[r] if not np.isnan(vals[r]) else "#ffffff"
            ax.add_patch(plt.Rectangle((xi, nrow - r - 1), 0.94, 0.9, fc=fc, ec="none"))
            if S[c].iloc[r]:
                ax.plot(xi + 0.47, nrow - r - 0.55, "o", ms=1.4, color=INK, zorder=3)
    ax.set_xlim(-0.2, x - gap + 0.2); ax.set_ylim(0, nrow)
    ax.set_yticks(np.arange(nrow) + 0.45)
    ax.set_yticklabels([label_of(m) for m in M.index][::-1], fontsize=5.2)
    ax.set_xticks(np.array(xs) + 0.47)
    ax.set_xticklabels([c.split("|")[-1] for c in allcols], fontsize=6, rotation=0)
    ax.xaxis.tick_top()
    ax.tick_params(length=0, pad=1)
    for sp in ax.spines.values(): sp.set_visible(False)
    i = 0
    for name, cs in col_groups:
        ax.text(xs[i] + (len(cs) * 1 - 0.06) / 2, nrow + 1.6, name, ha="center", fontsize=7.5, fontweight="bold")
        i += len(cs)
    cax = fig.add_axes([0.86, 0.6 / (row_h * nrow + 1.2) + (row_h * nrow) / (row_h * nrow + 1.2) * 0.75, 0.025, 0.1])
    sm = plt.cm.ScalarMappable(cmap=DIVERGING, norm=plt.Normalize(-clip, clip))
    cb = fig.colorbar(sm, cax=cax); cb.outline.set_visible(False); cb.ax.tick_params(labelsize=5.5)
    cb.set_label(cbar_label, fontsize=6)
    return svg(fig)


def count_bars(table, colors, ylabel="metabolites (q < 0.05)"):
    """table: DataFrame index=x categories, columns=series."""
    fig, ax = plt.subplots(figsize=(7.0, 2.1))
    n = len(table.columns); w = 0.8 / n
    for i, c in enumerate(table.columns):
        xs = np.arange(len(table)) + (i - (n - 1) / 2) * w
        ax.bar(xs, table[c].values, width=w * 0.9, color=colors[i], label=c, zorder=2)
        for xx, v in zip(xs, table[c].values):
            if v > 0: ax.text(xx, v, f"{int(v)}", ha="center", va="bottom", fontsize=5.5, color=INK2)
    ax.set_xticks(range(len(table))); ax.set_xticklabels(table.index)
    ax.set_ylabel(ylabel); ax.grid(axis="y", color=GRID, lw=0.4); ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=6.5, ncol=n, loc="upper left", bbox_to_anchor=(0, 1.18))
    return svg(fig)
