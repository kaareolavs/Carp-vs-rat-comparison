"""Statistics shared by the reports."""
import numpy as np
import pandas as pd
from scipy import stats
from data import SHARED_TIMES, TIMES

MIN_DETECTED = 3   # metabolite must be detected in >= this many samples of the tissue to be tested


def transformed(tissue, rows):
    """Return (Y, is_log) where Y is the matrix used for statistics.
    Compounds: log2 with N.D. -> half min. Params: log2 (energy charges kept linear), N.A. dropped."""
    out, islog = [], {}
    comp = [r for r in rows if r in tissue.compounds]
    par = [r for r in rows if r in tissue.params]
    if comp:
        out.append(tissue.log2_imputed(comp))
        islog.update({r: True for r in comp})
    if par:
        P = tissue.data.loc[par].copy()
        for r in par:
            if "Energy Charge" in r:
                islog[r] = False
            else:
                P.loc[r] = np.log2(P.loc[r].where(P.loc[r] > 0))
                islog[r] = True
        out.append(P)
    return pd.concat(out).loc[rows], pd.Series(islog)


def detected(tissue, rows):
    return tissue.data.loc[rows].notna()


def bh(p):
    p = np.asarray(p, float)
    q = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    if ok.sum():
        pv = p[ok]; n = len(pv); o = np.argsort(pv)
        r = pv[o] * n / np.arange(1, n + 1)
        r = np.minimum.accumulate(r[::-1])[::-1]
        qq = np.empty(n); qq[o] = np.minimum(r, 1)
        q[ok] = qq
    return q


def welch(Y, D, a, b, min_det=2):
    """Welch t-test per row of Y between column sets a (numerator) and b.
    Returns DataFrame effect (mean a - mean b, log2 ratio for log rows) and p."""
    eff, p = [], []
    for r in Y.index:
        ya, yb = Y.loc[r, a].dropna(), Y.loc[r, b].dropna()
        da, db = D.loc[r, a].sum(), D.loc[r, b].sum()
        if len(ya) < 2 or len(yb) < 2 or max(da, db) < min_det:
            eff.append(np.nan); p.append(np.nan); continue
        e = ya.mean() - yb.mean()
        if ya.var() + yb.var() == 0:
            eff.append(e); p.append(np.nan if e == 0 else 0.0); continue
        eff.append(e); p.append(stats.ttest_ind(ya, yb, equal_var=False).pvalue)
    return pd.DataFrame({"effect": eff, "p": p}, index=Y.index)


def ols_contrasts(Y, design, contrasts, ftest=None):
    """Per-row OLS. design: DataFrame (samples x terms) incl. intercept.
    contrasts: dict name -> vector over design columns. ftest: list of column names tested jointly.
    Returns dict name -> DataFrame(effect, se, t, p) and F-test p series."""
    X = design.values.astype(float)
    cols = list(design.columns)
    res = {k: [] for k in contrasts}
    fp = []
    for r in Y.index:
        y = Y.loc[r, design.index].values.astype(float)
        ok = ~np.isnan(y)
        Xo, yo = X[ok], y[ok]
        n, k = Xo.shape
        if n - k < 2 or np.linalg.matrix_rank(Xo) < k:
            for c in contrasts: res[c].append((np.nan,) * 4)
            fp.append(np.nan); continue
        XtXi = np.linalg.inv(Xo.T @ Xo)
        beta = XtXi @ Xo.T @ yo
        rss = ((yo - Xo @ beta) ** 2).sum(); dfr = n - k; s2 = rss / dfr
        for c, v in contrasts.items():
            v = np.asarray(v, float); e = v @ beta; se = np.sqrt(v @ XtXi @ v * s2)
            t = e / se if se > 0 else np.nan
            res[c].append((e, se, t, 2 * stats.t.sf(abs(t), dfr) if se > 0 else np.nan))
        if ftest:
            keep = [i for i, c in enumerate(cols) if c not in ftest]
            Xr = Xo[:, keep]
            br = np.linalg.lstsq(Xr, yo, rcond=None)[0]
            rss0 = ((yo - Xr @ br) ** 2).sum(); q = len(ftest)
            F = ((rss0 - rss) / q) / s2 if s2 > 0 else np.nan
            fp.append(stats.f.sf(F, q, dfr) if s2 > 0 else np.nan)
        else:
            fp.append(np.nan)
    out = {c: pd.DataFrame(v, index=Y.index, columns=["effect", "se", "t", "p"]) for c, v in res.items()}
    return out, pd.Series(fp, index=Y.index)


def summarise(tissue, rows, groups):
    """mean, sd, n detected on the raw (nmol/g or unitless) scale for each (label -> columns)."""
    D = tissue.data.loc[rows]
    out = {}
    for lab, cols in groups.items():
        out[(lab, "mean")] = D[cols].mean(axis=1)
        out[(lab, "sd")] = D[cols].std(axis=1)
        out[(lab, "n")] = D[cols].notna().sum(axis=1)
    return pd.DataFrame(out)


# ---------------------------------------------------------------- fish vs rat
def fish_vs_rat(tissue):
    rows = list(tissue.compounds) + list(tissue.params)
    Y, islog = transformed(tissue, rows)
    D = detected(tissue, rows)
    keep = D.sum(axis=1) >= MIN_DETECTED
    groups = {}
    for sp in ["Fish", "Rat"]:
        for t in TIMES:
            c = tissue.cols(sp, t)
            if c: groups[(sp, t)] = c
    summ = summarise(tissue, rows, {f"{sp}|{t}": c for (sp, t), c in groups.items()})
    R = {"Y": Y, "islog": islog, "detected": D, "keep": keep, "groups": groups, "summary": summ}
    comp_rows = [r for r in rows if r in tissue.compounds]
    par_rows = [r for r in rows if r in tissue.params]

    def adjust(df):
        df["q"] = np.nan
        for rr in (comp_rows, par_rows):
            idx = [r for r in rr if keep[r]]
            df.loc[idx, "q"] = bh(df.loc[idx, "p"])
        df.loc[~keep, ["effect", "p", "q"]] = np.nan
        return df

    # 1. fish vs rat at each shared time point
    R["species"] = {t: adjust(welch(Y, D, groups[("Fish", t)], groups[("Rat", t)])) for t in SHARED_TIMES}
    # 2. response vs own control
    R["response"] = {}
    for sp in ["Fish", "Rat"]:
        for t in TIMES[1:]:
            if (sp, t) in groups:
                R["response"][(sp, t)] = adjust(welch(Y, D, groups[(sp, t)], groups[(sp, "Ctr")]))
    # 3. species x time model on shared time points
    cols = [c for t in SHARED_TIMES for sp in ["Fish", "Rat"] for c in groups[(sp, t)]]
    S = tissue.samples.loc[cols]
    des = pd.DataFrame(index=cols)
    des["Intercept"] = 1.0
    des["Fish"] = (S.species == "Fish").astype(float)
    for t in SHARED_TIMES[1:]:
        des[t] = (S.group == t).astype(float)
    for t in SHARED_TIMES[1:]:
        des["Fish:" + t] = des["Fish"] * des[t]
    contr = {}
    for t in SHARED_TIMES[1:]:
        v = np.zeros(des.shape[1]); v[list(des.columns).index("Fish:" + t)] = 1
        contr[t] = v
    inter, fp = ols_contrasts(Y, des, contr, ftest=["Fish:" + t for t in SHARED_TIMES[1:]])
    R["interaction"] = {t: adjust(df) for t, df in inter.items()}
    fdf = adjust(pd.DataFrame({"effect": np.nan, "p": fp}))
    R["interaction_F"] = fdf
    return R


# ---------------------------------------------------------------- male vs female rats
def sex_comparison(tissue):
    rows = list(tissue.compounds) + list(tissue.params)
    Y, islog = transformed(tissue, rows)
    D = detected(tissue, rows)
    keep = D.sum(axis=1) >= MIN_DETECTED
    g = {(sx, t): tissue.cols("Rat", t, sx) for sx in ["M", "F"] for t in ["Ctr", "30min"]}
    summ = summarise(tissue, rows, {f"{sx}|{t}": c for (sx, t), c in g.items()})
    comp_rows = [r for r in rows if r in tissue.compounds]
    par_rows = [r for r in rows if r in tissue.params]

    def adjust(df):
        df["q"] = np.nan
        for rr in (comp_rows, par_rows):
            idx = [r for r in rr if keep[r]]
            df.loc[idx, "q"] = bh(df.loc[idx, "p"])
        df.loc[~keep, ["effect", "p", "q"]] = np.nan
        return df

    R = {"Y": Y, "islog": islog, "detected": D, "keep": keep, "groups": g, "summary": summ}
    R["sex_at"] = {t: adjust(welch(Y, D, g[("F", t)], g[("M", t)])) for t in ["Ctr", "30min"]}
    R["response"] = {sx: adjust(welch(Y, D, g[(sx, "30min")], g[(sx, "Ctr")])) for sx in ["M", "F"]}
    cols = [c for v in g.values() for c in v]
    S = tissue.samples.loc[cols]
    des = pd.DataFrame(index=cols)
    des["Intercept"] = 1.0
    des["F"] = (S.sex == "F").astype(float)
    des["30min"] = (S.group == "30min").astype(float)
    des["F:30min"] = des["F"] * des["30min"]
    out, _ = ols_contrasts(Y, des, {"sex": [0, 1, 0, 0.5], "time": [0, 0, 1, 0.5], "interaction": [0, 0, 0, 1]})
    R["anova"] = {k: adjust(v) for k, v in out.items()}
    return R
