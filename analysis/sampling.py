"""Read Rat_sampling_times.xlsx into a tidy table (one row per rat x organ)."""
import datetime as dt
from pathlib import Path
import numpy as np
import pandas as pd
import openpyxl

ROOT = Path(__file__).resolve().parent.parent
FILE = ROOT / "Rat_sampling_times.xlsx"


def _secs(t):
    return t.hour * 3600 + t.minute * 60 + t.second


def load():
    ws = openpyxl.load_workbook(FILE)["Sampling times"]
    rows = {r[0].value: [c.value for c in r[1:]] for r in ws.iter_rows() if r[0].value}
    out = []
    for j, rat in enumerate(rows["Rat"]):
        if not rat:
            continue
        arrest = rows["Circulatory arrest (clock)"][j]
        base = dict(rat=rat, sex=rows["Sex"][j], weight=rows["Weight (g)"][j], day=rows["Sampling day"][j],
                    date=rows["Date"][j], temp=rows["Room / buffer temp (°C)"][j], comment=rows["Comments"][j])
        for o in ["Heart", "Liver", "Brain"]:
            clock = rows[f"{o} – removed (clock)"][j]
            s = rows[f"{o} – seconds after arrest"][j]
            if isinstance(s, str):           # formula -> compute
                s = _secs(clock) - _secs(arrest) if (clock and arrest) else None
            out.append(dict(base, organ=o, group=rows[f"{o} – ischemia group (min)"][j],
                            removed=clock, seconds=s, hmt=rows[f"{o} – HMT sample"][j]))
    df = pd.DataFrame(out)
    df["weight"] = pd.to_numeric(df["weight"], errors="coerce")
    df["date"] = pd.to_datetime(df["date"], errors="coerce")
    df["num"] = df.rat.str.extract(r"(\d+)").astype(int)
    return df


def for_hmt(organ):
    """Rows for samples sent to HMT for one organ, indexed by HMT sample id."""
    df = load()
    return df[(df.organ == organ) & df.hmt.notna()].set_index("hmt")
