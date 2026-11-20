"""Build Rat_sampling_times.xlsx: one column per rat, organs in rows, comments last.

Times are transcribed from 'Copy of Owerview Rat_Fish sampled_weight_sex_ischemic time_22062026.xlsx'
(sheets 'Sampling times male rats', 'Female rats', 'Room temperatures', 'Male rats', 'Sendt metabolomics').
"Circulatory arrest" = time of death (the heart was cut out). Ischemia groups for FR1/FR2 follow the
lab notes and the metabolite data (they were sampled in the reverse of the planned order).
"""
import datetime as dt
import re
from pathlib import Path
import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.comments import Comment
from openpyxl.utils import get_column_letter

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "Copy of Owerview Rat_Fish sampled_weight_sex_ischemic time_22062026.xlsx"
OUT = ROOT / "Rat_sampling_times.xlsx"
ORGANS = ["Heart", "Liver", "Brain", "Kidney 1", "Kidney 2"]

# ---------------------------------------------------------------- transcription
# per rat: day, arrest clock, organ -> clock time "hh:mm:ss" or ("s", seconds) when only seconds were noted
D1, D2, D3, D5 = "2026-04-07", "2026-04-08", "2026-04-09", "2026-05-07"
TEMP = {
    "Day 1": "room 21.5", "Day 2": "room 20.5", "Day 3": "room 20.9 / buffer 19.8 (21.2 at 14:48)",
    "Day 4": "buffer/water 22.3–22.6", "Day 5": "room 21.0 / buffer 21.0", "Day 6": "room 22.0",
    "18 June": "room 23 / buffer 22", "22 June": "room 19.5 / buffer 21.5",
}
R = {}
for n in (1, 2):
    R[n] = dict(day="Day 1 (pilot)", date=D1, tkey="Day 1",
                comment="Pilot day (rat number assumed). Direct procurement slow (~2 min for heart and liver); "
                        "not to be used for metabolomics or direct procurement (DR).")
R[3] = dict(day="Day 2", date=D2, tkey="Day 2", Liver=("s", 5), comment="Liver direct 5 s.")
R[4] = dict(day="Day 2", date=D2, tkey="Day 2", Heart=("s", 5), comment="Heart direct 5 s.")
R[5] = dict(day="Day 2", date=D2, tkey="Day 2", Brain=("s", 54), comment="Brain direct 54 s, then heart, then liver.")
R[6] = dict(day="Day 2", date=D2, tkey="Day 2",
            comment="Liver direct, then heart, then brain, then both kidneys. Heart not dried before freezing.")
R[7] = dict(day="Day 2", date=D2, tkey="Day 2", Heart=("s", 5), comment="Heart direct 5 s, then liver, kidney, brain.")
R[8] = dict(day="Day 2", date=D2, tkey="Day 2", Brain=("s", 65),
            comment="Brain direct 1 min 5 s, then heart, liver, kidney. Kidney stayed ~3 min in the animal before "
                    "cooling. Both kidneys taken at 60 min.")
R[9] = dict(day="Day 3", date=D3, tkey="Day 3", arrest="12:46:00", Liver=("s", 5), Brain="12:49:00",
            comment="Arrest noted as 12:46 (minute only). Liver direct 5 s; kidney direct 57 s after liver was cut; "
                    "an organ put in buffer ~12:48; brain 12:49.")
R[10] = dict(day="Day 3", date=D3, tkey="Day 3", arrest="13:06:00", Heart=("s", 4), Liver="13:07:00",
             Brain="13:08:00", **{"Kidney 1": "13:08:00"},
             comment="Arrest noted as 13:06 (minute only). Heart direct 4 s; liver, brain and kidney to buffer at "
                     "13:07–13:08.")
R[11] = dict(day="Day 3", date=D3, tkey="Day 3", Brain=("s", 45), comment="Brain direct 45 s, then heart, kidney, liver.")
R[12] = dict(day="Day 3", date=D3, tkey="Day 3", arrest="13:44:50", Liver=("s", 4), Heart="13:45:40",
             Brain="13:47:45", **{"Kidney 1": ("s", 12), "Kidney 2": "13:47:45"},
             comment="Liver direct 4 s, kidney direct 12 s. Livers up to this rat were cut in two; later livers in "
                     "more pieces (possibly faster cooling).")
R[13] = dict(day="Day 3", date=D3, tkey="Day 3", arrest="14:02:40", Heart=("s", 8), **{"Kidney 1": "14:03:10"},
             Liver="14:04:00", Brain="14:06:00",
             comment="Heart direct 8 s (cut-to-freeze 7 s). Only one kidney taken.")
R[14] = dict(day="Day 3", date=D3, tkey="Day 3", Brain="14:37:40", Liver="14:39:00", Heart="14:39:00",
             **{"Kidney 1": "14:40:00"},
             comment="Brain cut 14:37:40, 47 s direct procurement (arrest ≈ 14:36:53, derived). Buffer with organs "
                     "21.2 °C at 14:48.")
R[14]["arrest"] = "14:36:53"
for n, a, org, c in [
    (15, "10:34:25", dict(Liver="10:34:32", Heart="10:34:58", **{"Kidney 1": "10:36:26"}, Brain="10:36:36"), ""),
    (16, "10:51:52", dict(Heart="10:51:57", Liver="10:52:30", **{"Kidney 1": "10:53:03"}, Brain="10:53:40"),
     "One brain half got 10 min extra ischemia at room temperature."),
    (17, "11:07:00", dict(Liver="11:07:03", **{"Kidney 1": ("s", 11), "Kidney 2": "11:12:20"}, Heart="11:08:30",
                          Brain="11:10:00"), "Kidney 1 at 11 s (clamp slightly warm)."),
    (18, "11:40:05", dict(Liver="11:40:11", **{"Kidney 1": ("s", 8), "Kidney 2": "11:42:05"}, Heart="11:40:51",
                          Brain="11:42:30"), "Good acute kidney."),
    (19, "12:04:43", dict(Liver="12:04:49", **{"Kidney 1": ("s", 8), "Kidney 2": "12:06:21"}, Heart="12:05:20",
                          Brain="12:07:05"), "Kidney 1 taken before the animal was killed (noted as 8 s)."),
    (20, "12:21:04", dict(Liver="12:21:20", **{"Kidney 1": ("s", 8), "Kidney 2": "12:22:52"}, Heart="12:21:46",
                          Brain="12:23:12"), ""),
]:
    R[n] = dict(day="Day 4", date=None, tkey="Day 4", arrest=a, comment=c, **org)
for n, a, org, c in [
    (21, "12:38:45", dict(Heart="12:38:52", Liver="12:39:07", **{"Kidney 1": "12:40:25"}, Brain="12:40:38"), ""),
    (22, "12:52:27", dict(Liver="12:52:33", Heart="12:52:45", **{"Kidney 1": "12:53:50"}, Brain="12:54:08"), ""),
    (23, "13:05:32", dict(Heart="13:05:37", Liver="13:05:50", **{"Kidney 1": "13:07:14"}, Brain="13:07:50"), ""),
    (24, "13:21:04", dict(Heart="13:21:08", Liver="13:21:15", **{"Kidney 1": "13:21:59"}, Brain="13:22:28"), ""),
    (25, "13:31:07", dict(Heart="13:31:12", Liver="13:31:22", **{"Kidney 1": "13:32:23"}, Brain="13:32:23"), ""),
    (26, None, {}, "Sampling times not written. Weight not recorded."),
    (27, "13:30:02", dict(Heart="13:30:10", Liver="13:30:50", **{"Kidney 1": "13:31:33"}, Brain="13:32:00"),
     "Arrest written as 1:30:02, read as 13:30:02. Noted 'a degree warmer'. Brain splashed onto the wall during "
     "removal, so an extra rat was taken (rat 31; see notes below the table)."),
]:
    R[n] = dict(day="Day 5", date=D5, tkey="Day 5", arrest=a, comment=c, **org)
for n, a, org, c in [
    (28, "13:41:56", dict(Brain="13:42:41", Heart="13:43:38", Liver="13:43:47", **{"Kidney 1": "13:44:25"}),
     "Kidney splashed onto the wall."),
    (29, "13:55:01", dict(Liver="13:55:09", Heart="13:55:45", Brain="13:57:15"), "No kidney from this rat."),
    (30, "14:14:17", dict(Liver="14:14:24", Heart="14:15:18", Brain="14:16:40"),
     "Heart ischemia 4 min too long (34 min)."),
]:
    R[n] = dict(day="Day 6", date=None, tkey="Day 6", arrest=a, comment=c, **org)
R[31] = dict(day=None, date=None, tkey=None, comment="No sampling times recorded.")

F = {
    "FR1": dict(day="18 June", date="2026-06-18", tkey="18 June", Brain="10:49:00",
                Heart=("s", 90), Liver=("s", 90),
                comment="Brain sampled first (within 1 min; clock noted '10>49'); other tissues <2 min (~90 s). "
                        "Both kidneys to 30 min. Planned as heart/liver control, brain 30 min – sampled in reverse; "
                        "groups here follow the actual order (confirmed by the metabolite data)."),
    "FR2": dict(day="18 June", date="2026-06-18", tkey="18 June", arrest="11:05:21", Heart="11:05:37",
                Liver="11:05:37", Brain="11:07:45", **{"Kidney 1": "11:07:45"},
                comment="Heart and liver frozen 16 s after death; brain and kidney 11:07:45. Both kidneys to 30 min. "
                        "Planned as brain control, heart/liver 30 min – sampled in reverse; groups here follow the "
                        "actual order (confirmed by the metabolite data)."),
    "FR3": dict(day=None, date=None, tkey=None, arrest="11:47:46", Heart="11:48:01", Liver="11:48:01",
                **{"Kidney 1": "11:48:11"}, Brain="11:49:45",
                comment="Heart and liver clamped 11:48:01. One kidney direct (11:48:11), one kidney to 30 min. "
                        "Date not noted."),
    "FR4": dict(day="22 June", date="2026-06-22", tkey="22 June", arrest="11:38:59", Brain="11:39:45",
                Heart="11:40:10", Liver="11:40:25", **{"Kidney 1": "11:41:00", "Kidney 2": "11:41:15"},
                comment="Kidney 1 direct, kidney 2 to 30 min ischemia."),
    "FR5": dict(day="22 June", date="2026-06-22", tkey="22 June", arrest="11:54:05", Brain="11:54:40",
                Heart="11:55:30", Liver="11:55:40", **{"Kidney 1": "11:56:00", "Kidney 2": "11:56:10"},
                comment="Kidney 1 direct, kidney 2 to 30 min ischemia."),
    "FR6": dict(day="22 June", date="2026-06-22", tkey="22 June", arrest="12:15:41", Liver="12:15:59",
                Heart="12:16:03", **{"Kidney 1": "12:16:30"}, Brain="12:18:30",
                comment="Kidney 1 direct."),
}
# corrected ischemia groups for FR1/FR2 (actual order)
FR_GROUPS = {"FR1": {"Heart": 30, "Liver": 30, "Brain": 0}, "FR2": {"Heart": 0, "Liver": 0, "Brain": 30}}
EXTRA_MALE = {
    32: "Rat 32 (sick, used for the oedema study) was originally numbered 31; the numbers were changed because an "
        "extra rat was taken after the brain splash on rat 27.",
}

# ---------------------------------------------------------------- groups & weights from the source workbook
wb = openpyxl.load_workbook(SRC, data_only=True)
blocks = {"Brain": 2, "Liver": 8, "Heart": 14, "Kidney": 20}
TIMES = [0, 30, 60, 180, 360]


def grid(sheet, prefix):
    out, wt = {}, {}
    for r in list(wb[sheet].iter_rows(values_only=True))[2:]:
        if not isinstance(r[0], int):
            continue
        key = f"{prefix}{r[0]}"
        wt[key] = r[1]
        g = {}
        for o, c0 in blocks.items():
            hits = [TIMES[i] for i in range(5) if r[c0 + i] == 1]
            g[o] = hits
        out[key] = g
    return out, wt


gm, wm = grid("Male rats", "")
gf, wf = grid("Female rats", "FR")
# HMT sample ids
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from data import Tissue
hmt = {}
for tn, suf in [("Heart", "h"), ("Brain", "b"), ("Liver", "L")]:
    x = Tissue(tn, correct_fr=False)
    for s in x.cols("Rat"):
        rid = s[:-1] if s.startswith("FR") else s[:-1]
        hmt[(rid, tn)] = s


def tval(s):
    h, m, sec = map(int, s.split(":"))
    return dt.time(h, m, sec)


# ---------------------------------------------------------------- write workbook
out = openpyxl.Workbook()
ws = out.active
ws.title = "Sampling times"
cols = [str(n) for n in range(1, 32)] + list(F)
rows = ["Rat", "Sex", "Weight (g)", "Sampling day", "Date", "Room / buffer temp (°C)", "Circulatory arrest (clock)"]
for o in ORGANS:
    rows += [f"{o} – ischemia group (min)", f"{o} – removed (clock)", f"{o} – seconds after arrest"]
    if o in ("Heart", "Liver", "Brain"):
        rows.append(f"{o} – HMT sample")
rows.append("Comments")
rowidx = {r: i + 1 for i, r in enumerate(rows)}

bold = Font(bold=True)
hdr_fill = PatternFill("solid", fgColor="DDE7F5")
org_fill = {"Heart": "FBE5D6", "Liver": "E2EFDA", "Brain": "EDE7F6", "Kidney 1": "FFF2CC", "Kidney 2": "FFF2CC"}
thin = Side(style="thin", color="D0D0D0")
for r, name in enumerate(rows, 1):
    c = ws.cell(r, 1, name); c.font = bold
    for o, col in org_fill.items():
        if name.startswith(o + " –"):
            c.fill = PatternFill("solid", fgColor=col)
ws.cell(1, 1).fill = hdr_fill

italic_grey = Font(italic=True, color="7F7F7F")
for j, rid in enumerate(cols, 2):
    L = get_column_letter(j)
    female = rid.startswith("FR")
    rec = F[rid] if female else R.get(int(rid), {})
    gsrc = (gf if female else gm).get(rid)
    wt = (wf if female else wm).get(rid)
    ws.cell(1, j, rid if female else f"Rat {rid}").font = bold
    ws.cell(1, j).fill = hdr_fill
    ws.cell(rowidx["Sex"], j, "F" if female else "M")
    ws.cell(rowidx["Weight (g)"], j, wt if isinstance(wt, (int, float)) else ("not recorded" if wt in (None, "?") else wt))
    ws.cell(rowidx["Sampling day"], j, rec.get("day"))
    if rec.get("date"):
        c = ws.cell(rowidx["Date"], j, dt.date.fromisoformat(rec["date"])); c.number_format = "yyyy-mm-dd"
    ws.cell(rowidx["Room / buffer temp (°C)"], j, TEMP.get(rec.get("tkey")))
    if rec.get("arrest"):
        c = ws.cell(rowidx["Circulatory arrest (clock)"], j, tval(rec["arrest"])); c.number_format = "hh:mm:ss"
    for o in ORGANS:
        # group
        base = "Kidney" if o.startswith("Kidney") else o
        gl = list((gsrc or {}).get(base, []))
        if female and rid in FR_GROUPS and o in FR_GROUPS[rid]:
            gl = [FR_GROUPS[rid][o]]
        if o.startswith("Kidney"):
            k = 0 if o == "Kidney 1" else 1
            g = gl[k] if len(gl) > k else None
            if female and rid in ("FR1", "FR2"):
                g = 30
            if female and rid == "FR3":
                g = 0 if o == "Kidney 1" else 30
        else:
            g = ", ".join(map(str, gl)) if gl else None
            if g and "," not in g:
                g = int(g)
        ws.cell(rowidx[f"{o} – ischemia group (min)"], j, g)
        v = rec.get(o)
        rr, rs = rowidx[f"{o} – removed (clock)"], rowidx[f"{o} – seconds after arrest"]
        if isinstance(v, str):
            c = ws.cell(rr, j, tval(v)); c.number_format = "hh:mm:ss"
            if rec.get("arrest"):
                ws.cell(rs, j, f'=ROUND(({L}{rr}-{L}{rowidx["Circulatory arrest (clock)"]})*86400,0)')
        elif isinstance(v, tuple):
            c = ws.cell(rs, j, v[1]); c.font = italic_grey
        if o in ("Heart", "Liver", "Brain"):
            ws.cell(rowidx[f"{o} – HMT sample"], j, hmt.get((rid if female else rid, o)))
    com = rec.get("comment", "")
    ws.cell(rowidx["Comments"], j, com or None)

# layout
ws.column_dimensions["A"].width = 30
for j in range(2, len(cols) + 2):
    ws.column_dimensions[get_column_letter(j)].width = 12
for r in range(1, len(rows) + 1):
    for j in range(1, len(cols) + 2):
        c = ws.cell(r, j)
        c.border = Border(top=thin, bottom=thin, left=thin, right=thin)
        c.alignment = Alignment(horizontal="left" if j == 1 else "center", vertical="top",
                                wrap_text=(rows[r - 1] == "Comments"))
ws.row_dimensions[rowidx["Comments"]].height = 150
ws.freeze_panes = "B2"
n = len(rows) + 2
notes = [
    "Notes",
    "Circulatory arrest = time of death (heart cut out). 'Seconds after arrest' is a formula where both clock times "
    "exist; grey italic values were noted only as seconds (e.g. 'heart direct 5 s').",
    "Ischemia group 0 = direct procurement (control). Groups from the 'Male rats' / 'Female rats' sheets; FR1 and FR2 "
    "follow the actual sampling order (the reverse of the plan), which the metabolite data confirm.",
    "HMT sample = sample ID in the HMT metabolomics data (h = heart, L = liver, b = brain).",
    "Days 4 and 6 dates and rat 31 times were not recorded in the source file. Rats 32–34 were used for the oedema "
    "study and are not listed. " + EXTRA_MALE[32],
]
for i, t in enumerate(notes):
    c = ws.cell(n + i, 1, t)
    c.font = bold if i == 0 else Font(italic=True, color="52514E")

# original notes sheet
ws2 = out.create_sheet("Original notes")
ws2.append(["Source sheet", "Text (verbatim)"])
for c in ws2[1]:
    c.font = bold
for sh in ["Sampling times male rats", "Female rats", "Room temperatures"]:
    for r in wb[sh].iter_rows(values_only=True):
        for v in r:
            if isinstance(v, str) and len(v) > 12 and not re.fullmatch(r"[A-Za-z #]+", v):
                ws2.append([sh, v])
ws2.column_dimensions["A"].width = 26
ws2.column_dimensions["B"].width = 160
out.save(OUT)
print("wrote", OUT)
