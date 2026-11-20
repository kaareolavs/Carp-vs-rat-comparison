"""Build all PDF reports into ../reports. Usage: python3 build_reports.py [Heart Brain Liver sex]"""
import sys, subprocess
from pathlib import Path
import pandas as pd
HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "reports"
OUT.mkdir(exist_ok=True)


def pdf(html, name, title):
    h = OUT / f"{name}.html"
    h.write_text(html, encoding="utf-8")
    subprocess.run(["node", str(HERE / "print_pdf.js"), str(h), str(OUT / f"{name}.pdf"), title], check=True)
    h.unlink()
    print("wrote", OUT / f"{name}.pdf")


def excel(sheets, name):
    with pd.ExcelWriter(OUT / f"{name}.xlsx") as w:
        for k, df in sheets.items():
            df.to_excel(w, sheet_name=k[:31])


if __name__ == "__main__":
    what = sys.argv[1:] or ["Heart", "Brain", "Liver", "sex"]
    for t in what:
        if t == "sex":
            import build_sex
            html, xl = build_sex.build()
            name = "Rat_male_vs_female_Ctr_30min"
            pdf(html, name, "Rat ischemia: male vs female")
        else:
            import build_species
            html, xl = build_species.build(t)
            name = f"{t}_crucian_carp_vs_rat"
            pdf(html, name, f"Crucian carp vs rat – {t}")
        excel(xl, name + "_statistics")
