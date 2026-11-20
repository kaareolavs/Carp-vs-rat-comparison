# Crucian carp vs rat – metabolome re-analysis

Builds the PDF reports in `../reports/` from the HMT workbooks in the repository root.

| Report | Content |
|---|---|
| `Heart_crucian_carp_vs_rat.pdf`, `Brain_…`, `Liver_…` | Crucian carp vs rat ischemia time-course for one tissue: PCA, response heat maps, hallmark metabolites, pathway maps, metabolic parameters, statistical tables |
| `Rat_male_vs_female_Ctr_30min.pdf` | Male vs female rats at control and 30 min, all three tissues |
| `*_statistics.xlsx` | All group means, effects and q-values behind each report |

## Rebuild

```bash
pip install numpy pandas scipy scikit-learn matplotlib openpyxl
cd analysis
python3 build_reports.py            # all reports
python3 build_reports.py Heart sex  # selected reports
```

PDF printing uses Node + Playwright (`print_pdf.js`) with the Inter font.

## Files

- `data.py` – loads `Table4 Absolute Concentration`; applies the FR1/FR2 label correction (`FR_CORRECTION`).
- `stats.py` – Welch t-tests, species × time and sex × time linear models, Benjamini–Hochberg FDR.
- `figures.py` – matplotlib figures (rendered as vector SVG).
- `pathways.py` – simplified pathway map layouts.
- `findings.py` – hand-written interpretation text shown in the summaries; update it if the data change.
- `build_species.py`, `build_sex.py`, `report_html.py` – report layout.

## Label correction

Female rats FR1 and FR2 are swapped relative to the HMT labels in all three tissues
(heart/liver: FR2 = control, FR1 = 30 min; brain: FR1 = control, FR2 = 30 min).
Set `Tissue(name, correct_fr=False)` to use the original labels.
