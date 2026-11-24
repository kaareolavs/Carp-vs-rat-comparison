const fs = require('fs');
const { Document, Packer, Paragraph, TextRun, LevelFormat, AlignmentType } = require('docx');
const F = 'Calibri', SZ = 20; // 10 pt
const p = (runs, opts = {}) => new Paragraph({ spacing: { after: 80, line: 260 }, ...opts,
  children: runs.map(r => typeof r === 'string' ? new TextRun({ text: r, font: F, size: SZ }) : new TextRun({ font: F, size: SZ, ...r })) });
const b = (t) => ({ text: t, bold: true });
const i = (t) => ({ text: t, italics: true });
const bullet = (runs) => p(runs, { numbering: { reference: 'bul', level: 0 }, spacing: { after: 50, line: 250 } });
const doc = new Document({
  styles: { default: { document: { run: { font: F, size: SZ } } },
    paragraphStyles: [{ id: 'Title', name: 'Title', run: { font: F, size: 28, bold: true }, paragraph: { spacing: { before: 0, after: 100 }, keepNext: true } }] },
  numbering: { config: [{ reference: 'bul', levels: [{ level: 0, format: LevelFormat.BULLET, text: '•', alignment: AlignmentType.LEFT,
    style: { paragraph: { indent: { left: 300, hanging: 200 } } } }] }] },
  sections: [{ properties: { page: { margin: { top: 1134, bottom: 1134, left: 1134, right: 1134 } } }, children: [
    new Paragraph({ style: 'Title', children: [new TextRun({ text: 'Sex differences in rat tissue metabolites: re-analysis with body weight and sampling date', font: F, size: 28, bold: true })] }),
    p(['All female rats were sampled on 18–22 June 2026, after all males (7–10 April, 7 May and 20 May 2026). All rats were born on 15 January 2026, so the males were ~12–18 weeks old at sampling and the females ~22 weeks, i.e. 4–11 weeks older; sex and age are fully confounded and cannot be separated statistically. Instead, we used the 14 males per organ, sampled in numbered order over 6 weeks while gaining weight (', i('ρ'), ' = 0.63 between rat number and weight), to estimate within-male trends with sampling order and weight after removing ischemia-time effects, and asked whether the female–male differences at control follow these trends.']),
    bullet([b('Free hydroxyproline – age. '), 'Declined with sampling date in males in all three organs (r = −0.54 to −0.69 with rat number; −4% to −17% per week of age within ischemia groups) and was 2–4× lower in females, who were sampled last and continued the same gradient (median heart 63 → 31 → 17 nmol/g in April males ~12 weeks, May males 16–18 weeks and females 22 weeks). Free hydroxyproline, a collagen-turnover product, is the most consistent age-declining metabolite across organs in a multi-organ mouse ageing atlas (11 of 12 organs; Pilley et al.), so the female difference fits the females being older. This is distinct from total (collagen-bound) hydroxyproline, which measures collagen content and increases with age in the heart. Because ischemia-time groups were unevenly spread over April and May, part of the differences in rat free hydroxyproline between time groups probably reflects age at sampling.']),
    bullet([b('Liver carnitine – age. '), 'Rose with sampling order in males (r = +0.73, p = 0.003) and was 2.2–2.6× higher in females. Ruoppolo et al. found no sex difference in liver carnitine in age-matched rats, so this is more likely age or cohort than sex.']),
    bullet([b('Amino acids in heart – sex. '), 'Val, Leu/Ile, Phe, Tyr and Gly were 20–40% lower (ornithine ~55% lower) in female control hearts but showed no significant age or weight trend in males. Ruoppolo et al. (age-matched Sprague-Dawley rats, protein-normalised) likewise found higher Glu, Gly, Ala, Leu/Ile, Val, Phe, ornithine, citrulline and Lys in male hearts, supporting a true sex difference. The gap narrowed after 30 min ischemia (e.g. Val 0.72 → 0.83 female/male).']),
    bullet([b('Overall profile and heart energy state – unresolved. '), 'The male sampling-order signature did not predict the global ~30% lower metabolite levels in female controls (Spearman ρ = −0.27 to +0.16 across organs). In male hearts, early ischemia markers (AMP, ADP, xanthine) fell with sampling order, suggesting the procedure improved over the sampling days; the better energy state of female control hearts (sampled last) may partly reflect this rather than sex.']),
    p([b('Conclusion. '), 'The ischemic response did not differ between sexes. Of the baseline differences, free hydroxyproline and liver carnitine are best explained by age, whereas lower amino acids in female hearts agree with published sex dimorphism. Future studies should use age-matched animals of both sexes sampled in alternating order.'], { spacing: { before: 60, after: 80, line: 260 } }),
    p([{ text: 'References: Ruoppolo M et al. Targeted metabolomic profiling in rat tissues reveals sex differences. Sci Rep 2018;8:4663. Pilley et al. A metabolic atlas of mouse aging. bioRxiv 2024, doi:10.1101/2024.05.04.592445.', size: 16, color: '52514E' }]),
  ] }],
});
Packer.toBuffer(doc).then(buf => fs.writeFileSync(process.argv[2], buf));
