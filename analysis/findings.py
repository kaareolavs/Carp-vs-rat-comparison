"""Curated interpretation text, written after reviewing the statistics of each report.
Numbers quoted are geometric-mean fold changes vs each species' own control unless stated otherwise."""

SPECIES = {
    "Heart": [
        "<b>The crucian carp heart keeps its ATP for the first hour; the rat heart does not.</b> Carp ATP and "
        "adenylate energy charge were essentially unchanged at 30 and 60 min (energy charge 0.98 → 0.97 → 0.93), whereas rat ATP fell "
        "~4-fold by 30 min and ~10-fold by 60 min (energy charge 0.94 → 0.53). ATP, phosphocreatine and energy charge "
        "all have significant species × time interactions.",
        "<b>Both hearts switch to anaerobic glycolysis at similar magnitude.</b> Lactate rose ~20–60-fold "
        "(not significant in carp because control values varied 10-fold), succinate ~8–16-fold and glycerol "
        "3-phosphate ~8–18-fold in both species, so the carp's preserved "
        "energy status is not explained by a smaller glycolytic response but by maintaining ATP while the rat heart "
        "cannot. The carp heart used phosphocreatine as a buffer: it fell gradually (~4-fold at 30 min, ~13-fold at 60 "
        "min) while ATP stayed constant.",
        "<b>Adenylate breakdown is delayed in carp.</b> In the rat heart AMP (~30-fold), adenosine, inosine and "
        "hypoxanthine rose already at 30 min. In carp, AMP did not rise until 180 min and inosine/xanthine "
        "accumulated later and more slowly; hypoxanthine rose in both.",
        "<b>Glycogen mobilisation in carp; loss of glycolytic intermediates in rat.</b> Glucose 1-phosphate rose "
        "progressively in the carp heart (~18–25-fold at 180–360 min) but not in rat, consistent with a large "
        "glycogen store feeding glycolysis, and carp kept UDP-glucose constant while the rat heart lost it (~8-fold "
        "by 180 min). The rat heart also lost 3-phosphoglycerate, 2,3-bisphosphoglycerate and PEP; several of these "
        "were below detection in carp, so their interaction results partly reflect detection. At 360 min carp ribose "
        "5-phosphate rose ~35-fold.",
        "<b>The carp heart eventually fails too.</b> By 180–360 min carp ATP had fallen to rat-like levels (energy "
        "charge 0.21 at 360 min), so the difference is a delay of hours rather than complete protection.",
        "<b>Baseline differences.</b> 48 of 106 metabolites differed between species at control (e.g. rat higher in "
        "2,3-BPG, fumarate, malate, urea, uric acid, asparagine; carp higher in SAH, choline, putrescine, "
        "phosphocreatine, PRPP), which is why the response-based comparison is used for the main conclusions.",
    ],
    "Brain": [
        "<b>The carp brain loses far less energy than the rat brain.</b> Carp brain ATP was not significantly reduced "
        "until 360 min (~2-fold lower at 30 min, n.s.) and energy charge fell from 0.97 to 0.75 at 30 min, whereas "
        "rat brain ATP fell ~14-fold by 30 min and energy charge from 0.61 to 0.14. ADP and AMP have the strongest species × time interactions: in carp ADP (~3-fold) and AMP (~20-fold) "
        "accumulate, whereas in rat ADP falls and AMP rises only slightly because the adenylates are degraded "
        "further to inosine and hypoxanthine. ATP and energy charge also differ significantly.",
        "<b>Note on rat brain controls.</b> The rat brain control samples already have a low energy charge (0.61 vs "
        "0.97 in carp control) and high AMP, consistent with the brain being removed 45–65 s after circulatory arrest (35–60 s in females). Rat brain "
        "responses at 30 min are therefore measured against a partly ischemic baseline, and the true rat response is, "
        "if anything, underestimated.",
        "<b>Glycolytic fuel is preserved in carp.</b> Rat brain glucose 6-phosphate (~6-fold), F6P, FBP (~9-fold) and "
        "pyruvate (~9-fold) dropped within 30 min, i.e. substrate depletion. In carp, G6P and F6P were maintained and glucose "
        "1-phosphate rose ~2-fold (glycogen breakdown) while it fell ~30-fold in rat; carp FBP fell only "
        "moderately (~2.5-fold).",
        "<b>Similar lactate and succinate accumulation.</b> Lactate rose ~2-fold and succinate ~2–3-fold in both "
        "species; the lactate/pyruvate ratio rose much more in rat (~13-fold vs ~2.5-fold) because pyruvate fell in "
        "rat only.",
        "<b>TCA intermediates.</b> Malate, fumarate and citrate fell rapidly in rat brain (malate ~5-fold by 30 min) "
        "but were stable in carp until 360 min. 2-Oxoglutarate fell sharply in carp brain (~25-fold at 30 min) and was "
        "mostly below detection in rat. Carp brain had markedly higher control levels of phosphocreatine, proline, betaine and 2-oxoglutarate.",
        "<b>Purine catabolism is smaller in carp.</b> Inosine rose ~8-fold in rat but only ~2.5-fold (30 min, then back "
        "to baseline) in carp; hypoxanthine increased in both, more in rat.",
    ],
    "Liver": [
        "<b>The carp liver is the most hypoxia-resistant tissue in this data set.</b> Carp liver ATP and energy charge "
        "were essentially unchanged through 60 min and only modestly reduced at 180–360 min (ATP ~2–3-fold), whereas "
        "rat liver ATP fell ~9-fold by 30 min and ~85-fold by 360 min (energy charge 0.75 → ~0.3).",
        "<b>The carp liver shows little sign of anaerobic stress.</b> Lactate rose only ~2–5-fold in carp vs ~10–14-fold "
        "in rat; succinate rose ~9-fold in rat but was nearly unchanged in carp (one of the strongest interactions); "
        "glycerol 3-phosphate rose ~9-fold in rat but not in carp.",
        "<b>Rat liver TCA cycle collapses; carp's does not.</b> Citrate fell >70-fold in rat liver (to near detection "
        "limit) but only ~3–8-fold in carp; malate and fumarate fell late in both species.",
        "<b>Carp liver maintains glycogen-linked and nucleotide pools.</b> UDP-glucose, GTP/GDP, ATP/ADP and "
        "6-phosphogluconate were maintained (sedoheptulose 7-phosphate even increased) in carp but fell in rat "
        "(UDP-glucose ~60-fold, GTP ~30-fold by 360 min) (strong interactions), "
        "consistent with continued glycogen-fuelled metabolism and pentose-phosphate activity in carp liver.",
        "<b>Purine breakdown occurs in both.</b> Hypoxanthine, xanthine and uric acid increased strongly in both "
        "species (significant in rat; no significant interaction), from a much lower baseline in carp.",
        "<b>Baseline.</b> Large species differences at control include higher arginine, proline, phosphocreatine and "
        "citrulline in carp liver and higher carnitine, adenosine, betaine and N-acetylglutamate in rat liver.",
    ],
}

SEX = [
    "<b>No sex difference after 30 min ischemia.</b> In none of the three tissues did any metabolite differ "
    "significantly between female and male rats at 30 min (q &lt; 0.05); the ischemic response (ATP loss, lactate, "
    "succinate and hypoxanthine accumulation) was of similar size in both sexes.",
    "<b>Heart: differences are at control and involve purine breakdown products.</b> Male control hearts had higher "
    "AMP (~20-fold), inosine (~10-fold) and hypoxanthine (N.D. in females) and a lower energy charge (0.90 vs 0.98) "
    "than female control hearts. This drives a sex main effect for 17 metabolites and a sex × ischemia interaction "
    "for 9 (inosine, hypoxanthine, adenosine, IMP, ribose phosphates) – the females simply start lower. The sampling records rule out "
    "faster handling of the females: female control hearts were removed 15–22 s after circulatory arrest, male "
    "control hearts 4–8 s. Female control hearts thus preserved their energy state better in the first seconds, which "
    "may reflect sex or cohort differences (all rats were born on 15 January 2026, so the females were 4–11 weeks older: ~22 vs 12–18 weeks).",
    "<b>Brain and liver: very few sex differences</b> (3–4 metabolites with a sex main effect, e.g. hydroxyproline in "
    "both tissues, ribulose 5-phosphate and ADP-ribose in brain; no significant interactions).",
    "<b>But the overall profiles separate by sex.</b> In PCA, female and male rats are partly separated along "
    "PC2 in all three tissues (clearest in heart and liver, and among brain controls), i.e. many small, consistent "
    "differences that individually do not survive FDR correction. Because the male and female rats were different "
    "cohorts (and the females were sampled with the paired tissue scheme), this may reflect sex or other cohort "
    "differences.",
    "<b>Statistical power is limited</b> (n = 3 per sex and time point), so moderate sex differences cannot be "
    "excluded. Since the ischemic response itself did not differ, pooling males and females in the fish-vs-rat "
    "comparison at control and 30 min appears reasonable, but sex/cohort should be kept in mind for "
    "baseline comparisons.",
]

SEX_TISSUE = {
    "Heart": ["Both sexes lose ~3/4 of ATP and accumulate lactate ~15–30-fold within 30 min. Female control hearts "
              "had a higher energy charge (0.98 vs 0.90) and much lower AMP and inosine than male controls, although they "
              "were removed later after circulatory arrest (15–22 s vs 4–8 s); at 30 min "
              "no single metabolite differed, although the PCA still separates "
              "the sexes along PC2."],
    "Brain": ["Few individual metabolites differ between male and female rat brains, although in the PCA the "
              "female control brains lie apart from the male controls along PC2. Both groups' control samples show a "
              "partly ischemic profile (energy charge ~0.6), reflecting brain removal time."],
    "Liver": ["Male and female rat livers responded alike to 30 min ischemia (ATP ~10-fold lower, lactate ~10-fold "
              "higher). Female control livers had lower inosine and undetectable hypoxanthine, a pattern similar to the "
              "heart but not significant after FDR correction."],
}
