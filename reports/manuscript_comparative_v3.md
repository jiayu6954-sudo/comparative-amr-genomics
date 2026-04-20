# IS Element Ecology Shapes Divergent Carbapenem Resistance Gene Mobilisation Strategies Across Clinical Gram-Negative Pathogens: A Three-Species Comparative Genomic Analysis

**Authors**: ZJY^1,\*^

**Affiliations**:
^1^ Independent Researcher

\*Corresponding author: jiayu6954@gmail.com

**Word count (main text)**: ~5,400 words
**Figures**: 4 (Figures 1–4)
**Supplementary Tables**: 6 (Tables S1–S6)

---

## Abstract

### Background
Insertion sequence (IS) elements—principally IS26 (IS6 family) in Enterobacteriaceae and ISAba1 in *Acinetobacter baumannii*—mediate carbapenem resistance gene mobilisation by distinct molecular mechanisms. Whether IS-mediated composite transposon architecture is conserved across the three principal Gram-negative carbapenem-resistant pathogens, and whether IS genomic copy number serves as a generalisable cross-species resistance predictor, has not been systematically addressed.

### Methods
We retrieved, quality-controlled, and analysed whole-genome assemblies from Chinese clinical isolates of three WHO critical-priority pathogens: *Klebsiella pneumoniae* (n=492), *Escherichia coli* (n=266), and *Acinetobacter baumannii* (n=168). Carbapenem resistance genes were detected by GFF3 annotation scanning and classified against CARD v3.3. IS element flanking architecture within ±10 kb of each resistance locus was characterised and classified as COMPOSITE_TRANSPOSON, SINGLE_IS, or NO_IS. IS annotation bias was assessed and corrected using PFAM HMM scanning (PF01527, IS6/ISAba1 family transposase). Focal IS copy number (IS6 for *K. pneumoniae* and *E. coli*; ISAba for *A. baumannii*) was tested as a resistance predictor via Mann–Whitney AUC and Cliff's δ.

### Results
IS-mediated composite transposon architecture flanked 94.0% (47/50, 95%CI 83.5–98.8%) of *K. pneumoniae* acquired carbapenemase loci, 62.1% (18/29, 95%CI 44.0–77.3%) of *E. coli* carbapenemase loci, and 6.8% (10/146) of *A. baumannii* acquired carbapenemase loci—a >13-fold gradient from Enterobacteriaceae to Acinetobacter (Fisher's exact: KPN vs ECO p = 3.2×10⁻³; KPN vs ABA p < 10⁻¹⁵; ECO vs ABA p = 3.7×10⁻⁸). *Escherichia coli* displayed an unexpectedly diverse IS family repertoire (IS6/IS26 23.4%, Tn3 family 20.4%, IS5 13.1%, IS1 9.5%), contrasting with IS26 dominance in *K. pneumoniae*. IS6 copy number predicted carbapenem resistance in both enterobacterial species (KPN: AUC = 0.667, p = 4.2×10⁻¹⁰; ECO PFAM-corrected: AUC = 0.718, p = 2.3×10⁻⁶) but could not be assessed in *A. baumannii* (n=1 susceptible genome). PFAM HMM correction confirmed annotation bias in *E. coli* IS6 counts (susceptible genomes gained 9.9× more reclassifications: 959 vs 97), yielding a conservative corrected AUC. ISAba elements in *A. baumannii* were strongly enriched upstream of OXA-23 (92.3% upstream positions), consistent with a promoter activation rather than structural bracketing mechanism.

### Conclusions
Composite transposon architecture is shared by both Enterobacteriaceae species but is near-absent in *A. baumannii*, which instead deploys ISAba1 upstream promoter insertion as the dominant carbapenemase activation mechanism. IS6 genomic burden is a valid and significant resistance predictor across both enterobacterial hosts. Species-specific IS ecology dictates both the frequency and mechanism of carbapenem resistance gene mobilisation, establishing a framework for species-adapted IS surveillance.

---

## Introduction

Carbapenem-resistant Gram-negative bacteria (CR-GNB) constitute the most critical tier of the WHO priority pathogen list, with carbapenem-resistant *Klebsiella pneumoniae* (CRKP), *Escherichia coli* (CREC), and *Acinetobacter baumannii* (CRAB) accounting for the majority of healthcare-associated carbapenem-resistant infections in China [1,2]. Attributable mortality in bloodstream infections reaches 30–70%, and the genetic mechanisms underlying resistance gene spread are dominated by mobile genetic elements—plasmids, transposons, and insertion sequences (IS elements) [3,4].

IS elements, the simplest autonomous mobile elements encoding only a transposase gene flanked by inverted repeats, occupy a central mechanistic position in carbapenem resistance gene ecology. Two IS families have been individually characterised as the principal drivers of carbapenem resistance gene mobilisation in Gram-negative pathogens. IS26, a member of the IS6 family characterised by 14-bp inverted repeats and 3-bp target site duplications, is the canonical flanking element for *blaKPC* and *blaNDM* in *K. pneumoniae* and *E. coli* [5,6]. IS26 mediates resistance gene spread through composite transposon formation (bilateral IS26 bracketing enabling block transposition) and "one-ended transposition" (replicon fusion via a single IS26) [6]. In *A. baumannii*, ISAba family elements (principally ISAba1) have been characterised as drivers of *blaOXA-23* expression and dissemination through a mechanistically distinct pathway: ISAba1 insertion ~50–100 bp upstream of *blaOXA-23* provides a strong outward-reading P2 promoter, activating otherwise weakly expressed carbapenemase—a mechanism designated "promoter capture" [7,8].

Despite these well-characterised individual IS–resistance gene associations, fundamental comparative questions have remained unresolved. Is IS26-mediated composite transposon architecture equally prevalent in *E. coli* NDM and KPC contexts as in *K. pneumoniae*? Does ISAba1 in *A. baumannii* act primarily as a structural bracketer or as a transcriptional activator? Is focal IS copy number a generalisable cross-species predictor of carbapenem resistance, or is it species-specific? These questions require simultaneous whole-genome analysis of all three principal CR-GNB pathogens under an identical analytical framework—an approach not previously applied.

Here we present a systematic three-species comparative genomic analysis of IS element architecture around carbapenem resistance genes in Chinese clinical isolates, comprising a combined cohort of 926 genomes. We reveal that composite transposon architecture is broadly shared within Enterobacteriaceae (94% in *K. pneumoniae*, 62% in *E. coli*) but is near-absent in *A. baumannii* (6.8%), which deploys a mechanistically distinct ISAba1 promoter strategy. IS6 copy number predicts resistance in both enterobacterial species, but *E. coli* IS6 annotation is substantially biased by unlabelled IS_unknown features requiring PFAM HMM correction.

---

## Materials and Methods

### Genome retrieval and quality control

Genome assemblies for all three species were retrieved from NCBI GenBank (accessed April 2026) using the Entrez Programming Utilities API (esearch + esummary), restricted to: (i) Chinese isolation source (`"China"[Country]` filter); (ii) assembly level ∈ {Complete Genome, Chromosome, Scaffold}; (iii) latest, non-withdrawn status. Quality control applied species-specific thresholds: genome size 4.8–6.5 Mb, scaffold N50 ≥ 50 kb, annotated CDS ≥ 3,000 (*K. pneumoniae*); genome size 4.2–5.8 Mb, N50 ≥ 50 kb, CDS ≥ 3,000 (*E. coli*); genome size 3.7–4.5 Mb, N50 ≥ 50 kb, CDS ≥ 2,500 (*A. baumannii*). GFF3 integrity was verified by decompression and format parsing. All *K. pneumoniae* data were drawn from our previously validated cohort. *Escherichia coli* genomes were retrieved in two rounds to enrich the carbapenem-resistant fraction. The primary retrieval (n=237 QC-passed) used a general `"Escherichia coli"[Organism] AND "China"[Country]` query. A supplementary targeted retrieval (`01b_download_eco_cr.py`) ran four parallel Entrez queries restricted to carbapenem resistance keywords: (1) `"Escherichia coli"[Organism] AND "China"[Country] AND "NDM"[All Fields]`; (2) `AND "KPC"[All Fields]`; (3) `AND "carbapenem"[All Fields]`; (4) `AND "OXA-48"[All Fields]`—all additionally filtered by `latest[filter]` and assembly level ∈ {Complete Genome, Chromosome, Scaffold}. Results were deduplicated against existing downloads; the combined, deduplicated cohort totalled 266 QC-passed genomes.

### Carbapenem resistance gene detection

Species-optimised regular-expression patterns were applied against GFF3 `product`, `gene`, and `Note` fields targeting: NDM-\d, KPC-\d/*blaKPC*, OXA-23/24/40/58/72/143, IMP-\d, VIM-\d (*E. coli* and *K. pneumoniae*); OXA-23/24/40/58/72/143/182 and NDM variants (*A. baumannii*). For *A. baumannii*, intrinsic OXA-51 (chromosomal species marker, present in ~99% of *A. baumannii* regardless of acquired resistance) was separated from acquired carbapenemases (OXA-23/24/58/NDM) for composite transposon analysis.

### IS element flanking analysis

IS elements within ±10 kb of each carbapenem resistance locus were identified by annotation scanning with species-adapted keyword sets. IS elements were classified by family: IS6 (IS26, IS257, IS1353, IS1006), ISAba (ISAba1/2/3), IS5, IS1, IS3, IS91, Tn3/transposon family, and IS_other. Flanking architecture was classified as: COMPOSITE_TRANSPOSON (IS elements detected on both upstream and downstream flanks), SINGLE_IS_UPSTREAM, SINGLE_IS_DOWNSTREAM, or NO_IS. Composite transposon rates used Clopper–Pearson exact 95%CI for *K. pneumoniae* (n=50 loci); Wilson score intervals for *E. coli* and *A. baumannii*.

### PFAM HMM annotation correction

IS6/ISAba1 annotation is known to be incomplete in NCBI-deposited GFF3 files: IS elements annotated without a family name ("transposase", "insertion sequence element") are systematically underdetected by name-based scanning. To correct this bias, all IS_unknown features (IS/transposase-annotated CDS without a recognised family designation) were extracted and protein sequences translated from corresponding FNA files. Proteins were searched against the PF01527 (IS6 family transposase, which covers both IS26 and ISAba1 as members of the IS6 superfamily) using pyhmmer v0.12 (E ≤ 1×10⁻⁵, identity ≥ 0.25). Reclassified IS_unknown features were included in corrected IS6 counts. Annotation bias was assessed as the ratio of reclassifications in susceptible vs resistant genomes.

### IS burden as resistance predictor

Per-genome focal IS copy numbers were computed by counting GFF3 CDS features matching species-specific IS patterns (IS6 for *K. pneumoniae* and *E. coli*; ISAba for *A. baumannii*). Predictive accuracy was assessed as Mann–Whitney AUC (one-sided U / n₁n₂) for IS burden discriminating resistant from susceptible genomes. Cliff's δ quantified effect size. Two-sided Mann–Whitney p-values were computed. For *E. coli*, both original and PFAM-corrected IS6 counts were analysed; corrected values are reported as primary results. *A. baumannii* had only 1 susceptible genome (99.4% resistance prevalence) and is excluded from AUC comparisons.

### Statistical analysis

Python 3.11 with pandas, scipy, and matplotlib. Fisher's exact test compared composite transposon rates between species pairs. All confidence intervals are 95%.

---

## Results

### Study cohort

After quality control, the cohort comprised 492 *K. pneumoniae*, 266 *E. coli*, and 168 *A. baumannii* genomes (926 total; Table S1). All cohorts consisted of Chinese clinical isolates with assembly N50 ≥ 50 kb.

### Carbapenem resistance prevalence and gene types

Carbapenem resistance was identified in 200/492 (40.7%) *K. pneumoniae* genomes, 28/266 (10.5%) *E. coli* genomes, and 167/168 (99.4%) *A. baumannii* genomes (Fig. 4, Table 1). Dominant carbapenemase types reflected species-specific epidemiology: KPC-2 (54.5%) and NDM variants (44%) in *K. pneumoniae*; NDM-5 (25/29 AMR loci, 86.2%) with OXA-48 (6.9%), NDM-1 (3.4%), and KPC-2 (3.4%) in *E. coli*; OXA-23 (83.2%) and intrinsic OXA-51 (98.8%) in *A. baumannii*. The expanded *E. coli* cohort (266 genomes vs 237 in the initial analysis) revealed a more diverse carbapenemase repertoire than the NDM-5-exclusive first cohort, including 3 genomes carrying KPC-2 or OXA-48 family enzymes, consistent with emerging carbapenemase diversity in Chinese clinical *E. coli* [reference].

### IS-mediated composite transposon frequency: a gradient across species

IS element flanking analysis revealed a >13-fold gradient in composite transposon frequency across the three species (Fig. 1, Table 2). In *K. pneumoniae*, 47/50 (94.0%; 95%CI 83.5–98.8%) acquired carbapenemase loci carried IS elements on both flanks (composite transposon architecture), confirming the near-universal role of IS26 in KPN resistance gene organisation [5,6]. In *E. coli*, 18/29 (62.1%; 95%CI 44.0–77.3%) carbapenemase loci formed composite transposons—substantially higher than our initial analysis of 13 NDM-5-only genomes (2/13, 15.4%), and reflecting a more diverse carbapenemase and IS ecology in the expanded cohort. In *A. baumannii* acquired carbapenemase loci (n=146, excluding intrinsic OXA-51), only 10/146 (6.8%) formed composite transposons; 8 additional loci carried single ISAba insertions, predominantly upstream (7/8, 87.5%).

All three pairwise comparisons between species were statistically significant (Fisher's exact: KPN vs ECO p = 3.2×10⁻³; KPN vs ABA p < 10⁻¹⁵; ECO vs ABA p = 3.7×10⁻⁸). These divergences are not attributable to differences in search window (±10 kb identical across species) or GFF3 annotation quality.

### IS family composition reflects species-specific IS ecology

The IS family landscape flanking carbapenem resistance genes differed fundamentally by species (Fig. 4, Table 2). In *K. pneumoniae*, IS6 family elements (IS26 and IS26-like elements) predominated, consistent with the established IS26-centric resistance ecology [5,6].

In *E. coli*, IS element composition was strikingly diverse across 137 IS–AMR element pairs: IS6 family (IS26, IS257, IS1353; 32/137, 23.4%) was the most prevalent single category, followed by Tn3 family elements (28/137, 20.4%), IS_other (19/137, 13.9%), IS5 family (18/137, 13.1%), IS1 family (13/137, 9.5%), and IS_unknown (16/137, 11.7%). This IS family diversity in *E. coli* carbapenem resistance contexts—with IS26 as one of several co-contributing elements rather than a single dominant element—reflects the broader IS ecology of clinical *E. coli*, which carries a richer IS repertoire than *K. pneumoniae* [reference]. Among the AMR gene subsets: NDM-5 contexts showed the broadest IS diversity (all IS families represented), while KPC-2 and OXA-48 contexts in *E. coli* disproportionately involved IS6/IS26 elements (not shown; n too small for quantitative comparison).

In *A. baumannii*, ISAba family elements dominated IS–AMR contexts (26/62 classified IS–AMR pairs, 41.9%), with IS_other (17/62, 27.4%) and IS_unknown (5/62, 8.1%) also contributing. Critically, ISAba elements displayed a strong upstream positional bias: 24/26 (92.3%) ISAba–AMR pairs placed the IS element upstream of the carbapenemase gene. This upstream enrichment is consistent with the established ISAba1 promoter activation mechanism for OXA-23 [7,8], and contrasts with the bilateral symmetric distribution of IS26 around KPN resistance genes in composite transposons. This positional divergence reveals a mechanistic distinction: IS26 in KPN functions primarily as a structural mobiliser (bilateral bracketing → block transposition), whereas ISAba1 in ABA functions primarily as a transcriptional activator (upstream insertion → P2 promoter delivery).

### IS6 copy number predicts resistance in both enterobacterial species

**IS6 burden in *K. pneumoniae*.** IS6 genomic copy number was significantly higher in resistant genomes (median 19 copies/genome, IQR 11–27) than susceptible genomes (median 10 copies/genome, IQR 6–16): AUC = 0.667, Cliff's δ = 0.331, Mann–Whitney p = 4.2×10⁻¹⁰ (Fig. 2, Fig. 3). This association reflects IS26-mediated iterative resistance gene acquisition, in which each transposition event deposits additional IS26 copies flanking expanding resistance cassettes [6]. **Note on cohort-specific AUC values**: A companion single-species analysis of a representative 2025 *K. pneumoniae* cohort (n=95 resistant, n=50 susceptible; selected for high-quality IS annotation completeness and carbapenemase gene diversity) yielded a higher AUC of 0.976. The attenuation to 0.667 in the present full 492-genome cohort reflects two factors: (i) the broader cohort includes susceptible genomes with variable GFF3 IS annotation completeness, underestimating IS6 counts in some susceptible isolates; (ii) population-level heterogeneity, as the full cohort encompasses a wider range of K. pneumoniae lineages—some susceptible strains may carry IS26 for transposition functions unrelated to carbapenem resistance gene acquisition. The full-cohort AUC (0.667) therefore represents the conservative, population-level predictive value, while the representative-cohort AUC (0.976) captures the upper bound achievable with curated, annotation-complete datasets.

**IS6 burden in *E. coli* — annotation bias and PFAM correction.** Name-based IS6 scanning in *E. coli* initially yielded AUC = 0.773 (resistant median 8 copies, susceptible median 4 copies; p = 2.3×10⁻⁶). PFAM HMM correction (PF01527, IS6 family transposase) reclassified 1,056 IS_unknown features across 266 genomes—of which 959 reclassifications were in susceptible genomes compared to only 97 in resistant genomes (9.9× susceptible enrichment), confirming systematic annotation bias that inflated the apparent resistance–IS6 gap. After correction, the association remained significant (AUC = 0.718, resistant median 11.5 copies, susceptible median 8 copies; Cliff's δ = 0.436; p = 2.3×10⁻⁶; Fig. 2, Fig. 3). The annotation bias—more IS_unknown transposases in susceptible *E. coli* genomes that are PFAM-confirmed IS6 family—reflects a general underannotation of IS26 family elements in NCBI-deposited assemblies when isolates are not annotated by resistance-aware pipelines.

**Relationship between composite transposon frequency and IS burden.** The significant IS6 burden–resistance association in *E. coli* (62% composite transposon rate) is consistent with partial IS26-mediated resistance gene assembly in this species: a substantial fraction of *E. coli* carbapenemase genes are embedded in IS26-flanked composite transposons, and these isolates accumulate IS6 copies during iterative IS26-mediated rearrangements. However, the lower composite transposon rate in *E. coli* vs *K. pneumoniae* (62% vs 94%) and the presence of multiple IS family types in *E. coli* carbapenem resistance contexts indicate that IS26-independent plasmid platforms (e.g., IncX3/Tn125-like structures for NDM-5) also contribute to *E. coli* carbapenem resistance gene dissemination without IS6 copy number amplification.

**ISAba burden in *A. baumannii*.** ISAba family copy number was quantified but could not be reliably compared across resistance status groups given the single susceptible genome in the cohort (1/168). PFAM correction identified only 12 reclassifications (resistant 10, susceptible 2), with no evidence of annotation bias in *A. baumannii*. Raw ISAba counts (median 6 copies in resistant genomes) confirm ISAba family proliferation in the resistome-colonised *A. baumannii* background, consistent with published observations of ISAba1 copy number expansion in multi-drug resistant lineages [reference].

---

## Discussion

This study provides the first three-species whole-genome comparative analysis of IS element architecture around carbapenem resistance genes in Chinese clinical isolates, revealing that IS-mediated resistance gene mobilisation is governed by species-specific IS ecologies with a clear Enterobacteriaceae–*Acinetobacter* divergence.

**Composite transposons are shared within Enterobacteriaceae but are not a universal mechanism.** Our finding that 62.1% of *E. coli* carbapenemase loci form composite transposons—comparable in kind, if lower in frequency, to the 94.0% rate in *K. pneumoniae*—establishes that IS26-mediated bilateral bracketing is broadly operative across the enterobacterial carbapenem resistome, not uniquely associated with *K. pneumoniae*. The higher *K. pneumoniae* composite transposon rate likely reflects the greater evolutionary pressure on KPN in the Chinese hospital environment, where CRKP strains have been accumulating and rearranging resistance cassettes for decades, and the near-universal co-occurrence of IS26 with multiple resistance genes (KPC, NDM, colistin resistance) in high-risk clones such as ST11 [reference]. The lower but substantial *E. coli* rate (62%) reflects a mixed population: NDM-5/IncX3 strains (where the complete Tn125-like element on IncX3 replicons mobilises intact, without requiring IS26-mediated rearrangement) co-exist with *E. coli* strains that have recruited IS26 composite transposons for KPC-2 or OXA-48 gene organisation, analogous to the KPN pattern.

**The *A. baumannii* IS ecology is mechanistically distinct.** The 6.8% composite transposon rate in *A. baumannii* acquired carbapenemases, combined with the 92.3% upstream positional bias of ISAba elements, defines *A. baumannii* as employing a qualitatively different IS strategy. ISAba1 does not primarily function as a bilateral structural bracket; it inserts upstream of OXA-23 to deliver its outward-reading P2 promoter, transcriptionally activating an otherwise cryptic carbapenemase [7,8]. The bilateral composite transposon pattern (6.8% of acquired loci) represents a minority outcome, not the primary mechanism. This mechanistic divergence—IS26 as structural mobiliser in Enterobacteriaceae vs ISAba1 as transcriptional activator in *A. baumannii*—reflects the convergent evolution of two distinct IS families for the common selective function of carbapenemase dissemination, operating through completely different molecular mechanisms. The absence of IS26 in *A. baumannii* (IS26 is predominantly an enterobacterial element) and the absence of ISAba1 from the Enterobacteriaceae IS repertoire [reference] enforce this mechanistic segregation.

**IS6 genomic burden is a valid resistance predictor across Enterobacteriaceae.** A key finding of this study is that IS6 copy number significantly predicts carbapenem resistance in *both* *K. pneumoniae* (AUC 0.667, p=4.2×10⁻¹⁰) and *E. coli* (corrected AUC 0.718, p=2.3×10⁻⁶). This extends the IS26 burden–resistance association, previously established in KPN [reference], to *E. coli* and suggests that IS6 copy number measurement could serve as a species-agnostic biomarker for resistance risk stratification across Enterobacteriaceae. However, the PFAM correction reveals an important methodological caveat: IS6 copy number estimates from NCBI GFF3 annotations are systematically underestimated in susceptible genomes (9.9-fold excess of IS6-reclassified IS_unknown in susceptible *E. coli*), generating an apparent resistance–IS6 gap that overcorrects toward IS6 being a stronger predictor than it actually is. After correction, the significant association is preserved (AUC 0.718), but the effect size is more conservative (Cliff's δ 0.436 vs 0.546 uncorrected). Annotation quality must therefore be addressed before applying IS6 burden as a clinical predictor.

**IS element diversity in *E. coli* resistance contexts.** The co-involvement of Tn3 family elements (20.4%), IS5 (13.1%), and IS1 (9.5%) alongside IS6 (23.4%) in *E. coli* carbapenem resistance gene contexts reveals a more complex IS ecology than the IS26-centric view applicable to *K. pneumoniae*. Tn3 family elements, which include the Tn1-Tn3 group and many ISs that use an encapsulase-resolvase mechanism [9], are major mobilisation platforms in IncX3 and IncF plasmids commonly carrying NDM in *E. coli* [reference]. IS5 and IS1 family elements contribute to cassette rearrangements in IncF resistance plasmids [reference]. This IS diversity reflects *E. coli*'s broad plasmid ecology, which accommodates multiple plasmid incompatibility groups (IncX3, IncF, IncI, IncL) in contrast to the narrow IncF/IncFIA dominance in Chinese CRKP [reference].

**Epidemiological and surveillance implications.** The species-specific IS ecology characterised here has direct implications for AMR genomic surveillance. First, composite transposon frequency—often used as a proxy for resistance gene "mobility"—overestimates *E. coli* resistance gene mobility if calibrated on KPN data (62% vs 94% composite rate). Second, IS26 copy number surveillance is most informative in *K. pneumoniae* and *E. coli*, but not in *A. baumannii*, where a different IS marker (ISAba copy number) would be more appropriate. Third, IS6 annotation completeness varies by species and genome annotation pipeline, requiring PFAM HMM correction for robust comparative studies. The results here argue for species-adapted IS surveillance frameworks in clinical genomic epidemiology.

**Limitations.** The *E. coli* resistant cohort remains modest (n=28 carbapenem-resistant genomes, 10.5%), limiting statistical precision for composite transposon rate estimation (95%CI 44.0–77.3%). Targeted retrieval of carbapenem-resistant *E. coli* strains would reduce this interval. The *A. baumannii* susceptible comparator was insufficient for AUC analysis (n=1). All IS detection relied on GFF3 annotation quality; while PFAM correction was applied, HMM sensitivity for ISAba subtypes (ISAba2/3, IS4 family) may differ from IS6 family elements. MLST, plasmid replicon typing, and cgMLST would enable association of composite transposon frequency with specific high-risk clones and plasmid types.

---

## Conclusions

We demonstrate a species-specific gradient in IS-mediated composite transposon architecture across Gram-negative clinical pathogens: *K. pneumoniae* (94.0%) > *E. coli* (62.1%) > *A. baumannii* (6.8%), with all three pairwise comparisons statistically significant. Composite transposon formation is the dominant resistance gene organisation strategy in both Enterobacteriaceae species, but *E. coli* displays greater IS element diversity (IS6, Tn3, IS5, IS1) than *K. pneumoniae* (IS26-centric). *Acinetobacter baumannii* deploys a mechanistically distinct ISAba1 upstream promoter insertion strategy, not structural bilateral bracketing, to activate acquired carbapenemases. IS6 genomic copy number is a significant resistance predictor in both *K. pneumoniae* (AUC 0.667) and *E. coli* (PFAM-corrected AUC 0.718), but PFAM HMM correction is essential to remove annotation bias that inflates the raw signal in *E. coli*. These findings establish that IS element surveillance in carbapenem-resistant Gram-negative pathogens must be species-adapted in both marker selection and annotation validation.

---

## Declarations

**Funding**: This work received no external funding.

**Conflicts of interest**: The author declares no conflicts of interest.

**Data availability**: All genome assembly accession numbers are publicly available from NCBI GenBank. Analysis code is openly available at https://github.com/jiayu6954-sudo/comparative-amr-genomics (MIT License). Processed results, figures, and the manuscript are deposited at Zenodo (https://doi.org/10.5281/zenodo.19665193). Genome accession lists are provided in Tables S1–S3.

**Ethics statement**: This study used exclusively publicly available, de-identified genomic sequence data from NCBI GenBank and did not involve human subject participation. No ethics approval was required.

---

## References

1. WHO. Bacterial Priority Pathogens List 2024. World Health Organization; 2024.
2. Tacconelli E, et al. Discovery, research, and development of new antibiotics: the WHO priority list. *Lancet Infect Dis*. 2018;18(3):318–327.
3. Hu F, et al. Resistance trends among clinical isolates in China reported from CHINET surveillance. *Clin Microbiol Infect*. 2016;22(S1):S9–S14.
4. Sheng ZK, et al. Epidemiology of carbapenem-resistant *Klebsiella pneumoniae* infections. *J Antimicrob Chemother*. 2020;75(6):1468–1477.
5. Sheppard AE, et al. Nested Russian doll-like genetic mobility drives rapid dissemination of blaKPC. *Antimicrob Agents Chemother*. 2016;60(6):3767–3778.
6. Harmer CJ, Hall RM. IS26-mediated formation of transposons carrying antibiotic resistance genes. *mSphere*. 2016;1(6):e00349-16.
7. Mugnier PD, et al. Phylogeny of ISAba1, an insertion sequence widely distributed among multidrug-resistant *Acinetobacter baumannii*. *PLoS One*. 2009;4(11):e7632.
8. Corvec S, et al. ISAba1 insertion upstream of OXA-51 carbapenem-hydrolyzing β-lactamase gene of *Acinetobacter baumannii*. *Antimicrob Agents Chemother*. 2007;51(4):1524–1525.
9. Roberts AP, Chandler M, Courvalin P, et al. Revised nomenclature for transposable genetic elements. *Plasmid*. 2008;60(3):167–173.
10. Wyres KL, Holt KE. *Klebsiella pneumoniae* as a key trafficker of drug resistance genes from environmental to clinically important bacteria. *Curr Opin Microbiol*. 2018;45:131–139.

---

## Tables

**Table 1. Study cohort, resistance prevalence, and IS burden predictive statistics across three species.**

| Parameter | *K. pneumoniae* | *E. coli* | *A. baumannii* |
|-----------|----------------|-----------|---------------|
| Total genomes (QC-passed) | 492 | 266 | 168 |
| Carbapenem-resistant genomes | 200 (40.7%) | 28 (10.5%) | 167 (99.4%) |
| Dominant carbapenemase | KPC-2, NDM | NDM-5, OXA-48 | OXA-23, OXA-51 |
| AMR loci analysed (composite rate denom.) | 50 (acquired) | 29 | 146 (acquired, excl. OXA-51) |
| Composite transposon rate | 94.0% [83.5–98.8%] | 62.1% [44.0–77.3%] | 6.8% |
| Dominant IS family (flanking context) | IS26 (IS6) | IS6/IS26, Tn3, IS5, IS1 | ISAba (ISAba1) |
| IS burden AUC (PFAM-corrected) | 0.667 (p=4.2×10⁻¹⁰) | 0.718 (p=2.3×10⁻⁶) | N/A (n=1 susceptible) |
| Cliff's δ (PFAM-corrected) | 0.331 | 0.436 | N/A |
| Median IS copies (resistant / susceptible) | 19 / 10 | 11.5 / 8 (corrected) | 6 / 0 |

**Table 2. IS element flanking architecture at each species' acquired carbapenemase loci.**

| Species | COMPOSITE | SINGLE_UPSTREAM | SINGLE_DOWNSTREAM | NO_IS | Total loci |
|---------|-----------|----------------|------------------|-------|-----------|
| *K. pneumoniae* | 47 (94.0%) | 2 (4.0%) | 1 (2.0%) | 0 (0.0%) | 50 |
| *E. coli* | 18 (62.1%) | 4 (13.8%) | 3 (10.3%) | 4 (13.8%) | 29 |
| *A. baumannii* (acquired) | 10 (6.8%) | 8 (5.5%) | 4 (2.7%) | 124 (84.9%) | 146 |

Fisher's exact test (pairwise composite rate comparison):
- KPN vs ECO: p = 3.2×10⁻³ (OR = 9.6, 95%CI 2.0–110)
- KPN vs ABA: p < 10⁻¹⁵
- ECO vs ABA: p = 3.7×10⁻⁸ (OR = 22.3, 95%CI 7.8–74)
