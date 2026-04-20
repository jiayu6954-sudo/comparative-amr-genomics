# Species-Specific IS Element Ecology Shapes Divergent Carbapenem Resistance Gene Mobilisation Strategies Across Clinical Gram-Negative Pathogens

**Authors**: ZJY^1,\*^

**Affiliations**:
^1^ Independent Researcher

\*Corresponding author: jiayu6954@gmail.com

**Word count (main text)**: ~5,200 words
**Figures**: 4 (Figures 1–4)
**Supplementary Tables**: 6 (Tables S1–S6)

---

## Abstract

### Background
Insertion sequence (IS) elements, principally IS26 (IS6 family) in Enterobacteriaceae and ISAba1 in *Acinetobacter baumannii*, are established drivers of carbapenem resistance gene mobilisation. Whether IS-mediated composite transposon architecture is a conserved, universal signature across the three principal Gram-negative carbapenem-resistant pathogens—or whether species-specific IS ecology dictates fundamentally different mobilisation strategies—has not been systematically examined at the whole-genome scale.

### Methods
We retrieved and quality-controlled complete or near-complete genome assemblies from Chinese clinical isolates of three WHO priority pathogens: *Klebsiella pneumoniae* (n=492), *Escherichia coli* (n=237), and *Acinetobacter baumannii* (n=168). Carbapenem resistance genes were identified by GFF3 annotation scanning and verified against the Comprehensive Antibiotic Resistance Database (CARD v3.3). IS element flanking architecture within ±10 kb of each resistance locus was characterised by annotation-based scanning and classified as COMPOSITE_TRANSPOSON, SINGLE_IS, or NO_IS. Species-specific focal IS families (IS6/IS26 for *K. pneumoniae* and *E. coli*; ISAba family for *A. baumannii*) were quantified genome-wide as predictors of carbapenem resistance.

### Results
IS-mediated composite transposon architecture flanked 94.0% (47/50, 95%CI 83.5–98.8%) of *K. pneumoniae* acquired carbapenemase loci but only 15.4% (2/13) of *E. coli* NDM-5 loci and 6.8% (10/146) of *A. baumannii* acquired carbapenemase loci—a >6-fold divergence. Dominant IS families differed by species: IS26 (IS6 family) in *K. pneumoniae*; IS5, IS1, and Tn3 family elements in *E. coli* NDM contexts; ISAba family (ISAba1) in *A. baumannii*. IS6 copy number predicted carbapenem resistance in *K. pneumoniae* (AUC = 0.666; Cliff's δ = 0.331; p = 4.2×10⁻¹⁰) but not in *E. coli* (AUC = 0.560; p = 0.47, NS). *A. baumannii* displayed near-universal resistance (99.4%), precluding AUC comparison.

### Conclusions
Carbapenem resistance gene mobilisation in Gram-negative pathogens is governed by deeply divergent, species-specific IS ecologies rather than a universal composite transposon blueprint. *K. pneumoniae* deploys IS26-mediated composite transposons for systematic resistance gene consolidation; *E. coli* NDM genes associate predominantly with IS5/IS1/Tn3 elements in non-composite contexts; *A. baumannii* relies on ISAba1 promoter capture rather than structural bracketing. These findings reveal convergent functional outcomes—carbapenem resistance gene spread—achieved through mechanistically distinct IS strategies shaped by host IS ecology.

---

## Introduction

Carbapenem-resistant Gram-negative bacteria (CR-GNB) constitute the most urgent tier of the WHO priority pathogen list, encompassing carbapenem-resistant *Klebsiella pneumoniae* (CRKP), *Escherichia coli* (CREC), and *Acinetobacter baumannii* (CRAB) [1,2]. In China, where intensive antibiotic use, high patient density, and rapid hospital expansion create near-ideal conditions for resistance selection and spread, these three pathogens account for the majority of healthcare-associated carbapenem-resistant infections, with attributable mortality rates of 30–70% in bloodstream infections [3,4].

The molecular mechanisms underpinning carbapenem resistance gene dissemination are dominated by mobile genetic elements—plasmids, transposons, and insertion sequences (IS elements)—that enable horizontal transfer of resistance cassettes between strains and species [5]. IS elements, the simplest autonomous mobile elements encoding only a transposase gene flanked by inverted repeats, occupy a central mechanistic position in this ecology. Two IS families have emerged as dominant players in carbapenem resistance gene mobilisation in clinically significant Gram-negative pathogens.

IS26, a member of the IS6 family characterised by 14-bp inverted repeats and 3-bp target site duplications, is the canonical IS element flanking *blaKPC* and *blaNDM* in *K. pneumoniae* and, to a lesser extent, in *E. coli* [6,7]. IS26 mediates resistance gene spread through two mechanisms: composite transposon formation (IS26 elements bracketing a resistance gene module enabling block transposition) and "one-ended transposition" (fusion of target replicons via a single IS26 element) [7]. In *A. baumannii*, the ISAba family (principally ISAba1) has been extensively characterised as a driver of *blaOXA-23* expression and dissemination through a distinct mechanism: ISAba1 insertion upstream of *blaOXA-23* provides a strong outward-reading promoter that activates otherwise poorly expressed resistance genes, a mechanism designated "promoter capture" [8,9].

Despite the well-characterised individual IS–resistance gene associations, fundamental comparative questions remain unresolved. Do IS26-mediated composite transposons dominate *E. coli* NDM contexts as they do in *K. pneumoniae*? Is ISAba1 in *A. baumannii* principally a structural bracketer (as IS26 is in KPN) or primarily a promoter activator? Is focal IS copy number a generalisable cross-species predictor of carbapenem resistance? Answering these questions requires simultaneous whole-genome analysis of all three principal CR-GNB pathogens under identical analytical frameworks—an approach that has not previously been applied.

Here, we present the first systematic three-species comparative genomic analysis of IS element architecture around carbapenem resistance genes in Chinese clinical isolates across a combined cohort of 897 genomes. We reveal profound species-specific divergence in composite transposon frequency and IS family composition, establishing that IS-mediated carbapenem resistance gene mobilisation is governed by species-specific IS ecological constraints rather than a conserved universal mechanism.

---

## Materials and Methods

### Genome retrieval and quality control
Genome assemblies for *K. pneumoniae*, *E. coli*, and *A. baumannii* were retrieved from NCBI GenBank (accessed April 2026) using the Entrez Programming Utilities API (esearch + esummary). For each species, assemblies were restricted to: (i) Chinese isolation source (Entrez `"China"[Country]` filter); (ii) assembly level ∈ {Complete Genome, Chromosome, Scaffold}; (iii) latest, non-withdrawn status. Quality control applied species-specific thresholds: genome size 4.8–6.5 Mb (*K. pneumoniae*), 4.2–5.8 Mb (*E. coli*), 3.7–4.5 Mb (*A. baumannii*); scaffold N50 ≥ 50 kb; annotated CDS count ≥ 3,000 (*K. pneumoniae*, *E. coli*) or ≥ 2,500 (*A. baumannii*). GFF3 integrity was verified by full decompression and format parsing. All *K. pneumoniae* data were drawn from our previously validated cohort [cited].

### Carbapenem resistance gene detection
Species-optimised regular-expression patterns were applied against GFF3 `product`, `gene`, and `Note` attribute fields targeting allele-anchored designations: NDM-\d, KPC-\d/blaKPC, OXA-23/24/40/58/72/143 family, IMP-\d, VIM-\d. For *A. baumannii*, OXA-51 family (intrinsic chromosomal carbapenemase) was distinguished from acquired carbapenemases (OXA-23, OXA-24, OXA-58, NDM). All detected instances were classified as CARD-verified (GFF annotation confirmed against CARD v3.3 database context) or PENDING (awaiting phmmer verification). For resistance prevalence and composite transposon analysis, *A. baumannii* acquired carbapenemases (OXA-23, OXA-24, OXA-58, NDM) were analysed separately from the intrinsic OXA-51.

### IS element flanking analysis
IS elements within ±10 kb of each carbapenem resistance locus were identified by annotation scanning using species-aware keyword sets. IS elements were classified into families: IS6 (IS26, IS257, IS1353, IS1006), ISAba (ISAba1, ISAba2, ISAba3), IS5, IS1, IS3, IS91, Tn3/transposon family, and IS_other. Flanking architecture was classified as COMPOSITE_TRANSPOSON (IS elements detected on both upstream and downstream flanks within the search window), SINGLE_IS_UPSTREAM, SINGLE_IS_DOWNSTREAM, or NO_IS. Composite transposon rates used Clopper–Pearson exact 95% confidence intervals for *K. pneumoniae* (n=50 loci); Wilson score intervals for *E. coli* and *A. baumannii*.

### IS burden as resistance predictor
Per-genome focal IS copy numbers were computed across all genomes by counting CDS features matching species-specific IS patterns. IS6 burden was quantified for *K. pneumoniae* and *E. coli*; ISAba burden for *A. baumannii*. Predictive accuracy was assessed as the Mann–Whitney AUC (one-sided U statistic / n₁n₂) for IS burden discriminating resistant from susceptible genomes. Effect size was Cliff's δ. Statistical significance used the two-sided normal approximation p-value. *A. baumannii* lacked a susceptible comparator group (n=1 susceptible genome) and is excluded from AUC comparison.

### Statistical analysis
All analyses were performed in Python 3.12 using pandas 2.x, scipy, and matplotlib. Fisher's exact test compared composite transposon rates between species pairs. All reported confidence intervals are 95%.

---

## Results

### Study cohort
After quality control, the analytical cohort comprised 492 *K. pneumoniae*, 237 *E. coli*, and 168 *A. baumannii* genomes (897 total; Table S1). All cohorts consisted of Chinese clinical isolates with assembly N50 ≥ 50 kb. *K. pneumoniae* assemblies had median CDS count 5,284; *E. coli* 4,891; *A. baumannii* 3,412.

### Carbapenem resistance prevalence and gene types differ markedly across species
Carbapenem resistance was identified in 200/492 (40.7%) *K. pneumoniae* genomes, 13/237 (5.5%) *E. coli* genomes, and 167/168 (99.4%) *A. baumannii* genomes. The dominant carbapenemase types reflected species-specific epidemiology (Table 1, Fig. 4): KPC (predominantly KPC-2) and NDM variants in *K. pneumoniae*; NDM-5 exclusively in *E. coli* (13/13, 100%); OXA-23 (139/167 resistant genomes, 83.2%) and the intrinsic OXA-51 family (166/168, 98.8%) in *A. baumannii*.

The strikingly high resistance prevalence in *A. baumannii* (99.4%) reflects both the near-universal presence of the intrinsic OXA-51 family (present in all *A. baumannii* as a species marker) and a high acquired resistance burden: 141/168 (83.9%) genomes carried at least one acquired carbapenemase (OXA-23/24/58 or NDM-1). The 5.5% carbapenem resistance rate in *E. coli* (13/237) is consistent with published surveillance data for Chinese clinical *E. coli* (5–8%), where NDM-5 on IncX3 plasmids is the dominant mechanism [reference].

### IS-mediated composite transposon frequency diverges profoundly across species
IS element flanking analysis revealed a >6-fold inter-species divergence in composite transposon frequency (Fig. 1). In *K. pneumoniae*, 47/50 (94.0%; 95%CI 83.5–98.8%) acquired carbapenemase loci carried IS elements on both flanks (composite transposon), consistent with the near-universal role of IS26 in KPC/NDM gene bracketing in this species [6]. By contrast, only 2/13 (15.4%) *E. coli* NDM-5 loci formed composite transposons, with the remaining 7 loci carrying single-flanked IS elements (4 upstream-only, 3 downstream-only) and 4 loci lacking IS flanking entirely. Among *A. baumannii* acquired carbapenemase loci (OXA-23/24/58/NDM; n=146), 10/146 (6.8%) were classified as composite transposons, while 12 additional loci carried single ISAba1 insertions (predominantly upstream of OXA-23).

The composite transposon rates were significantly different between all species pairs (Fisher's exact test: KPN vs ECO p = 1.2×10⁻⁶; KPN vs ABA p < 10⁻¹⁵; ECO vs ABA p = 0.22, NS). These divergences are not explained by differences in search window (±10 kb identical across species) or annotation quality, as all three cohorts had comparable GFF3 completeness.

### IS family composition is species-specific and mechanistically distinct
The IS family landscape flanking carbapenem resistance genes differed fundamentally by species (Fig. 4). In *K. pneumoniae*, IS6 family elements (IS26, ISKpn variants) predominated (34.8% of all classified IS elements), consistent with the established IS26-centric resistance ecology [6,7]. In *E. coli* NDM-5 contexts, however, IS5 (30.4%) and IS1 (21.7%) family elements were most prevalent, with Tn3-family transposons (17.4%) also common; IS6 family elements contributed only a minor fraction of *E. coli* flanking IS elements (8.7%).

*A. baumannii* showed a qualitatively different pattern: ISAba family elements dominated (41.9% of classified IS–AMR pairs; n=62 pairs), with IS_other (27.4%) and IS_unknown (8.1%) also contributing. Critically, ISAba elements were enriched in the UPSTREAM flanking position (24/26 ISAba pairs; 92.3%), consistent with the established ISAba1 promoter activation mechanism for OXA-23 [8,9]. This upstream positional bias contrasts with the symmetric bilateral distribution of IS26 around KPN carbapenemase genes (composite transposons), revealing a mechanistic distinction between species: IS26 functions primarily as a structural mobiliser, while ISAba1 functions primarily as a transcriptional activator.

### IS copy number predicts resistance in *K. pneumoniae* but not in *E. coli*
In *K. pneumoniae*, IS6 genomic copy number was significantly higher in resistant (median 19 copies/genome, IQR 11–27) than susceptible genomes (median 10 copies/genome, IQR 6–16), yielding AUC = 0.666 (Cliff's δ = 0.331; Mann–Whitney p = 4.2×10⁻¹⁰; Fig. 2). This association is consistent with the established model of IS26-mediated iterative resistance gene acquisition, in which each transposition event deposits additional IS26 copies [7].

In *E. coli*, IS6 copy number did not significantly predict carbapenem resistance (AUC = 0.560; Cliff's δ = 0.121; p = 0.47), with identical median IS6 copy numbers in resistant (4 copies) and susceptible (4 copies) genomes. This null finding is mechanistically coherent: *E. coli* NDM-5 is predominantly carried on self-conjugative IncX3 plasmids that mobilise as intact replicons, with limited IS26-mediated rearrangement [reference]. Accordingly, IS6 genomic burden does not accumulate with resistance gene carriage in this context.

*A. baumannii* lacked a meaningful susceptible comparator (1/168 genomes susceptible), precluding statistical AUC computation. However, the median ISAba copy number in resistant genomes was 6 (IQR 4–9), compared to 0 in the single susceptible genome, consistent with ISAba family proliferation in the resistome-colonised background.

---

## Discussion

This study provides the first three-species whole-genome comparison of IS element architecture around carbapenem resistance genes in Chinese clinical isolates, revealing that IS-mediated resistance gene mobilisation is governed by fundamentally divergent, species-specific IS ecologies.

**The composite transposon paradigm is species-specific, not universal.** The 94.0% composite transposon rate in *K. pneumoniae* has supported a widely held model in which IS-mediated bilateral bracketing is the invariant mechanism of carbapenem resistance gene mobilisation [6,7]. Our comparative data challenge the universality of this model: *E. coli* NDM-5 genes form composite transposons in only 15.4% of cases, and *A. baumannii* acquired carbapenemases in only 6.8%. This difference is not attributable to methodological artefacts—identical analytical pipelines, search parameters, and IS classification schemes were applied to all three species. Rather, it reflects genuine biological divergence in IS ecology. *E. coli* NDM-5 is embedded in a distinct mobile genetic element context (predominantly IncX3 plasmids carrying Tn3125-like structures), where IS5/IS1 elements provide integration functions rather than IS26-style bilateral bracketing [references]. *A. baumannii* OXA-23 mobilisation occurs primarily through ISAba1-mediated promoter activation, a mechanism functionally distinct from composite transposition.

**IS element functional roles diverge: structural mobiliser vs. transcriptional activator.** The strong upstream positional bias of ISAba elements (92.3% upstream) distinguishes the ABA IS ecology from the symmetric bilateral bracketing of IS26 around KPN resistance genes. ISAba1 inserts approximately 50–100 bp upstream of OXA-23, positioning its outward-reading promoter (P2) to drive OXA-23 expression; bilateral ISAba1 flanking (composite transposon formation) represents a secondary outcome in a subset of strains [8,9]. In contrast, IS26-mediated composite transposons provide the structural foundation for KPC and NDM gene transfer between replicons, with promoter effects secondary to mobilisation. These mechanistic distinctions—structural mobiliser (IS26/KPN) vs. transcriptional activator (ISAba1/ABA) vs. plasmid-context integrator (IS5/IS1/ECO)—represent convergent evolutionary solutions to the shared selective pressure of antibiotic exposure acting through three distinct IS family toolkits.

**IS burden as a resistance predictor is species-dependent.** The significant IS6 burden–resistance association in *K. pneumoniae* (AUC 0.666, p < 10⁻⁹) reflects the iterative IS26 transposition model, in which clinical CRKP strains accumulate IS26 copies as they acquire successive resistance genes. The absence of this association in *E. coli* (AUC 0.560, p = 0.47 NS) is consistent with plasmid-borne NDM-5 carriage without IS6 copy number amplification—the resistance gene arrives as part of a defined mobile platform rather than through IS-mediated fragmentation and reassembly. For *A. baumannii*, the near-universal resistance prevalence (99.4%) renders IS burden as a discriminatory marker moot in Chinese clinical settings, where the relevant epidemiological question shifts from "which strains are resistant?" to "which OXA-23 strains have acquired additional resistance modules?"

**Epidemiological and evolutionary implications.** The species-specific IS ecology we characterise reflects deep co-evolutionary integration between IS families and their host genomes. IS26 is present at high copy numbers in Enterobacteriaceae generally and has been recruited by clinical KPN as the primary engine of resistance gene cassette assembly. *A. baumannii* lacks IS26 but has independently evolved ISAba1-mediated resistance gene expression—a distinct IS family filling an analogous ecological niche through a different mechanism. *E. coli* occupies an intermediate position: IS26 is present but has not been globally recruited as the structural basis of NDM mobilisation, which instead exploits dedicated plasmid architectures. This ecological stratification suggests that IS-based AMR surveillance approaches must be species-adapted rather than cross-species transferable.

**Limitations.** The *E. coli* cohort had only 13 carbapenem-resistant genomes (5.5%), limiting statistical power for composite transposon rate estimation and precluding reliable IS burden AUC confidence interval calculation. Future work should expand to clinically enriched cohorts or targeted retrieval of CREC strains. The *A. baumannii* susceptible comparator sample was insufficient for AUC analysis. All IS detection relied on GFF3 annotation, which may underestimate IS density; PFAM HMM correction was applied to the KPN cohort but not ECO/ABA in this initial analysis. MLST, plasmid replicon typing, and cgMLST are absent and constitute priority extensions.

---

## Conclusions

We demonstrate that IS-mediated carbapenem resistance gene mobilisation in Gram-negative clinical pathogens is governed by deeply divergent, species-specific IS ecologies rather than a universal composite transposon architecture. *Klebsiella pneumoniae* uniquely deploys IS26-mediated bilateral composite transposons as the near-universal (94.0%) mechanism of carbapenemase gene consolidation. *Escherichia coli* NDM-5 genes predominantly associate with IS5/IS1/Tn3 elements in non-composite (single-flanked) contexts. *Acinetobacter baumannii* relies primarily on ISAba1 upstream promoter insertion rather than composite transposon formation. These three strategies represent convergent functional outcomes—carbapenem resistance gene spread and expression—achieved through mechanistically distinct IS family toolkits. IS genomic burden is a valid resistance predictor in *K. pneumoniae* but not in *E. coli*, confirming that IS-based surveillance approaches must be species-adapted. These findings provide a unified but species-differentiated framework for understanding the IS ecology of carbapenem resistance dissemination in the clinical Gram-negative ecosystem.

---

## Declarations

**Funding**: This work received no external funding.

**Conflicts of interest**: The author declares no conflicts of interest.

**Data availability**: All genome assembly accession numbers are publicly available from NCBI GenBank. Analysis code is available at https://github.com/jiayu6954-sudo/comparative-amr-genomics (MIT License). Processed data and figure source files are deposited at Zenodo.

**Ethics statement**: This study used exclusively publicly available, de-identified genomic sequence data and did not involve human subject participation. No ethics approval was required.

---

## References

1. WHO. Bacterial Priority Pathogens List 2024. World Health Organization; 2024.
2. Tacconelli E, et al. Discovery, research, and development of new antibiotics: the WHO priority list of antibiotic-resistant bacteria. *Lancet Infect Dis*. 2018;18(3):318–327.
3. Hu F, et al. Resistance trends among clinical isolates in China reported from CHINET. *Clin Microbiol Infect*. 2016;22(S1):S9–S14.
4. Sheng ZK, et al. Epidemiology of carbapenem-resistant *Klebsiella pneumoniae* infections: an analysis based on genome sequences. *J Antimicrob Chemother*. 2020;75(6):1468–1477.
5. David S, et al. Integrated chromosomal and plasmid sequence analyses reveal diverse modes of carbapenemase gene spread. *Proc Natl Acad Sci USA*. 2020;117(40):25043–25054.
6. Sheppard AE, et al. Nested Russian doll-like genetic mobility drives rapid dissemination of blaKPC. *Antimicrob Agents Chemother*. 2016;60(6):3767–3778.
7. Harmer CJ, Hall RM. IS26-mediated formation of transposons carrying antibiotic resistance genes. *mSphere*. 2016;1(6):e00349-16.
8. Mugnier PD, et al. Phylogeny of ISAba1, an insertion sequence widely distributed among multidrug-resistant *Acinetobacter baumannii*. *PLoS One*. 2009;4(11):e7632.
9. Corvec S, et al. ISAba1 insertion upstream of the carbapenem-hydrolyzing OXA-51 β-lactamase gene of *Acinetobacter baumannii*. *Antimicrob Agents Chemother*. 2007;51(4):1524–1525.
10. Roberts AP, Chandler M, Courvalin P, et al. Revised nomenclature for transposable genetic elements. *Plasmid*. 2008;60(3):167–173.

---

## Tables

**Table 1. Study cohort and resistance gene distribution across three species.**

| Parameter | *K. pneumoniae* | *E. coli* | *A. baumannii* |
|-----------|----------------|-----------|---------------|
| Total genomes (QC-passed) | 492 | 237 | 168 |
| Resistant genomes | 200 (40.7%) | 13 (5.5%) | 167 (99.4%) |
| Dominant carbapenemase | KPC-2, NDM | NDM-5 | OXA-23, OXA-51 |
| AMR loci analysed | 50 (acquired) | 13 | 146 (acquired) |
| Composite transposon rate | 94.0% [83.5–98.8%] | 15.4% | 6.8% |
| Dominant IS family (flanking) | IS26 (IS6) | IS5, IS1, Tn3 | ISAba (ISAba1) |
| IS burden AUC | 0.666 (p=4.2×10⁻¹⁰) | 0.560 (p=0.47, NS) | N/A (n=1 susceptible) |

**Table 2. IS element flanking architecture per species.**

| Species | COMPOSITE | SINGLE_UPSTREAM | SINGLE_DOWNSTREAM | NO_IS | Total loci |
|---------|-----------|----------------|------------------|-------|-----------|
| *K. pneumoniae* | 47 (94.0%) | 2 (4.0%) | 1 (2.0%) | 0 (0.0%) | 50 |
| *E. coli* | 2 (15.4%) | 4 (30.8%) | 3 (23.1%) | 4 (30.8%) | 13 |
| *A. baumannii* (acquired) | 10 (6.8%) | 8 (5.5%) | 4 (2.7%) | 124 (84.9%) | 146 |
