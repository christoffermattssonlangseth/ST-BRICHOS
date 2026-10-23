# ST-BRICHOS run catalogue

One entry per analysis run, built 2026-10-23 from notebooks, scripts, result files, git history and project notes. Numbers are fixed once assigned; cite runs as §N. Corrections are added as dated notes.

### 1. SpaceRanger outputs are converted into ten per-section AnnData files
*Notebook:* `notebooks/00/00_create_h5ad.ipynb`

**Rationale.** Infrastructure. This was the first step of the project (initial commit 2025-06-17; logbook 2025-06-17). It brings the 10x Visium SpaceRanger output for every section into Scanpy.

**Technical approach.** For each SpaceRanger sample folder, the notebook reads `outs/filtered_feature_bc_matrix` with `sc.read_10x_mtx` (gene symbols, made unique). It then joins `tissue_positions.csv` or `tissue_positions_list.csv`, whichever exists (the logbook notes that this metadata varied between runs). Full-resolution pixel coordinates go into `.obsm['spatial']`, and the scalefactors and hires image go into `.uns['spatial'][sample]`. One `<sample>.h5ad` is written per section. Unit: Visium spot.

**Results.** The output lists 10 sections saved: P24215_101, P24215_201, P24215_301 and P28052_102, 103, 201, 202, 203, 301, 401. The files went to an absolute path on the author's other machine (`/Users/christoffer/work/karolinska/development/ST_BRICHOS/results/`), not to this repo.

**Caveats.** Input and output paths are hard-coded and live outside the repo. No QC is done at this step. Spots get no treatment labels here; those are added in §3.

### 2. All sections are merged into one combined AnnData with harmonised spatial coordinates
*Notebook:* `notebooks/00/01_concatenation_of_samples.ipynb`

**Rationale.** Infrastructure. The analyses need a single object holding all 10 sections, with per-section image metadata kept and the X/Y axes in the right orientation (logbook 2025-06-17: "corrected for axis flipping").

**Technical approach.** The notebook does an outer-join `ad.concat` of the per-section h5ad files from §1, with `sample_id` labels and barcode suffixes `-<sample>`. `.uns['spatial']` is rebuilt per sample. `.obsm['spatial']` is reset to `[pxl_col, pxl_row]` (X, Y) to fix the flip. The result is written as `combined_samples.h5ad`. As QC, it draws spatial plots of Snap25 (P24215_101) and Plp1 (every section).

**Results.** "Spatial metadata harmonized" and "Combined AnnData written". The Plp1 spatial plots were rendered for all 10 sections. Memory notes give the size of `combined_samples.h5ad` as 25,671 spots × 19,465 genes of raw counts. They also say the file now lives only on the remote analysis machine.

**Caveats.** The output went to an absolute path outside the repo, and the object holds raw counts only. The notebook leaves the axis convention undocumented, beyond the code comment.

### 3. Filtering, normalisation and Leiden clustering give the analysis object, but the treatment map lists a section that does not exist
*Notebook:* `notebooks/00/02_preprocessing.ipynb`

**Rationale.** Infrastructure. This step produces the normalised, clustered object with treatment labels that every later notebook in `notebooks/00/` reads (`ST_BRICHOS_new.h5ad`).

**Technical approach.** The notebook keeps spots with ≥200 genes and genes seen in ≥3 spots, and stores raw counts in `layers['counts']`. It then runs normalize_total (1e4), log1p, PCA (arpack) on all genes (the HVG selection and scaling are commented out), a 15-NN graph on 30 PCs, UMAP and Leiden at resolution 0.5. Treatment is mapped from `sample_id` through a hard-coded dictionary.

**Results.** A UMAP coloured by Leiden cluster and by sample was produced (image only). The object read downstream (§8) has 25,660 spots × 18,751 genes. Set against the 25,671 × 19,465 combined object, this step therefore removed about 11 spots and 714 genes. Treatment labels as stored in the object (printed in nb 04): WT = P24215_101; BRICHOS = P24215_201, P24215_301, P28052_401; PBS = P28052_102, 103, 201, 202, 203, 301. That is 1 WT, 6 PBS and 3 BRICHOS sections.

**Caveats.** The `sample_groups` dictionary also contains `"P24215_401": "PBS"`, a section that §1 never created. It is harmless, since nothing maps to it, but it is a wrong record of the design. The README gives n = 1/3/3 *mice* for WT/PBS/BRICHOS, but the code has 1/6/3 *sections*. No section-to-mouse mapping exists anywhere in the code, so every downstream "replicate" is a section. There is no batch correction. Treatment is fully confounded with Visium slide: all PBS sections are on P28052 and 2 of 3 BRICHOS sections are on P24215 (memory: batch-confound). The global UMAP therefore mixes slide and treatment.

### 4. Collaborator (Cloé) manual region annotations are mapped onto 17,613 of the 25,660 spots
*Notebook:* `notebooks/00/03_regional_annotation.ipynb`

**Rationale.** Infrastructure. The project asks *in which brain regions* BRICHOS acts, so spots need anatomical labels. Cloé made these by hand (logbook 2025-06-18; memory: user profile).

**Technical approach.** The notebook reads every CSV in `data/anatomical_annotations/` (barcode → region) and suffixes barcodes with the sample ID parsed from the file name. It maps them to `obs['region_annotation']`, then harmonises spelling variants into one set of names (for example "Cortex Layers 2/3/4" → "Cortex L2/3/4"; "Olfactory Area/region" → "Olfactory Region"; "Ophtalmic region" → "Optical Region"). It plots the regions on the UMAP and per section, and writes `ST_BRICHOS_region_new.h5ad`.

**Results.** "Mapped region annotation to 17613 spots", so about 8,050 spots have no manual label. Labels include cortical layer groups, Caudoputamen, Striatum, Putamen, Corpus Callosum, Hippocampus, Hypothalamus, Thalamus, Amygdala, Septal Nucleus, Ventricular / Periventricular, Olfactory and Optical regions, and "Excluded".

**Caveats.** The annotation CSVs now in `data/anatomical_annotations/` cover only the seven P28052 sections (102, 103, 201, 202, 203, 301, 401). That would leave WT and both P24215 BRICHOS sections unannotated, which fits the partial coverage. It may be why the labels were later replaced by cluster-derived regions (§8). The region vocabulary is inconsistent across sections (for example cortex layers split differently). An unused `mapper` variable is built from the last file only.

### 5. The plaque-induced gene (PIG) module score is higher in PBS than in WT at spot level
*Notebook:* `notebooks/00/04_region_independent_analysis.ipynb`

**Rationale.** Before asking about BRICHOS, the analysis first checks that the APP-NLGF model shows its expected plaque-associated glial signature: README objectives 1–2 and the logbook from 2025-06-18.

**Technical approach.** Spots in `ST_BRICHOS_new.h5ad`. The notebook shows stacked bars of Leiden-cluster composition per treatment, Wilcoxon markers per cluster and per treatment, and a dotplot of 25 PIG genes by treatment. It then computes `sc.tl.score_genes` scores for the PIG list (Apoe, Clu, Ctsd, C1qa/b/c, C4b, Tyrobp, Trem2, Cst7, …, C3) and a 14-gene oligodendrocyte list (Plp1, Mbp, Mobp, …). PBS and WT spots are compared with a Mann–Whitney U test.

**Results.** Mann–Whitney p for PIG score, PBS vs WT, printed as 0.0000e+00 (below floating-point precision). There are violin plots of PIG and OLIG scores by treatment, but no numbers are printed for BRICHOS or for the OLIG score.

**Caveats.** Spots are treated as independent replicates, so the p-value is pseudoreplicated: there is 1 WT section against 6 PBS. WT and PBS are on different slides. No test of BRICHOS vs PBS was done here.

### 6. Spot-level Wilcoxon volcanoes contrast PBS vs WT and BRICHOS vs WT across the whole tissue
*Notebook:* `notebooks/00/04_region_independent_analysis.ipynb` (§1–2)

**Rationale.** README objective 1 (baseline AD signature, WT vs PBS) and objective 2 (cross-reference with known AD/PIG genes).

**Technical approach.** Spots are subset to WT+PBS or WT+BRICHOS. The notebook runs `sc.tl.rank_genes_groups` (Wilcoxon) on log-normalised `.X`, removes genes with |LFC| ≥ 5, and calls significance at FDR < 0.05 and |log2FC| > 1. Volcano plots label the top 20 (PBS) or top 10 (BRICHOS) genes. A flag marks DE genes that are in the PIG list.

**Results.** Volcano images and styled tables only. Neither the tables nor the gene counts are saved as text. Downstream (§7), 143 genes pass FDR < 0.05 and log2FC > 1 for PBS vs WT.

**Caveats.** Spot-level test with a single WT section, so p-values are strongly inflated. The contrast is also confounded with slide (WT is on P24215 and all PBS on P28052). Nothing is written to `results/`.

### 7. No disease-up gene is "strictly rescued" by BRICHOS, but a sign error in the contrast makes that zero meaningless
*Notebook:* `notebooks/00/04_region_independent_analysis.ipynb` (§3–4)

**Rationale.** README objectives 3–4: find genes that BRICHOS shifts back toward WT, and check whether they overlap AD/PIG genes.

**Technical approach.** Spot-level Wilcoxon of PBS against BRICHOS (`group="PBS", reference="BRICHOS"`). The table is stored as `de_table_BRICHOS_PBS`. `compute_brichos_rescue` takes disease-up genes (PBS vs WT log2FC > 1, FDR < 0.05) and computes delta_vs_WT = LFC(PBS vs WT) − LFC(BRICHOS vs WT). It calls "rescued (strict)" when delta ≥ 0.5, lfc_b_pbs ≤ −0.5 and FDR(b_pbs) < 0.05, and "partial" when delta ≥ 0.5 without the strict criteria. A separate check intersects the PIG list with the BRICHOS-vs-PBS table.

**Results.** Disease-up genes: 143. Rescued (strict): 0. Partial rescue: 24. Not rescued: 119. Top partial-rescue genes by delta_vs_WT: Pmch 4.99, Pomc 1.94, Fmo2 1.59, H2-Aa 1.09, Gbp5 0.96, Cd74 0.93, Col6a4 0.89, H2-Ab1 0.86, H2-Eb1 0.86, Cybb 0.62, Nlrc5 0.59. All of these have lfc_b_pbs > 0 (for example Pmch 4.99), with padj ≈ 0. The "PIG genes modulated by BRICHOS" print-out lists 23 of the 25 PIG genes in the list used (Apoe, Trem2, C1qa/b/c, Tyrobp, …).

**Reading.** The table called `BRICHOS_PBS` is really **PBS vs BRICHOS**: a positive lfc_b_pbs means PBS > BRICHOS. The strict criterion requires lfc_b_pbs ≤ −0.5, so it can only select genes that BRICHOS *raises*. "0 strict rescues" is therefore an artefact of the sign. With the sign corrected, the 20 printed "partial" genes (all lfc_b_pbs ≥ 0.58, all padj_b_pbs ≤ 0.04) and likely most of the other 4 would count as strict rescues. That claim is inferred from the printed values and was not re-run. Several top genes (Pmch, Pomc) are hypothalamic neuropeptides. The memory notes (classifier thread) attribute these genes to cutting-plane/anatomy differences, not treatment. The MHC-II genes (H2-Aa/Ab1/Eb1, Cd74) are a more biologically plausible lead.

**Caveats.** Sign bug as above. The PIG-overlap print-out applies no significance filter: it intersects PIG names present in the full DE table with names in the PBS-vs-WT table, so it is not evidence of modulation. Spot-level pseudoreplication. Fully confounded with slide.

**Next.** Re-run with the contrast sign fixed and at pseudobulk level. The later notebook 12 (§17) does pseudobulk; the notes describe the attenuation/rescue tests in notebooks/03 as handling this properly.

### 8. Leiden clustering at resolution 0.3 is hand-mapped to ten brain "regions" that every later regional analysis uses
*Notebook:* `notebooks/00/06_further_clustering.ipynb`

**Rationale.** Infrastructure. The manual annotation (§4) is partial and inconsistent between sections. A uniform, data-driven region label was needed for all 10 sections (logbook 2025-06-19: "Refined the clustering").

**Technical approach.** On `ST_BRICHOS_region_new.h5ad` (25,660 × 18,751), the notebook runs Leiden at resolutions 0.1, 0.2, 0.3, 0.5, 0.75, 1, 1.5, 2 and 2.5 on the existing graph. It shows t-test markers for leiden_0.3 and Wilcoxon markers for leiden_0.75 (top 20 genes per cluster printed). The 13 leiden_0.3 clusters are mapped by hand to `re_annotation_regions`: 0/11 Caudoputamen, 1 Supragranular layers (spelled "Supragranualar" in code), 2 Infragranular layers, 3/6 Hypothalamus, 4 Corpus callosum, 5 Olfactory areas, 7 Meninges, 8 Hippocampal formation, 9/10 Thalamus, 12 Ventricular area. PIG (52-gene extended list) and OLIG module scores are dotplotted by treatment and by leiden_0.75 within WT and within BRICHOS. The result is written as `ST_BRICHOS_region_subcluster.h5ad`.

**Results.** Ten region labels covering all spots. Leiden_0.75 markers are anatomically coherent: cluster 1 is striatal (Rasd2, Ppp1r1b, Drd2), 9 astro/Apoe/Cst3/Clu/Gfap, 10 oligodendrocyte (Mbp, Plp1, Mobp), 12 Pmch/Hcrt hypothalamic, 13/14 thalamic (Prkcd, Tcf7l2, Slc17a6). Module dotplots are images only.

**Caveats.** The "regions" are expression clusters, not atlas registrations. A cluster such as "Meninges" or "Ventricular area" can mix tissue types, and its composition can differ by section. Clusters were built on the joint, uncorrected embedding, so slide effects could shape the boundaries. The mapping was made by eye.

### 9. Regional PIG scores rise in every region in PBS vs WT, and BRICHOS lowers them only slightly in 6 of 10 regions
*Notebook:* `notebooks/00/05_comparing_regions.ipynb`

**Rationale.** This is the first test of the project's core question: where is the AD signature strongest, and does BRICHOS lower it region by region?

**Technical approach.** Spots in `ST_BRICHOS_region_subcluster.h5ad`, grouped by `re_annotation_regions` × treatment. PIG score (`score_genes`) and OLIG score are averaged per region × treatment over spots, and the differences PBS−WT and BRICHOS−PBS are computed. The notebook also has PIG dotplots by region within PBS and within BRICHOS, a ranking of regions by mean PIG score within PBS, a within-PBS t-test of each region against the rest (Meninges volcano), spatial maps of PIG and OLIG scores and Cst7, and spot-level t-test volcanoes of PBS vs BRICHOS for Meninges and Corpus callosum.

**Results.** Mean PIG score (BRICHOS / PBS / WT; PBS−WT; BRICHOS−PBS):

| Region | BRICHOS | PBS | WT | PBS−WT | BRICHOS−PBS |
|---|---|---|---|---|---|
| Meninges | 0.474 | 0.508 | 0.302 | +0.207 | −0.034 |
| Olfactory areas | 0.215 | 0.259 | 0.053 | +0.206 | −0.044 |
| Supragranular layers | 0.179 | 0.214 | 0.022 | +0.193 | −0.035 |
| Infragranular layers | 0.155 | 0.182 | 0.019 | +0.163 | −0.027 |
| Thalamus | 0.197 | 0.232 | 0.074 | +0.157 | −0.035 |
| Caudoputamen | 0.200 | 0.216 | 0.070 | +0.146 | −0.016 |
| Hypothalamus | 0.283 | 0.282 | 0.137 | +0.145 | +0.001 |
| Hippocampal formation | 0.249 | 0.245 | 0.121 | +0.125 | +0.003 |
| Corpus callosum | 0.435 | 0.385 | 0.268 | +0.118 | +0.050 |
| Ventricular area | 0.481 | 0.406 | 0.360 | +0.047 | +0.075 |

OLIG score, BRICHOS−PBS: positive in Corpus callosum (+0.115), Hippocampal formation (+0.108), Ventricular area (+0.086), Meninges (+0.066), Hypothalamus (+0.060), Olfactory (+0.059) and Thalamus (+0.032); negative in Supragranular (−0.040), Infragranular (−0.018) and Caudoputamen (−0.013). OLIG PBS−WT is largest in Ventricular area (+0.135) and lowest in Corpus callosum (−0.233).

**Reading.** BRICHOS lowers the PIG score by about 11–22 % of the disease increase in 6 regions: supra/infragranular cortex, olfactory areas, meninges, thalamus and caudoputamen. It has no effect in hypothalamus and hippocampus, and the score is higher in CC and ventricular area. This regional profile (cortex/meninges/olfactory most "rescued", CC/hypothalamus least) recurs in the later section-centred analyses (memory: slide-removal thread). Those analyses found it not significant and reversed in the only within-slide BRICHOS section (401).

**Caveats.** Means over spots with no statistical test and no unit of replication. The PIG score definition is ambiguous: the notebook defines both a 52-gene and a 25-gene PIG list, and its execution counts show out-of-order runs. The table may therefore use either version. The BRICHOS dotplot subsets on `region_annotation.notna()` (the manual labels), unlike the PBS one, which may restrict BRICHOS to annotated sections. Slide-confounded. WT is one section.

### 10. In meninges, PBS and BRICHOS show the same DAM/PIG gene induction vs WT
*Notebook:* `notebooks/00/05_comparing_regions.ipynb`

**Rationale.** Meninges had the highest PIG score (§9). This run checks whether BRICHOS dampens individual plaque/glial markers there.

**Technical approach.** Meninges spots only. Scanpy t-tests of PBS vs WT and BRICHOS vs WT are run on log-normalised `.X`, looking at 18 markers (Gfap, B2m, Cst7, Tyrobp, Hexb, Ttr, Cst3, C1qa/b/c, Fcer1g, Ppp1r1b, Arc, Trem2, Cxcl1, Laptm5, Tmsb4x; "Cab" is not a gene). Wilcoxon PBS-vs-BRICHOS and vs-WT rankings are plotted, along with a barplot of LFCs for the markers that are PIG genes. A pseudobulk count matrix per section is built (`pb`) but never tested.

**Results.** PBS vs WT log2FC (padj): Cst7 5.11 (3.5e-34), Tyrobp 2.35 (6.4e-15), Gfap 1.60 (1.4e-12), Trem2 1.76 (5.9e-8), Fcer1g 1.71 (1.1e-8), B2m 1.24 (1.8e-9). BRICHOS vs WT: Cst7 4.87 (3.4e-18), Tyrobp 2.28 (6.6e-11), Gfap 1.57 (1.2e-11), Trem2 1.62 (3.2e-4), Fcer1g 1.39 (1.9e-3), B2m 0.80 (8.6e-3). Cst3 is unchanged in both (padj ≈ 1).

**Reading.** BRICHOS-treated meninges keep almost the full DAM/astro-reactive induction. Only B2m, Fcer1g and Trem2 are slightly lower.

**Caveats.** Spot-level t-test, so p-values are pseudoreplicated. WT is a single section. The "Meninges" label is a cluster (§8). The pseudobulk aggregation was set up but no model was fitted.

### 11. Genes reported dysregulated in Chen et al. 2020 are inspected by treatment (dotplot only)
*Notebook:* `notebooks/00/07_random_from_chen_et_al.ipynb`

**Rationale.** README objective 2 (validate AD signatures against published data). Chen et al. 2020 (spatial transcriptomics in an AD model) reported Cst7 (LFC 1.91), Cd68 (1.70), C1qa (1.04), Slc1a3 (0.57), Clu (0.41) and Mbp (0.37) as dysregulated.

**Technical approach.** A dotplot of these 6 genes by treatment on `ST_BRICHOS_new.h5ad` (spots, variance-scaled).

**Results.** Dotplot image only. No statistics were computed.

**Caveats.** Qualitative only. No LFCs were computed for a direct comparison with the published values.

### 12. Itm2c is unchanged across conditions, while tight-junction genes shift in opposite directions in PBS and BRICHOS
*Notebook:* `notebooks/00/08_expression_of_Bri_tighjunctions_related_genes.ipynb`

**Rationale.** BRICHOS is the chaperone domain of the Bri2/Bri3 proteins (Itm2b/Itm2c). The README TODO "Verify expression levels of Itm2c (WT vs BRICHOS; WT vs PBS) [DONE]" points here. Tight-junction genes (Cldn5, Ocln, F11r, Tjp1) were added as blood–brain-barrier readouts. First committed 2025-09-29, after a meeting that day.

**Technical approach.** Spots in `ST_BRICHOS_region_subcluster.h5ad`. Dotplot by treatment, then Scanpy t-tests on log-normalised `.X` for BRICHOS vs WT, PBS vs WT and BRICHOS vs PBS, reported for the 6 genes.

**Results.** log2FC (padj):
- Itm2c: BRICHOS vs WT 0.0003 (1.0); PBS vs WT −0.015 (0.35); BRICHOS vs PBS 0.016 (0.16). Not changed.
- Itm2b: BRICHOS vs WT 0.069 (3.3e-5); PBS vs WT 0.053 (2.8e-4); BRICHOS vs PBS 0.017 (0.11).
- Cldn5: PBS vs WT −0.377 (7.6e-15); BRICHOS vs WT −0.255 (3.6e-6); BRICHOS vs PBS +0.122 (4.2e-4).
- Tjp1: PBS vs WT +0.184 (9.8e-10); BRICHOS vs PBS −0.185 (4.2e-19).
- F11r: PBS vs WT +0.385 (2.8e-10); BRICHOS vs PBS −0.191 (9.4e-6).
- Ocln: PBS vs WT +0.356 (3.5e-4); BRICHOS vs PBS −0.241 (5.3e-4).

**Reading.** Endogenous Itm2c is not affected by disease or treatment. The four tight-junction genes move partway back toward WT in BRICHOS. The effect sizes are small (|log2FC| < 0.4).

**Caveats.** Spot-level t-tests, so even tiny LFCs come out "significant". Slide-confounded, and WT is one section. No regional breakdown.

### 13. GO and endocytosis-focused enrichment finds no significant endocytosis term in the PBS-vs-WT ranking, and the BRICHOS-vs-PBS top hits are hypothalamic neuropeptides
*Notebook:* `notebooks/00/09_GOterm_analysis.ipynb`

**Rationale.** README TODO 5 asks for DE followed by GO enrichment, with a focus on endocytosis/endocytic cycling, noting possible microglial effects seen in bulk data. Committed 2025-09-29 after a meeting that day.

**Technical approach.** Spots in `ST_BRICHOS_region_subcluster.h5ad`. Scanpy t-tests (log-normalised `.X`) for PBS vs WT, BRICHOS vs WT and BRICHOS vs PBS, written to `de_bulk_treatment.csv` in the notebook folder (not present on disk now). Three enrichment analyses follow:
- g:Profiler GO:BP over-representation on significant up/down genes (FDR ≤ 0.05, |log2FC| ≥ 0.25), filtered to endocytosis/trafficking keywords.
- gseapy prerank of sign(LFC)·−log10(padj) against GO_Biological_Process_2023 (1,000 permutations) for PBS vs WT.
- g:Profiler on the top 500 genes by |LFC| at FDR ≤ 0.05, for each contrast.

**Results.** Endocytosis-filtered prerank for PBS vs WT: no term reaches FDR < 0.3. The best are Positive Regulation of Endocytosis (NES +1.46, FDR 0.31; lead genes B2m, Trem2, Mertk, Fcer1g, Cyba), Presynaptic Endocytosis (NES −1.73, FDR 0.32), Protein Targeting to Lysosome (NES +1.51, FDR 0.33), Synaptic Vesicle Endocytosis (NES −1.79, FDR 0.35) and Synaptic Vesicle Recycling (NES −1.74, FDR 0.35). The direction is microglial/lysosomal up and synaptic-vesicle down. 9.49 % of ranks were tied. In the BRICHOS vs PBS table, the most down-regulated genes are Pmch (log2FC −4.99, padj 0), Hcrt (−4.46), Cox6c and Stmn1. The most "up-regulated" are zero-count artefacts (Vmn2r88, Olfr genes, LFC ≈ 16–18, padj > 0.1). Over-representation results are bubble plots only, with no text tables.

**Caveats.** Spot-level p-values feed the ranking, so they are inflated and tied. Slide-confounded. The BRICHOS-vs-PBS signal is dominated by hypothalamic neuropeptides (Pmch, Hcrt). The memory notes (classifier thread) attribute these to differences in section cutting plane, not treatment.

### 14. Pseudobulk DESeq2 plus decoupler: PBS vs WT shows NF-κB, TNFα, hypoxia and p53 pathway activation and suppressed oxidative phosphorylation
*Notebook:* `notebooks/00/11_decoupler_WT_PBS.ipynb`

**Rationale.** The logic is to move from spot-level to pseudobulk DE, which the notebook's own text calls more statistically robust, then infer transcription-factor and pathway activity for the disease contrast. First committed 2025-10-14.

**Technical approach.** Raw counts from `ST_BRICHOS_region_subcluster.h5ad` are pseudobulked by section × leiden_0.3 cluster with `dc.pp.pseudobulk` (sum), giving 130 profiles. Of these, 91 are WT+PBS, including zero-count cluster×section combinations. The notebook plots PCA of the pseudobulks. The `filter_by_expr` and `filter_by_prop` calls are plotting functions only, so no genes are actually filtered. pydeseq2 is run with `design_factors=["treatment"]` (no cluster/region term), contrast PBS vs WT, on all 18,751 genes. The Wald `stat` vector is scored with ULM against CollecTRI (mouse TFs), PROGENy (mouse pathways) and MSigDB Hallmark, keeping results with padj < 0.05.

**Results.** DESeq2 warned "Iterative size factor fitting did not converge". The number of significant genes was not printed. Low-count Olfr genes appear as highly significant (for example Olfr881, baseMean 0.005, padj 1.2e-12). There are 87 significant TFs, including Cebpb (+3.82), Cebpd (+4.36), Stat5a (+4.41), Cebpa (+4.51), Cebpg (+5.08), Vdr (+3.78) and Trp53 (−3.78). Network plots cover Stat1, Elk1, Stat5b, Myc, Id2, Stat3 and Hnf1a. PROGENy pathways: NFkB +17.35, Hypoxia +12.20, p53 +12.17, TNFa +12.06, EGFR +9.96, Trail +6.98, JAK-STAT +5.10, VEGF −5.66, PI3K −3.97, TGFb −3.07. Hallmark: 30 sets significant, including Oxidative phosphorylation −10.45, DNA repair −8.03, Allograft rejection +7.22, Protein secretion −7.07, Unfolded protein response −6.17, Adipogenesis −5.95 and Coagulation +5.59.

**Caveats.** The pseudobulks are cluster × section profiles pooled into one treatment term. Region/cluster is not modelled, and the 13 WT profiles all come from **one** section, so this is still pseudoreplication. No expression filter was applied, and size-factor fitting did not converge. Slide-confounded. The contrast is labelled "disease.vs.normal", copied from a tutorial.

### 15. BRICHOS vs WT gives almost the same pathway profile as PBS vs WT, so pathway activity shows no reversal
*Notebook:* `notebooks/00/11_decoupler_WT_BRICHOS.ipynb`

**Rationale.** This mirrors §14 for BRICHOS-treated mice. If BRICHOS reversed the disease programme, the inflammatory pathway activation should be weaker here than in PBS vs WT.

**Technical approach.** Same pipeline as §14, restricted to WT + BRICHOS sections: 48 section × cluster pseudobulks from 4 sections. pydeseq2 is run with `~treatment` and contrast BRICHOS vs WT. ULM scoring against CollecTRI, PROGENy and Hallmark, padj < 0.05.

**Results.** Size-factor fitting again did not converge. There are 103 significant TFs, for example Cebpa +4.72, Cebpd +4.80, Cdx2 +5.19 and Trp53 −3.83. PROGENy: NFkB +18.53, Hypoxia +13.13, p53 +12.94, TNFa +12.82, EGFR +11.02, Trail +7.49, JAK-STAT +4.12, VEGF −6.43, PI3K −3.67, TGFb −2.99. Hallmark: 33 sets significant, for example DNA repair −8.21, Protein secretion −6.85, Adipogenesis −6.78, UPR −6.38, Allograft rejection +6.22, Coagulation +6.00.

**Reading.** Every PROGENy pathway has the same sign and similar or larger magnitude than in PBS vs WT (for example NFkB 18.5 vs 17.3). At this resolution, BRICHOS-treated tissue looks as "diseased" as PBS tissue relative to WT.

**Caveats.** Same design flaws as §14: pseudoreplicated cluster profiles, one WT section, no filtering, non-converged size factors. The comparison between §14 and §15 is descriptive, with no interaction test. WT and 2 of 3 BRICHOS sections share slide P24215, so slide confounding here differs from §14.

### 16. The direct pseudobulk BRICHOS-vs-PBS decoupler analysis was never completed
*Notebook:* `notebooks/00/11_decoupler_PBS_BRICHOS.ipynb`

**Rationale.** This was meant to complete the trio in §14–15 with the treatment contrast itself.

**Technical approach.** Same pipeline: 130 section × cluster pseudobulks, then the 117 profiles that are not WT. pydeseq2 `~treatment` with contrast BRICHOS vs PBS, followed by CollecTRI, PROGENy and Hallmark ULM.

**Results.** Not executed to completion. The DESeq2 cell shows partial log output (size factors and dispersions) and has no execution count. Every later cell (results, volcano, TF and pathway scoring) is unexecuted, with no saved outputs.

**Caveats.** No results exist. Had it run, it would carry the same pseudoreplication issues as §14 and would be fully confounded with slide.

### 17. Region-wise pseudobulk DESeq2 finds zero BRICHOS-vs-PBS DE genes in every region, while the disease signal concentrates in ventricular area and thalamus
*Notebook:* `notebooks/00/12_DGE_regions_across_tissue.ipynb`

**Rationale.** This is the core regional question done with the section as the unit (pseudobulk). Does BRICHOS change expression relative to PBS in any brain region, and where is the disease signal? First committed 2025-10-14 as `13_DGE_regions_across_tissue`.

**Technical approach.** Raw counts from `ST_BRICHOS_region_subcluster.h5ad` are summed per section within each `re_annotation_regions` region, plus a tissue-wide `ALL_TISSUE` pseudobulk. That gives 10 section pseudobulks per region (WT 1, PBS 6, BRICHOS 3), and all 3 conditions are present in every region. For each region, pydeseq2 fits `~condition` and reports Wald contrasts PBS vs WT, BRICHOS vs WT and BRICHOS vs PBS. Per-condition mean CPM is joined in. CSVs were written to `notebooks/00/visium_region_de/`, which is not present in the repo or on disk now.

**Results.** Genes at FDR < 0.05 (PBS vs WT / BRICHOS vs WT / BRICHOS vs PBS):
- ALL_TISSUE: 19 / 6 / 0
- Ventricular area: 127 / 50 / 0
- Thalamus: 24 / 29 / 0
- Corpus callosum: 3 / 0 / 0
- Hypothalamus: 1 / 2 / 0
- Supragranular layers: 1 / 0 / 0
- Infragranular layers: 1 / 0 / 0
- Olfactory areas, Caudoputamen, Meninges, Hippocampal formation: 0 / 0 / 0

**Reading.** At section level, no region shows a single BRICHOS-vs-PBS DE gene. The disease contrasts are dominated by the ventricular area and thalamus (see §18 for what drives them). The partial "rescue" seen in the spot-level analyses (§7, §9) does not reach significance with sections as replicates.

**Caveats.** WT n = 1 section, so every "vs WT" contrast rests on one sample's dispersion. BRICHOS n = 3 and PBS n = 6 sections, and the README's mouse counts (3/3) imply sections are not independent mice. Fully confounded with slide (memory: batch-confound). Regions are expression clusters (§8), so region composition can differ by section.

### 18. Per-region scatter of PBS-vs-WT against BRICHOS-vs-WT log fold changes shows BRICHOS tracking the disease axis, with the thalamic DAM genes fully retained
*Notebook:* `notebooks/00/12_DGE_regions_across_tissue.ipynb` · *Figure:* `results/figures/ALL_TISSUE_PBSvsWT_vs_BRICHOSvsWT.png`, `results/figures/Thalamus_PBSvsWT_vs_BRICHOSvsWT.png`, `results/figures/Ventricular_area_PBSvsWT_vs_BRICHOSvsWT.png`, `results/figures/Meninges_PBSvsWT_vs_BRICHOSvsWT.png`, `results/figures/Corpus_callosum_PBSvsWT_vs_BRICHOSvsWT.png`, `results/figures/Hippocampal_formation_PBSvsWT_vs_BRICHOSvsWT.png`, `results/figures/Hypothalamus_PBSvsWT_vs_BRICHOSvsWT.png`, `results/figures/Caudoputamen_PBSvsWT_vs_BRICHOSvsWT.png`, `results/figures/Olfactory_areas_PBSvsWT_vs_BRICHOSvsWT.png`, `results/figures/Supragranualar_layers_PBSvsWT_vs_BRICHOSvsWT.png`, `results/figures/Infragranular_layers_PBSvsWT_vs_BRICHOSvsWT.png`

**Rationale.** This visualises rescue directly. Genes that BRICHOS normalises should fall below the diagonal: up in PBS vs WT, closer to 0 in BRICHOS vs WT.

**Technical approach.** Uses the §17 DESeq2 results per region. Genes are filtered to mean CPM ≥ 3. The plot shows log2FC(PBS vs WT) on x and log2FC(BRICHOS vs WT) on y, coloured by ΔLFC (BRICHOS − PBS), sized by mean CPM, with points outlined at padj < 0.05. Labels go to the top 20 genes per contrast by |LFC| among genes with padj < 0.5 (top 30 for ALL_TISSUE). The figures are saved as png/svg/pdf to `results/figures/` (file mtimes 2025-12-05). Top-gene tables are printed. The notebook ends with a dotplot of curated Thalamus and Ventricular-area gene panels (complement, interferon, stress and mitokine genes) by region, and a Cd74 spatial map.

**Results.** Highlights from the printed top-gene tables:
- **Thalamus.** DAM/microglial genes are up similarly in both contrasts, for example Cst7 log2FC 7.88 (PBS) vs 8.19 (BRICHOS), Tyrobp 1.93 vs 1.90, Cd52 2.35 vs 2.45, C3ar1 2.44 vs 2.13. Cbln1 is down in both (−1.72 vs −2.23).
- **Ventricular area.** The top "PBS vs WT" genes are neuronal (Chgb, Npy, Scn1b, Kcnq2), and BRICHOS vs WT likewise (Slc17a7, Nrgn, Cck, Gabra1). This points to differences in neuron content of the cluster, not pathology.
- **ALL_TISSUE PBS vs WT.** The top genes are mostly choroid-plexus/ependymal genes *down* in PBS (Kcnj13 −1.95, F5 −2.11, Col8a2 −1.78, Kl −1.59, Aqp1 −1.82, Igfbp2 −1.05). The top BRICHOS vs WT genes include Prelp −0.88, Emb −0.69 and Clu −0.30.
- **Corpus callosum.** Vim, A2m, Ccl12 and S1pr3 are down in both contrasts.
- **Hypothalamus.** Microglial genes are lower in BRICHOS vs WT: P2ry12 −0.75 and Ctss −0.99 (both padj 0.139), and Cx3cr1 −1.00. P2ry12 is also lower in PBS vs WT (−0.68).

**Reading.** In the regions with signal, BRICHOS largely reproduces the PBS-vs-WT changes (points on the diagonal). Much of the "disease" signal in ventricular area and tissue-wide looks like differences in tissue composition (choroid plexus, neurons) between the single WT section and the others, not plaque biology.

**Caveats.** The label threshold of padj < 0.5 is loose, although the code comment says "only label significant genes". These are descriptive scatters with no formal test of slope or attenuation; the formal attenuation tests were done later in `notebooks/03/`, per the memory notes. WT n = 1. Slide-confounded.

### 19. Single-cell label transfer onto spots (scVI/scANVI from Linnarsson and ABC references) failed, calling 99 % of spots "Neurons"
*Notebook:* `notebooks/00/_scratch/Untitled.ipynb`

**Rationale.** This was an early attempt to get cell-type composition per spot for interpreting regional effects (README TODO 2: "Run stereoscope to integrate single cell and spatial data"). It was later superseded by the cell2location/stereoscope deconvolution in notebooks/03 (memory: deconvolution thread). Committed 2025-11-05.

**Technical approach.** The notebook builds a class-balanced reference from Allen ABC whole-mouse-brain atlas parts (saved as `ABC_WBA_subset.h5ad`) and/or the Linnarsson adolescent mouse atlas. The balanced reference has Neurons 12,000, Oligos 12,000, Astrocytes 8,000, Vascular 6,000, Immune 4,000, PeripheralGlia 4,000 and Ependymal 1,826 cells. It trains scVI and then scANVI on MPS and predicts a hard `pred_celltype` per spot of `ST_BRICHOS_region_subcluster.h5ad`. A further cell (multinomial logistic regression on 3,002 HVGs with temperature calibration) was started but has no prediction output. The sibling `notebooks/00/_scratch/Untitled1.ipynb` tried a cell2location RegressionModel on the Linnarsson reference (160,796 cells, 7 classes, 12,692 genes). Training crashed with NaN parameters.

**Results.** scANVI spot labels: Neurons 25,333; Vascular 120; Ependymal 108; Astrocytes 74; Oligos 24; Immune 1 (of 25,660 spots). The cell2location attempt raised `ValueError` with NaN loc and produced no output.

**Caveats.** Hard single-label classification of multi-cell Visium spots is the wrong model. The near-all-Neuron output shows it failed. Only a few training epochs were run (max_epochs = 3 in the scVI smoke test). Exploratory scratch code with no saved results.

### 20. APP-NLGF (PBS) vs WT whole-section pseudobulk gives 803 DEGs led by microglial/DAM and plaque-induced genes (Cst7, Itgax, Ccl3, Ccl6, Tyrobp, Gfap)
*Notebook:* `notebooks/01/00-NLGF-vs-WT-figure1.ipynb` · *Figure:* `results/html/nb_figures/nb00_volcano_pbs_vs_wt.png`

**Rationale.** First step of the project: define the baseline Alzheimer's-model transcriptional signature (project objective 1, "WT vs APP-NLGF non-treated") before asking whether BRICHOS reverses it. The notebook was meant to be the source of manuscript Figure 1. Its title cell, "mtDSB x Age DE analysis (pydeseq2, LRT-free)", was copied from another project.

**Technical approach.** Input: `ST_BRICHOS_region_subcluster.h5ad`, raw counts layer, subset to PBS and WT spots. Counts were summed per `sample_id`, giving one pseudobulk per Visium section: 6 PBS sections (P28052_102/103/201/202/203/301) and 1 WT section (P24215_101). pyDESeq2 was fitted with `design_factors=["condition"]` and Cook's refitting, and the Wald contrast PBS vs WT was computed. Significance was padj < 0.05. A custom volcano plot flags |LFC| ≥ 1 and padj < 0.05.

**Results.** 18,751 genes were tested and 803 were significant at padj < 0.05. The top up-regulated genes by LFC were Cst7 (LFC 6.85, padj 1.2e-85), Itgax (5.18), Ccl3 (5.14), Fgr (3.46), Ccl6 (3.42), Cxcl5 (2.77), Cd74 (2.55), Ccl4, Pdcd1, H2-Aa (2.26). Tyrobp (1.84, p 9e-76), Gfap (1.81) and Cd52 (2.07) were also strongly significant. The most down-regulated significant genes included Calca (−3.27), Crb1 (−2.02), Ngfr (−1.93) and Anxa8 (−1.86). The code saved the tables to `results/DGE/PBS-wt-sig.csv` and `PBS-wt-all.csv`.

**Caveats.** There is only one WT section, so dispersion comes almost entirely from the PBS group and the WT-anchored p-values are fragile. All PBS sections are on slide P28052 and the WT is on P24215, so disease and batch are confounded. The replicate unit is the section, not the mouse (README: 3 PBS mice, but 6 PBS sections). The volcano was only shown inline. The PNG is a later export, and `results/DGE/` does not exist on disk now. The paths are relative to the old `notebooks/` location, from before the 2026-01-27 refactor.

**Next.** §21 repeats the same contrast for BRICHOS, and §23 refits all three groups jointly.

### 21. BRICHOS-treated APP-NLGF vs WT still shows the same DAM/PIG up-regulation, but with fewer DEGs (46)
*Notebook:* `notebooks/01/01-NLGF-BRI-vs-WT-figure2.ipynb` · *Figure:* `results/html/nb_figures/nb01_volcano_bri_vs_wt.png`

**Rationale.** This is the treatment arm of the first figure pair: does the disease signature from §20 persist when the mice are treated with BRICHOS?

**Technical approach.** The pipeline is the same as §20: section-level pseudobulk of raw counts, then pyDESeq2 `~condition` with a Wald contrast of BRICHOS vs WT. There are 3 BRICHOS sections (P24215_201, P24215_301, P28052_401) and 1 WT (P24215_101). Significance was padj < 0.05.

**Results.** 46 genes were significant at padj < 0.05. The top up-regulated genes were again microglial/plaque genes: Cst7 (LFC 6.72, p 6.2e-70), Itgax (5.39), Ccl3 (4.88), Ccl6 (3.30), Cxcl5 (2.85), Cd52 (1.98), Tyrobp (1.88), Gfap (1.67), C3ar1, Cd68, Trem2 (1.54), the C1q genes (≈1.0–1.1) and Ctss (1.08). The tables were written to `results/DGE/BRICHOS-wt-sig.csv` and `BRICHOS-wt-all.csv`.

**Reading.** BRICHOS does not abolish the plaque/DAM response. The smaller DEG count (46 vs 803) mostly reflects n = 3 vs n = 6 replicates and is not evidence of rescue.

**Caveats.** n = 1 WT. Two of the three BRICHOS sections share slide P24215 with the WT, so this contrast is less batch-confounded than §20, which makes the two DEG counts not directly comparable. The tables are no longer on disk.

### 22. A gene-wise Δlog2FC comparison of separately fitted PBS-vs-WT and BRICHOS-vs-WT tables puts almost every shared DEG on the diagonal
*Notebook:* `notebooks/01/02-comparing-bri-pbs.ipynb` · *Figure:* `results/html/nb_figures/nb02_volcano_bri_vs_pbs.png`

**Rationale.** This was a quick way to find genes that BRICHOS modulates differently from PBS. It compares the two WT-anchored significant-gene tables from §20 and §21 directly, without a BRICHOS-vs-PBS model.

**Technical approach.** The notebook reads `results/DGE/PBS-wt-sig.csv` (803 genes) and `BRICHOS-wt-sig.csv` (46 genes). It scatters LFC(PBS vs WT) against LFC(BRICHOS vs WT) and colours genes with |Δ| = |LFC_BRI − LFC_PBS| ≥ 0.2. The table that builds `df` (`log2FC_pbs`, `log2FC_brichos`, `delta`) is not created anywhere in the saved notebook, so the plot depended on hidden kernel state.

**Results.** No numeric output was printed. In the saved figure all plotted genes lie close to the identity line. Only four genes pass |Δ| ≥ 0.2: Itgax and Thbs4 are higher under BRICHOS, and Ccl3 and C4b are lower. These were read off the figure, not from printed output.

**Caveats.** The notebook cannot be reproduced as saved because `df` is undefined. It compares two independently fitted models that share the single WT. The exported PNG is labelled "Volcano: BRICHOS vs PBS" in `results/html/nb_figures/manifest.json`, but it is actually this Δ scatter.

### 23. A joint three-group pseudobulk model finds only 27 BRICHOS-vs-PBS DEGs, and the core disease markers are not reversed
*Notebook:* `notebooks/01/03-NLGF-vs-WT-vs-BRI.ipynb` · *Figure:* `results/html/nb_figures/nb03_volcano_bri_vs_pbs.png`, `results/html/nb_figures/nb03_lfc_correlation.png`, `results/html/nb_figures/nb03_heatmap.png`, `results/html/nb_figures/nb03_pig_dotplot.png`

**Rationale.** §20–22 fitted the groups separately. This run fits WT, PBS and BRICHOS in one DESeq2 model so the key contrast, BRICHOS vs PBS ("does it do anything?"), is estimated directly with shared dispersions. It also asks whether plaque-induced genes (PIGs) are attenuated.

**Technical approach.** All 10 sections were pseudobulked per `sample_id` and fitted with pyDESeq2 `~condition` (WT reference). Three Wald contrasts were computed: PBS vs WT, BRICHOS vs WT and BRICHOS vs PBS. A global LFC correlation was computed between PBS-vs-WT and BRICHOS-vs-WT. Marker genes were checked, and a "rescue scatter" plots LFC(BRI vs PBS) against LFC(PBS vs WT). For the 53 of 58 PIGs found in the data, the PBS-vs-WT LFC was plotted against the BRICHOS-vs-WT LFC, classing genes as attenuated (Δ < −0.1) or enhanced (Δ > 0.1). Dot plots by treatment were also drawn.

**Results.**
- DEG counts at padj < 0.05: PBS vs WT 382, BRICHOS vs WT 170, BRICHOS vs PBS 27.
- Global LFC correlation between PBS-vs-WT and BRICHOS-vs-WT: r = 0.81.
- Disease markers (LFC PBS vs WT → LFC BRI vs PBS): Gfap 1.81 → −0.13, Cst7 6.85 → −0.13, C1qa 1.03 → −0.04, Trem2 1.63 → −0.09, Cd68 1.57 → +0.03.
- Top BRICHOS-vs-PBS genes by padj: Insig2 (+0.21, padj 6.7e-4), Hdgfl2 (−0.23, padj 3.7e-3), Gh (−1.53, padj 0.014), Ift122, Ppp2r3c, Pomc (−2.02), Gbp5 (−0.80), Fosl2 (+0.65), Ddx3x, Setdb2, H2-D1 (−0.32), Cd93 (+0.39). All have padj ≤ 0.04.
- Apoe was not changed (PBS vs WT LFC 0.13, padj 0.50).
- PIG LFCs for PBS vs WT included Cst7 6.85, Msr1 5.61, Itgax 5.18, Ccl6 3.42 and Lyz2 2.71.

**Reading.** At whole-section level BRICHOS barely moves the disease signature. The few BRICHOS-vs-PBS hits are small-effect or pituitary/hypothalamic hormone genes (Gh, Pomc), not the PIG/DAM programme.

**Caveats.** The counts of attenuated and enhanced PIGs were not printed. The BRICHOS-vs-PBS contrast is confounded with slide: all PBS sections are on P28052. Gh and Pomc probably reflect how much pituitary or hypothalamic tissue each section captured, not treatment. The exported panel names in `results/html/nb_figures/manifest.json` do not always match the content: the PIG scatter is called "Three-way heatmap".

### 24. A section-level regional plaque-gene score is lower in BRICHOS than PBS in most regions, but no region passes FDR
*Notebook:* `notebooks/01/04-regional-analysis.ipynb` · *Figure:* `results/html/nb_figures/nb04_spatial_pig.png`, `results/html/nb_figures/nb04_pig_by_region.png`, `results/html/nb_figures/nb04_spatial_regions.png`, `results/html/nb_figures/nb04_region_composition.png`

**Rationale.** The project's central question is whether BRICHOS acts in specific brain regions. This run asks whether a plaque-induced-gene score differs between BRICHOS and PBS within each annotated region.

**Technical approach.**
- **Regions.** Leiden-0.3 clusters were relabelled into 10 compartments (`re_annotation_regions`): Caudoputamen, Supragranular layers, Infragranular layers, Hypothalamus, Corpus callosum, Olfactory areas, Meninges, Hippocampal formation, Thalamus, Ventricular area. These are expression clusters mapped to anatomy, not atlas registration.
- **Score.** The plaque score is the mean of the stored `adata.X` over 53 PIGs, computed per spot (25,660 spots).
- **Tests.**
  - (i) Spot-level mixed model `plaque_score ~ treatment * region` with a random intercept per section.
  - (ii) Spot-level Mann–Whitney test per region.
  - (iii) Spot-level two-way OLS ANOVA.
  - (iv) Section-level means per region, compared by one-way ANOVA across WT/PBS/BRICHOS.
  - (v) Section-level Mann–Whitney test of PBS vs BRICHOS per region with BH FDR. This is the appropriate test.

**Results.**
- (v) Section-level, PBS n = 6 vs BRICHOS n = 3:
  - Infragranular layers: Δ = −0.029, p = 0.024, FDR 0.107.
  - Supragranular layers: Δ = −0.039, p = 0.024, FDR 0.107.
  - Caudoputamen: Δ = −0.028, p = 0.095, FDR 0.29.
  - Meninges: Δ = −0.043, p = 0.17.
  - Hypothalamus, Ventricular area (Δ +0.020), Corpus callosum (Δ +0.0001), Hippocampal formation (PBS 5 vs BRICHOS 2) and Olfactory areas (Δ −0.035): p ≥ 0.55.
  - Thalamus dropped out of this test.
- (iv) Section-level ANOVA including WT: infragranular FDR 2.4e-5, supragranular 3.5e-4, caudoputamen 3.9e-3, thalamus 7.3e-3, meninges 0.051.
- (ii) Spot-level Mann–Whitney: every region had p < 0.05 (e.g. supragranular p = 8e-45). This is pseudoreplicated.
- (i) Mixed model: WT intercept shift −0.096 (p = 0.024). PBS vs BRICHOS main effect +0.033 (p = 0.20). Significant PBS × region interactions for corpus callosum, hippocampus, hypothalamus and ventricular area. The fit reported convergence warnings and a non-positive-definite Hessian.

**Reading.** The direction is consistent: BRICHOS is lower in 7 of 9 regions, most clearly in cortex and caudoputamen. With n = 6 vs 3 sections none of this survives FDR. These are the same regions that the later rescue tests in notebooks/03 single out.

**Caveats.**
- An earlier per-region Mann–Whitney cell (cell 21) printed n = 10 vs 10 and NaN p-values because of a groupby bug. The bar chart in cell 23, exported as panel C of manuscript Figure 3, therefore shows "FDR=nan".
- Treatment is confounded with slide.
- The replicate unit is the section, not the mouse.
- The score is the mean of normalised expression, so it is sensitive to spot depth.

### 25. Region-wise pseudobulk DESeq2 finds 0–43 BRICHOS-vs-PBS DEGs per region, and almost none are plaque-induced genes
*Notebook:* `notebooks/01/04-regional-analysis.ipynb` · *Results:* `results/tables/DGE_BRICHOS_vs_PBS_Caudoputamen.csv` (one CSV per region) · *Figure:* `results/figures/regional_BRICHOS_vs_PBS_summary.png`, `results/html/nb_figures/nb04_de_summary_barplot.png`, `results/html/nb_figures/nb04_heatmap_regions.png`

**Rationale.** This moves from a single score (§24) to genome-wide regional differential expression, to see where BRICHOS changes expression and whether the changed genes are PIGs.

**Technical approach.**
- **Pseudobulk.** Raw counts were summed per section × region, giving up to 10 sections × 10 regions.
- **Per-region DESeq2.** For each region, pyDESeq2 `~treatment` was fitted on all three groups, so WT is included in the dispersion fit, with contrasts PBS vs WT, BRICHOS vs WT and BRICHOS vs PBS.
- **Summaries.** DEG counts at padj < 0.05, up and down. PIG overlap and mean PIG LFC. A heatmap of BRICHOS-vs-PBS LFCs of PIGs per region, restricted to baseMean > 50.
- **Marker ranking.** Spot-level scanpy t-test markers per treatment within each region.
- **Interaction model.** An attempted region × treatment interaction DESeq2.

**Results.**
- BRICHOS vs PBS DEGs at padj < 0.05:

  | Region | DEGs | Up | Down |
  |---|---:|---:|---:|
  | Infragranular layers | 43 | 40 | 3 |
  | Supragranular layers | 40 | 31 | 9 |
  | Hypothalamus | 38 | 24 | 14 |
  | Thalamus | 32 | 14 | 18 |
  | Caudoputamen | 21 | 19 | 2 |
  | Corpus callosum, Hippocampal formation, Meninges, Olfactory areas, Ventricular area | 0 | 0 | 0 |

- Significant PIGs: 1 each in caudoputamen, thalamus and supragranular layers, and 0 elsewhere.
- Mean PIG LFC (BRICHOS vs PBS) ranged from −0.32 (meninges), −0.26 (olfactory), −0.14 (hippocampus), −0.10 (infragranular), −0.08 (supragranular) to +0.10 (hypothalamus) and +0.20 (corpus callosum).
- Spot-level top markers: BRICHOS spots in most regions ranked Egr1, Faim2, Syt1, Pfn2 and Prkar1b highest. PBS spots ranked Pmch, Gfap, Cd74, Mt2, B2m and C4b highest.
- The interaction model was fitted with additive `design_factors=["treatment","region"]` only, so no interaction coefficients existed and the cell printed "No interaction results available".

**Reading.** Regional DEGs are mostly up-regulated under BRICHOS in cortex and striatum and are not PIGs. The PIG programme shifts only weakly downward, most in meninges and olfactory areas. Pmch, a hypothalamic gene, appearing as a PBS marker in many regions points to dissection or section-level differences, not regional biology.

**Caveats.** n = 6 PBS vs 3 BRICHOS sections per region, and all PBS sections are on slide P28052, so the DEG counts may be largely batch. Memory notes that regional DGE needs re-checking for the slide confound. The interaction analysis did not run. The single WT makes the WT contrasts unreliable.

**Next.** The notebooks/03 attenuation series (scaffolded in §33) replaced DEG counting with disease-anchored rescue statistics.

### 26. Among disease genes, the ones BRICHOS attenuates most are interferon/MHC-II genes (H2-Aa, Gbp5, Cd74, Zbp1), not plaque-induced genes
*Notebook:* `notebooks/01/05-brichos-attenuation.ipynb` · *Figure:* `results/html/nb_figures/nb05_pig_barplot.png`, `results/html/nb_figures/nb05_lfc_comparison.png`, `results/html/nb_figures/nb05_attenuation_scatter.png`, `results/html/nb_figures/nb05_volcano_bri_vs_pbs.png`

**Rationale.** Instead of asking which genes are DE between BRICHOS and PBS, this asks which disease genes (significant in PBS vs WT) move back toward WT under BRICHOS. This is an "attenuation" framing.

**Technical approach.** Section-level pseudobulk was fitted with the same three-group pyDESeq2 `~treatment` model as §23 (WT reference). For genes with padj(PBS vs WT) < 0.05, Δ = LFC(BRI vs WT) − LFC(PBS vs WT) was computed. The genes were ranked to give the top 20 attenuated (most negative Δ) and top 20 "increased" (most positive Δ, i.e. disease-down genes moving back up). Bars were coloured by membership of the 53-gene PIG list. A BRICHOS-vs-PBS volcano was drawn with baseMean > 150, |LFC| ≥ 0.5 and padj < 0.1.

**Results.**
- Top attenuated genes (Δ): H2-Aa −0.81, Gbp5 −0.80, Cd74 −0.72, Zbp1 −0.60, Galnt15 −0.60, Nlrc5 −0.45, Cybb −0.44, Sprtn −0.41, Stk38 −0.41, Tap1 −0.40, Gpr65 −0.39, Psmb9 −0.38. Further down the list were Irgm1, Gsn, Xdh and Trim14 (−0.33).
- Only H2-Aa among the top 20 is in the PIG list.
- In the direct BRICHOS-vs-PBS test only Gbp5 is significant (padj 0.020). Adgrl4 has padj 0.11, and the rest have padj 0.19–0.77.
- Disease-down genes restored most were Calca (Δ +0.82), Ngfr (+0.70), Anxa8 (+0.58), Slc39a14 (+0.57) and Grin2d (+0.54).

**Reading.** The genes BRICHOS attenuates point to an interferon / antigen-presentation axis (Gbp5, Zbp1, Nlrc5, Tap1, Psmb9, Irgm1, MHC-II) rather than the core DAM/PIG genes. This anticipates the later interferon-response-microglia (IRM) rescue finding recorded in `results/manuscript_findings.md`.

**Caveats.** The ranking is on point estimates; almost no gene is individually significant for BRICHOS vs PBS. The analysis inherits the n = 1 WT problem and the slide confound. Several cells errored in the saved notebook (`top_atten` used before it was defined; `marker_genes` undefined; `sc.pl.he` does not exist). One gene in a hard-coded list is misspelt as "Hp gds". The executed and non-executed copies have the same cells.

### 27. Tangram cell-type-stratified regional analysis was written but never run, so manuscript Figure 5 is placeholders
*Notebook:* `notebooks/01/06-celltype-regional-analysis.ipynb` · *Code:* `scripts/assemble_figures.py` · *Figure:* `results/figures/manuscript/figure5_celltype.png`

**Rationale.** The aim was to test whether BRICHOS effects are carried by particular cell types, e.g. microglia/immune. The plan was to combine Tangram deconvolution (`notebooks/00/tangram/12_tangram.ipynb`) with the regions from §24–25.

**Technical approach (code only).** Load `data/ST_BRICHOS_with_tangram.h5ad`. Then, using the top Tangram label per spot:
1. Compute the cell-type composition per section × region and test PBS vs BRICHOS with Mann–Whitney and BH, writing `results/tables/celltype_composition_tests.csv`.
2. Compare a PIG score (`sc.tl.score_genes` on a different 43-gene PIG list) by cell type × treatment × region, at section level.
3. Run immune-enriched pseudobulk DESeq2 per region on the top quartile of summed immune/microglia Tangram probability, writing `results/tables/DGE_immune_BRICHOS_vs_PBS_*.csv`.

**Results.** None. The notebook has no executed outputs. The output tables and PDFs it would write are absent from `results/`, and `data/ST_BRICHOS_with_tangram.h5ad` is not present in the repo. In `figure5_celltype.png` all six panels (A–F: cell-type proportions, cell-type × region, microglia rescue, endocytosis module, Tangram overlay, cell-type DGE) are grey "missing nb06_*.png" placeholders.

**Caveats.** The code imports a non-existent module (`notebooks_04_funcs`) and would fail at step 3. Its PIG list differs from the 53-gene list used elsewhere. The cell-type question was later answered with cell2location/stereoscope deconvolution in notebooks/03 (09–09e), not here.

### 28. The PIG score shows no gradient with distance to the tissue border/meninges and no treatment × distance interaction
*Notebook:* `notebooks/01/07-PIG_spatial_gradients_meninges.ipynb` · *Figure:* `results/figures/manuscript/figure4_spatial_gradient.png`, `results/html/nb_figures/nb07_pig_gradient_pbs.png`, `results/html/nb_figures/nb07_pig_gradient_bri.png`, `results/html/nb_figures/nb07_tissue_rim.png`, `results/html/nb_figures/nb07_spatial_dist_meninges.png`

**Rationale.** Meninges had the highest regional plaque score (§24: intercept shift +0.141). The idea was that a peripherally delivered agent such as BRICHOS might act first near the brain surface, giving a surface-to-depth gradient.

**Technical approach.**
- **Data.** Spot-level, using the precomputed `PIG_score` stored in the h5ad.
- **Distances.** Per section, (a) the distance to the nearest Meninges-labelled spot (KD-tree on `array_row`/`array_col`; range 0–50.6 array units), and (b) the distance to a tissue-border rim. The rim is spots within 1.2 array units of the per-section convex hull of spot positions.
- **Profiles.** Binned PIG-vs-distance profiles (bin 2, section means then mean ± SD across sections, per treatment). PIG by region per treatment.
- **Model.** Linear mixed model `PIG_score ~ treatment * log1p(dist_to_border)` with a random intercept per section, PBS reference, on 22,703 spots in 9 sections (6 PBS + 3 BRICHOS).
- **Zone means.** Means in a border zone (distance < 5) vs the rest.

**Results.**
- Mixed model:

  | Term | Coefficient | p |
  |---|---:|---:|
  | Intercept | 0.260 | — |
  | BRICHOS | −0.018 | 0.229 |
  | log distance | −0.003 | 0.101 |
  | BRICHOS × log distance | +0.004 | 0.201 |

  Group variance was ≈ 0.
- Zone means:

  | | Interior | Border zone |
  |---|---:|---:|
  | PBS | 0.255 | 0.251 |
  | BRICHOS | 0.248 | 0.230 |

**Reading.** There is no detectable surface gradient and no treatment-dependent gradient. BRICHOS is slightly lower near the border in raw means, but this is not supported statistically.

**Caveats.** Distances are in Visium array-grid units, but the axes and comments label them µm. The convex hull is a poor stand-in for the real tissue edge in coronal sections with ventricles. The "border" rim and the meninges cluster overlap only partially. The analysis is spot-level with only 9 sections. The data path is hard-coded to `/Volumes/processing2/...`, which is not mounted now.

### 29. Canonical section orientation tooling exists, but every rotation is still 0° and no oriented output is on disk
*Notebook:* `notebooks/01/08-canonical-orientation.ipynb` · *Code:* `scripts/apply_canonical_orientation.py`, `scripts/preview_orientations.py`, `scripts/_build_orientation_notebook.py`, `utils/canonical_orientation.py` · *Results:* `data/rotations.yaml`

**Rationale.** Infrastructure, so that all coronal sections face the same way in figure panels. The oriented h5ad is the declared input (`H5AD_ORIENTED`) of the notebooks/03 attenuation series.

**Technical approach.**
- `utils/canonical_orientation.py` rotates each library's spot coordinates (`obsm['spatial']`, full-resolution pixels) and its hires/lowres images by a per-sample counter-clockwise angle, with an optional horizontal flip, using `scipy.ndimage.rotate`.
- `data/rotations.yaml` holds the per-sample table.
- `preview_orientations.py` renders a grid of each section with treatment-coloured spots and a centre crosshair. It writes `results/figures/manuscript/orientation_preview.{png,svg}`.
- `apply_canonical_orientation.py` writes `ST_BRICHOS_region_subcluster_oriented.h5ad`.
- Notebook 08 is the interactive version of the same loop. It is generated by `_build_orientation_notebook.py`.

**Results.** All entries in `rotations.yaml` are `angle_deg: 0, flip_h: false`, so no reorientation has been set. Notebook 08 has no executed outputs, and `orientation_preview.png` is not in `results/figures/manuscript/`.

**Caveats.** `rotations.yaml` lists 11 libraries, including "P24215_401" under PBS. That section does not appear in the analysed h5ad, which has 10 sections and no PBS on P24215 (see the batch-confound memory note). It conflicts with the project-wide section list. The default `sample_key="library_id"` is not an obs column in the h5ad used elsewhere (`sample_id`). Paths point to an unmounted volume.

### 30. Six Cell-style manuscript composites (Figures 0–5) are tiled from exported notebook panels, and Figure 5 is entirely missing panels
*Code:* `scripts/assemble_figures.py` · *Figure:* `results/figures/manuscript/figure0_schematic.png`, `results/figures/manuscript/figure1_disease_baseline.png`, `results/figures/manuscript/figure2_global_rescue.png`, `results/figures/manuscript/figure3_regional_heterogeneity.png`, `results/figures/manuscript/figure4_spatial_gradient.png`, `results/figures/manuscript/figure5_celltype.png`

**Rationale.** Infrastructure: a first manuscript figure set (status 2026-05-05) built from the notebooks/01 outputs.

**Technical approach.** Cell journal style: Helvetica 7–9 pt, 174 mm double column, uppercase panel labels, SVG plus 300-dpi PNG. Figure 0 is a drawn schematic: WT n = 1, APP-NLGF + PBS n = 3, APP-NLGF + BRICHOS n = 3, with mouse, syringe, brain and capture-chip icons. Figures 1–5 tile PNGs from `results/html/nb_figures/`, which were exported from the notebook HTMLs in `results/html/` (the exporting script is not in the repo), plus `results/figures/regional_BRICHOS_vs_PBS_summary.png`.

| Figure | Content | Source |
|---|---|---|
| 1 | Disease baseline | §20, §24 |
| 2 | Global rescue | §21–23, §26 |
| 3 | Regional heterogeneity | §24–26, §23 |
| 4 | Spatial gradient | §28 |
| 5 | Cell type | §27 |

Missing inputs are drawn as grey "missing" boxes.

**Results.** All six composites exist. Figure 5 consists only of placeholders. Several panel titles in Figure 3 do not match their content: "Significant gene ranks" is the plaque-score Δ bar chart with FDR=nan, "PIG response per region" is the top-attenuated-genes bar chart, and "Gene-set heatmap" is the PIG LFC scatter.

**Caveats.** These figures reuse inline notebook renders, so they inherit the caveats of §20–28: n = 1 WT, the slide confound, and the FDR=nan panel. They were superseded by the story figures (§31).

### 31. Three story figures present convergent but partial BRICHOS rescue in cortex and caudoputamen (sign concordance 0.74–0.82; infragranular slope permutation p = 0.012)
*Code:* `scripts/assemble_story_figures.py` · *Results:* `results/tables/attenuation/` · *Figure:* `results/figures/manuscript/figure_story1_design_disease.png`, `results/figures/manuscript/figure_story2_regional_rescue.png`, `results/figures/manuscript/figure_story3_mechanistic_axis.png`

**Rationale.** These replace Figures 0–5 (§30) with a three-part narrative: design and disease, regional rescue, and mechanistic axis. They are drawn directly from the notebooks/03 attenuation tables (scaffolded in §33), not from notebook screenshots.

**Technical approach.**
- **Story 1.** Reuses the schematic, the top 15 up and down genes of `disease_signature.tsv` (whole-tissue PBS vs WT), and PBS-vs-WT LFCs of `PIG_CHEN2020` genes.
- **Story 2.**
  - (A) Sign concordance per region with binomial stars, from `regional_rescue_stats.tsv`.
  - (B) Rescue slope (LFC_BRI−PBS regressed on LFC_disease) with bootstrap CI and permutation p, from `rescue_permutation_null.tsv`; 84 permutations.
  - (C) Per-section PIG-rescued spot fraction, PBS vs BRICHOS, with Mann–Whitney p, from `responder_fraction_by_region.tsv`.
- **Story 3.** A hand-drawn cascade diagram, a heatmap of 21 signatures from `broader_signature_rescue.tsv`, and GSEA top terms for 6 regions (`gsea_*.tsv`).

**Results.** The values below are from the TSVs plotted, for 7 regions with n_BRI = 3 and n_PBS = 6.

| Region | Concordance | Binomial p | Slope | Perm p | Responder fraction, PBS → BRI | MWU p |
|---|---:|---:|---:|---:|---|---:|
| Infragranular | 0.82 | 4e-23 | −0.070 | 0.012 | 0.066 → 0.105 | 0.024 |
| Supragranular | 0.81 | 5e-22 | −0.011 | 0.26 | 0.059 → 0.170 | 0.048 |
| Olfactory | 0.79 | — | −0.091 | 0.12 | — | 0.27 |
| Meninges | 0.76 | — | −0.104 | 0.12 | BRI < PBS | 0.87 |
| Caudoputamen | 0.74 | — | −0.004 | 0.52 | 0.102 → 0.182 | 0.048 |
| Hypothalamus | 0.47 | 0.84 | — | — | — | — |
| Corpus callosum | 0.38 | — | +0.061 | — | — | — |

The corpus callosum slope CI is 0.004–0.179, i.e. in the wrong direction.

**Caveats.** The underlying statistics belong to the notebooks/03 runs. This entry only documents their presentation. All slopes are far from the full-normalisation value of −1, so rescue is partial. Every contrast is BRICHOS vs PBS and therefore confounded with slide. The later slide-removal work (memory: slide-removal thread) found that the cortex rescue survives RUVg but that section 401 disagrees. The permutation p-value label has a formatting bug: it prints "p≤.012" whenever p ≤ 0.001.

### 32. Manuscript Results and Discussion text was drafted claiming convergent cortical/striatal rescue and an interferon-driven microglia–astrocyte cascade
*Code:* `scripts/update_manuscript_docx.py`

**Rationale.** Infrastructure: write the attenuation findings into the manuscript draft `brichos-st.docx`.

**Technical approach.** python-docx opens `~/Downloads/brichos-st.docx`, strips trailing empty paragraphs and appends new text. It adds a Results subsection, "Three orthogonal tests converge on cortical and striatal BRICHOS-mediated rescue", a subsection "BRICHOS attenuates an interferon-driven microglial-astrocyte cascade and concurrently restores neuronal signaling activity", and a DISCUSSION section. The result is written to `~/Downloads/brichos-st_updated.docx`, leaving the original untouched. The prose mirrors `results/manuscript_findings.md` (status 2026-05-05).

**Results.** Text only. It states sign concordance 0.78–0.85 (binomial p < 1e-3), per-mouse Mann–Whitney p = 0.048 / 0.024 / 0.048 and infragranular slope permutation p = 0.012. Olfactory and meninges reach only p = 0.12. Corpus callosum and hypothalamus fail all three tests. Hippocampus, thalamus and ventricular area are excluded (≤ 2 PBS sections). Whether the script was run cannot be verified, because its outputs live outside the repo.

**Caveats.** The drafted text says "n = 3 vs n = 7 sections" and "84 unique reassignments … across the seven non-WT mice". The TSVs have n_PBS = 6, and 84 = C(9,3), which corresponds to 9 non-WT sections, not 7 mice. It describes the Mann–Whitney test as "per-mouse", although the unit is the section. It predates the discovery of the slide confound (2026-09-25), so its causal claims need hedging.

### 33. The disease-anchored attenuation statistics roadmap was scaffolded as seven generated notebooks (pseudobulk → rescue stats → mixed model → responder spots → permutation null → confounders → broader signatures)
*Code:* `scripts/_build_attenuation_notebooks.py`

**Rationale.** Infrastructure. Because there is only one WT, per-region WT-anchored DE is not valid (see §25). The plan was to anchor the disease signature once, globally, and test rescue per region on BRICHOS vs PBS. This script turned that roadmap into runnable notebooks.

**Technical approach.** nbformat code generation writes `notebooks/03/01-pseudobulk.ipynb`, `02-attenuation-stats.ipynb`, `03-signature-mixed-model.ipynb`, `04-responder-spots.ipynb`, `05-permutation-null.ipynb`, `06-confounder-checks.ipynb` and `07-broader-signatures.ipynb`. They share a header that reads the oriented h5ad (`ST_BRICHOS_region_subcluster_oriented.h5ad`, from §29) from the processing volume.

**Results.** The seven notebooks exist in `notebooks/03/`. Their outputs (`results/tables/attenuation/*.tsv`) feed §31 and §32. The scientific results belong to the notebooks/03 entries.

**Caveats.** The generated notebooks depend on the oriented h5ad, which is identical to the unoriented one because all rotations are 0 (§29). Regenerating would overwrite any manual edits made in `notebooks/03/01–07`.

### 34. Spot counts collapse into 80 section × region pseudobulk profiles, with only 7 regions present in at least 2 BRICHOS and 2 PBS sections
*Notebook:* `notebooks/03/01-pseudobulk.ipynb` · *Code:* `utils/attenuation.py`, `scripts/_build_attenuation_notebooks.py` · *Results:* `results/tables/attenuation/pseudobulk_counts.tsv`, `results/tables/attenuation/pseudobulk_meta.tsv`

**Rationale.** The attenuation pipeline (notebooks 03/01–07) needs a unit of replication above the spot. Treating thousands of spots as independent would inflate significance, so counts are summed per section and region. Everything downstream (disease signature, regional rescue tests, permutation null, broader signatures) reads this table.

**Technical approach.** Input is the region-annotated AnnData on the processing volume (`ST_BRICHOS_region_subcluster[_oriented].h5ad`), with raw integer counts in `layers['counts']`. `make_pseudobulk` sums the counts for each (sample_id, region) pair and keeps pairs with at least 10 spots (the earlier threshold was 20; the comment says it was relaxed to keep small regions). The sums are rounded to integers. Each sample_id is one Visium section, and downstream notebooks call it a "mouse". The notebook skeleton was generated by `scripts/_build_attenuation_notebooks.py`.

**Results.** The notebook has no saved outputs. The tracked tables (committed 2026-05-05) hold 80 section × region rows over about 18,750 genes. The groups are 3 BRICHOS sections (P24215_201, P24215_301, P28052_401), 6 PBS sections (all on P28052) and 1 WT section (P24215_101). Coverage by region (BRICHOS/PBS/WT):

| Region | BRICHOS | PBS | WT |
|---|---|---|---|
| Caudoputamen, Corpus callosum, Hypothalamus, Infragranular layers, Meninges, Olfactory areas, Supragranular layers | 3 | 6 | 1 |
| Hippocampal formation, Thalamus | 1 | 2 | 1 |
| Ventricular area | 1 | 1 | 0 |

**Caveats.**
- The code sets `REGION_KEY='anatomical_region'`, but the saved region names match the `re_annotation_regions` column that notebooks 04, 07 and 08 use. The committed notebook is probably not the exact version that produced the table.
- Each section is treated as an independent mouse. The repo README lists n=3 PBS and n=3 BRICHOS animals, while the data hold 6 PBS sections. Whether some sections come from the same animal is not documented.
- Treatment is confounded with Visium slide: all PBS sections are on P28052, and 2 of 3 BRICHOS sections are on P24215. This affects every downstream contrast.

### 35. The global PBS-vs-WT "disease signature" is led by microglial and plaque-induced genes (Cst7, Ccl6, Cd74, Ccl3) and gives 167 up and 32 down genes
*Notebook:* `notebooks/03/02-attenuation-stats.ipynb` (§1) · *Code:* `utils/attenuation.py`, `utils/gene_sets.py` · *Results:* `results/tables/attenuation/disease_signature.tsv` · *Figure:* `results/figures/manuscript/figure_story1_design_disease.png` (assembled by `scripts/assemble_story_figures.py`)

**Rationale.** With a single WT animal, WT-anchored tests within each region are not valid. The design therefore defines "what disease does" once, on whole tissue, and then asks region by region whether BRICHOS opposes it (§36).

**Technical approach.**
- **Data:** the pseudobulk from §34, summed over regions so each section is one sample.
- **Statistic:** `pseudobulk_lfc` computes log2(CPM+1) and a t-test of PBS (6 sections) against WT (1 section). With `min_n=1`, the PBS replicates supply a pooled variance. P-values are BH-adjusted, and the notebook itself states they are inflated and descriptive only.
- **Gene sets:** "disease-up" is the union of data-driven genes (LFC > 0.5 and padj < 0.1) and the Chen 2020 plaque-induced gene list (PIG, 52 genes, all present in the table). "Disease-down" is the data-driven LFC < −0.5 set.

**Results.** The notebook has no saved outputs; the numbers below come from the tracked table.
- 16,531 genes were tested.
- The data-driven set has 167 up genes and 32 down genes.
- Union disease-up = 190 genes (the 167 data-driven plus PIG; 29 PIG genes are also in the data-driven set).
- Top LFCs: Cst7 5.37 (padj 7.8e-5), Ccl6 3.06, Cd74 2.80, Ccl3 2.36, Cd52 1.93, Itgax 1.87 and Gfap 1.81. These are microglial, plaque-associated and astrocytic genes.
- Hypothalamic neuropeptides also rank high (Pmch 2.19, Hcrt 2.13), but with padj around 0.8.

**Caveats.**
- WT is a single section. The disease axis is therefore one animal against six sections, and padj values are not valid inference.
- WT (P24215) and PBS (P28052) sit on different slides, so the disease axis also contains the slide effect.
- **Corrected 2026-10-22:** the slide-removal work in notebooks 03/14–16 found that the disease axis (PBS on P28052 − WT on P24215) and the rescue axis (BRICHOS, mostly on P24215 − PBS) share the slide effect with opposite signs. The slide alone can therefore produce an apparent rescue slope (see §36).
- The high-ranked neuropeptides (Pmch, Hcrt) probably reflect the cutting plane and anatomy rather than disease. Notebook 03/12 later found such genes driving the treatment classifier.

### 36. BRICHOS moves PBS-vs-WT disease genes back toward WT in cortex, meninges and olfactory areas, but not in corpus callosum or hypothalamus
*Notebook:* `notebooks/03/02-attenuation-stats.ipynb` (§2–4) · *Code:* `utils/attenuation.py` · *Results:* `results/tables/attenuation/regional_rescue_stats.tsv` · *Figure:* `results/figures/manuscript/regional_rescue_slope.png`, `results/figures/manuscript/regional_sign_concordance.png`, `results/figures/manuscript/figure_story2_regional_rescue.png`

**Rationale.** This is the project's central question: does BRICHOS reverse Alzheimer's-related expression, and in which regions? It is tested with BRICHOS vs PBS per region, where there are replicates, measured against the global disease axis from §35.

**Technical approach.**
- **Unit:** section × region pseudobulk.
- **Per-region contrast:** `pseudobulk_lfc` runs a Welch t-test on log2(CPM+1), BRICHOS vs PBS. Regions need at least 2 sections per group, so 7 regions qualify (3 vs 6 in each).
- **Sign concordance (T2′):** the fraction of disease-signature genes (union-up plus data-driven down) that BRICHOS moves in the rescue direction, tested with a binomial test against 0.5 and a one-sided Wilcoxon test on the signed LFCs.
- **Rescue slope (T3′):** an ordinary least-squares slope of each region's BRICHOS-vs-PBS LFC on the global disease LFC, over the 190 disease-up genes. A negative slope means rescue and −1 means full rescue. The 95% CI comes from 1,000 bootstrap resamples of genes, not of mice.

**Results.** The notebook has no saved outputs; the numbers below come from the tracked table.

| Region | Concordant fraction | Binomial p | Slope | 95% CI |
|---|---|---|---|---|
| Meninges | 0.76 | 8.9e-15 | −0.104 | −0.24 to −0.05 |
| Olfactory areas | 0.79 | 3.1e-18 | −0.091 | −0.21 to −0.04 |
| Infragranular layers | 0.82 | 3.8e-23 | −0.070 | −0.15 to −0.04 |
| Supragranular layers | 0.81 | 4.7e-22 | −0.011 | −0.08 to 0.03 |
| Caudoputamen | 0.74 | 2.5e-13 | −0.004 | −0.07 to 0.04 |
| Hypothalamus | 0.47 | 0.84 | −0.027 | −0.10 to 0.01 |
| Corpus callosum | 0.38 | 1.00 | +0.061 | 0.004 to 0.18 |

Each region has 212–220 signature genes. In Corpus callosum the slope is positive, meaning BRICHOS moves disease genes further from WT.

**Reading.** The rescue is partial: all slopes are far from −1. Its regional pattern is cortex, meninges and olfactory areas rescued, and white matter and hypothalamus not rescued.

**Caveats.**
- The gene-level binomial and bootstrap statistics treat genes as independent, with only 3 vs 6 sections behind them. §39 gives the section-level permutation test.
- **Corrected 2026-10-22:** notebooks 03/14–16 found that the slide confound alone mimics a negative rescue slope. A pure slide contrast (P24215 BRICHOS vs P28052_401) gives a PIG rescue slope of −0.10. The within-slide contrast (401 vs P28052 PBS) gives +0.03. 401 ranks only 2/7 within P28052 on the PIG rescue. Section-centred profiles (03/14) keep the regional ranking, with no permutation p below 0.13. With RUVg slide removal (03/15), the cortical rescue survives: infragranular p = 0.012 under every method. That rescue rests on the P24215 BRICHOS sections, and 401 alone shows the opposite regional profile.

### 37. In the spot-level mixed model, PIG, DAM, ARM, homeostatic and endocytosis scores are all higher in both PBS and BRICHOS than in WT, and BRICHOS sits only slightly below PBS
*Notebook:* `notebooks/03/03-signature-mixed-model.ipynb` · *Code:* `utils/gene_sets.py` · *Results:* `results/tables/attenuation/signature_mixed_model.tsv`

**Rationale.** The aim was to replace thousands of gene-level tests with one test per signature per region, at spot resolution, with a random effect absorbing between-animal variation. The README describes this as robust to the single WT animal.

**Technical approach.**
- **Scores:** `sc.tl.score_genes` on log-normalised `.X` for PIG (Chen 2020, 52 genes), "microglia" (`MICROGLIA_HOMEOSTATIC`, 9 genes), endocytosis (13), DAM (Keren-Shaul 2017, 20) and ARM (Sala Frigerio 2019, 14).
- **Model:** `statsmodels.MixedLM` with `score ~ C(treatment, ref=WT) * C(region)` and a random intercept per sample_id, fitted by ML with lbfgs. The markdown describes a (1|mouse) + (1|section) model, but only one grouping is fitted, and mouse equals section equals sample_id.
- **Reference levels:** WT for treatment and Caudoputamen for region.

**Results.** The notebook has no saved outputs; the numbers below come from the tracked table. Main effects relative to WT in Caudoputamen:

| Signature | BRICHOS coefficient | PBS coefficient |
|---|---|---|
| PIG | +0.128 | +0.154 |
| Homeostatic ("microglia") | +0.119 | +0.120 |
| Endocytosis | +0.077 | +0.091 |
| DAM | +0.082 | +0.089 |
| ARM | +0.075 | +0.089 |

All of these have p < 1e-5. The interaction terms are tested against the WT reference only. For PIG, BRICHOS × region interactions reach p < 0.05 in Corpus callosum (0.017), Meninges (0.023), Olfactory areas (0.011) and Supragranular layers (0.036). No BRICHOS-vs-PBS contrast was computed or saved.

Combining the coefficients by hand (my arithmetic, not a model output) gives BRICHOS − PBS point estimates for PIG:
- below zero in Caudoputamen, Infragranular and Supragranular layers, Meninges and Olfactory areas (about −0.03 to −0.04);
- above zero in Corpus callosum, Hippocampal formation, Hypothalamus, Thalamus and Ventricular area.

This is the same split as §36.

**Caveats.**
- Each "AD vs WT" coefficient compares against one WT section, which is the only member of its random-effect group. The model is not truly robust to n=1.
- No slide term was fitted. The notebook's own notes say to add `C(slide_batch)`, but this was never done.
- The "microglia" score is the homeostatic gene set, yet it rises in AD here. That is consistent with it partly tracking microglial abundance at spot resolution.
- Slide confound: the later slide-adjusted microglia-state model in notebook 03/11 reduced the BRICHOS state effects to non-significance (per memory notes).

### 38. More spots move to low plaque-induced gene scores under BRICHOS than under PBS in cortex and caudoputamen (section-level p = 0.02–0.05), but not in meninges
*Notebook:* `notebooks/03/04-responder-spots.ipynb` · *Code:* `utils/attenuation.py` · *Results:* `results/tables/attenuation/responder_fraction_by_region.tsv` · *Figure:* `results/figures/manuscript/responders_WT.png`, `results/figures/manuscript/responders_PBS.png`, `results/figures/manuscript/responders_BRICHOS.png`

**Rationale.** This is a spot-level view that does not depend on the WT animal. It asks where in the tissue BRICHOS spots look less diseased than PBS, and whether "responder" spots cluster by region or along a meningeal-to-parenchymal axis.

**Technical approach.**
- **Score:** the PIG score (`sc.tl.score_genes`, Chen 2020) per spot.
- **z-shift:** `responder_zshift` gives each spot's z relative to the reference spots' mean and SD within the same region. Two versions are computed, one against WT and one against PBS.
- **Responder:** a spot with z below −1 relative to the regional PBS distribution.
- **Test:** the per-section responder fraction in each region, compared with a one-sided Mann-Whitney test (BRICHOS greater than PBS). Sections are the unit.
- **Region key:** `re_annotation_regions`.
- **Maps:** spatial maps per treatment, coloured by the z relative to PBS.

**Results.** These are saved outputs (first committed with outputs 2026-09-16; the table dates from 2026-05-05). Mean responder fraction per section, PBS vs BRICHOS:

| Region | PBS | BRICHOS | Mann-Whitney p |
|---|---|---|---|
| Infragranular layers | 0.066 | 0.105 | 0.024 |
| Supragranular layers | 0.059 | 0.170 | 0.048 |
| Caudoputamen | 0.102 | 0.182 | 0.048 |
| Hypothalamus | 0.140 | 0.163 | 0.19 |
| Olfactory areas | 0.139 | 0.177 | 0.27 |
| Corpus callosum | 0.109 | 0.146 | 0.36 |
| Meninges | 0.107 | 0.080 | 0.87 |

Hippocampal formation, Thalamus and Ventricular area have n=1 BRICHOS, so they were not tested or the test is meaningless. As a sanity check, WT spots score as responders at 0.25–0.70 per region (from the spot-level output; e.g. Supragranular layers 0.704, Ventricular area 0.250).

**Caveats.**
- The README says the reference is the "PBS median", but the code uses the PBS mean and SD.
- The spatial-figure titles say "relative to regional WT", but the plotted colour is the z relative to PBS.
- Nine one-sided tests are not corrected for multiplicity.
- The meninges result contradicts §36, where meninges had the steepest slope.
- Responder status partly tracks lower PIG overall and therefore slide. P24215 sections are shallower, and per the classifier thread (03/12), BRICHOS-ness is roughly low depth. Nothing here separates the BRICHOS effect from the slide.

### 39. A section-label permutation test of the rescue slope reaches p = 0.012 only in infragranular cortex, which is the smallest p-value 84 permutations allow
*Notebook:* `notebooks/03/05-permutation-null.ipynb` · *Code:* `utils/attenuation.py` · *Results:* `results/tables/attenuation/rescue_permutation_null.tsv` · *Figure:* `results/figures/manuscript/figure_story2_regional_rescue.png` (rescue slope + permutation null panel)

**Rationale.** The bootstrap CIs and binomial tests in §36 resample genes and so ignore that only 9 sections carry the treatment contrast. Permuting treatment labels at section level gives an honest null with exact resolution.

**Technical approach.**
- WT is held fixed.
- Every assignment of 3 of the 9 non-WT sections to "BRICHOS" is enumerated: C(9,3) = 84 permutations.
- For each permutation and region, the BRICHOS-vs-PBS LFC is recomputed and the rescue slope refitted without bootstrap.
- The one-sided p is the fraction of null slopes at or below the observed slope.
- The observed slope is read from `regional_rescue_stats.tsv` (§36).

**Results.** The notebook has no saved outputs; the numbers below come from the tracked table.

| Region | Observed slope | Null median | p (one-sided) |
|---|---|---|---|
| Infragranular layers | −0.070 | 0.002 | 0.012 (= 1/84, the floor) |
| Meninges | −0.104 | — | 0.12 |
| Olfactory areas | −0.091 | — | 0.12 |
| Hypothalamus | — | — | 0.25 |
| Supragranular layers | — | — | 0.26 |
| Caudoputamen | — | — | 0.52 |
| Corpus callosum | — | — | 0.82 |

The steepest observed slopes (meninges, olfactory areas) are not significant at section level.

**Caveats.**
- The markdown says "3 BRICHOS vs 7 PBS, C(10,3)=120 permutations", and the README repeats it. The data have 6 PBS sections, so there are 84 permutations.
- The observed slope (from §36) was fitted over the 190-gene union of data-driven and PIG genes. The null slopes are fitted over the 167 data-driven up genes only, so the observed value and the null use different gene sets.
- Permuting labels does not break the slide structure: two of the three real BRICHOS sections are the only non-WT sections on P24215.
- **Corrected 2026-10-22:** notebook 03/15 reports that after RUVg slide removal the infragranular p stays at the floor of 0.012 under every method. Supragranular layers also reach 0.024 with RUVg k1, but that is weakened with k2. P28052_401 alone does not replicate the cortical rescue.

### 40. The planned confounder checks (section QC, slide batch, cell composition, meningeal rim) were scaffolded but never run
*Notebook:* `notebooks/03/06-confounder-checks.ipynb` · *Code:* `utils/attenuation.py`

**Rationale.** These are robustness checks a reviewer would request:
1. dropping low-quality sections;
2. slide batch as a fixed effect;
3. regressing out deconvolved cell-type fractions before scoring;
4. sensitivity of hippocampal attenuation to spots near the meningeal rim.

**Technical approach.**
- **(1) QC:** per-section QC table (spot count, median UMI, %mt), with `BAD_SAMPLES` left empty. It then re-runs pseudobulk (min_spots=20) and fits a PBS-vs-WT and BRICHOS-vs-WT attenuation slope per region (`attenuation_slope`, genes with padj < 0.05).
- **(2) Slide batch:** listed in the markdown, but there is no code for it.
- **(3) Composition:** requires `obsm['cell_type_fractions']` and skips otherwise.
- **(4) Rim:** requires `obs['distance_to_rim_px']` and skips otherwise.

**Results.** No outputs were saved, and `results/tables/attenuation/qc_filtered_attenuation.tsv` does not exist, so the notebook was never executed. As written, step (1) would produce an empty table: `pseudobulk_lfc` is called with its default `min_n=2`, and WT has one section per region, so every region would raise and be skipped.

**Caveats.** Nothing to interpret. The slide-batch question was later addressed outside this notebook:
- notebook 03/11 (slide-adjusted microglia-state model);
- notebooks 03/12–17 (classifier, within-section contrasts, RUV slide removal and overview).

The composition question was later addressed by the cell2location deconvolution in notebook 03/09.

### 41. Across 21 cell-type and pathway signatures, BRICHOS pushes reactive-astrocyte, oligodendrocyte-stress, immediate-early and microglial disease programs toward WT in cortex, meninges and olfactory areas, but not in corpus callosum or hypothalamus
*Notebook:* `notebooks/03/07-broader-signatures.ipynb` (§1) · *Code:* `utils/gene_sets.py`, `utils/attenuation.py` · *Results:* `results/tables/attenuation/broader_signature_rescue.tsv` · *Figure:* `results/figures/manuscript/broader_signatures_heatmap.png`, `results/figures/manuscript/figure_story3_mechanistic_axis.png`

**Rationale.** PIG is microglia-centric. BRICHOS might also act on astrocyte reactivity, oligodendrocyte stress, barrier integrity, synapses or chaperone networks. A region × signature map shows which biology moves where.

**Technical approach.**
- **Contrast:** for the 7 regions with at least 2 sections per group, the BRICHOS-vs-PBS LFC from §36.
- **Signatures:** every entry in `SIGNATURE_REGISTRY`, 21 at run time:
  - microglial: PIG, DAM, ARM, IRM, MGnD, homeostatic;
  - astrocyte: A1, A2, pan-reactive;
  - oligodendrocyte: mature, OPC, stress;
  - vascular: tight junction, endothelial, pericyte;
  - neuronal and synaptic: neuronal, synapse, immediate-early plasticity;
  - cellular stress and trafficking: endocytosis, heat shock, UPR.
- **Signed score:** each gene's expected disease sign comes from the data-driven disease set (§35) when the gene is in it, otherwise from the signature-level prior (`SIGNATURE_DIRECTION`). The score is the median of signed LFC (positive = rescue) and the fraction of genes moved in the rescue direction.
- No significance test is applied.

**Results.** These are saved outputs (run 2026-05-05). Positive median signed LFC means rescue.

| Signature | Infragranular | Supragranular | Meninges | Olfactory | Caudoputamen | Corpus callosum | Hypothalamus |
|---|---|---|---|---|---|---|---|
| PIG | +0.114 | +0.110 | +0.169 | +0.120 | +0.059 | −0.163 | −0.052 |
| DAM | +0.083 | +0.050 | +0.227 | +0.093 | +0.084 | −0.145 | −0.074 |
| MGnD | +0.089 | +0.086 | +0.274 | +0.081 | +0.088 | −0.217 | −0.017 |

Other notable signatures:
- **A1 astrocyte:** +0.27 to +0.41 in cortex, meninges, olfactory areas and caudoputamen, with concordance of 0.75–1.00.
- **Oligodendrocyte stress:** +0.22 to +0.43 in the same regions.
- **Immediate-early plasticity:** +0.44 in infragranular and +0.34 in supragranular layers.
- **Programs that move against rescue in most regions:** endocytosis is negative in all 7 regions (−0.05 to −0.12). Tight junction, pericyte and mature oligodendrocyte scores are mostly negative.
- **Corpus callosum and hypothalamus:** largely negative for the microglial programs. IRM is −0.30 in hypothalamus.

**Caveats.**
- The method is descriptive: medians over 7–52 genes, with no null or multiple-testing control.
- The 3 vs 6 section contrast carries the slide confound (see §36 corrections).
- The immediate-early gene "rescue" may be a depth or handling artifact. Notebook 03/12 found Egr1, Arc and Nr4a1 to be the top treatment-classifier genes, highest in 401 and P28052_301.
- `SIGNATURE_REGISTRY` now holds 23 signatures (BA4_MIC7 and BA4_MIC7_core were added 2026-09-16), but the saved table and heatmap have only 21. The notebook was not re-run after the Rexach sets were added, despite the memory note saying 07 "picks them up".
- A copy of the heatmap also sits at `notebooks/03/broader_signatures_heatmap.png`.

### 42. Few genes outside the PIG list pass FDR 0.1 for BRICHOS vs PBS in any region, and GSEA ranks interferon and antiviral responses as lower under BRICHOS in cortex, olfactory areas and meninges
*Notebook:* `notebooks/03/07-broader-signatures.ipynb` (§2–3) · *Results:* `results/tables/attenuation/discovery_top_DEGs_per_region.tsv`, `results/tables/attenuation/gsea_*.tsv` · *Figure:* `results/figures/manuscript/figure_story3_mechanistic_axis.png` (GSEA panels)

**Rationale.** This is discovery mode. It asks which non-PIG genes and pathways BRICHOS modulates in each region, to suggest mechanisms beyond microglial plaque-induced genes.

**Technical approach.**
- **DEG discovery:** per region (the 7 testable regions), BRICHOS-vs-PBS Welch t-test on log-CPM, padj < 0.1, PIG genes excluded, top 30 by |LFC|, annotated with signature membership.
- **GSEA:** `gseapy.prerank` on the ranked BRICHOS-vs-PBS LFCs against GO_Biological_Process_2021 (mouse), keeping the top 20 terms per region.

**Results.** These are saved outputs (GSEA run 2026-05-05).
- **Discovery:** only 6 genes pass padj < 0.1 in total, and none belongs to any signature:

| Region | Genes passing padj < 0.1 (LFC, padj) |
|---|---|
| Caudoputamen | Dvl3 (−0.23), Canx (+0.13), padj 0.065 |
| Corpus callosum | Abra (−0.63, padj 0.087) |
| Hypothalamus | Slc25a47 (−0.44, padj 0.034) |
| Supragranular layers | Kif1bp (−0.14), Cab39 (−0.13), padj 0.047 |
| Infragranular layers, Meninges, Olfactory areas | none |

- **GSEA, top 5 terms by NES per region** (negative NES = lower in BRICHOS):
  - **Supragranular layers:** response to interferon-gamma (NES −2.37), defense response to virus (−2.25), innate immune response (−2.19).
  - **Infragranular layers:** response to interferon-gamma (−2.22). Regulated secretory pathway is positive (+2.01).
  - **Olfactory areas:** viral entry (−2.11), interferon-gamma (−2.06), interferon-beta (−2.01).
  - **Meninges:** interferon-gamma (−1.85), T cell proliferation (−1.84).
  - **Corpus callosum** (opposite direction): positive regulation of inflammatory response (+2.09), defense response (+2.07) and neutrophil migration (+2.03) are higher in BRICHOS.
  - **Caudoputamen:** positive chemotaxis (+2.24) and neuropeptide signalling (+2.22).
  - **Hypothalamus:** cholinergic transmission (+2.14).
  - Nominal p is 0.001, the permutation floor, for most top terms.

**Caveats.**
- With 3 vs 6 sections, gene-level FDR is essentially powerless. The few hits are small effects with no clear biology.
- The GSEA output reports only nominal p-values; FDR q-values were not reported.
- Both analyses carry the slide confound. Notebook 03/12 found that hypothalamic neuropeptides (Avp, Oxt, Pmch, Hcrt) and the slide drive treatment separation, which likely contaminates the caudoputamen and hypothalamus terms.
- The interferon "calming" in cortex agrees in direction with the inflammation thread (03/10–11). That signal also attenuated to non-significance once slide was adjusted.

### 43. The test of whether BRICHOS induces the Rexach 2024 BA4_MIC-7 (ITM2B⁺) resilience microglial program was built but never executed, so there is no result
*Notebook:* `notebooks/03/08-rexach-mic7.ipynb` · *Code:* `utils/gene_sets.py`, `utils/orthologs.py`, `utils/attenuation.py`

**Rationale.** A collaborator pointed to Rexach et al., *Cell* 2024. That paper describes BA4_MIC-7, an AD-specific microglial state in motor cortex (a region spared in AD) that the authors read as protective. Its top up-marker is ITM2B, which encodes Bri2, the BRICHOS-domain precursor dosed into these mice. The question is whether BRICHOS shifts cortical spots toward this program. The direction is left open, because BA4_MIC-7 is up in AD yet read as protective.

**Technical approach.**
1. **Gene list:** the BA4_MIC-7 AD up-genes from Rexach Table S5 (padj < 0.05), mapped human→mouse with `utils.orthologs`. An assertion checks the mapped list against `BA4_MIC7_UP_REXACH2024`.
2. **Five signatures:**
   - full (18 mouse genes);
   - "core" (13 genes, dropping the probable ambient contaminants Plp1, Grid2, Dpysl2, Cntn2 and Ppdpf);
   - amyloid processing (4 genes);
   - chaperone-mediated autophagy (4 genes);
   - oxidative (5 genes).
   The notebook also scores a microglia-content proxy (Csf1r, C1qa/b/c, Fcer1g, Ctss, Tyrobp, Aif1, Hexb).
3. **Itm2b itself:** pseudobulk LFC per region for PBS vs WT and BRICHOS vs PBS, plus spatial maps.
4. **Regional shift:** median BRICHOS-vs-PBS LFC per signature and region, with a Wilcoxon test.
5. **Mixed model:** `score ~ treatment * region + microglia_content + (1|section)`, with PBS as reference.
6. **Responder map:** z relative to regional PBS above +1 counts as "induced"; Mann-Whitney test per region.
7. **Specificity:** Jaccard overlap and spot correlation with PIG and DAM, and a residual after regressing out microglia content and PIG.

**Results.**
- **None from this notebook.** It has no saved outputs, and none of its outputs exist on disk: `itm2b_regional_lfc.tsv`, `ba4_mic7_*.tsv`, `itm2b_*.png` and `mic7_zshift_*.png`. The memory note confirms it was "not yet executed" because `/Volumes/processing2` was unmounted.
- **Provenance check (my own read of `data/external/rexach2024_BA4_MIC7_AD_DE.tsv`):** 3,134 genes were tested in Rexach; 17 are up and 33 down at padj < 0.05. ITM2B has the top LFC at 0.74.
- **Indirect (from §35's `disease_signature.tsv`):** Itm2b is not changed globally between PBS and WT (LFC −0.012, padj 0.97).

**Caveats.**
- `notebooks/03/09b-deconvolution-results.ipynb` says its result "dovetails with the BA4_MIC-7 resilience result in `03/08`", but no such result exists.
- The Rexach sets were added to `SIGNATURE_REGISTRY`, but notebook 07 was not re-run, so BA4_MIC-7 is absent from §41 as well.
- If run, the analysis would be a hypothesis-generating port of a human snRNA-seq state onto 55 µm mouse Visium spots, not a validation.
- It would inherit the treatment/slide confound.
- `data/` is gitignored, so the Rexach table is not version-controlled.

### 44. Early Tangram mapping gave near-uniform, anatomically implausible spot calls and left no outputs that later analyses use
*Notebook:* `notebooks/00/tangram/12_tangram_take2.ipynb`

**Rationale.** This was the first attempt at per-spot cell-type composition for the BRICHOS Visium data. It was meant to give cell-type context, and fractions to regress out, for regional analyses. `notebooks/01/06-celltype-regional-analysis.ipynb` expects a Tangram-annotated `data/ST_BRICHOS_with_tangram.h5ad`. `notebooks/03/06-confounder-checks.ipynb` also plans to regress "tangram cell-type fractions" out of signature scores.

**Technical approach.** Tangram (`tg.map_cells_to_space`, `mode='clusters'`) was used with two references. (a) Linnarsson/Zeisel adolescent mouse CNS (`linnarsson_adolescence_full.h5ad`, 160,796 cells, `Class` = 7 types: Neurons, Oligos, Astrocytes, Vascular, Ependymal, PeripheralGlia, Immune). Run in `12_tangram_take2.ipynb` with 16,562 training genes and 300 epochs, and again in `Untitled.ipynb` on MPS. (b) Allen ABC atlas WMB-10Xv2 **midbrain (MB)** subset (29,878 cells, 26 `class` labels). Run in `12_tangram-Copy1.ipynb` / `12_tangram.ipynb` with 1,000 seurat_v3 HVGs (857 training genes) and 500 epochs. Both data sets were normalize_total + log1p. Per-spot `tangram_ct_pred` was row-normalised and argmax-called. The Linnarsson reference-prep notebook is `notebooks/00/tangram/14_linnarsson_mouse_adolenscence.ipynb`. The notebooks entered git in a "refactor" commit on 2026-01-27; the actual run date is not recorded.

**Results.** Linnarsson run: final mapping score 0.470 (0.473 in the MPS rerun). Raw `tangram_ct_pred` values were ~1e-5 for every class. After row-normalisation each class took roughly 0.05–0.23 of a spot (near-uniform). Argmax calls across 25,660 spots were Neurons 5,353, **Immune 4,555**, Vascular 3,810, Astrocytes 3,732, Oligos 3,205, Ependymal 2,880 and PeripheralGlia 2,125. ABC-midbrain run: top labels were dominated by GABA/Glut neuronal classes (07 CTX-MGE GABA 2,770; 03 OB-CR Glut 2,674; 04 DG-IMN Glut 2,426 …). Only 128 spots were OPC-Oligo, 95 Astro-Epen, 79 Vascular and 21 Immune.

**Caveats.** Immune as the top call in ~18 % of spots is not plausible for brain tissue. The midbrain-only ABC subset does not match coronal forebrain sections. The spatial input was already log-transformed before re-normalisation (scanpy warning "adata.X seems to be already log-transformed"). No Tangram output file exists in `data/`, and `01/06-celltype-regional-analysis.ipynb` has no executed outputs. No memory note documents this work.

**Next.** Per-spot deconvolution moved to cell2location with a combined reference (§45–§46).

### 45. A combined Linnarsson + Chhatbar reference (25,113 cells, 18 types) supplied the brain background and fine immune identities for deconvolution
*Notebook:* `notebooks/03/09-microglia-deconvolution.executed.ipynb` (§0–3) · *Code:* `utils/deconvolution.py`, `utils/export_seurat_atlas.R`

**Rationale.** The Chhatbar et al. 2026 mouse CNS-immune atlas (Zenodo 10.5281/zenodo.16938034) contains immune cells only. On its own it would force neuronal signal into immune states, so it was grafted onto a whole-brain background reference.

**Technical approach.** The immune atlas was loaded from a pre-subsampled h5ad (13,287 cells × 33,040 genes). The alternative route exports the Seurat RDS to a Matrix-Market bundle via `utils/export_seurat_atlas.R`. Labels came from `Subclass_label`, which has 12 identities; `Module_label` (8 functional modules) was deliberately not used. The background was Linnarsson adolescent CNS, `Class`, capped at 2,000 cells per class. Its own `Immune` class was stripped (`dv.strip_myeloid`; dropped labels printed `['Immune']`). Immune identities were grafted with a `MYE:` prefix (`dv.build_combined_reference`) and re-capped at 2,000 per type. `ref_source` was kept as the batch key. Helpers live in `utils/deconvolution.py`: mtx loader, reference builder, `subsample_per_type`, and q05 → fraction conversion. The work ran on the Tailscale peer in conda env `c2l`, CPU only.

**Results.** The combined reference had 25,113 cells × 19,917 genes and 18 types: 6 brain classes plus 12 `MYE:` states. Neurons, Astrocytes, Oligos, PeripheralGlia, Vascular, MYE:Microglia, MDM, CAM and CD8 TRM have 2,000 cells each. Ependymal has 1,826, and the rare states run down to NK at 336 and Treg at 353. After `filter_genes` (cell_count_cutoff 5, cell_percentage_cutoff2 0.03, nonz_mean_cutoff 1.12), the signature matrix was 13,070 genes × 18 types. 11,844 of those genes were shared with the Visium data.

**Caveats.** There is a cross-dataset batch effect between the two references. The background is coarse (neurons are not sub-typed). Several immune states have few cells (<400). The plain notebook's interpretation text and the memory note say "15,428 genes shared with Visium", which conflicts with the 11,844 printed in the executed output.

### 46. First cell2location run anatomised the brain correctly but over-assigned immune cells (~24 %/spot); its BRICHOS < PBS microglia trend later proved to be slide
*Notebook:* `notebooks/03/09-microglia-deconvolution.executed.ipynb` (§4–8) · *Code:* `utils/deconvolution.py` · *Results:* `results/tables/deconvolution/` · *Figure:* `results/figures/manuscript/deconv_MYE_Microglia_PBS.png`, `results/figures/manuscript/deconv_myeloid_frac_PBS.png`, `results/figures/manuscript/deconv_myeloid_region_heatmap.png`

**Rationale.** The aim was to estimate per-spot cell-type composition, resolve microglia and other myeloid states, and test whether BRICHOS shifts them relative to PBS region by region.

**Technical approach.** cell2location 0.1.5 on CPU (no CUDA, and MPS is impossible because of float64). Reference: §45. `RegressionModel` used batch `ref_source`, 250 epochs and 1,000 posterior samples. Spatial `Cell2location` used `batch_key=sample_id`, `N_cells_per_location=30`, `detection_alpha=20` and 5,000 epochs (30,000 recommended on GPU), all 25,660 spots and 11,844 genes. Fractions are q05 abundances (`q05_cell_abundance_w_sf`) normalised per spot. A region × treatment composition table was produced. BRICHOS vs PBS was tested per region × state on sample-level means with Mann–Whitney U. The run completed 2026-09-16; per the memory note, spatial mapping took ~3 h.

**Results.** Mean composition was Neurons 0.375, Astrocytes 0.172, Oligos 0.082, Ependymal 0.055, Vascular 0.043, PeripheralGlia 0.034 and MYE:Microglia 0.030. Every other MYE: state sat at ~0.019–0.023 (Plasmablasts 0.023, CAM 0.021, MAIT/Treg/Th17/cDC1/MDM 0.020, NK 0.019). Total immune was 0.239 per spot. The notebook interpretation says regional composition tracks anatomy: corpus callosum oligo-rich (24 %), meninges astro- and vascular-rich, ventricular area ependymal-rich. In spot-level region means, microglia were lowest in WT in all 10 regions, and BRICHOS was below PBS in 9/10. Examples: Caudoputamen 0.0256 vs 0.0284, Meninges 0.0336 vs 0.0370, Olfactory 0.0326 vs 0.0387; corpus callosum was the exception. In the sample-level test (9 regions × 12 states = 108 rows), microglia were lower in BRICHOS in 7/9 regions. Five rows hit the discrete floor p = 0.0238: Infragranular Microglia, MDM in Olfactory and in Meninges, Ventricular MAIT and Ventricular Plasmablasts.

**Caveats.** The 24 % immune share is inflated. The eight adaptive and peripheral states sit at a flat ~1.5–2 % floor with no spatial structure, i.e. signature bleed (notebook part 8 / memory note). p-values sit at the 3-vs-6 MWU floor, with no FDR. Treatment is confounded with slide.

**Corrected 2026-10-07:** The BRICHOS < PBS microglia trend is explained by Visium slide P24215, not treatment (see §51).

### 47. Pruning the reference to innate myeloid types still left a ~17 % immune floor; post-hoc renormalisation showed only MYE:Microglia is separable
*Code:* `utils/run_deconv_pruned.py`

**Rationale.** Run 1 spread immune mass across 12 similar rare signatures (§46). Dropping the adaptive lymphoid states was expected to sharpen the myeloid estimates.

**Technical approach.** `run_deconv_pruned.py` is a headless script that reuses the cached combined reference. It keeps the 6 brain classes plus `MYE:Microglia`, `CAM`, `MDM`, `cDC1` and `Granulocyte`, dropping CD4 Th17, CD8 TRM, MAIT, NK, Treg, TEx and Plasmablasts. The c2l settings match run 1: regression 250 epochs; spatial `N_cells_per_location=30`, `detection_alpha=20`, 5,000 epochs on CPU. It writes to `data/external/deconv_reference_pruned/`, `results/tables/deconvolution_pruned/` and `deconv_pruned_*` figures. It ran on the peer on 2026-09-17. The script was committed 2026-09-17 (8d1ff84). On 2026-09-18 the saved q05 posterior was renormalised post hoc after dropping further types (memory note).

**Results.** All figures come from the memory note. Total immune was 17.4 % per spot. A uniform ~3 % floor from cDC1/MDM/Granulocyte remained, and CAM did not resolve to border/ventricular niches (flat ~0.03 everywhere). Post-hoc renormalisation gave a total immune of 8.8 % without cDC1/MDM/Granulocyte, and a microglia-only fraction of 5.2 %. The microglia fraction itself barely moved (4.5 → 5.2 %) across these scenarios. BRICHOS < PBS microglia in 9/10 regions (Olfactory −0.73 pp, Meninges −0.49 pp), with WT lowest. The diagnosis: at Visium resolution with this atlas, only `MYE:Microglia` is separable, and rare peripheral myeloid signatures act as ambient collectors.

**Caveats.** None of this can be checked from the repo. `results/tables/deconvolution_pruned/` and the `deconv_pruned_*` figures are not on local disk, and no notebook executes this run. The 5.2 % microglia-only estimate does not match the 10.1 % the dedicated microglia-only run produced (§48).

**Corrected 2026-10-07:** The BRICHOS < PBS microglia direction is a slide P24215 offset (§51).

### 48. The microglia-only cell2location run (run 3) gives ~10 % microglia per spot, with WT lowest in every region
*Notebook:* `notebooks/03/09b-deconvolution-results.ipynb` (§1–2) · *Results:* `results/tables/deconvolution_microglia/` · *Figure:* `results/figures/manuscript/deconv_micro_MYE_Microglia_PBS.png`, `results/figures/manuscript/deconv_micro_myeloid_region_heatmap.png`

**Rationale.** Runs 1–2 showed that only microglia are separable (§46–§47). This run collapses the immune compartment to a single `MYE:Microglia` class, so mass cannot bleed into junk-collector signatures.

**Technical approach.** The run used `utils/run_deconv_micro.py`, a copy of the pruned script with `KEEP_MYELOID={'MYE:Microglia'}` and `N_cells_per_location=12`, on CPU with 250 + 5,000 epochs. It was launched on the peer 2026-09-18, and its output tables are dated 2026-09-19. The reference is 6 Linnarsson classes plus Microglia. The c2l signature gene set has 11,264 genes (stated in 09e/09d and reused by stereoscope). Outputs are `per_spot_fractions.csv` and `region_treatment_composition.csv`, plus `deconv_micro_*` spatial maps. A compact `spatial_bundle_micro.h5ad` was pulled back for 09c. 09b was repointed to this run on 2026-09-25 (commit 4a06b33).

**Results.** 25,660 spots × 7 types. Mean composition: Neurons 0.461, Astrocytes 0.176, **MYE:Microglia 0.101**, Oligos 0.091, Ependymal 0.061, Vascular 0.059, PeripheralGlia 0.050. Total structural 0.899. Regional microglia (spot-level means) ranged from 0.060 (WT, ventricular area) to 0.116 (PBS, caudoputamen). WT was lowest of the three groups in all 10 regions (e.g. Meninges WT 0.078 / BRICHOS 0.096 / PBS 0.106).

**Caveats.** `run_deconv_micro.py` is not in the repo; its settings are known only from the memory note and the stereoscope docstring. ~10 % is about double the "~5 % expected" mentioned in 09e. WT is a single section on slide P24215, which runs low on microglia (§50–§51).

### 49. Per-region review of run 3 showed BRICHOS below PBS in all 9 testable regions (−0.3 to −0.9 pp), an apparent partial normalisation later traced to slide
*Notebook:* `notebooks/03/09b-deconvolution-results.ipynb` (§3–5) · *Results:* `results/tables/deconvolution_microglia/myeloid_state_BRIvPBS_shift.tsv`

**Rationale.** 09b is a lightweight review of run 3's tables that needs no compute. It checks whether BRICHOS moves microglia relative to PBS in each region.

**Technical approach.** Per region, sample-level mean microglia fractions were compared between BRICHOS and PBS with Mann–Whitney U. Regions with < 2 sections per arm are dropped, which excludes thalamus. The bar chart is of BRICHOS − PBS.

**Results.** BRICHOS < PBS in 9/9 regions. Meninges 0.1055 → 0.0966 (−0.0089, p = 0.167); Hippocampal formation −0.0085 (p = 0.571); Caudoputamen −0.0082 (p = 0.095); Ventricular area −0.0069 (p = 0.048); Olfactory −0.0067 (p = 0.262); Supragranular −0.0046; Corpus callosum −0.0039; Hypothalamus −0.0034; Infragranular −0.0030 (p 0.26–0.55). At the time this was read as BRICHOS partly normalising injection-induced microgliosis toward WT.

**Caveats.** 09b's caveats cell (part 5) still describes run 1 (24 % immune, CAM, p floor 0.0238). Commit 4a06b33 claimed that caveat was removed. The header still says "completed 2026-09-16". There are few replicates and no FDR correction.

**Corrected 2026-10-07:** Region-matched and slide-adjusted analysis in 09c (§51) attributes the uniform BRICHOS deficit to slide P24215. BRICHOS 401, the only BRICHOS section on the PBS slide, sits inside the PBS range.

### 50. Per-spot pie maps and a z-scored "most-enriched" map recover expected anatomy (oligos in white matter, ependymal on ventricles)
*Notebook:* `notebooks/03/09c-spatial-pie-maps.ipynb` (§1–5) · *Figure:* `results/figures/manuscript/deconv_pie_coarse_PBS.png`, `results/figures/manuscript/deconv_pie_zoom_regions_P24215_101.png`, `results/figures/manuscript/deconv_dominant_enriched_P24215_101.png`

**Rationale.** These maps check visually that the deconvolved compositions are anatomically sensible. Pies blend colours and neurons are the plurality almost everywhere, so extra views were added for readability.

**Technical approach.** The notebook reads a compact spatial bundle (coords, H&E, `prop_*`). It was first built on run 1 (2026-09-17, 8 coarse categories plus a myeloid-compartment split) and repointed to run 3's `spatial_bundle_micro.h5ad` on 2026-09-25. The myeloid-split pies were then dropped, which left part 4 empty. Views:
- per-treatment scatterpie maps (6 structural types + Microglia, cropped to tissue);
- a single-section zoom;
- "Zoom-in QC" windows of ±6 spots around region centroids with known expected composition;
- a "Dominant-type segregation map", where each type is z-scored across spots and each spot is assigned its most-enriched type.

Sections are numbered part 1–5; the QC and segregation views are unnumbered sections that follow them.

**Results.** 25,660 spots across 10 samples (PBS 15,107, BRICHOS 7,596, WT 2,957). Most-enriched counts over all spots: Neurons 7,817, Astrocytes 4,340, Microglia 3,892, Oligos 3,200, PeripheralGlia 2,919, Vascular 2,216, Ependymal 1,276. Per the notebook, the zoom-QC targets were Corpus callosum → Oligos, Ventricular area → Ependymal, Supragranular → Neurons and Caudoputamen → Neurons + Microglia; the result is judged visually from the figure. The segregation map places ependymal on the ventricles, oligos in white matter and vascular/peripheral glia at the meningeal rim.

**Caveats.** Pies show proportions, not counts. Z-scored argmax shows relative enrichment, not abundance (§51). The 09c intro text still says "Pass-1 fractions (18 types)". `deconv_pie_myeloid_{WT,PBS,BRICHOS}.png` are stale run-1 figures dated 2026-09-17.

### 51. Section means show no microglia-abundance effect of BRICHOS: the argmax grid exaggerated small shifts, and region-matched −5 % becomes +1 % once slide is modelled
*Notebook:* `notebooks/03/09c-spatial-pie-maps.ipynb` · *Results:* `results/tables/deconvolution_microglia/microglia_abundance_region_matched.csv`, `results/tables/deconvolution_microglia/microglia_section_region_means.csv` · *Figure:* `results/figures/manuscript/deconv_microglia_section_means_by_slide.png`, `results/figures/manuscript/deconv_microglia_region_matched.png`, `results/figures/manuscript/deconv_microglia_fraction_grid_shared_scale.png`, `results/figures/manuscript/deconv_dominant_enriched_grid_by_treatment.png`

**Rationale.** In the z-scored microglia grid (`deconv_microglia_zgrid_by_treatment.png`) and the all-sections most-enriched grid, PBS sections looked microglia-rich. The question was whether this reflects abundance. Sections were also cut at different levels: hippocampus and thalamus appear only in 101, 102, 103 and P24215_201.

**Technical approach.** These are the unnumbered 09c sections "Microglia fraction — shared absolute scale", "Per-section microglia abundance, by treatment and slide" and "Region-matched comparison", added 2026-10-07 (commits 60b9eed, 707be59, e02eb75). The steps were:
1. Plot microglia fraction on one shared scale (1st–99th percentile).
2. Compute per-section mean microglia and the share of argmax-"Microglia" spots.
3. Compare PBS vs BRICHOS with MWU (6 vs 3 sections).
4. Restrict to regions with ≥ 20 spots in every section, standardise section means with equal region weights, and fit the mixed model `mg ~ treatment + region (+ slide) + (1|section)` (statsmodels, REML).
5. Describe 401 vs PBS on the same slide, n = 1.

**Results.**
- **Section means.** PBS 0.1055 vs BRICHOS 0.1002 (−5.0 %, MWU p = 0.38). Section means range 0.086 (WT 101) to 0.112 (PBS 203), about 1.3×. The argmax red-spot share ranges 4.0 % to 26.8 %, about 7×.
- **Region-matched.** Seven common regions: Caudoputamen, Corpus callosum, Hypothalamus, Infragranular, Meninges, Olfactory, Supragranular. BRICHOS vs PBS was −5.2 % region-standardised (MWU p = 0.26). The region-adjusted mixed model gave −5.2 % (p = 0.071) and **region + slide-adjusted gave +1.2 % (p = 0.70)**.
- **Same-slide comparison.** BRICHOS 401 (P28052), region-standardised, is +1.2 % vs the PBS mean (0.1077 vs PBS 0.1025–0.1105).
- **Slide P24215.** BRICHOS 201 (0.0942), 301 (0.1006) and WT 101 (0.0864) are low in every region.

Conclusion: there is no evidence that BRICHOS changes microglia abundance; the −5 % is the P24215 slide offset.

**Caveats.** The slide term is estimable only through 401 (one BRICHOS section on P28052) and WT. Mixed models gave boundary convergence warnings. Treatment and slide cannot be fully separated.

### 52. Stereoscope (an independent model) puts microglia at 0.17/spot behind a 31 % noise term, and its reference model was still converging at 400 epochs
*Notebook:* `notebooks/03/09e-stereoscope-results.ipynb` (§1–3) · *Code:* `utils/run_deconv_stereoscope.py` · *Results:* `results/tables/deconvolution_stereoscope/` · *Figure:* `results/figures/manuscript/deconv_stereo_MYE_Microglia_PBS.png`

**Rationale.** This is a cross-method check that the run-3 microglia signal is not a cell2location artifact.

**Technical approach.** The run used `scvi.external` RNAStereoscope and SpatialStereoscope (scvi 1.4.2, peer env `c2l`). Reference and gene set match run 3: 6 brain classes + `MYE:Microglia` from the cached combined reference, and the 11,264 c2l run-3 signature genes. The models were trained for 400 reference + 2,000 spatial epochs with batch 1,024 on CPU (`nice 19`, 4 threads). A first launch on 2026-10-07 was stopped at reference epoch 196. The 2026-10-22 run trained fully but crashed in `get_proportions(keep_noise=True)` (scvi `.append` on an ndarray). That was fixed by building the frame from `module.get_proportions`, and proportions were extracted 2026-10-23 with `STEREO_RESUME=1`. The noise "cell type" is saved separately (`noise_fraction.csv`) and fractions are renormalised over the real types. The run used the base h5ad, because the oriented one was not on the peer. ELBO histories are saved as `loss_rna.csv` and `loss_spatial.csv`.

**Results.** Mean composition (noise excluded): Neurons 0.323, Astrocytes 0.176, **Microglia 0.169**, Oligos 0.144, PeripheralGlia 0.097, Vascular 0.067, Ependymal 0.024. The noise term averaged 0.314 per spot before renormalisation. Reference (RNA) ELBO fell from 14.87 M to 6,615 and was still dropping **−40.6 % over the last 100 epochs**. The spatial ELBO had plateaued (−0.01 %).

**Caveats.** The under-trained reference means the signatures are not fully fitted. That may explain part of the gap to c2l (higher microglia, lower neurons, large noise term). A rerun with `STEREO_SC_EPOCHS=2000` needs full retraining; it was not launched and is pending the user's decision. Absolute fractions are not cell counts.

### 53. Stereoscope reproduces where microglia are only moderately (per-spot ρ 0.48, hotspot Jaccard 0.55) but reproduces the slide pattern of section means
*Notebook:* `notebooks/03/09d-stereoscope-vs-cell2location.ipynb` (§1–3) · *Code:* `utils/run_deconv_stereoscope.py` · *Results:* `results/tables/deconvolution_stereoscope/c2l_vs_stereo_agreement.csv` · *Figure:* `results/figures/manuscript/deconv_stereo_vs_c2l_agreement.png`, `results/figures/manuscript/deconv_stereo_vs_c2l_sections.png`

**Rationale.** This tests whether the microglia map and the P24215-low section pattern hold under a second model.

**Technical approach.** 25,660 shared spots and 7 shared types. Steps: per-type Pearson and Spearman correlation of c2l vs stereo fractions; section means by slide; hotspots defined as global z > 1.5 of each method's microglia fraction, compared by overlap. Run locally from CSVs on 2026-10-23 (kernel `brichos-deconv`).

**Results.**
- **Per-spot agreement.** Microglia r = 0.541, ρ = 0.481, the weakest type. The other types range r = 0.63–0.75; Ependymal has ρ 0.337. Mean microglia was 0.101 under c2l vs 0.169 under stereo, and Neurons 0.461 vs 0.323.
- **Section means.** Section-level Spearman between methods was 0.503. Under stereo, PBS 0.1764 vs BRICHOS 0.1623 (MWU p = 0.167); under c2l, 0.1055 vs 0.1002 (p = 0.381). Stereo also puts P24215_301 (0.145) and WT (0.144) low, and 401 (0.176) inside the PBS range (0.165–0.190). However, P24215_201 (0.166) is not low under stereo.
- **Hotspots.** c2l 1,644, stereo 1,300, shared 1,050, Jaccard 0.55. By treatment (c2l / stereo): WT 45/58, PBS 1,123/872, BRICHOS 476/370.

**Caveats.** Stereo's reference model was under-converged (§52). The noise term varies by region (§56).

### 54. The homeostatic-up microglial state shift in BRICHOS reproduces with stereoscope weights, and collapses with slide under both methods
*Notebook:* `notebooks/03/09d-stereoscope-vs-cell2location.ipynb` (§4) · *Results:* `results/tables/deconvolution_stereoscope/c2l_vs_stereo_weighted.csv`, `results/tables/deconvolution_stereoscope/c2l_vs_stereo_hotspot_mixedmodel.csv` · *Figure:* `results/figures/manuscript/deconv_stereo_vs_c2l_state_effect.png`

**Rationale.** Notebooks 10/10c/11 (outside this catalogue part) found that BRICHOS shifts microglia toward a homeostatic state, using c2l microglia weights and hotspots. This entry checks whether that result depends on the deconvolution method.

**Technical approach.** Per-spot signature scores (`sig_Homeostatic`, `sig_DAM`, `sig_DAMratio`) come from `results/inflammation/tables/per_spot_scores.csv`. (A) Section means weighted by each method's microglia fraction, compared PBS vs BRICHOS with MWU (6 vs 3). (B) Hotspot spots (z > 1.5, PBS + BRICHOS) fitted with `score ~ treatment + (1|section)` and with `+ C(slide)`. Scores are z-scored over all spots, so β is in SD units (REML).

**Results.**
- **(A) Weighted section means.** DAMratio Δ was −0.0258 (c2l) and −0.0271 (stereo), both MWU p = 0.0238, the discrete floor. Homeostatic +0.0145 (p = 0.26) and +0.0326 (p = 0.095). DAM −0.0113 (p = 0.38) and +0.0054 (p = 0.71).
- **(B) Hotspot mixed models.** Homeostatic: c2l β = 0.385 SD (p = 0.026, n = 1,599) and stereo β = 0.495 SD (p = 0.003, n = 1,242). With slide: c2l 0.222 (p = 0.41) and stereo 0.287 (p = 0.22). DAMratio: c2l −0.251 (p = 0.089) → −0.090 (p = 0.67) with slide; stereo −0.262 (p = 0.121) → −0.088 (p = 0.72).

**Caveats.** Slide-adjusted models emit convergence warnings; the notebook says betas are stable across optimizers. The state-shift result is not method-dependent but remains slide-confounded.

### 55. Region-matched BRICHOS vs PBS microglia is −5.2 % (c2l) and −7.8 % (stereo, p = 0.03), and both go to ~+1 % once slide is modelled
*Notebook:* `notebooks/03/09d-stereoscope-vs-cell2location.ipynb` (§5) · *Results:* `results/tables/deconvolution_stereoscope/c2l_vs_stereo_region_delta.csv`, `results/tables/deconvolution_stereoscope/c2l_vs_stereo_region_mixedmodel.csv` · *Figure:* `results/figures/manuscript/deconv_stereo_vs_c2l_regions.png`

**Rationale.** This repeats the 09c region-matched abundance test (§51) under both methods.

**Technical approach.** The test uses regions with ≥ 30 spots in all 10 sections; these are the same 7 as §51. It computes per-region BRICHOS vs PBS as a % of each method's PBS mean (MWU), then a section × region mixed model with the response scaled to % of PBS mean, `y ~ treatment + region (+ slide) + (1|section)`.

**Results.** Per region, c2l / stereo (%): Caudoputamen −7.1 / −14.7; Corpus callosum −3.8 / +0.7; Hypothalamus −3.5 / −0.6; Infragranular −2.7 / −9.0; Meninges −8.4 (p = 0.167) / −9.9 (p = 0.048); Olfactory −6.1 / −10.9; Supragranular −4.5 / −15.7. The sign agrees in 6/7 regions, and the region-delta Spearman is 0.536. Mixed model: c2l −5.21 % (p = 0.071) → **+1.23 % (p = 0.70) with slide**; stereo −7.76 % (p = 0.034) → **+0.81 % (p = 0.88) with slide**.

**Reading.** The apparent microglia reduction comes from slide P24215 under both methods, not from a c2l artifact and not from treatment.

**Caveats.** Stereo % effects in cortex and striatum are inflated by noise-term renormalisation (§56). There are 6 vs 3 sections.

### 56. Stereoscope's noise term is regional (4 % ventricular → 55 % infragranular cortex), not slide or treatment, and stretches stereo effects in cortex and striatum
*Notebook:* `notebooks/03/09d-stereoscope-vs-cell2location.ipynb` (§6) · *Results:* `results/tables/deconvolution_stereoscope/noise_fraction.csv`

**Rationale.** Stereo fractions are renormalised after the noise type is dropped. If noise differed by slide or treatment, the renormalisation could create or hide a microglia difference.

**Technical approach.** Mean noise share per section, by slide (MWU) and by region. Per-spot Spearman correlation of noise with each method's microglia fraction.

**Results.** Section means were 0.233–0.384 and split into ~0.24 (P28052_102, 103 and P24215_201) vs ~0.35 (the rest). By slide, P28052 was 0.328 vs P24215 0.295 (MWU p = 0.517). By region: Ventricular 0.037, Hypothalamus 0.076, Thalamus 0.100, Corpus callosum 0.105, Meninges 0.200, Hippocampal 0.336, Olfactory 0.338, Supragranular 0.394, Caudoputamen 0.465, Infragranular 0.549. Spearman of noise vs stereo microglia was −0.521; vs c2l microglia it was +0.307.

**Caveats.** Noise is unrelated to slide or treatment, but regions with ~50 % noise (caudoputamen, infragranular cortex) are where stereo shows its largest BRICHOS–PBS gaps, so read those as amplified. The reference under-convergence (§52) may also inflate the noise term.

### 57. Stereoscope review: BRICHOS < PBS microglia in 5/7 regions (up to −16 %), carried by the two P24215 BRICHOS sections, with only Meninges at p < 0.05
*Notebook:* `notebooks/03/09e-stereoscope-results.ipynb` (§4–7) · *Results:* `results/tables/deconvolution_stereoscope/microglia_BRIvPBS_shift_sections.csv` · *Figure:* `results/figures/manuscript/deconv_stereo_microglia_per_section_regions.png`

**Rationale.** 09e reviews the stereoscope run in the same way 09b reviews run 3 (§49), but its test uses section means.

**Technical approach.** Per section × region mean microglia fraction (≥ 30 spots), BRICHOS (3) vs PBS (6) with MWU; regions need ≥ 2 sections per arm. Per-section dot plots by region, with slide as marker shape. Noise term by region and section, spatial maps, and a caveats list. Written 2026-10-23; plots revised the same day (0863be0).

**Results.** BRICHOS − PBS per region: Caudoputamen 0.1854 → 0.1582 (−14.7 %, p = 0.26); Meninges −9.9 % (**p = 0.048**); Supragranular −15.7 % (p = 0.26); Olfactory −10.9 % (p = 0.71); Infragranular −9.0 % (p = 0.26); Hypothalamus −0.6 % (p = 1.0); Corpus callosum +0.7 % (p = 0.90). Per the notebook, the shift is carried by BRICHOS 201 and 301 on P24215, while 401 sits within the PBS range. With slide in the model, the effect goes to ~0 (09d part 5, −7.8 % → +0.8 %; see §55).

**Caveats.** Listed in 09e part 7:
- the reference model is under-converged (§52);
- the absolute level (0.17) is not a cell count;
- noise is regional and stretches cortex and striatum effects (§56);
- the result is slide-confounded;
- the smallest possible MWU p is 0.024 (6 vs 3);
- per-spot agreement with c2l is only moderate (§53).

### 58. In microglia-enriched spots, BRICHOS raises the homeostatic signature without changing microglia abundance
*Notebook:* `notebooks/03/11-inflammation-results.ipynb` (§1–2, 5) · *Code:* `notebooks/03/10-microglia-inflammation-brichos.py` · *Results:* `results/inflammation/` · *Figure:* `results/inflammation/figs/inflam_state_by_treatment.png`

**Rationale.** Follow-up to the cell2location deconvolution: having found microglia-enriched spots, does BRICHOS reduce neuroinflammation, and if so by recruiting fewer microglia (abundance) or by moving the microglia that are there toward a less inflammatory state? Run 2026-09-25 (commit bf649d1).

**Technical approach.** Full Visium expression (`combined_samples.h5ad`, raw counts; only available on the remote analysis machine, env `sc_py312`), joined by barcode to the deconvolution annotations: 25,660 spots across 10 sections. Counts normalised to 10⁴ and log1p. Five curated signatures scored per spot with `scanpy.tl.score_genes`: Homeostatic (P2ry12, Tmem119, Cx3cr1, Hexb… 12/12 genes), DAM (Trem2, Tyrobp, Cst7, Itgax, Apoe, Cd68… 16/17 genes; Clec7a absent), Neuroinflammation (15), Complement (7) and ReactiveAstro (7). DAM:homeostatic "ratio" = DAM − Homeostatic. Microglia hotspots are spots where the global z-score of `prop_MYE:Microglia` is above 1.5: 1,644 spots (PBS 1,123 / BRICHOS 476 / WT 45). Primary test: a spot-level linear mixed model (score ~ treatment + random intercept per section) on the 1,599 PBS+BRICHOS hotspot spots. Sensitivity test: Mann-Whitney on section pseudobulk means (PBS n=6 vs BRICHOS n=3). WT (n=1) is descriptive only. Differential expression: pydeseq2 on hotspot pseudobulk counts per section (BRICHOS vs PBS, 13,079 genes after filtering), then gseapy prerank of the five signatures (1,000 permutations). Also plotted: a spatial map of the DAM:homeostatic ratio for every section (`inflam_spatial_DAMratio_grid.png`), regional bar means (`inflam_regional_DAMratio.png`) and a volcano plot (`inflam_de_volcano.png`).

**Results.**
- **Abundance did not change.** Microglia fraction in hotspots: PBS 0.151 vs BRICHOS 0.152 (MWU p=0.26).
- **Homeostatic signature was higher under BRICHOS.** Mixed-model β=+0.077, p=0.013; pseudobulk 0.445 → 0.521, MWU p=0.095.
- **DAM:homeostatic ratio was lower, as a trend.** β=−0.040, p=0.056; pseudobulk MWU p=0.17.
- **No change in the other signatures.** DAM β=+0.033 (p=0.40), Neuroinflammation p=0.38, Complement p=0.72, ReactiveAstro p=0.39.
- **GSEA:**
  - Homeostatic NES +2.06 (FDR q<0.001; leading edge P2ry12, Cx3cr1, Hexb, Csf1r…).
  - ReactiveAstro NES −1.81 (q=0.005; Cp, Vim, Gfap).
  - DAM NES +1.46 (q=0.14).
- **No single gene passed FDR.** The best was Ganab, padj=0.18; 0 genes had padj<0.1.
- **Regional means (spot-level, untested).** The largest DAM:homeostatic reductions were in thalamus (PBS 0.54 vs BRICHOS 0.27), olfactory areas (0.29 vs 0.15) and meninges (0.33 vs 0.23).

**Caveats.** Only 6 vs 3 sections. Signatures are scored on the whole spot, not on pure microglia. WT has n=1. The regional bars are spot-level means with no test. Treatment is confounded with slide.

**Corrected 2026-09-25:** with slide added as a covariate, the homeostatic effect falls to p=0.30 and the DAM:homeostatic effect to p=0.61 (§62).

**Corrected 2026-10-14:** within-section niche contrasts show no treatment effect on any microglial state (§67).

### 59. The DAM:homeostatic reduction holds across hotspot thresholds and gets stronger after adjusting for microglia fraction
*Notebook:* `notebooks/03/11-inflammation-results.ipynb` (§3) · *Code:* `notebooks/03/10b-inflammation-sensitivity.py` · *Results:* `results/inflammation/tables/sensitivity_threshold_adjust.csv` · *Figure:* `results/inflammation/figs/inflam_sensitivity.png`

**Rationale.** This checks that the state shift in §58 is not an artefact of residual abundance differences or of the arbitrary z>1.5 hotspot cutoff. Run 2026-09-25 (commit b5faa3e).

**Technical approach.** Same data assembly and signature scoring as §58 (Homeostatic, DAM, DAMratio only). Hotspot thresholds were z > 1.0, 1.5 and 2.0. Each was fitted with the same mixed model (score ~ treatment + 1|section), both without and with microglia fraction as a fixed covariate. Output is the BRICHOS−PBS β with 95% CI and p.

**Results.** Hotspot spot counts (PBS/BRICHOS) were 2427/884 at z>1.0, 1123/476 at z>1.5 and 550/249 at z>2.0.
- **DAM:homeostatic ratio.** Abundance-adjusted p was below 0.05 at every threshold.

  | Threshold | Unadjusted p | Adjusted p | Adjusted β |
  |---|---|---|---|
  | z>1.0 | 4.5e-5 | 2.9e-6 | −0.047 |
  | z>1.5 | 0.056 | 0.022 | −0.044 |
  | z>2.0 | 0.054 | 0.033 | −0.070 |

- **Homeostatic.** Adjusted p was 0.0024 at z>1.0 and 0.016 at z>1.5, but only 0.090 at z>2.0.
- **DAM.** Never significant (adjusted p 0.21–0.74).

**Caveats.** These are spot-level p-values with section as a random effect. With 9 sections, the random-effect variance is poorly estimated. The check does not address the slide confound.

**Corrected 2026-09-25:** robust to threshold and to abundance adjustment, but not to slide (§62).

### 60. Weighting every spot by its microglia fraction (no hotspot cutoff) reproduces a lower DAM:homeostatic ratio under BRICHOS
*Notebook:* `notebooks/03/11-inflammation-results.ipynb` (§4) · *Code:* `notebooks/03/10c-inflammation-weighted-regional.py` · *Results:* `results/inflammation/tables/weighted_summary.csv` · *Figure:* `results/inflammation/figs/inflam_weighted_regional.png`

**Rationale.** This removes the hard hotspot threshold altogether, letting microglia-rich spots dominate continuously so the result is closer to a "pure" microglia state. Part A of 10c, run 2026-09-25 (commit cafef43).

**Technical approach.** Signatures scored as in §58 on all spots. For each section, the mean of each signature is weighted by the deconvolution microglia fraction (clipped at ≥0). PBS (n=6) vs BRICHOS (n=3) is compared with Mann-Whitney on the section values. WT is descriptive.

**Results.**
- **DAM:homeostatic ratio.** PBS 0.324 vs BRICHOS 0.298 (Δ −0.026, MWU p=0.024). WT was 0.253.
- **Homeostatic.** Δ +0.015, p=0.26.
- **DAM.** Δ −0.011, p=0.38.
- **Per section.** The two P24215 BRICHOS sections had the lowest weighted ratios among AD sections (0.285, 0.301). P28052_401 had 0.309. PBS sections ranged 0.313–0.332.

**Caveats.** The lower values come mostly from the two P24215 BRICHOS sections, and P24215 is also the WT slide. 401 sits just below the PBS range. This is the same slide confound described in §62.

### 61. The microglial state shift does not differ significantly by brain region
*Notebook:* `notebooks/03/11-inflammation-results.ipynb` (§4) · *Code:* `notebooks/03/10c-inflammation-weighted-regional.py` · *Results:* `results/inflammation/tables/region_stratified.csv` · *Figure:* `results/inflammation/figs/inflam_weighted_regional.png`

**Rationale.** The regional bars in §58 suggested the strongest calming in thalamus, olfactory areas and meninges. This tests whether the BRICHOS effect is region-specific. Part B of 10c, run 2026-09-25.

**Technical approach.** The unit is hotspot spots (z>1.5) in regions with ≥20 PBS and ≥20 BRICHOS hotspot spots. Eight regions qualified; thalamus was dropped for too few BRICHOS spots.
- **Interaction test.** Likelihood-ratio test of two mixed models with a section random intercept: score ~ treatment × region vs score ~ treatment + region.
- **Per-region deltas.** BRICHOS−PBS differences from section-level means, tested with MWU.

**Results.**
- **No significant treatment × region interaction.** LR p = 0.19 for DAMratio, 0.42 for Homeostatic and 0.60 for DAM. These values come from the 10c script printout as recorded in notebook 11 and in memory; they are not saved in a table.
- **DAM:homeostatic deltas (BRICHOS−PBS).** Six of eight regions were negative and none was significant (all MWU p ≥ 0.17).

  | Region | Delta | MWU p |
  |---|---|---|
  | Meninges | −0.159 | 0.26 |
  | Olfactory areas | −0.134 | 0.17 |
  | Hippocampal formation (only 2 vs 1 sections) | −0.057 | — |
  | Corpus callosum | −0.041 | — |
  | Infragranular layers | −0.034 | — |
  | Supragranular layers | −0.007 | — |
  | Caudoputamen | +0.027 | — |
  | Hypothalamus | +0.029 | — |

**Caveats.** Per-region power is very low (3 vs 6 sections; the minimum MWU p is about 0.024). The earlier "strongest in thalamus" impression was descriptive only, and thalamus could not be tested here. The analysis is not slide-adjusted.

### 62. Slide alone explains the apparent BRICHOS drop in microglia abundance and halves the state effect
*Notebook:* `notebooks/03/11-inflammation-results.ipynb` (§6) · *Results:* `results/inflammation/tables/state_slide_adjusted.csv` · *Figure:* `results/inflammation/figs/abundance_slide_confound.png`

**Rationale.** A spatial z-score map of microglia abundance seemed to show fewer microglia under BRICHOS. Inspecting sections by Visium slide revealed that all six PBS sections are on slide P28052, while two of the three BRICHOS sections (and the WT) are on P24215. This analysis was added 2026-09-25 (commit ec95965) as a caveat to the whole inflammation thread.

**Technical approach.**
- **Abundance.** Section-mean microglia z-score, coloured by slide. Whole-section abundance compared PBS vs BRICHOS with a section-level test.
- **State.** The hotspot mixed model of §58 refitted with slide as an additional fixed effect, using `per_spot_scores.csv`. Within-slide, the only treatment contrast left is 6 PBS vs 1 BRICHOS (P28052_401).
- No script for this step is in the repo. It was computed ad hoc, and only the figure and table were committed.

**Results.**
- **Whole-section abundance.** PBS 0.106 vs BRICHOS 0.100, section-level p=0.38. Hotspot burden was 7.5% vs 6.3% (p=0.55), as stated in the nb 11 text.
- **Slide P24215 runs low in microglia whatever the treatment.** Section-mean z-scores: WT −0.72, BRICHOS 201 −0.44, BRICHOS 301 −0.05. BRICHOS 401, on P28052, sits within the PBS range (+0.32; PBS range about −0.19 to +0.48). These values are read off the figure.
- **The state effects lose significance with slide in the model:**

  | Signature | β unadjusted | p unadjusted | β slide-adjusted | p slide-adjusted |
  |---|---|---|---|---|
  | Homeostatic | 0.077 | 0.013 | 0.045 | 0.30 |
  | DAM:homeostatic | −0.040 | 0.056 | −0.013 | 0.61 |
  | DAM | — | — | — | 0.63 |

**Reading.** Direction is preserved, but the effect roughly halves and is not significant. The homeostatic-restoration signal is internally consistent but hypothesis-generating, not confirmatory.

**Caveats.**
- Slide and treatment are not separately identifiable. The slide-adjusted estimate rests on n=1 (401).
- The hotspot-burden p-value and the abundance test have no code in the repo.
- Notebook 17 (part 2) later shows a region-matched version of this abundance comparison: −5% becomes +1% with slide in the model, p=0.70. That table comes from `results/tables/deconvolution_microglia/` (notebook 09c), not from this run.

### 63. A treatment classifier separates PBS from BRICHOS sections mainly by sequencing depth and slide
*Notebook:* `notebooks/03/17-slide-confound-overview.ipynb` (§4) · *Code:* `notebooks/03/12-treatment-classifier.py` · *Results:* `results/classifier/`

**Rationale.** Question: is there any treatment signal in the spots that a model can recognise in an unseen section, and is that signal anything other than slide? Started 2026-10-14 after the slide confound (§62). Committed 2026-10-23 (a412756).

**Technical approach.**
- **Validation.** Leave-one-section-out (LOSO), because random spot splits leak section identity. 800 training spots are sampled per section. The section score is the mean p(BRICHOS) over the held-out section's spots.
- **Preprocessing.** HVGs, scaling and PCA are refitted per fold on training sections only.
- **Weighting.** Weights give every section, and each class, equal total weight.
- **Model.** L2 logistic regression, C=0.05.
- **Feature sets.**
  - qc: log counts, log genes, %mt, %ribo. This is a technical negative control.
  - sig: 23 curated signature scores.
  - pca: 30 PCs of 2,000 HVGs.
  - hvg: the 2,000 HVGs themselves.
- **Tasks.**
  - Treatment: 9 AD sections.
  - Slide: P24215 vs P28052, 10 sections, as a batch baseline.
  - WT projection: the treatment model applied to the WT section.
- **Nulls.**
  - All C(9,3)=84 BRICHOS relabelings, giving a section-AUC null (floor p=0.012).
  - Slide-conditional test: the P24215 BRICHOS sections stay labelled BRICHOS, and each of the 7 P28052 sections in turn plays the third BRICHOS. This gives 401's rank (minimum p=1/7).
- **Model inspection.** Treatment-vs-slide HVG coefficient correlation from the all-fold models.

**Results.** This was the run with no confound removal (no flags).

| Feature set | Treatment AUC | Perm p | Slide AUC | 401 rank in P28052 |
|---|---|---|---|---|
| qc | 1.00 | 0.012 | 1.00 | 1/7 |
| sig | 0.61 | 0.32 | 0.67 | 1/7 |
| pca | 0.89 | 0.048 | 0.90 | 1/7 |
| hvg | 0.78 | 0.083 | 0.95 | 2/7 |

- **Where 401 ranks.** In hvg, P28052_301 (PBS) ranked first. In pca, 401 (0.444) and 301 (0.434) were nearly tied.
- **Treatment and slide models use the same genes.** Their coefficients correlated at r=0.57.
- **Top BRICHOS-ward genes** were immediate-early genes (Egr1, Arc, Nr4a1, Klf2), hypothalamic neuropeptides (Avp, Oxt) and Hexb.
- **Top PBS-ward genes** were Pmch, Hcrt, Ttr, Pomc and Slc17a7.

**Reading.** The best treatment model is the QC-only model. According to memory and notebook 17, that reflects sequencing depth: P24215 sections have about 9.5–12k median UMI/spot, P28052 PBS about 13.5–20k, and 401 is the shallowest P28052 section (about 12.2k). The remaining signal comes from immediate-early genes and the cutting plane.

**Caveats.** The depth medians are not printed in any saved output of this run. The pca perm p=0.048 is at the edge of significance, and the slide-aligned relabelings also score well. With slide perfectly nested in this design, a slide-free treatment signal cannot be shown.

### 64. Downsampling to equal depth removes the QC-feature "treatment" signal, but slide stays predictable
*Notebook:* `notebooks/03/17-slide-confound-overview.ipynb` (§4) · *Code:* `notebooks/03/12-treatment-classifier.py` · *Results:* `results/classifier_depth/`

**Rationale.** This tests whether removing sequencing depth (the driver found in §63) leaves any treatment signal. Run 2026-10-14.

**Technical approach.** Same pipeline as §63, run with `--downsample 5000 --tag depth`, so every spot is binomially downsampled to exactly 5,000 UMI. Spots below 5,000 UMI were dropped; 24,653/25,660 were kept. Per the run log, immediate-early genes were not dropped and regions were not matched.

**Results.**

| Feature set | Treatment AUC | Perm p | Slide AUC | 401 rank |
|---|---|---|---|---|
| qc | 0.67 | 0.26 | 0.95 | 7/7 |
| sig | 0.61 | 0.27 | 0.67 | 1/7 |
| pca | 0.67 | 0.095 | 0.76 | 2/7 (behind PBS 301) |
| hvg | 0.72 | 0.071 | 0.81 | 2/7 (behind PBS 301) |

- **Treatment and slide models still use the same genes.** Coefficient correlation was r=0.64.
- **Top BRICHOS-ward genes** were still Egr1, Nr4a1 and Arc, plus Oxt and Avp.
- **PBS-ward genes now include immune genes:** Cst3, Igkc, Cd74, C3 and Apoe.

**Caveats.** Downsampling equalises total counts, but library complexity and %ribo still differ by slide, so the QC model still predicts slide. No expression model beats the permutation null.

### 65. Also removing immediate-early genes and matching regions still leaves only slide signal, no treatment signal
*Code:* `notebooks/03/12-treatment-classifier.py` · *Results:* `results/classifier_matched/`

**Rationale.** After depth (§64), the remaining top genes were immediate-early (handling-stress) genes and region-mix / cutting-plane markers. This run removes both. Run 2026-10-14.

**Technical approach.** Flags were `--downsample 5000 --drop-ieg --match-regions --tag matched`.
- **Genes.** 36 immediate-early and stress genes were dropped.
- **Regions.** Only the 7 regions with ≥30 spots in every section were kept: Caudoputamen, Corpus callosum, Hypothalamus, Infragranular, Meninges, Olfactory and Supragranular; 21,867 spots remained.
- **Region balancing.** Each region gets equal weight within a section, in both the training weights and the held-out section score.

**Results.**

| Feature set | Treatment AUC | Perm p | Slide AUC | 401 rank |
|---|---|---|---|---|
| qc | 0.67 | 0.21 | 1.00 | 6/7 |
| sig | 0.50 | 0.30 | 0.71 | 6/7 |
| pca | 0.83 | 0.060 | 0.86 | 2/7 (behind PBS 301) |
| hvg | 0.67 | 0.11 | 0.86 | 2/7 (behind PBS 301) |

- **Treatment and slide models still use the same genes.** Coefficient correlation was r=0.64.
- **Top genes still include hypothalamic neuropeptides and TFs:** Oxt, Sim1, Avp and Shox2 BRICHOS-ward; Pmch, Hcrt and Prl PBS-ward. Region matching therefore does not remove the sub-region cutting plane.

**Caveats.** This run is not shown in notebook 17. Its "depth-matched" row actually loads `classifier_depth` (§64) but describes these settings. Thalamus, hippocampus and ventricular area are excluded.

### 66. Removing hypothalamus and 40 neuropeptide/choroid genes changes nothing: treatment is still not predictable beyond slide
*Notebook:* `notebooks/03/17-slide-confound-overview.ipynb` (§4) · *Code:* `notebooks/03/12-treatment-classifier.py` · *Results:* `results/classifier_noanat/`

**Rationale.** Hypothalamic neuropeptides survived region matching (§65), which pointed to cutting-plane anatomy. This run removes that anatomy as well. Run 2026-10-14.

**Technical approach.** Flags were `--downsample 5000 --drop-ieg --match-regions --drop-anatomy --tag noanat`. This drops the Hypothalamus region and 40 genes: hypothalamic neuropeptides and hormones (Oxt, Avp, Pmch, Hcrt, Pomc…), hypothalamic TFs (Sim1, Otp, Shox2…) and choroid-plexus genes (Ttr, Kl, Folr1…). Six shared regions remained, 17,809 spots.

**Results.**

| Feature set | Treatment AUC | Perm p | Slide AUC | 401 rank |
|---|---|---|---|---|
| qc | 0.67 | 0.30 | 1.00 | 6/7 |
| sig | 0.72 | 0.16 | 0.95 | 4/7 |
| pca | 0.83 | 0.060 | 0.86 | 2/7 (behind PBS 301) |
| hvg | 0.67 | 0.095 | 0.81 | 2/7 |

- **Treatment and slide models still use the same genes.** Coefficient correlation was r=0.62.
- **Top PBS-ward genes** include Slc17a7, Cd74, Igkc, Apoe and Cst3, which are plaque-induced / immune-like genes.
- **Top BRICHOS-ward genes** include Syt1, S100a9 and Prkar1b.

**Reading.** Depth, immediate-early genes, region mix and anatomy are all removable nuisances. Slide is not removable in this design. The PBS-ward immune genes fit the "BRICHOS calms" story, but they cannot be separated from slide.

**Caveats.** n=9 sections. The pca perm p of 0.06 is the closest any model gets.

### 67. Within sections, the plaque-niche DAM/PIG microglial response is fully present under BRICHOS, with no treatment effect
*Notebook:* `notebooks/03/17-slide-confound-overview.ipynb` (§3) · *Code:* `notebooks/03/13-within-section-contrasts.py` · *Results:* `results/within_section/` · *Figure:* `results/within_section/figs/within_section_contrasts.png`

**Rationale.** Slide and depth act mostly as section-wide shifts, so a contrast between two groups of spots inside the same section should largely cancel them. Question: does BRICHOS change the microglial response relative to each section's own background? Started 2026-10-14 as a follow-up to the classifier (§63–66). Committed 2026-10-23.

**Technical approach.** Expression was downsampled to 5,000 UMI/spot. The 7 regions with ≥30 spots in every section were kept (21,874 spots).
- **Signatures.** Homeostatic, DAM (Keren-Shaul 2017), PIG (Chen 2020, 52 genes), MGnD, Neuroinflammation, Complement and ReactiveAstro, plus negative controls: neuronal, synapse, mature oligodendrocyte and endothelial.
- **Contrast.** Within each section × region, spots are ranked by deconvolution microglia fraction. The contrast is hot (top 20%) minus cold (bottom 50%) mean score, plus the OLS slope of score on rank. A section's value is the mean over regions.
- **Tests.**
  - MWU on section values (PBS 6 vs BRICHOS 3).
  - Mixed model on section × region values: diff ~ treatment + region, with a section random intercept.
  - 401's rank among the 7 P28052 sections.
  - Same-slide P24215 BRICHOS − WT as a reference.

**Results.**
- **The plaque-niche response is clear.** Hot−cold contrasts in AD vs WT: DAM 0.17 vs 0.02, PIG 0.22 vs 0.08, MGnD 0.14 vs 0.02, Complement 0.29 vs 0.16.
- **It is present under BRICHOS on WT's own slide.** P24215 BRICHOS − WT: DAM +0.14, PIG +0.15, MGnD +0.11.
- **No treatment effect on any microglial niche contrast.** For example DAM:homeostatic MWU p=0.55 / mixed-model p=0.36, Homeostatic 0.26/0.26, PIG 0.38/0.45.
- **401 shows the opposite of "calming".** It has the highest DAM, PIG and MGnD contrast on P28052 (rank 1/7).
- **Homeostatic slope is weakly higher under BRICHOS** (Δ +0.035, MWU p=0.095, mixed-model p=0.068, 401 rank 1/7). This is the same direction as §58.
- **Negative control: the oligodendrocyte contrast failed.** BRICHOS hot spots are less oligo-depleted (Δ +0.075, MWU p=0.024, mixed-model p=0.003, 401 rank 1/7). The difference is largest in meninges (BRICHOS +0.04 vs PBS −0.22), hypothalamus and caudoputamen.

**Reading.** Cancelling slide removes the "BRICHOS calms microglia" contrast at the niche level, while the AD plaque-niche response persists under BRICHOS.

**Caveats.**
- The oligodendrocyte negative control is unexplained. It is either slide leaking into the contrast or a biological lead.
- Hot/cold quantiles are arbitrary.
- The deconvolution fraction is a proxy for microglia, and only 3 BRICHOS sections are available.

### 68. Within-section centring keeps the regional rescue ranking (cortex, meninges, olfactory) but no region reaches significance, and 401 alone shows the opposite profile
*Notebook:* `notebooks/03/16-slide-confound-rescue.ipynb` (§2) · *Code:* `notebooks/03/14-within-section-rescue.py` · *Results:* `results/within_section_rescue/` · *Figure:* `results/within_section_rescue/rescue_within_section.png`

**Rationale.** In the regional rescue tests of 03/02, the disease axis (PBS on P28052 − WT on P24215) and the rescue axis (BRICHOS mostly on P24215 − PBS) contain the slide effect with opposite signs. Slide alone can therefore produce "rescue". This run is a difference-in-differences version that cancels region-uniform slide effects. Run locally 2026-10-22 (commit 3ef42c8).

**Technical approach.** Input is the (section, region) pseudobulk from 03/01 (`results/tables/attenuation/pseudobulk_*.tsv`): 10 sections × 7 shared regions, 14,927 genes, log2 CPM.
- **Centring.** Each section is centred on its own mean over the 7 regions. Any region-uniform section, slide or depth effect cancels, and so does any global rescue; only the regional profile remains.
- **Disease genes.** PIG (52, slide-free selection) and a "data" set (03/02 PBS-vs-WT genes: 190 up / 29 down).
- **Per-region statistics.** Concordance (fraction of disease genes moved toward rescue) and rescue slope (rescue LFC on disease LFC).
- **Designs.** Uncentred (03/02 analogue), centred, and within-slide (401 vs P28052 PBS, raw and centred).
- **Nulls.** 84 BRICHOS relabelings, or 7 for the within-slide designs.
- **Assumption check.** Variance of the slide contrast (P24215 BRICHOS − 401) is split into region-uniform and region-specific parts and compared with same-slide PBS 2-vs-1 noise splits.

**Results.**
- **Centred PIG design, regional ranking.** Concordance is high in olfactory (0.85), infragranular (0.83), meninges and supragranular (0.79). It is low in corpus callosum (0.06) and hypothalamus (0.10).
- **No region reaches significance.** The best centred permutation p values are about 0.13. The lowest is infragranular PIG concordance (p=0.131); the lowest slope p values are 0.167 (data-set olfactory) and 0.214 (PIG supragranular).
- **401 alone (within P28052, PIG) shows the opposite profile.** Corpus callosum (0.85) and caudoputamen (0.83) are the most rescued; supragranular is the least (0.50).
- **The slide effect is only roughly region-uniform.**
  - The region-uniform slide variance is 0.029, about 1.7× the noise median of 0.017.
  - The region-specific slide variance is 0.070, about 1.3× the noise median of 0.054. It exceeds 88% of the same-slide noise splits.
  - Correlations of each region's slide effect with the region-average slide effect are 0.54–0.62.

**Caveats.** The centring assumption (no slide × region interaction) holds only approximately. The global rescue is lost by construction. The 401-only analysis is n=1.

### 69. After control-gene slide removal (RUVg, k=1), the infragranular cortex rescue survives at the permutation floor
*Notebook:* `notebooks/03/16-slide-confound-rescue.ipynb` (§3) · *Code:* `notebooks/03/15-slide-removal.py` · *Results:* `results/slide_removal/` · *Figure:* `results/slide_removal/rescue_after_slide_removal.png`

**Rationale.** Within-section centring (§68) also removes global rescue. Model-based removal instead tries to subtract slide while keeping treatment, with each method relying on a different assumption. Run locally 2026-10-22 (commit 3ef42c8).

**Technical approach.** Same (section, region) pseudobulk as §68. Four methods:
- **tech (CQN-like).** A gradient-boosting model predicts the slide contrast (P24215 BRICHOS − 401) from gene GC content, gene/transcript length (Ensembl BioMart, `data/external/gene_features/`), mean expression and mito/ribo membership, with gene-wise 5-fold CV. The prediction is subtracted from P24215 samples.
- **RUVs, k=1 and 2.** SVD of residuals within each treatment × region group.
- **RUVg, k=1 and 2.** SVD of region-centred expression of 1,000 control genes: expressed, stable across regions and across PBS sections, and excluding PIG and disease genes.

After correction, factors are regressed out jointly with the treatment × region design. The per-region rescue tests are then rerun, with the disease axis recomputed on corrected data. The null is 84 relabelings, and RUV factors are refitted per permutation. Diagnostics: residual slide variance relative to same-slide noise, and correlation of each factor with slide and treatment.

**Results.**
- **The slide effect is not technical.** The tech model's CV R² for the slide contrast is 0.013, versus 0.029 for same-slide noise, so it removes nothing.
- **What each method removed:**

  | Method | Factor–slide r | Region-uniform slide vs noise | Region-specific slide vs noise |
  |---|---|---|---|
  | No correction | — | 1.69× | 1.30× |
  | RUVs k1 | 0.13 (learns animal noise) | 1.47× | — |
  | RUVs k2 | 0.75 | 1.13× | — |
  | RUVg k1 | 0.97 | 0.42× | 1.27× |
  | RUVg k2 | 0.97 | 0.30× | 1.22× |

  No method reduces the region-specific part, which stays at 1.13–1.27×.
- **PIG rescue-slope permutation p after RUVg k1.** Infragranular 0.024 and supragranular 0.024, versus 0.060/0.060 without correction.
- **Data gene set after RUVg k1.** Infragranular 0.012, supragranular 0.012, meninges 0.036, olfactory 0.048.
- **Data-set infragranular slope p = 0.012 (the floor) under every method:** none, tech, RUVs k1/k2 and RUVg k1/k2.
- **Two factors weaken the PIG result.** With k=2, PIG cortex p rises to 0.12–0.18.
- **Non-responders.** Corpus callosum (p ≥ 0.62) and hypothalamus (p ≥ 0.26) do not respond under any method.

**Caveats.**
- RUVg assumes BRICHOS does not change the control genes.
- Results depend on k.
- The cortex rescue rests on the two P24215 BRICHOS sections; 401 does not reproduce it (§68).
- The "data" gene set was selected on a slide-contaminated PBS-vs-WT contrast.

### 70. A pure slide contrast between BRICHOS sections produces a steeper "rescue slope" than the real treatment comparison
*Notebook:* `notebooks/03/16-slide-confound-rescue.ipynb` (§1)

**Rationale.** This is a direct demonstration that the 03/02 rescue-slope statistic is contaminated by slide. It is a small live computation in notebook 16, built 2026-10-22, and repeated in notebook 17 part 5a.

**Technical approach.** Whole-section pseudobulk (log2 CPM) for the 52 PIG genes. Each contrast's LFC is regressed on the 03/02 disease LFC (PBS − WT). Three contrasts:

| Contrast | Definition | What it isolates |
|---|---|---|
| Slide only | P24215 BRICHOS (201, 301) − P28052_401 | Slide; both arms are BRICHOS |
| Treatment only | 401 − the 6 P28052 PBS | Treatment; same slide, n=1 |
| Original | All BRICHOS − all PBS | Treatment plus slide |

Within-slide reference: each of the 7 P28052 sections in turn is compared with the other six.

**Results.**

| Contrast | Rescue slope | Fraction of PIG genes down |
|---|---|---|
| Slide only | −0.098 | 0.37 |
| Treatment only | +0.026 | 0.67 |
| Original | −0.040 | 0.60 |

- **401 is not unusual on its own slide.** It ranks 2/7 for fraction of PIG genes down; PBS section P28052_202 is higher at 0.71.

**Reading.** The rescue slope is driven by slide, so it is not interpretable as treatment without correction. Concordance is not driven by slide (slide pushes PIG genes up), but there is no within-slide evidence of a global PIG rescue.

**Caveats.** The treatment-only contrast is n=1, and the slide-only contrast includes animal-to-animal noise.

### 71. The slide confound is the study's largest limitation: only the infragranular cortex rescue partly survives
*Notebook:* `notebooks/03/17-slide-confound-overview.ipynb`

**Rationale.** A single overview of every analysis touching the treatment × slide confound, intended for manuscript framing. Built 2026-10-23 (commit a412756). It only loads saved outputs from 09c, 10–11, 12, 13, 14–16, plus the live slide-contrast computation of §70.

**Technical approach.** The notebook brings together:
- the design table and section spot counts;
- the region-matched microglia abundance table (`results/tables/deconvolution_microglia/microglia_abundance_region_matched.csv`, from 09c);
- the slide-adjusted state model (§62);
- the niche contrasts (§67);
- a classifier summary over the original, depth and noanat runs (§63, §64, §66);
- the slide-contrast demo (§70);
- rescue permutation p-values across corrections (§68–69).

**Results (status table).**

| Claim | Status after slide checks |
|---|---|
| AD disease response (PIG/DAM up vs WT) | Robust: present within sections and on WT's own slide |
| BRICHOS lowers microglia abundance | Gone: −5.0% (whole section, p=0.38; region-adjusted mixed model p=0.071) becomes +1.2% with slide in the model (p=0.70) |
| BRICHOS shifts microglia to a homeostatic state | Weakened: p 0.013 → 0.30 with slide |
| BRICHOS calms plaque-niche microglia | Not supported |
| Classifier tells PBS from BRICHOS | Learns slide and depth |
| Regional rescue | Partly survives: infragranular cortex p=0.012–0.024 after RUVg k1, but rests on the P24215 BRICHOS sections |
| Original rescue-slope statistics | Not usable as-is |

The notebook includes a draft manuscript caveat paragraph.

**Reading.**
- **Defensible core.** A clear AD response that BRICHOS does not abolish; a modest cortical rescue that survives control-gene slide removal; no change in microglia abundance.
- **Drop or recast.** "BRICHOS lowers microglia", "BRICHOS calms plaque-niche microglia", and the uncorrected rescue slopes.

**Next.** In order of value:
1. Remaining tissue, so PBS can be placed on the same slide as BRICHOS.
2. Shared samples from the same processing batches.
3. Orthogonal validation independent of the slides: IBA1, CD68 and P2RY12 staining, plaque staining, or qPCR in cortex.

**Caveats.**
- The classifier section labels `classifier_depth` as "downsampled, IEG dropped, regions matched". Per its run log, that run was downsample-only; those settings belong to `classifier_matched`, which is not shown.
- The depth medians quoted (about 9.5–12k vs 13.5–20k UMI) are not in any saved output.
