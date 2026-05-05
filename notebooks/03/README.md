# notebooks/03 — Statistical attenuation analysis

Run order:

1. `01-pseudobulk.ipynb` — collapse spot counts to (sample, region) (min_spots=10 to keep small regions).
2. `02-attenuation-stats.ipynb` — global PBS-vs-WT disease signature (n=1 WT pooled-variance fallback), then per-region BRI-vs-PBS rescue tests: sign concordance + rescue slope with bootstrap CI. Outputs `disease_signature.tsv`, `regional_rescue_stats.tsv`.
3. `03-signature-mixed-model.ipynb` — mixed-effects on PIG / microglia / endocytosis scores at spot level (mouse random effect, robust to n=1 WT).
4. `04-responder-spots.ipynb` — per-spot z-shift vs PBS median (treatment-blind, avoids n=1 WT bottleneck), responder maps, Mann-Whitney on per-mouse responder fraction.
5. `05-permutation-null.ipynb` — mouse-level permutation of BRI/PBS labels (WT held out). C(10,3)=120 unique permutations, empirical p per region for the rescue slope.
6. `06-confounder-checks.ipynb` — QC filter, slide-batch, composition adjustment, meninges-rim sensitivity.
7. `07-broader-signatures.ipynb` — beyond PIG: rescue heatmap across ~20 cell-type and pathway signatures (astrocyte, oligo, barrier, synapse, heat-shock, etc.); discovery of top non-PIG DEGs per region; optional GSEA via gseapy.

Design note — the WT cohort is a single mouse, which blocks per-region WT-anchored tests. The pipeline anchors disease once globally and then tests rescue per region using BRI vs PBS, where replication is fine.

Outputs land in `results/tables/attenuation/` and `results/figures/manuscript/`.

Top-of-notebook variables to set on first run:

* `SAMPLE_KEY` — obs column for mouse/section id (try `sample_id` or `library_id`).
* `REGION_KEY` — obs column for region label (`anatomical_region`, `region`, etc.).
* `TREATMENT_KEY` — usually `treatment` with values WT / PBS / BRICHOS.
