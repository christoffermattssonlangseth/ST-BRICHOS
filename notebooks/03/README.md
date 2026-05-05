# notebooks/03 — Statistical attenuation analysis

Run order:

1. `01-pseudobulk.ipynb` — collapse spot counts to (mouse, region).
2. `02-attenuation-stats.ipynb` — T1 (LFC), T2 (sign concordance), T3 (slope + bootstrap CI). Produces `regional_attenuation_stats.tsv`.
3. `03-signature-mixed-model.ipynb` — T4 mixed-effects on PIG / microglia / endocytosis scores.
4. `04-responder-spots.ipynb` — per-spot z-shift, responder maps, regional fraction test.
5. `05-permutation-null.ipynb` — mouse-level label permutation, empirical p per region.
6. `06-confounder-checks.ipynb` — QC filter, slide-batch, composition adjustment, meninges-rim sensitivity.

Outputs land in `results/tables/attenuation/` and `results/figures/manuscript/`.

Top-of-notebook variables to set on first run:

* `SAMPLE_KEY` — obs column for mouse/section id (try `sample_id` or `library_id`).
* `REGION_KEY` — obs column for region label (`anatomical_region`, `region`, etc.).
* `TREATMENT_KEY` — usually `treatment` with values WT / PBS / BRICHOS.
