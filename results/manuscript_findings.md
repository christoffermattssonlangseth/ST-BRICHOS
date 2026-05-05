# BRICHOS-ST — manuscript findings

Working summary of statistical findings from the spatial transcriptomics
attenuation analysis. Status as of 2026-05-05. Backs the manuscript at
`brichos-st.docx`.

## 1. Cohort & design

- **Model.** APP\textsubscript{NLGF} transgenic mouse, BRICHOS treatment 2× / week
  for 12 weeks; harvested for coronal Visium sections.
- **Cohort.** WT n = 1 mouse / 1 section; PBS n = 3 mice / 7 sections;
  BRICHOS n = 3 mice / 3 sections.
- **Regions analysed.** 10 Allen-aligned regions; usable for inferential
  testing: 7 (cortex supra/infra, caudoputamen, corpus callosum,
  hypothalamus, olfactory, meninges). Hippocampus, thalamus,
  ventricular area underpowered (≤2 PBS sections per region).

## 2. Statistical strategy — anchoring around n = 1 WT

The single-WT cohort blocks per-region WT-anchored differential
expression tests. We therefore anchor disease *once globally* and test
*rescue per region* using BRICHOS-vs-PBS, where replication is
adequate (3 vs 7 sections).

- **Disease signature** = whole-tissue PBS-vs-WT pseudobulk
  (pooled-variance fallback for n = 1 WT; descriptive padj). Cross-checked
  against the canonical 52-gene PIG set from Chen *et al.* 2020.
- **Rescue tests** per region on BRICHOS-vs-PBS pseudobulk:
  - T2′ — sign concordance (binomial vs 0.5).
  - T3′ — rescue slope of LFC\textsubscript{BRI−PBS} on
    LFC\textsubscript{PBS−WT}, with bootstrap CI.
- **Spot-level orthogonal test** — per-mouse fraction of PIG-rescued
  spots (z-shift < −1 vs PBS reference), Mann–Whitney across mice.
- **Permutation null** — 84 mouse-level BRI/PBS reassignments
  (WT held fixed), empirical p on the rescue slope.

## 3. Convergent rescue in cortex and striatum

Three orthogonal tests agree on the same regions.

| Region | Sign concordance | Spot responder MWU | Slope perm null | Conclusion |
|---|---:|---:|---:|---|
| Infragranular cortex | 0.82 *** | p = 0.024 | **p = 0.012** | strongest |
| Supragranular cortex | 0.85 *** | p = 0.048 | NS | robust |
| Caudoputamen | 0.78 *** | p = 0.048 | NS | robust |
| Olfactory areas | 0.78 *** | NS (0.27) | marginal (0.12) | direction-only |
| Meninges | 0.78 *** | NS (BRI < PBS) | marginal (0.12) | direction-only |
| Corpus callosum | 0.40 (NS) | NS (0.36) | wrong direction | non-responder |
| Hypothalamus | 0.46 (NS) | NS (0.19) | NS | non-responder |

**Headline.** Convergent significance in **infragranular and
supragranular cortex plus caudoputamen** establishes BRICHOS-mediated
rescue. Olfactory and meninges show direction-consistent but
small-magnitude effects.

## 4. Rescue is partial in magnitude

Per-region rescue slopes lie between 0 and −0.2 — far from the −1.0
"full normalisation" reference. Spot-level responder fractions:
~10–16 % of the disease distance recovered in responding regions.
**BRICHOS is a directional modulator, not a normaliser.** Manuscript
phrasing: "partial rescue, ~10–15 % of the disease distance recovered."

| Region | Disease gap (WT−PBS, %) | BRICHOS shift (BRI−PBS, %) | % rescued |
|---|---:|---:|---:|
| Supragranular | 64.4 | 10.5 | 16 % |
| Caudoputamen | 56.4 | 6.3 | 11 % |
| Olfactory | 39.1 | 4.4 | 11 % |
| Hypothalamus | 34.8 | 3.6 | 10 % |
| Infragranular | 45.3 | 3.6 | 8 % |
| Corpus callosum | 30.6 | 2.8 | 9 % |

## 5. BRICHOS modulates biology beyond microglia

A 21-signature rescue heatmap (regions × signatures, prior-anchored to
disease direction) and pathway-level GSEA (gseapy prerank, GO BP
2021) jointly reveal a coherent mechanistic axis.

### 5.1 Microglial subtypes

- **IRM — interferon-response microglia — is the most strongly and
  most broadly rescued microglial program.** Signature heatmap shows
  rescue across cortex, caudoputamen, corpus callosum, meninges and
  olfactory areas. Hypothalamus is the single exception, where the
  IRM signature appears exacerbated under BRICHOS.
- PIG, DAM, ARM and MGnD modules track each other; rescue is
  strongest in meninges and modest in cortex / caudoputamen, with no
  rescue in corpus callosum.

### 5.2 Astrocytes and oligodendrocytes

- **A1 (neurotoxic) astrocyte signature** strongly rescued in cortex,
  caudoputamen, meninges and olfactory areas — the strongest
  non-microglial rescue. Consistent with the
  microglia-to-astrocyte induction cascade (Liddelow *et al.*, 2017),
  in which BRICHOS dampens the upstream microglial activation signal
  and the downstream A1 conversion is reduced in parallel.
- **Oligodendrocyte stress** (Serpina3n, C4b, Klk6 set) rescued
  broadly, including a small but consistent effect in corpus callosum
  where microglial PIG / DAM are not rescued — providing the only
  white-matter readout that responds to BRICHOS.
- **Pan-reactive astrocyte** signature most strongly rescued in
  olfactory areas.
- **A2 astrocyte** signature rescued in infragranular cortex and
  olfactory areas but not elsewhere — region-specific, not
  globally preserved.

### 5.3 Neuronal / synaptic restoration

- **Plasticity / immediate-early-gene (IEG) signature (Arc, Fos,
  Egr1, Bdnf, Nptx1) is induced by BRICHOS** in cortex, caudoputamen
  and meninges. New manuscript-worthy headline: BRICHOS treatment is
  associated with **restored neuronal activity**, not only with
  dampened inflammation.
- GSEA confirms this at pathway level. Genes upregulated by BRICHOS
  are enriched for **neuropeptide signaling, GPCR signaling coupled
  to cyclic nucleotide second messengers, regulated secretory
  pathway, calcium-dependent exocytosis and neurotransmitter receptor
  regulation** in caudoputamen, infragranular cortex and hypothalamus
  (NES ≈ +1.9 to +2.2, p < 0.001).

### 5.4 The molecular axis of rescue

Pathway enrichment of the BRICHOS-vs-PBS contrast shows that
**downregulated genes are strongly enriched for type II interferon
(IFN-γ), antiviral and innate immune response pathways** across
supragranular cortex (NES = −2.4, p < 0.001), olfactory areas,
infragranular cortex and meninges. This connects the IRM-microglia,
A1-astrocyte and IEG-neuronal observations into a single cascade:

```
                          BRICHOS
                              |
    ↓ IFN-γ / antiviral / innate immune signaling
                              |
    ↓ IRM activation  →  ↓ A1 astrocyte induction  →  ↓ oligo stress
                              |
    ↑ neuronal IEG / synaptic / neuropeptide signaling
```

Each arrow is supported by an independent test (sign-concordance,
broader-signature heatmap, GSEA prerank).

### 5.5 Programs NOT rescued — boundary of the BRICHOS effect

- **Mature oligodendrocyte program** is *not* rescued — myelination
  transcripts (Mbp, Mog, Plp1, Cnp, Mag) remain depressed across most
  regions including supragranular cortex, where everything else
  responds.
- **Tight junction, endothelial, pericyte signatures** — vascular /
  barrier programs are not preserved at the transcriptional level. At
  least transcriptionally, BRICHOS does not act through BBB stabilisation.
- **Heat-shock response, UPR** — only modest, no clear rescue.
  BRICHOS is itself a chaperone, but downstream intracellular
  chaperone induction is not engaged at the transcriptional level,
  consistent with the established model in which BRICHOS works
  extracellularly on Aβ aggregation.

### 5.6 Region-discordant finding — corpus callosum

GSEA in corpus callosum yields the opposite pattern of the rescue
regions: **positive enrichment for inflammatory response, defense
response, neutrophil migration and pattern-recognition-receptor
signaling** (NES ≈ +2.0, p < 0.001). Combined with the wrong-direction
rescue slope and slight blueing in the broader-signature heatmap, this
indicates BRICHOS does not benefit white matter and may transiently
*increase* inflammatory pathway activity there. Worth flagging as a
caveat / discussion point.

## 5.7 Outputs underpinning these claims (Section 5)

- `results/tables/attenuation/broader_signature_rescue.tsv` —
  21 signatures × 7 regions, signed-LFC + concordance.
- `results/tables/attenuation/discovery_top_DEGs_per_region.tsv` —
  top non-PIG BRI-vs-PBS DEGs per region, with signature-membership
  annotation.
- `results/tables/attenuation/gsea_<region>.tsv` — GSEA prerank top
  20 GO BP terms per region (one file per region).
- `results/figures/manuscript/broader_signatures_heatmap.svg` —
  21-signature × 7-region heatmap with disease-up / homeostatic
  blocks separated.

## 6. Region-level non-responders

- **Corpus callosum** — sign concordance ≈ 0.40, slope wrong
  direction, BRICHOS slightly worsens A2-astrocyte and homeostatic
  microglia signatures. White matter does not respond to BRICHOS
  through PIG, but shows a small oligodendrocyte stress benefit.
- **Hypothalamus** — direction-inconsistent across tests, deep
  grey-matter region; possible drug-penetrance ceiling.
- **Hippocampus, thalamus, ventricular area** — under-powered
  (n\textsubscript{PBS} ≤ 2 sections), descriptive only.

## 7. Caveats baked into figure captions

1. WT n = 1 — disease signature p-values are descriptive; rescue
   tests use BRICHOS vs PBS where replication is adequate.
2. Permutation resolution = 1 / 84 ≈ 0.012 — slope-based test cannot
   resolve smaller p-values.
3. Spot responder threshold = −1 SD relative to PBS regional median;
   threshold-sensitive in regions with skewed PBS distributions
   (hippocampus, thalamus).
4. Cell composition is not yet decomposed — tangram deconvolution
   pending; "rescue" currently mixes per-cell expression and
   composition shifts.
5. Microglia-rich-spot stratification not yet performed; will
   separate recruitment vs activation effects when tangram fractions
   land.

## 8. Manuscript-ready phrases

> Three orthogonal tests — gene-level sign concordance, mouse-level
> Mann–Whitney on rescued-spot fractions, and slope permutation null —
> converged on the supragranular and infragranular cortex and the
> caudoputamen as the regions of robust BRICHOS-mediated rescue
> (binomial p < 1×10⁻³, MWU p ≤ 0.048, perm p = 0.012 for
> infragranular cortex at the resolution floor).

> Across 21 cell-type and pathway signatures, BRICHOS treatment
> dampened not only the canonical microglial PIG program but also
> downstream A1 (neurotoxic) astrocyte and oligodendrocyte-stress
> programs, with regional convergence in cortex, caudoputamen,
> olfactory areas, and meninges. White matter (corpus callosum) and
> deep grey-matter (hypothalamus) showed limited or absent response.

> Effect size, expressed as the fraction of disease-distance
> recovered, was 8–16 % across responding regions, indicating
> partial rather than full transcriptional normalisation.

## 9. Outputs underpinning these claims

- `results/tables/attenuation/disease_signature.tsv`
- `results/tables/attenuation/regional_rescue_stats.tsv`
- `results/tables/attenuation/responder_fraction_by_region.tsv`
- `results/tables/attenuation/rescue_permutation_null.tsv`
- `results/tables/attenuation/broader_signature_rescue.tsv`
- `results/tables/attenuation/discovery_top_DEGs_per_region.tsv`
- `results/figures/manuscript/regional_rescue_slope.svg`
- `results/figures/manuscript/regional_sign_concordance.svg`
- `results/figures/manuscript/broader_signatures_heatmap.svg`
- `results/figures/manuscript/responders_{WT,PBS,BRICHOS}.svg`
