"""Curated gene sets used across notebooks.

Single source of truth so figures, signatures, and tests all reference
the same lists. Update here, regenerate notebooks if applicable.
"""
from __future__ import annotations

# ----- Plaque-induced genes (Chen, Lu, Craessaerts et al., Cell 2020) -------
# DOI: 10.1016/j.cell.2020.06.038, PMID: 32702314
# Paper reports a 57-gene PIG WGCNA module. The list below mirrors the
# 52-gene working set used in `notebooks/00/07_random_from_chen_et_al.ipynb`.
# Confirm against the paper's Supp Table S2 if a final manuscript-ready
# list is needed.
PIG_CHEN2020 = [
    "Apoe", "Arpc1b", "Axl", "B2m", "C1qa", "C1qb", "C1qc", "C4b",
    "Cd63", "Cd9", "Clu", "Csf1r", "Cst3", "Ctsa", "Ctsb", "Ctsd",
    "Ctsh", "Ctsl", "Ctss", "Ctsz", "Cx3cr1", "Cyba", "Fcer1g", "Fcgr3",
    "Fcrls", "Gfap", "Gns", "Grn", "Gusb", "H2-D1", "H2-K1", "Hexa",
    "Hexb", "Igfbp5", "Itgb5", "Itm2b", "Laptm5", "Lgals3bp", "Lgmn", "Ly86",
    "Lyz2", "Man2b1", "Mpeg1", "Npc2", "Olfml3", "Plek", "Prdx6",
    "S100a6", "Serpina3n", "Trem2", "Tyrobp", "Vsir",
]

# Disease-associated microglia (Keren-Shaul et al., Cell 2017) — top markers.
# Useful for cross-checking or DAM-vs-PIG comparisons.
DAM_KEREN_SHAUL2017 = [
    "Cst7", "Itgax", "Apoe", "Trem2", "Tyrobp", "Lpl", "Csf1", "Spp1",
    "Clec7a", "Ccl3", "Ctsb", "Ctsd", "Cd68", "Tlr2", "Gpnmb", "Fabp5",
    "Hif1a", "Lilrb4a", "Cd9", "Cd63",
]

# Activated response microglia (Sala Frigerio et al., Cell Rep 2019) — APP_NLGF
ARM_SALA_FRIGERIO2019 = [
    "Cst7", "Cstb", "Apoe", "Trem2", "Tyrobp", "Lpl", "Csf1", "Spp1",
    "Clec7a", "Cd9", "Cd63", "Itgax", "Axl", "Ank",
]

# Interferon-response microglia (Sala Frigerio 2019)
IRM_SALA_FRIGERIO2019 = [
    "Ifit1", "Ifit2", "Ifit3", "Isg15", "Irf7", "Stat1", "Oasl1",
    "Mx1", "Rsad2", "Usp18",
]

# Generic homeostatic microglia markers
MICROGLIA_HOMEOSTATIC = [
    "Tmem119", "P2ry12", "Cx3cr1", "Csf1r", "Hexb", "Sall1", "Mertk",
    "Selplg", "Siglech",
]

# Endocytosis / endosomal-lysosomal pathway
ENDOCYTOSIS = [
    "Rab5a", "Rab5b", "Rab5c", "Rab7", "Rab7a", "Eea1", "Lamp1",
    "Lamp2", "Cd63", "Ctsd", "Ctsb", "Hexb", "Tfeb",
]

# ----- Reactive astrocytes (Liddelow et al, Nature 2017) --------------------
# A1 = neurotoxic / pan-reactive inflammatory; A2 = neuroprotective
A1_ASTROCYTE_LIDDELOW2017 = [
    "Gbp2", "Serping1", "Iigp1", "Srgn", "H2-T23", "Ggta1", "Ligp1",
    "Gbp4", "Fbln5", "Ugt1a1", "Fkbp5", "Psmb8", "H2-D1", "Cxcl10",
]
A2_ASTROCYTE_LIDDELOW2017 = [
    "Clcf1", "Tgm1", "Ptx3", "S100a10", "Sphk1", "Cd109", "Ptgs2",
    "Emp1", "Slc10a6", "Tm4sf1", "B3gnt5", "Cd14",
]
PAN_REACTIVE_ASTROCYTE = [
    "Lcn2", "Steap4", "S1pr3", "Timp1", "Hspb1", "Cxcl10", "Cd44",
    "Osmr", "Cp", "Serpina3n", "Aspg", "Vim", "Gfap",
]

# ----- Oligodendrocyte programs --------------------------------------------
OLIGODENDROCYTE_MATURE = [
    "Plp1", "Mbp", "Mobp", "Cnp", "Cldn11", "Mal", "Apod", "Trf",
    "Fth1", "Plekhb1", "Ppp1r14a", "Ttyh2", "Fa2h", "Aspa", "Mog", "Mag",
]
OLIGODENDROCYTE_PRECURSOR = [
    "Pdgfra", "Cspg4", "Sox10", "Olig1", "Olig2", "Nkx2-2", "C1ql1",
]
# Stress / disease-associated oligodendrocytes (Falcao 2018; Pandey 2022)
OLIGODENDROCYTE_STRESS = [
    "Serpina3n", "C4b", "B2m", "Klk6", "Il33", "Mhc1", "H2-D1", "H2-K1",
]

# ----- Vascular / barrier ---------------------------------------------------
TIGHT_JUNCTION = [
    "Cldn1", "Cldn3", "Cldn5", "Cldn11", "Cldn12", "Ocln", "Tjp1",
    "Tjp2", "Tjp3", "Jam2", "Jam3", "F11r",
]
ENDOTHELIAL = [
    "Pecam1", "Cdh5", "Cldn5", "Vwf", "Flt1", "Kdr", "Tek", "Cd34",
    "Esam", "Slc2a1", "Mfsd2a",
]
PERICYTE = [
    "Pdgfrb", "Cspg4", "Notch3", "Acta2", "Rgs5", "Mcam", "Kcnj8",
    "Abcc9", "Anpep",
]

# ----- Neuronal & synaptic --------------------------------------------------
NEURONAL_MARKERS = [
    "Rbfox3", "Tubb3", "Syn1", "Snap25", "Syp", "Camk2a", "Slc17a7",
    "Gad1", "Gad2", "Vglut1", "Map2",
]
SYNAPSE = [
    "Snap25", "Syp", "Syn1", "Syn2", "Syn3", "Sv2a", "Sv2b", "Bsn",
    "Pclo", "Dlg4", "Gphn", "Nlgn1", "Nlgn3",
]
PLASTICITY_IEG = [
    "Arc", "Fos", "Egr1", "Egr2", "Junb", "Npas4", "Bdnf", "Homer1",
    "Nr4a1", "Nptx1", "Nptx2",
]

# ----- Stress / chaperones (BRICHOS itself is a chaperone) -----------------
HEAT_SHOCK = [
    "Hspa1a", "Hspa1b", "Hsph1", "Bag3", "Dnajb1", "Dnajb4", "Dnaja1",
    "Hspb1", "Hspa8", "Hsp90aa1", "Hsp90ab1",
]
UNFOLDED_PROTEIN_RESPONSE = [
    "Atf4", "Atf6", "Hspa5", "Ddit3", "Xbp1", "Ern1", "Eif2ak3",
    "Pdia3", "Calr",
]

# ----- Microglia neurodegenerative (Krasemann et al, Immunity 2017) --------
MGND_KRASEMANN2017 = [
    "Apoe", "Itgax", "Axl", "Lilrb4a", "Clec7a", "Csf1", "Spp1",
    "Lpl", "Cst7", "Igf1", "Ccl6", "Cd9", "Tyrobp", "Ccl2", "Tgfb1",
]


# Expected direction in disease (PBS vs WT). Used by the rescue
# heatmap so signatures with no overlap to the data-driven disease
# set still get scored against a biological prior.
#   "up"   : gene set is elevated in disease (rescue = BRICHOS pulls down)
#   "down" : gene set is reduced in disease  (rescue = BRICHOS pulls up)
SIGNATURE_DIRECTION = {
    # disease-up programs
    "PIG":            "up",
    "DAM":            "up",
    "ARM":            "up",
    "IRM":            "up",
    "MGnD":           "up",
    "A1_astro":       "up",
    "panreactive_astro": "up",
    "oligo_stress":   "up",
    "UPR":            "up",
    "heat_shock":     "up",
    # homeostatic / neuronal programs that go DOWN in disease
    "microglia_homeo": "down",
    "A2_astro":       "down",
    "oligo_mature":   "down",
    "oligo_OPC":      "down",
    "tight_junction": "down",
    "endothelial":    "down",
    "pericyte":       "down",
    "neuronal":       "down",
    "synapse":        "down",
    "plasticity_IEG": "down",
    "endocytosis":    "down",
}


# Convenience registry — keeps the master notebook compact
SIGNATURE_REGISTRY = {
    "PIG":            PIG_CHEN2020,
    "DAM":            DAM_KEREN_SHAUL2017,
    "ARM":            ARM_SALA_FRIGERIO2019,
    "IRM":            IRM_SALA_FRIGERIO2019,
    "MGnD":           MGND_KRASEMANN2017,
    "microglia_homeo": MICROGLIA_HOMEOSTATIC,
    "A1_astro":       A1_ASTROCYTE_LIDDELOW2017,
    "A2_astro":       A2_ASTROCYTE_LIDDELOW2017,
    "panreactive_astro": PAN_REACTIVE_ASTROCYTE,
    "oligo_mature":   OLIGODENDROCYTE_MATURE,
    "oligo_OPC":      OLIGODENDROCYTE_PRECURSOR,
    "oligo_stress":   OLIGODENDROCYTE_STRESS,
    "tight_junction": TIGHT_JUNCTION,
    "endothelial":    ENDOTHELIAL,
    "pericyte":       PERICYTE,
    "neuronal":       NEURONAL_MARKERS,
    "synapse":        SYNAPSE,
    "plasticity_IEG": PLASTICITY_IEG,
    "endocytosis":    ENDOCYTOSIS,
    "heat_shock":     HEAT_SHOCK,
    "UPR":            UNFOLDED_PROTEIN_RESPONSE,
}


def filter_to_var(genes, var_names) -> list[str]:
    """Drop genes not present in the AnnData var_names index."""
    s = set(var_names)
    return [g for g in genes if g in s]
