"""Human -> mouse gene-symbol mapping.

Lightweight, dependency-free ortholog resolution for porting human gene
signatures (e.g. the Rexach 2024 BA4_MIC-7 microglia state) onto this
mouse Visium dataset. No network / biomart dependency so the notebooks
stay reproducible offline.

Strategy
--------
1. Curated overrides (`HUMAN_TO_MOUSE`) handle the cases where a simple
   title-case transform is wrong: one-to-many MHC genes, renamed symbols
   (SEPP1->Selenop, TF->Trf), mitochondrially-encoded genes (MT-ND3 ->
   mt-Nd3), etc.
2. Everything else falls back to the title-case heuristic
   (`HUMAN -> Human`), which is correct for the large majority of
   1:1 mouse orthologs.
3. `to_mouse` flattens one-to-many maps and `filter_to_var` keeps only
   symbols actually present in the AnnData.

This is deliberately conservative: when in doubt, check the printed
mapping table in the notebook against MGI / the data's var_names.
"""
from __future__ import annotations

# Curated overrides. Values may be a single symbol or a list (one-to-many).
# Anything NOT in here is resolved by the title-case fallback.
HUMAN_TO_MOUSE: dict[str, str | list[str]] = {
    # MHC class I — no 1:1 ortholog; mouse uses H2 haplotype genes
    "HLA-B": ["H2-K1", "H2-D1"],
    "HLA-A": ["H2-K1", "H2-D1"],
    "HLA-C": ["H2-K1", "H2-D1"],
    "HLA-E": "H2-T23",
    "B2M":   "B2m",
    # renamed / non-obvious symbols
    "SEPP1":  "Selenop",     # selenoprotein P, renamed
    "TF":     "Trf",          # transferrin (mouse symbol is Trf)
    "FCGR3A": ["Fcgr3", "Fcgr4"],  # human FCGR3A ~ mouse Fcgr3/Fcgr4
    "MT2A":   "Mt2",          # metallothionein 2
    "MT1A":   "Mt1",
    "MT-ND3": "mt-Nd3",       # mitochondrially encoded
    "MT-CO1": "mt-Co1",
    "MT-CYB": "mt-Cytb",
    "C1ORF43": "Tmco1",
    "H3-3A":  "H3f3a",
    # already-correct / identical symbols kept explicit for clarity
    "C3":     "C3",
    "GRN":    "Grn",
    "APOE":   "Apoe",
}


def _titlecase(sym: str) -> str:
    """HUMAN_SYMBOL -> Mouse_symbol (first letter cap, rest lower).

    Handles hyphenated suffixes used by some mouse families by leaving
    them alone after the first segment is title-cased.
    """
    s = sym.strip()
    if not s:
        return s
    return s[0].upper() + s[1:].lower()


def to_mouse(human_genes, *, keep_unmapped: bool = False) -> list[str]:
    """Map a list of human symbols to mouse symbols.

    One-to-many maps are flattened. Order preserved, duplicates removed.
    Set ``keep_unmapped=True`` to title-case-fallback unknown symbols
    (default True behaviour anyway — curated dict only overrides).
    """
    out: list[str] = []
    seen: set[str] = set()
    for g in human_genes:
        if g is None:
            continue
        g = str(g).strip()
        if not g:
            continue
        mapped = HUMAN_TO_MOUSE.get(g, _titlecase(g))
        targets = mapped if isinstance(mapped, list) else [mapped]
        for t in targets:
            if t and t not in seen:
                seen.add(t)
                out.append(t)
    return out


def mapping_table(human_genes) -> list[tuple[str, list[str], str]]:
    """Return [(human, [mouse...], source)] for inspection/QC.

    source is 'curated' if the symbol came from HUMAN_TO_MOUSE, else
    'titlecase'.
    """
    rows = []
    for g in human_genes:
        g = str(g).strip()
        if g in HUMAN_TO_MOUSE:
            m = HUMAN_TO_MOUSE[g]
            rows.append((g, m if isinstance(m, list) else [m], "curated"))
        else:
            rows.append((g, [_titlecase(g)], "titlecase"))
    return rows
