import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def volcano_plot_region_within_group(adata, region, logfc_thresh=1, padj_thresh=0.05, top_n=10):
    de = sc.get.rank_genes_groups_df(adata, group=region)
    de = de.replace([np.inf, -np.inf], np.nan).dropna(subset=['logfoldchanges', 'pvals_adj'])
    de['-log10(padj)'] = -np.log10(de['pvals_adj'])

    plt.figure(figsize=(10, 6))
    plt.scatter(de['logfoldchanges'], de['-log10(padj)'], c='lightgrey', alpha=0.7)

    sig = (de['pvals_adj'] < padj_thresh) & (np.abs(de['logfoldchanges']) > logfc_thresh)
    plt.scatter(
        de.loc[sig, 'logfoldchanges'],
        de.loc[sig, '-log10(padj)'],
        c='crimson'
    )

    top_genes = de.loc[sig].sort_values('-log10(padj)', ascending=False).head(top_n)
    for _, row in top_genes.iterrows():
        plt.text(
            row['logfoldchanges'],
            row['-log10(padj)'],
            row['names'],
            fontsize=8,
            ha='right' if row['logfoldchanges'] < 0 else 'left'
        )

    plt.axhline(-np.log10(padj_thresh), linestyle='--', color='gray')
    plt.axvline(-logfc_thresh, linestyle='--', color='gray')
    plt.axvline(logfc_thresh, linestyle='--', color='gray')
    plt.xlabel("log2 fold change (region vs rest)")
    plt.ylabel("-log10 adjusted p-value")
    plt.title(f"Volcano plot — Region: {region} (within {adata.obs['treatment'][0]})")
    plt.tight_layout()
    plt.show()

def volcano_plot_br(adata, groupby="treatment", group1="BRICHOS", group2="PBS", top_n=10):
    """
    Generate a volcano plot comparing group1 vs group2 in a Scanpy AnnData object.

    Parameters:
    - adata: AnnData object
    - groupby: column in adata.obs to group by (default: 'treatment')
    - group1: group to compare (e.g. 'BRICHOS')
    - group2: reference group (e.g. 'PBS')
    - top_n: number of top DE genes to label
    """
    # Subset to relevant groups
    ad_sub = adata[adata.obs[groupby].isin([group1, group2])].copy()

    # Run DE
    sc.tl.rank_genes_groups(ad_sub, groupby=groupby, method="t-test_overestim_var", reference=group2)
    de_genes = sc.get.rank_genes_groups_df(ad_sub, group=group1)

    # Clean and compute -log10(padj)
    de_genes = de_genes.replace([np.inf, -np.inf], np.nan).dropna(subset=['logfoldchanges', 'pvals_adj'])
    de_genes['-log10(padj)'] = -np.log10(de_genes['pvals_adj'])

    # Plot
    plt.figure(figsize=(10, 7))
    plt.scatter(
        de_genes['logfoldchanges'],
        de_genes['-log10(padj)'],
        c='gray',
        alpha=0.7
    )

    # Highlight significant genes
    sig = (de_genes['pvals_adj'] < 0.05) & (np.abs(de_genes['logfoldchanges']) > 1)
    plt.scatter(
        de_genes.loc[sig, 'logfoldchanges'],
        de_genes.loc[sig, '-log10(padj)'],
        c='crimson',
        label='FDR<0.05 & |log2FC|>1'
    )

    # Label top DE genes
    top_genes = de_genes.loc[sig].sort_values('-log10(padj)', ascending=False).head(top_n)
    for _, row in top_genes.iterrows():
        if np.isfinite(row['logfoldchanges']) and np.isfinite(row['-log10(padj)']):
            plt.text(
                row['logfoldchanges'],
                row['-log10(padj)'],
                row['names'],
                fontsize=9,
                ha='right' if row['logfoldchanges'] < 0 else 'left'
            )

    # Decorations
    plt.axhline(-np.log10(0.05), linestyle='--', color='gray')
    plt.axvline(-1, linestyle='--', color='gray')
    plt.axvline(1, linestyle='--', color='gray')
    plt.xlabel(f'log2 fold change ({group1} vs {group2})')
    plt.ylabel('-log10 adjusted p-value')
    plt.title(f'Volcano plot: {group1} vs {group2}')
    plt.legend()
    plt.tight_layout()
    plt.show()

import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np

def plot_spatial_clusters_per_sample(
    adata, 
    color='leiden', 
    spot_size=500, 
    figsize=(20, 10), 
    layout='grid'  # 'single' or 'grid'
):
    """
    Plots spatial clusters for each sample in a multi-sample AnnData object.

    Parameters:
    - adata: AnnData object
    - color: column in .obs to color by (default: 'leiden')
    - spot_size: size of the spatial spots
    - figsize: tuple indicating figure size
    - layout: 'single' to plot each sample one by one, 'grid' to show side-by-side in a single figure
    """
    unique_samples = adata.obs['sample_id'].unique()
    n_samples = len(unique_samples)

    if layout == 'single':
        # Plot each sample individually
        for sample_id in unique_samples:
            print(f"📍 Plotting sample: {sample_id}")
            ad_int = adata[adata.obs['sample_id'] == sample_id].copy()

            # Patch spatial
            if "spatial_metadata_per_sample" in adata.uns and sample_id in adata.uns["spatial_metadata_per_sample"]:
                ad_int.uns["spatial"] = {sample_id: adata.uns["spatial_metadata_per_sample"][sample_id]}
            else:
                print(f"⚠️  No spatial metadata found for {sample_id}")
                continue

            with plt.rc_context({'figure.figsize': figsize}):
                sc.pl.spatial(ad_int, spot_size=spot_size, color=color, title=sample_id)
                plt.show()

    elif layout == 'grid':
        # Create subplots in a grid
        cols = min(3, n_samples)
        rows = int(np.ceil(n_samples / cols))
        fig, axs = plt.subplots(rows, cols, figsize=figsize, squeeze=False)

        for i, sample_id in enumerate(unique_samples):
            row, col = divmod(i, cols)
            ad_int = adata[adata.obs['sample_id'] == sample_id].copy()

            if "spatial_metadata_per_sample" in adata.uns and sample_id in adata.uns["spatial_metadata_per_sample"]:
                ad_int.uns["spatial"] = {sample_id: adata.uns["spatial_metadata_per_sample"][sample_id]}
            else:
                print(f"⚠️  No spatial metadata found for {sample_id}")
                continue

            ax = axs[row][col]
            sc.pl.spatial(
                ad_int,
                spot_size=spot_size,
                color=color,
                title=sample_id,
                show=False,
                ax=ax
            )

        # Hide any unused axes
        for j in range(i + 1, rows * cols):
            r, c = divmod(j, cols)
            fig.delaxes(axs[r][c])

        plt.tight_layout()
        plt.show()

def plot_dotplot_by_treatment(
    adata,
    genes,
    groupby='treatment',
    standard_scale='var',  # optional: standardizes each gene across groups
    title='Dotplot of selected genes by treatment'
):
    """
    Plots a dotplot of selected genes grouped by treatment conditions.

    Parameters:
    - adata: AnnData object
    - genes: list of gene names
    - groupby: column in adata.obs to group cells (e.g. "treatment")
    - standard_scale: None, 'var' (genes) or 'group' — whether to standardize
    - title: plot title
    """
    sc.pl.dotplot(
        adata,
        var_names=genes,
        groupby=groupby,
        standard_scale=standard_scale,
        show=False
    )
    plt.title(title)
    plt.tight_layout()
    plt.show()

def plot_relative_cluster_composition(adata, cluster_key="leiden", groupby="treatment"):
    """
    Plot the relative cluster composition across different treatment groups using scanpy-generated cluster colors.

    Parameters:
    - adata: AnnData object
    - cluster_key: str, column in .obs indicating cluster labels (default: "leiden")
    - groupby: str, column in .obs indicating group (e.g. treatment, sample) (default: "treatment")
    """
    import pandas as pd
    import matplotlib.pyplot as plt

    # Compute relative distribution
    cluster_distribution = pd.crosstab(adata.obs[groupby], adata.obs[cluster_key])
    relative_distribution = cluster_distribution.div(cluster_distribution.sum(axis=1), axis=0)

    # Get scanpy-assigned cluster colors
    cluster_colors = adata.uns.get(f"{cluster_key}_colors")
    cluster_labels = cluster_distribution.columns.astype(str)

    if cluster_colors is not None and len(cluster_colors) >= len(cluster_labels):
        color_map = dict(zip(cluster_labels, cluster_colors))
        colors = [color_map[label] for label in cluster_labels]
    else:
        colors = None  # fallback to default colormap

    # Plot
    ax = relative_distribution[cluster_labels].plot(
        kind="bar",
        stacked=True,
        figsize=(10, 6),
        color=colors,
        edgecolor="black"
    )

    plt.title(f"Relative Cluster Composition per {groupby.capitalize()}")
    plt.ylabel("Proportion of Cells")
    plt.xlabel(groupby.capitalize())
    plt.legend(title="Cluster", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.show()

def score_and_plot_modules(adata, modules, groupby="treatment", dendrogram=True):
    """
    Score gene modules and plot dotplot with optional dendrogram.

    Parameters:
    - adata: AnnData object
    - modules: dict of {module_name: [genes]}
    - groupby: str, column in .obs to group by
    - dendrogram: bool, whether to compute and show hierarchical clustering
    """

    # Score each module
    for name, genes in modules.items():
        genes_present = [g for g in genes if g in adata.var_names]
        if len(genes_present) < 3:
            print(f"⚠️ Skipping {name} (only {len(genes_present)} genes found)")
            continue
        sc.tl.score_genes(adata, gene_list=genes_present, score_name=f"{name}_score")
        print(f"✅ Scored {name} module with {len(genes_present)} genes")

    # Optional dendrogram
    if dendrogram:
        try:
            sc.tl.dendrogram(adata, groupby=groupby)
        except ValueError as e:
            print(f"⚠️ Could not compute dendrogram: {e}")
            dendrogram = False

    # Plot
    sc.pl.dotplot(
        adata,
        var_names=[f"{k}_score" for k in modules],
        groupby=groupby,
        dendrogram=dendrogram,
        standard_scale="var",
        color_map="coolwarm"
    )