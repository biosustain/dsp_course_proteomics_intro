# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: tags,-all
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Data Analysis PXD040621
#
# Plan
# - read data and log2 transform intensity values
# - aggregate peptide intensities to protein intensities
# - format data from long to wide format
# - remove contaminant proteins
# - check for missing values
# - Clustermap of sample and proteins
# - differential analysis (Volcano Plots)
# - Enrichment Analysis
# - check for maltose update pathway (Fig. 3 in paper)

# %%
from pathlib import Path

import acore.differential_regulation
import acore.enrichment_analysis
import acore.normalization
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import seaborn as sns
import vuecore
from acore.io.uniprot import fetch_annotations, process_annotations
from vuecore.viz import get_enrichment_plots

# %% [markdown]
# ## Read in the data
# - `file_in`: input file with the quantified peptide data in MSstats format as provided by quantms
#
# The file can be downloaded from [Google Drive](https://drive.google.com/drive/folders/1Nm5Ha-tCvjU-B323BLhna1GwHdNpK_lU?usp=drive_link)

# %% tags=["parameters"]
file_in = "data/PXD040621/processed/PXD040621.sdrf_openms_design_msstats_in.csv"
df = pd.read_csv(file_in, sep=",", header=0)  # .set_index([])
df.head()


# %% [markdown]
# define the output folder for our VueGen report which we will create later

# %%
out_dir = "data/PXD040621/report/"
out_dir = Path(out_dir)
out_dir.mkdir(parents=True, exist_ok=True)

# %% [markdown]
# We have the following columns in the data:
#
# ```python
# cols = [
#     "ProteinName",
#     "PeptideSequence",
#     "PrecursorCharge",
#     "FragmentIon",
#     "ProductCharge",
#     "IsotopeLabelType",
#     "Condition",
#     "BioReplicate",
#     "Run",
#     "Intensity",
#     "Reference",
# ]

# %% [markdown]
# ## Log2 transform the intensity values
# - log2 transformations are common for lognormal distributed data

# %%
df["Intensity"] = np.log2(df["Intensity"].astype(float))

# %% [markdown]
# ## Aggregate the peptide intensities to protein intensities
# - we use the median of the peptide intensities for each protein
#
# There are more sophisticated ways to do this, e.g. using MaxLFQ, iBAQ, FlashLFQ, DirectLFQ, etc.

# %%
proteins = (
    df.groupby(["ProteinName", "Reference"])["Intensity"].median().unstack(level=0)
)
proteins

# %% [markdown]
# ## Remove contaminant proteins
# Remove the contaminant proteins which were added to the fasta file used in the data processing.
# Contaminant proteins are e.g. creation which gets into the sample from the human skin or hair
# when the sample is prepared.
#
# These are filtered out as they are most of the time not relevant, but a contamination.

# %%
decoy_proteins = proteins.filter(like="CON_", axis=1)
proteins = proteins.drop(decoy_proteins.columns, axis=1)
proteins

# %% [markdown]
# Create a label for each sample based on the metadata.
# - we will use a string in the sample name, but you can see how the metadata is organized

# %%
meta = df[["Condition", "BioReplicate", "Run", "Reference"]].drop_duplicates()
meta

# %%
label_encoding = {0: "control", 1: "10 µm sulforaphane"}
label_suf = pd.Series(
    proteins.index.str.contains("_Suf_").astype(int),
    index=proteins.index,
    name="label_suf",
    dtype=np.int8,
).map(label_encoding)
label_suf

# %% [markdown]
# ## Plot the data completeness for each protein.

# %%
view_name = "Protein"
out_dir_subsection = out_dir / "1_data" / "completeness"
out_dir_subsection.mkdir(parents=True, exist_ok=True)

# %%
view_name = "Protein"
ax = (
    proteins.notna()
    .sum()
    .sort_values()
    .plot(
        rot=45,
        ylabel=f"Number of Samples {view_name.lower()} was observed in",
    )
)
ax.get_figure().savefig(
    out_dir_subsection / f"data_completeness_step_plot.png",
    bbox_inches="tight",
    dpi=300,
)

# %%
view_name = "Protein"
ax = (
    proteins.notna()
    .sum()
    .value_counts()
    .sort_index(ascending=False)
    .plot(
        kind="bar",
        title=f"Data Completeness per {view_name}",
        xlabel=f"Number of Samples {view_name.lower()} was observed in",
        ylabel=f"Number of {view_name}s",
        color="steelblue",
        figsize=(10, 6),
    )
)
ax.get_figure().savefig(
    out_dir_subsection / f"data_completeness_bar_plot.png",
    bbox_inches="tight",
    dpi=300,
)

# %%
# Explode column names to examine split by '|'
proteins_meta = (
    proteins.columns.str.split("|", expand=True)
    .to_frame()
    .dropna(how="any", axis=1)
    .reset_index(drop=True)
)
proteins_meta.columns = ["Source", "ProteinName", "GeneName"]
proteins_meta["GeneName"] = proteins_meta["GeneName"].str.split("_").str[0]
proteins_meta.index = proteins.columns
proteins_meta.index.name = "identifier"
proteins_meta


# %% [markdown]
# For later convenience let's replace the identifier with the UNIPROT ID

# %%
proteins.columns = proteins_meta["ProteinName"]
proteins

# %% [markdown]
# And let's save a table with the data for inspection

# %%
proteins_meta.to_csv(out_dir_subsection / "proteins_identifiers.csv")
proteins.to_csv(out_dir_subsection / "proteins.csv")


# %% [markdown]
# ## Hierarchical Clustering of data
# - using completely observed data only
# Find correlations in data

# %%
out_dir_subsection = out_dir / "1_data" / "clustermap"

# %%
_group_labels = label_encoding.values()
lut = dict(zip(_group_labels, [f"C{i}" for i in range(len(_group_labels))]))
row_colors = label_suf.map(lut).rename("group color")
row_colors

# %%
vuecore.set_font_sizes(7)
cg = sns.clustermap(
    proteins.dropna(how="any", axis=1),
    method="ward",
    row_colors=row_colors,
    figsize=(11, 6),
    robust=True,
    xticklabels=True,
    yticklabels=True,
)
fig = cg.figure
cg.ax_heatmap.set_xlabel("Proteins")
cg.ax_heatmap.set_ylabel("Sample ID")
vuecore.select_xticks(cg.ax_heatmap)
handles = [
    plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=lut[name], markersize=8)
    for name in lut
]
cg.ax_cbar.legend(
    handles, _group_labels, title="Groups", loc="lower left", bbox_to_anchor=(2, 0.5)
)
fname = out_dir_subsection / "clustermap_ward.png"
# vuecore.savefig(fig, fname, pdf=True, dpi=600, tight_layout=False)
fig.savefig(
    out_dir_subsection / "clustermap_ward.png",
    bbox_inches="tight",
    dpi=300,
)

# %% [markdown]
# ### Hierarchical Clustering of normalized data
# - using completely observed data only
# Checkout the [recipe on normalization methods](https://analytics-core.readthedocs.io/latest/api_examples/normalization_analysis.html).

# %%
normalization_method = "median"
X = acore.normalization.normalize_data(
    proteins.dropna(how="any", axis=1), normalization_method
)
X

# %%
X.median(axis="columns")

# %%
vuecore.set_font_sizes(7)
cg = sns.clustermap(
    X,
    method="ward",
    row_colors=row_colors,
    figsize=(11, 6),
    robust=True,
    xticklabels=True,
    yticklabels=True,
)
fig = cg.figure
cg.ax_heatmap.set_xlabel("Proteins")
cg.ax_heatmap.set_ylabel("Sample ID")
vuecore.select_xticks(cg.ax_heatmap)
handles = [
    plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=lut[name], markersize=8)
    for name in lut
]
cg.ax_cbar.legend(
    handles, _group_labels, title="Groups", loc="lower left", bbox_to_anchor=(2, 0.5)
)
fname = out_dir_subsection / "clustermap_ward.png"
# vuecore.savefig(fig, fname, pdf=True, dpi=600, tight_layout=False)
fig.savefig(
    out_dir_subsection / f"clustermap_ward_{normalization_method}.png",
    bbox_inches="tight",
    dpi=300,
)

# %% [markdown]
# ## Differential Regulation

# %%
out_dir_subsection = out_dir / "2_differential_regulation"
out_dir_subsection.mkdir(parents=True, exist_ok=True)

# %%
view = proteins
group = "label_suf"
diff_reg = acore.differential_regulation.run_anova(
    view.dropna(how="any", axis=1).join(label_suf),
    alpha=0.15,
    drop_cols=[],
    subject=None,
    group=group,
).sort_values("pvalue", ascending=True)
diff_reg["rejected"] = diff_reg["rejected"].astype(bool)
diff_reg.sort_values("pvalue")

# %%
diff_reg.plot(x="log2FC", y="-log10 pvalue", kind="scatter", title=group)

# %% [markdown]
# ## Interactive Volcano Plot

# %%
str_cols = diff_reg.dtypes[diff_reg.dtypes == "object"].index.tolist()
hover_data = {
    "rejected": ":.0f",
    **{
        c: ":.4f"
        for c in [
            "padj",
            "FC",
        ]
    },
    **{c: True for c in str_cols},
}
fig = px.scatter(
    diff_reg,
    x="log2FC",
    y="-log10 pvalue",
    color="rejected",
    hover_data=hover_data,
    width=1200,
    height=800,
    title=f"Volcano plot for {view_name}s",
)
fig

# %% [markdown]
# Save result to subsection folder

# %%
fig.write_json(
    out_dir_subsection / "0_volcano_plot.json",
    pretty=False,
)
diff_reg.to_csv(out_dir_subsection / "1_differential_regulation.csv")

# %% [markdown]
# # Enrichment Analysis

# %%
out_dir_subsection = out_dir / "uniprot_annotations"
out_dir_subsection.mkdir(parents=True, exist_ok=True)

# %% [markdown]
# Fetch the annotations from UniProt API.

# %%
fname_annotations = out_dir_subsection / "annotations.csv"
try:
    annotations = pd.read_csv(fname_annotations, index_col=0)
    print(f"Loaded annotations from {fname_annotations}")
except FileNotFoundError:
    print(f"Fetching annotations for {proteins.columns.size} UniProt IDs.")
    FIELDS = "go_p,go_c,go_f"
    annotations = fetch_annotations(proteins.columns, fields=FIELDS)
    annotations = process_annotations(annotations, fields=FIELDS)
    # cache the annotations
    fname_annotations.parent.mkdir(exist_ok=True, parents=True)
    annotations.to_csv(fname_annotations, index=True)

annotations

# %%
enriched = acore.enrichment_analysis.run_up_down_regulation_enrichment(
    regulation_data=diff_reg,
    annotation=annotations,
    min_detected_in_set=2,
    lfc_cutoff=1,
    pval_col="pvalue",
    correction_alpha=0.05,  # adjust the p-value to see more or less results
)
enriched

# %%
from vuecore.viz import get_enrichment_plots

fig = get_enrichment_plots(
    enriched,
    identifier="anything",  # ToDo: figure out what this does
    args=dict(title="Enrichment Analysis"),
)
fig = fig[0]
fig.write_json(
    out_dir_subsection / "enrichment_analysis.json",
    pretty=True,
)
fig


# %% [markdown]
# ## Check for Maltose Uptake

# %%
out_dir_subsection = out_dir / "3_maltose_uptake"
out_dir_subsection.mkdir(parents=True, exist_ok=True)

# %% [markdown]
# apply filtering of 'differentially abundant proteins' as described in the paper
# > Differentially abundant proteins were determined as those with log2 fold-change
# > > 1 and < -1, and p < 0.05
# This means not multiple testing correction was applied.

# %%
view = diff_reg.query("pvalue < 0.05 and FC > 1")  # .shape[0]
view.to_csv(
    out_dir_subsection / "1_differently_regulated_as_in_paper.csv",
    index=False,
)
view

# %% [markdown]
# Let's find the proteins highlighted in the volcano plot in Figure 3.

# %%
highlighted_genes = ["LamB", "MalE", "Malk", "CitF", "CitT", "CitE", "Frd"]
highlighted_genes = "|".join([p.upper() for p in highlighted_genes])
highlighted_genes = proteins_meta.query(
    f"`GeneName`.str.contains('{highlighted_genes}')"
)
highlighted_genes

# %%
highlighted_proteins = "|".join([p.upper() for p in highlighted_genes["ProteinName"]])
view = diff_reg.query(f"`identifier`.str.contains('{highlighted_proteins}')")
view.to_csv(
    out_dir_subsection / "2_highlighted_proteins_in_figure3.csv",
    index=False,
)
view

# %% [markdown]
# How to explain the differences?
