# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: tags,-all
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
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

# %% tags=["hide-output"]
# %pip install acore vuecore "pingouin<0.6.0"

# %%
from pathlib import Path

import acore.differential_regulation
import acore.enrichment_analysis
import acore.normalization
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import scipy.stats
import seaborn as sns
import vuecore
from acore.io.uniprot import fetch_annotations, process_annotations
from vuecore.viz import get_enrichment_plots

# %% [markdown]
# # Paramters
# - `file_in`: input file with the quantified peptide data in MSstats format
#    as provided by quantms
# - `out_dir`: output directory for the results of the data analysis, 
#    which will be used later for the report generation with VueGen.
#
# The file will be loaded from the online repository if it is not present.

# %% tags=["parameters"]
file_in: str = Path(
    "data/PXD040621/processed/PXD040621.sdrf_openms_design_msstats_in.csv"
) # input file with the quantified peptide data in MSstats format as provided by quantms
out_dir = "data/PXD040621/report/" # output directory for the results of the data analysis, which will be used later for the report generation with VueGen.
min_obs_per_group: int = 3 # minimum number of observations per group for a protein to be included in the differential regulation analysis

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


# %%
if not file_in.exists():
    file_in = (
        "https://raw.githubusercontent.com/biosustain/dsp_course_proteomics_intro/HEAD"
        "/data/PXD040621/processed/PXD040621.sdrf_openms_design_msstats_in.csv"
    )
df = pd.read_csv(file_in, sep=",", header=0)  # .set_index([])
df.head()


# %% [markdown]
# # VueGen report output folder
# define the output folder for our `VueGen` report which we will create later

# %%
out_dir = "data/PXD040621/report/"
out_dir = Path(out_dir)
out_dir.mkdir(parents=True, exist_ok=True)


# %% [markdown]
# # Log2 transform the intensity values
# - log2 transformations are common for lognormal distributed data

# %%
df["Intensity"] = np.log2(df["Intensity"].astype(float))
df.head()

# %% [markdown]
# # Exploratory and Data Quality Plots (peptide level)
# Shows the overal distribution of the log2-normal intensities per sample

# %%
df["BioReplicate"] = df["BioReplicate"].replace({5: 1, 6: 2, 7: 3, 8: 4})
fg = sns.displot(
    data=df.rename(columns={"BioReplicate": "Rep", "Condition": "C."}),
    x="Intensity",
    col="C.",
    row="Rep",
    # hue="Reactor_ID",
    kind="kde",
    height=2,
    aspect=1.1,
)

# %% [markdown]
# # Aggregate the peptide intensities to protein intensities
# - we use the median of the peptide intensities for each protein
#
# There are more sophisticated ways to do this, e.g. using MaxLFQ, iBAQ, FlashLFQ, 
# DirectLFQ, etc.
#
# - shorten sample name for readability

# %%
proteins = (
    df.groupby(["ProteinName", "Reference"])["Intensity"].median().unstack(level=0)
)
proteins.index = proteins.index.str.split("_").str[4:6].str.join("_")
proteins

# %% [markdown]
# # Remove contaminant proteins
# Remove the contaminant proteins which were added to the fasta file used in the data
# processing. Contaminant proteins are e.g. creation which gets into the sample from the
# human skin or hair when the sample is prepared.
#
# These are filtered out as they are most of the time not relevant, but a contamination.

# %%
decoy_proteins = proteins.filter(like="CON_", axis=1)
proteins = proteins.drop(decoy_proteins.columns, axis=1)
proteins

# %% [markdown]
# # Factor variable
# Create a label for each sample based on the metadata.
# - we will use a string in the sample name, but you can see how the metadata is organized
#
# > This could be provided in a separate metadata file.

# %%
meta = df[["Condition", "BioReplicate", "Run", "Reference"]].drop_duplicates()
meta

# %%
label_encoding = {0: "control", 1: "10 µm sulforaphane"}
label_suf = pd.Series(
    proteins.index.str.contains("Suf_").astype(int),
    index=proteins.index,
    name="label_suf",
    dtype=np.int8,
).map(label_encoding)
label_suf

# %% [markdown]
# # Plot the data completeness for each protein.
# - see the how many proteins are observed across how many samples.
#
# > We create a subfolder in our results folder to keep it organized for the 
#   report generation later.

# %% tags=["hide-input"]
view_name = "Protein"
out_dir_subsection = out_dir / "1_data" / "completeness"
out_dir_subsection.mkdir(parents=True, exist_ok=True)


# %% [markdown]
# From the line plot (in staircase shape) you can see that over half of the proteins
# are observed in all 8 samples.

# %% tags=["hide-input"]
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

# %% [markdown]
# The bar plot shows the same information, but it is easier to see the exact number of
# proteins observed in how many samples.

# %% tags=["hide-input"]
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

# %% [markdown]
# # Protein identifiers metadata
# Keep the protein identifiers metadata. Will be used to query `UNIProt`.

# %% tags=["hide-input"]
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
# For later in the enrichment analysis let's replace the protein identifier from the Fasta
# file with the UNIPROT ID

# %%
proteins.columns = proteins_meta["ProteinName"].rename("UniprotID")
proteins

# %% [markdown]
# And let's save a table with the data for inspection

# %%
proteins_meta.to_csv(out_dir_subsection / "proteins_identifiers.csv")
proteins.to_csv(out_dir_subsection / "proteins.csv")

# %% [markdown]
# # Hierarchical Clustering of data
# - using completely observed data only find correlations in data

# %%
out_dir_subsection = out_dir / "1_data" / "clustermap"
out_dir_subsection.mkdir(parents=True, exist_ok=True)

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
# # Analytical Plots
# - data distribution (e.g. histogram)
# - coefficient of variation (CV)
# - number of identified proteins per sample

# %%
_min_int, _max_int = proteins.min().min(), proteins.max().max()
bins = range(int(_min_int), int(_max_int) + 1, 1)
ax = proteins.T.hist(layout=(2, 4), bins=bins, sharex=True, sharey=True, figsize=(8, 4))

min_int, max_int = int(proteins.min().min()), int(proteins.max().max())
bins = range(min_int, max_int+1, 1)
ax = proteins.T.hist(layout=(2, 4), bins=bins, sharex=True, sharey=True, figsize=(8, 4))

# %%
fig, ax = plt.subplots(figsize=(8, 4))
proteins.T.boxplot(ax=ax)
ax.set_xlabel("Sample ID")
ax.set_ylabel("log2 intensity")
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
fname = out_dir_subsection / "boxplot_log2_unnormalized.png"
vuecore.savefig(fig, fname, pdf=True, dpi=600, tight_layout=False)
print(f"Saved boxplot to {fname}")

# %% [markdown]
# # Coefficient of Variation (CV)
# - CV = standard deviation / mean
#   $$ CV = \frac{\sigma}{\mu} $$
# - per group
#
# Use
# [`scipy.stats.variation`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.variation.html)
# to calculate the CV.

# %%
df_cvs = (
    proteins.groupby(label_suf)  # .join(metadata[grouping])
    # .agg(scipy.stats.variation)
    .agg([scipy.stats.variation, "mean"])  # .rename_axis(["feat", "stat"], axis=1)
)
df_cvs

# %%
df_cvs = df_cvs.stack(0, future_stack=True).reset_index().dropna()
df_cvs

# %% [markdown]
# Visualize the relationship between the coefficient of variation (CV) and
# the mean intensity per group: x-axis: CV ('variation'), y-axis: mean intensity.
#
# > Here no further filtering is needed, but many researchers filter for proteins with a
# > mean intensity above a certain threshold (aroung 0.3). This is best done on quality
# > control samples, e.g. a pooled sample, for large datasets.

# %% tags=["hide-input"]
default_args = dict(
    facet_col="label_suf",
    # facet_row="Time",
    labels={
        "label_suf": "group",
        "variation": "CV",
    },
)
fig = px.scatter(
    data_frame=df_cvs,
    x="variation",
    y="mean",
    trendline="ols",
    hover_data={"UniprotID": True, "mean": True},  # hide the index in the hover
    **default_args,
)
fname = "cv_vs_mean"
fig.write_json(
    out_dir_subsection / f"{fname}.json",
    pretty=True,
)
fig

# %% [markdown]
# # Hierarchical Clustering of normalized data
# - using completely observed data only
# Checkout the [recipe on normalization methods](https://analytics-core.readthedocs.io/latest/api_examples/normalization_analysis.html).

# %%
normalization_method = "median"
X = acore.normalization.normalize_data(
    proteins.dropna(how="any", axis=1), normalization_method
)
# X = acore.normalization.normalize_data(
#     X, "zscore"
# )  # to make it more comparable across proteins, but not necessary for the clustering
X

# %%
X.median(axis="columns")

# %% tags=["hide-input"]
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
# Exercise: 
# 1. Try out different normalization methods and see how the clustering changes.
# 2. Maybe you want to combine a sample based with a protein based normalization method.

# %% [markdown]
# # Differential Regulation

# %%
out_dir_subsection = out_dir / "2_differential_regulation"
out_dir_subsection.mkdir(parents=True, exist_ok=True)

# %% [markdown]
# Retain all proteins with at least 3 observations in each group
# - this is a requirement for a standard t-test
# - you could look into imputation methods to fill in missing values)
#   - protein in at least two samples per group?
#   - missing all in one condition?
#
# Let's not impute, but filter for proteins with at least 3 observations in each group

# %%
group_counts = proteins.groupby(label_suf).count()
group_counts

# %% [markdown]
# Let's see how many observations are present by condition (group)
# - The order is dermined by the treatment group here
# - Some proteins are absent in the control group that are observed in the treatment
#   group, which is interesting for a separate analysis (not shown).

# %% tags=["hide-input"]
frac_non_na_by_group = (
    proteins.join(label_suf, how="inner")
    .groupby(label_suf.name)
    .apply(lambda x: x.notna().sum(), include_groups=False)
    .T
)
frac_non_na_by_group = frac_non_na_by_group.sort_values(by=frac_non_na_by_group.columns[0])
fig, ax = plt.subplots(figsize=(7, 3))
n_groups = len(frac_non_na_by_group.columns)
jitter_span = 0.1
offsets = np.linspace(-jitter_span / 2, jitter_span / 2, n_groups) if n_groups > 1 else [0.0]
for i, g in enumerate(frac_non_na_by_group.columns):
    (frac_non_na_by_group[g] + offsets[i]).plot(
        kind="line",
        ax=ax,
        style=".",
        alpha=0.5,
    )
ax.legend(title="Group")
ax.set_ylabel("Ratio of non-missing values\n(completeness)")
ax.set_title("Number of observations")

# %% [markdown]
# For now, we filter the proteins to only those with at least 3 observations
# in each group.

# %%
mask = group_counts.groupby("label_suf").transform(
    lambda x: x >= min_obs_per_group
).all(axis=0)
mask

# %% [markdown]
# We use a simple t-test for the differential regulation analysis.
#
# > Experiment with the normalization method and it's effect on the results.

# %%
view = X.loc[:, mask].join(label_suf)
group = label_suf.name
diff_reg = acore.differential_regulation.run_anova(
    view,
    alpha=0.15,
    drop_cols=[],
    subject=None,
    group=group,
).sort_values("pvalue", ascending=True)
diff_reg["rejected"] = diff_reg["rejected"].astype(bool)
diff_reg.sort_values("pvalue")

# %%
diff_reg.sort_values("pvalue").head(20)

# %% [markdown]
# # Interactive Volcano Plot
# - y-axis: -log10(p-value)
# - x-axis: log2(fold-change)
# - color: rejected (significant or not after multiple testing correction)
#
# > Hover over the points to see additional information

# %% tags=["hide-input"]
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
#
# - You can reload the annotations, but for latency, I added a cached version
#   to the repository

# %%
out_dir_subsection = out_dir / "uniprot_annotations"
out_dir_subsection.mkdir(parents=True, exist_ok=True)

# %% [markdown]
# ## Fetch the annotations from UniProt API.

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

# %% [markdown]
# ## Run the enrichment analysis
# - background is the set of identified proteins in the experiment (not the whole proteome
#   of the organisim, here E. coli)
# - The enrichment is performed separately for the up- and down-regulated proteins ('rejected'),
#   which are few in our example where we had to set the adjusted p-value to 0.15.
#
# In the enrichment we set the cutoff for the adjusted p-value to 0.2, which is
# a bit arbitrary to see some results.

# %%
enriched = acore.enrichment_analysis.run_up_down_regulation_enrichment(
    regulation_data=diff_reg,
    annotation=annotations,
    min_detected_in_set=1,
    lfc_cutoff=1,
    pval_col="padj",  # toggle if it does not work
    correction_alpha=0.2,  # adjust the p-value to see more or less results
)
enriched

# %% [markdown]
# Plot the enrichment scores for the up- and down-regulated proteins separately.
# - y-axis: GO term
# - x-axis: enrichment score (e.g. -log10(p-pvalue adjusted))

# %% tags=["hide-input"]
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
# # Check for Maltose Uptake

# %%
out_dir_subsection = out_dir / "3_maltose_uptake"
out_dir_subsection.mkdir(parents=True, exist_ok=True)

# %% [markdown]
# Apply filtering of 'differentially abundant proteins' as described in the paper
# > "Differentially abundant proteins were determined as those with log2 fold-change
# > $log2FC > 1$  and $log2FC < -1$, and $p < 0.05$"
#
# This means no multiple testing correction was applied.

# %%
view = diff_reg.query("pvalue < 0.05 and FC > 1")  # .shape[0]
view.to_csv(
    out_dir_subsection / "1_differently_regulated_as_in_paper.csv",
    index=False,
)
view

# %% [markdown]
# Let's find the proteins highlighted in the volcano plot in Figure 3 in the
# [paper](https://www.sciencedirect.com/science/article/pii/S1756464623002451):
# ![Figure 3](https://ars.els-cdn.com/content/image/1-s2.0-S1756464623002451-gr3_lrg.jpg)

# %% tags=["hide-input"]
highlighted_genes = ["LamB", "MalE", "Malk", "CitF", "CitT", "CitE", "Frd"]
highlighted_genes = "|".join([p.upper() for p in highlighted_genes])
highlighted_genes = proteins_meta.query(
    f"`GeneName`.str.contains('{highlighted_genes}')"
)
highlighted_genes

# %% tags=["hide-input"]
highlighted_proteins = "|".join([p.upper() for p in highlighted_genes["ProteinName"]])
view = diff_reg.query(f"`identifier`.str.contains('{highlighted_proteins}')")
view = view.set_index("identifier").join(proteins_meta.set_index("ProteinName"))
view.to_csv(
    out_dir_subsection / "2_highlighted_proteins_in_figure3.csv",
    index=False,
)
sel_cols = [
    "identifier",
    "GeneName",
    "log2FC",
    "pvalue",
    "padj",
    "rejected",
    "group1",
    "group2",
    "Method",
]
view.reset_index()[sel_cols].sort_values("log2FC", ascending=False)

# %% [markdown]
# - See normalized data.

# %% tags=["hide-input"]
view_X = (
    (
        X.loc[:, X.columns.isin(highlighted_genes["ProteinName"].to_list())].T.join(
            proteins_meta.set_index("ProteinName")["GeneName"]
        )
    )
    .set_index("GeneName", append=True)
    .T
)  # to check]
view_X.to_csv(
    out_dir_subsection / "3_highlighted_proteins_in_figure3_intensities_normalized.csv",
    index=True,
)
view_X

# %% [markdown]
# Let's see their original data:
# - See unnormalized data.
# - note that proteins can be discarded based on your filtering criteria.
#
# > How would the analysis change with imputation?

# %% tags=["hide-input"]
view_proteins = (
    (
        proteins[highlighted_genes["ProteinName"].to_list()].T.join(
            proteins_meta.set_index("ProteinName")["GeneName"]
        )
    )
    .set_index("GeneName", append=True)
    .T
)  # to check]
view_proteins.to_csv(
    out_dir_subsection / "3_highlighted_proteins_in_figure3_intensities.csv",
    index=True,
)
view_proteins


# %% [markdown]
# How to explain the differences?
