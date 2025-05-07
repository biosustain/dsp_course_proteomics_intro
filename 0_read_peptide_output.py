# %% [markdown]
# # Read Peptide Output
# Formatted for input into MsStats

# %%
import acore.differential_regulation
import numpy as np
import pandas as pd
import plotly.express as px

# %% tags=["parameters"]
file_in = "data/PXD040621/processed/PXD040621.sdrf_openms_design_msstats_in.csv"
df = pd.read_csv(file_in, sep=",", header=0)  # .set_index([])
df.head()

# %%
df['Intensity'] = np.log2(df['Intensity'].astype(float))

# %%
cols = [
    "ProteinName",
    "PeptideSequence",
    "PrecursorCharge",
    "FragmentIon",
    "ProductCharge",
    "IsotopeLabelType",
    "Condition",
    "BioReplicate",
    "Run",
    "Intensity",
    "Reference",
]


# %%
proteins = (
    df.groupby(["ProteinName", "Reference"])["Intensity"].median().unstack(level=0)
)
proteins

# %% [markdown]
# Remove the contaminant proteins which were added to the fasta file used in the data processing.
# - contaminant proteins are e.g. creation which gets into the sample from the human skin or hair
# when the sample is prepared.
# These are filtered out as they are most of the time not relevant, but a contamination.

# %%
decoy_proteins = proteins.filter(like='CON_', axis=1)
proteins = proteins.drop(decoy_proteins.columns, axis=1)
proteins

# %% [markdown]
# Create a label for each sample based on the metadata.
# - we use a string in the sample name, but you can see how the metadata is organized

# %%
meta = df[["Condition", "BioReplicate", "Run", "Reference"]].drop_duplicates()
meta

# %%
label_suf = pd.Series(
    proteins.index.str.contains("_Suf_").astype(int),
    index=proteins.index,
    name="label_suf",
    dtype=np.int8,
)
label_suf

# %% [markdown]
# Plot the data completeness for each protein.

# %%
ax = proteins.notna().sum().sort_values().plot(rot=45)

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

# %%
view = proteins
group = "label_suf"
diff_reg = acore.differential_regulation.run_anova(
    view.dropna(how='any', axis=1).join(label_suf),
    alpha=0.15,
    drop_cols=[],
    subject=None,
    group=group,
).sort_values(
    "pvalue", ascending=True
)
diff_reg["rejected"] = diff_reg["rejected"].astype(bool) 
diff_reg.sort_values("pvalue")

# %%
diff_reg.plot(x='log2FC', y='-log10 pvalue', kind='scatter', title=group)

# %%  [markdown]
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
# apply filtering of 'differentially abundant proteins' as described in the paper
# > Differentially abundatn proteins were determined as those with log2 fold-change
# > > 1 and < -1, and p < 0.05
# This means not multiple testing correction was applied.

# %%
diff_reg.query('pvalue < 0.05 and FC > 1') #.shape[0]

# %% [markdown]
# Let's find the proteins hightlighted in the volcano plot in Figure 3.

# %%
highlighted_proteins = ["LamB", "MalE", "Malk", "CitF", 'CitT', 'CitE', 'Frd']
highlighted_proteins = '|'.join([p.upper() for p in highlighted_proteins])

diff_reg.query(
    f"`identifier`.str.contains('{highlighted_proteins}')"
)

# %% [markdown]
# How to explain the differences?
