# %% [markdown]
# # Explore processed data for PXD041301

import matplotlib.pyplot as plt
# %%
import pandas as pd

# %%
fname_in = "PXD041301.sdrf_openms_design_msstats_in.csv"

# %%
df = pd.read_csv(fname_in, sep=",", header=0)  # .set_index([])
df = df.query('not ProteinName.str.contains("CON_")')
df.head()

# %%
# df.columns.to_list()
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
df['Condition'].unique()

# %%
proteins = df.groupby(["ProteinName", "Condition", "BioReplicate"])["Intensity"].median().unstack(
    level=0
)
proteins

# %%
proteins.columns = proteins.columns.str.split("_").str[0].str.split("|").str[-1]
proteins['BioRep'] = proteins.groupby("Condition").cumcount() + 1
proteins = proteins.reset_index(level=1,drop=True).set_index('BioRep', append=True)
proteins

# %%
protein_counts = proteins.notna().sum()
ax = protein_counts.sort_values(ascending=True).plot(rot=90)
ax.set_ylabel("Number of samples")
ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))

# %%
protein_counts.value_counts().sort_index().plot.bar()

# %%
