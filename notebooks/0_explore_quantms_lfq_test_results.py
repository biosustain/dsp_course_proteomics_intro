# %% [markdown]
# # Explore quantMS LFQ test results
# using the tools used in the quantMS LFQ test workflow for DDA example data
# - OpenMS: raw spectra processing data analysis software
# - Triqler: tool differential error estimation on the protein level based on precursor-level quantifications
# - MSstats: downstream statistics software
#


# %%
import pandas as pd

# %% [markdown]
# ## Precursor quantification results

# %%
fpath_precursors_triqler = (
    "../example_outputs/results_lfq/proteomicslfq/BSA_design_urls_openms_design_triqler_in.tsv"
)
fpath_precursors_openms = (
    "../example_outputs/results_lfq/proteomicslfq/BSA_design_urls_openms_design_msstats_in.csv"
)

# %%
precursors_openms = pd.read_csv(fpath_precursors_openms)
precursors_openms.head()

# %%
_index = [
    "ProteinName",
    "PeptideSequence",
    "PrecursorCharge",
    "Condition",
    "BioReplicate",
]
precursors_openms.set_index(_index, inplace=True)
# %%
precursors_openms.sort_index().loc["P02769|ALBU_BOVIN"]

# %%
precursor_triqler = pd.read_csv(fpath_precursors_triqler, sep="\t")
precursor_triqler.head()

# %%
precursor_triqler.describe()

# %%
precursor_triqler.describe(exclude="number")

# %% [markdown]
# How many peptides are quantified for each protein?

# %%
precursor_triqler.groupby("proteins").size()

# %%
