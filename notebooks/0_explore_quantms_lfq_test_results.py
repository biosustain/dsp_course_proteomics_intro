# %% [markdown]
# # Explore quantMS LFQ test results
# using the tools used in the quantMS LFQ test workflow for DDA example data
# - OpenMS: raw spectra processing data analysis software
# - Triqler: tool differential error estimation on the protein level based on precursor-level quantifications
# - MSstats: downstream statistics software
#


# %%
import numpy as np
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
    "Run",
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

# %% [markdown]
# ## Protein Quantification
# protein quantification results are not save as outputs directly. We can use the sum 
# or median to summarize the peptide quantifications for each protein.
#
# Other tools are available to make the aggregation.

# %%
proteins_intensities = precursors_openms.groupby(["ProteinName", 'Run'])['Intensity'].sum().to_frame()
proteins_intensities

# %%
np.log2(proteins_intensities)

# %%
ax = np.log2(proteins_intensities).hist()

# %%
