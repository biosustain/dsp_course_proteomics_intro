# %% [markdown]
# # Differential analysis

# %%
import pandas as pd

# %% [markdown]
# ## MSstats results
# obtained using quantms DDA LFQ test example data

# %%
fpath_msstats_results = "../example_outputs/results_lfq/msstats/BSA_design_urls_openms_design_msstats_in_comparisons.csv"
msstats_results = pd.read_csv(fpath_msstats_results, sep="\t")  # tsv file, not csv...
msstats_results

# %% [markdown]
# Results are not between conditions it seems

# %% [markdown]
# ## Run t-test for `ALBU_BOVIN` protein manually?
# - conditions and biological replicates are the same?
