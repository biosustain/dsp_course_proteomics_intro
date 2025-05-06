# %%
import gdown

# %%
# Fasta File
url = "https://drive.google.com/drive/folders/1jN70x2eoFoRe-tiPkAcT1Qii1IH6cUCG"
gdown.download_folder(url, output="data/fasta")

# %%
# spectra files
# PXD040621
url = 'https://drive.google.com/drive/folders/150mvfggFJZjkvw2jImPHuchkfHXvOSmF'
gdown.download_folder(url, output="data/PXD040621/mzML")
