# %%
import gdown

# processed results file (for convenience)
url = "https://drive.google.com/drive/folders/1Nm5Ha-tCvjU-B323BLhna1GwHdNpK_lU?usp=drive_link"
gdown.download_folder(url, output="data/PXD040621/mzML", use_cookies=False)
