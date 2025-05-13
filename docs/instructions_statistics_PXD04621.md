# Data Analysis Hands-On for PXD04621

In the GitHub workspace, we will use the base conda environment which you can activate
with the following command:

```bash
conda activate base
```

## Download the data

> If you copied it from the quantms results before, you can skip this step.

```bash
# pip install gdown (if not already installed)
python 1_download_PXD040621_results.py
```


## Install the required packages

```bash
pip install acore 'numpy<2.1.0'
```

## Open and run the notebook

[2_read_peptide_output.ipynb](../2_read_peptide_output.ipynb)
is a Jupyter notebook that will guide you through the analysis of the data.

> We will go through it an discuss it together in the class.


## VueGen Report

We wrote some files in the notebook to `data/PXD04621/report` folder. 

You can use the following command to generate a VueGen report from them:

```bash
vuegen -dir data/PXD040621/report
```

This will create a streamlit app in the `streamlit_report` folder. Check the options 
availble for the command line tool vuegen with the following command:

```
vuegen -h
```

Create and start the streamlit app with the following command:

```bash
vuegen -dir data/PXD040621/report -st_autorun
```

### VueGen GUI (app)

We also ship a GUI for VueGen. Check out our 
[tutorials on youtube](https://www.youtube.com/playlist?list=PLTbkQyef1c2S3qGzzva_JLlgdwsXjHCHH).
