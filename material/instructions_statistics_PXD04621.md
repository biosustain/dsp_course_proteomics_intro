# Data Analysis Hands-On for PXD04621

In the GitHub workspace, we will use the base conda environment which you can activate
with the following command (it should be active by default when you open the codespace):

```bash
conda activate base
```

## Results data from QuantMS

Are stored in this repository in the
[`data/PXD04621/processed/`](https://github.com/biosustain/dsp_course_proteomics_intro/tree/main/data/PXD040621/processed)
folder.

## Required packages

are mainly:

- [acore](https://analytics-core.readthedocs.io) for some analytical functions (developed
  by DSP and MoNA at DTU)
- [vuegen](https://vuegen.readthedocs.io) for the report generation (developed by DSP and
  MoNA at DTU)
- [vuecore](https://vuecore.readthedocs.io) for the visualization (developed by DSP and
  MoNA at DTU)

The packages should have been installed automatically in the GitHub codespace. In case you
need to manually install them, you can use the following command (for local use):

```bash
pip install acore 'numpy<2.1.0' nbformat vuecore vuegen
```

## Open and run the notebook

[2_data_analysis.ipynb](../2_data_analysis.ipynb)
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

### Inspect streamlit report example online

A deployed version of the report can be accessed online at the following link:

- [dsp-course-proteomics-intro-report.streamlit.app/](https://dsp-course-proteomics-intro-report.streamlit.app/)

### Download the VueGen HTML report

- [html report](pmultiqc_report/multiqc_report_pmultiqc.html)

See instructions on [biosustain/dsp_course_proteomics_intro_report](https://github.com/biosustain/dsp_course_proteomics_intro_report).

### VueGen GUI (app)

We also ship a GUI for VueGen. Check out our
[tutorials on youtube](https://www.youtube.com/playlist?list=PLTbkQyef1c2S3qGzzva_JLlgdwsXjHCHH).
