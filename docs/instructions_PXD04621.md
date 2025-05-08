# Instructions for PXD04621

## Download the data

Currently the data is stored on a 
[Google Drive folder](https://drive.google.com/drive/folders/1gxUh9nMx9icFLrI0vn3zAB9dDjZf-1Nh).
A more persistent location would be great - probably on [Zenodo](https://zenodo.org/). 

```bash
conda activate base
pip install gdown
python download.py
```

## Run the analysis

```bash
# export NXF_VER=24.10.6
nextflow run bigbio/quantms \
         -revision 1.4.0 \
         -params-file PXD040621_w_contaminants-params.yaml \
         -profile docker,gitpod \
         -resume
```

If you run locally on a Mac with Apple Silicion (M-ships), you need to addtionally the `arm` profile:

```bash
# export NXF_VER=24.10.6
nextflow run bigbio/quantms \
         -revision 1.4.0 \
         -params-file PXD040621_w_contaminants-params.yaml \
         -profile docker,arm,gitpod \
         -resume
```


### Run quantms 1.3.0

```bash
nextflow run bigbio/quantms \
         -revision 1.3.0 \
         -params-file PXD040621_w_contaminants-params.yaml \
         -profile docker \
         -resume
```

