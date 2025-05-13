# QuantMS Hands-On for PXD04621

## Download the data

Currently the data is stored on a
[Google Drive folder](https://drive.google.com/drive/folders/1gxUh9nMx9icFLrI0vn3zAB9dDjZf-1Nh).
A more persistent location would be great - probably on [Zenodo](https://zenodo.org/).

```bash
conda activate base
pip install gdown
python download_PXD040621_data.py
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

## Copy files for further analysis

```bash
cp -aL results/PXD040621/proteomicslfq/. data/PXD040621/processed/
# And maybe save the parameters to reproduce the analysis (for 1.4.0 and above):
cp -aL results/PXD040621/pipeline_info/. data/PXD040621/processed/pipeline_info/
```

Now it is safe to delete the `results/PXD040621` and `work` folder.

```bash
rm -r results/PXD040621 work
```

## Clean up unused docker images

As we are running in a GitHub codespace, we have limited storage. Therefore let's
[clean up our docker images](https://docs.docker.com/engine/manage-resources/pruning/)
store after the analysis is done.

```bash
docker images # see all images
docker image prune -a
```

## Run a different analysis

- feel free to run the project `PXD041301` in quantms if you want an exercise

> [!WARNING]
> not tested in codespace and probably too big for the codespace storage
