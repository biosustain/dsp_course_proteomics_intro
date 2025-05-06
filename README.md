# Massspectrometry-based Proteomcis introduction

The course will introduce MS-based proteomics data, it's basic processing and
the downstream analysis. The course will be a mix of lectures and hands-on session 
where we will perform the basic processing steps using quantms, a nextflow workflow, and
acore for the downstream analysis.


## Run nextflow

```bash
nextflow run bigbio/quantms \
         -revision 1.3.0 \
         -params-file PXD040621_w_contaminants-params.yaml \
         -profile docker \
         -resume
```

## help commands


Check the local storage usage (you have maximum of 32GB in a GitHub codespace)

```bash
du -hd 1 
```

