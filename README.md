# Massspectrometry-based Proteomcis introduction

The course will introduce MS-based proteomics data, it's basic processing and
the downstream analysis. The course will be a mix of lectures and hands-on session 
where we will perform the basic processing steps using quantms, a nextflow workflow, and
acore for the downstream analysis.


## Run nextflow

```
nextflow run bigbio-quantms_1.3.0/1_3_0/main.nf \
         -params-file PXD040621_w_contaminants-params.yaml \
         -profile docker \
         -resume
```



