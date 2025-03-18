# Massspectrometry-based Proteomcis introduction

The course will introduce MS-based proteomics data, it's basic processing and
the downstream analysis. The course will be a mix of lectures and hands-on session 
where we will perform the basic processing steps using quantms, a nextflow workflow, and
acore for the downstream analysis.


## Test runs

### LFQ test run

```bash
nextflow run 'https://github.com/bigbio/quantms' -profile docker,test_lfq
```

will create a [folder `results_lfq`]() with the results of the test run.

### Troubleshooting

Check the 
[ci workflow on the bigbio/quantms repository](https://github.com/bigbio/quantms/blob/master/.github/workflows/ci.yml)
for seeing more test workflows.


## Added to devcontainer postCreateCommand

```bash
conda install matplotlib
```
