# Container info

The devcontainer used in this project is based on the one used for nf-core training on
nextflow, which you can find in the
[Nextflow training](https://github.com/nextflow-io/training/tree/master/.devcontainer)
`.devcontainer` directory. Their latest version can be different from the currently used
one in this project, so if things break, have a look there.

The devcontainer uses the
[`nextflow-io/training` container](https://github.com/nextflow-io/training/pkgs/container/training)
from GitHub Container Registry, as provided by the nextflow team. I copied the
[Dockerfile](https://github.com/biosustain/dsp_course_proteomics_intro/blob/HEAD/.devcontainer/Dockerfile)
which was used to build the container image `2.1.7` for reference of what the container contains.

> ⚠️ Dockerfile is to be verified

## Additional software

Using an
[`setup.sh` script](https://github.com/biosustain/dsp_course_proteomics_intro/blob/HEAD/.devcontainer/setup.sh)
which is linked to the `.devcontainer.json` file,
additional software is installed in the container. This includes:

- `gdown` for downloading on devcontainer startup the relevant `mzML` files from a private Google Drive folder (:warning: to be changed!)
- `acore` and `vuecore` to run the analysis notebook for the downstream data analysis
- set's up conda so the notebook kernel is recognized by VSCode

## Regular checks

Once a month I check that the website builds, including execution of the downstream analysis notebook
[2_data_analysis.ipynb](../2_data_analysis.ipynb).

That the devcontainer is working is not yet automatically checked, but means manually running one.

## How to update

- check for changes in
  [`nextflow-io/training`](https://github.com/nextflow-io/training/tree/master/.devcontainer)
  repo and then set a new image tag in
  [`.devcontainer.json`](https://github.com/biosustain/dsp_course_proteomics_intro/blob/HEAD/.devcontainer/devcontainer.json)
  file
- update the
  [`setup.sh` script](https://github.com/biosustain/dsp_course_proteomics_intro/blob/HEAD/.devcontainer/setup.sh)
  script if needed, e.g. to install new software or new data sources
