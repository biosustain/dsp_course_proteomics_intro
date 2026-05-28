# only for reference, we do not use this file directly
FROM ghcr.io/nextflow-io/training:3.5

# Default shell for following RUN steps
SHELL ["/bin/bash", "-lc"]

RUN conda init bash

RUN source /opt/conda/etc/profile.d/conda.sh && \
    conda activate base && \
    pip install --no-cache-dir -r requirements.txt
