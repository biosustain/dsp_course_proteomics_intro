# only for reference, we do not use this file directly
FROM ghcr.io/nextflow-io/training:3.5

# Default shell for following RUN steps
SHELL ["/bin/bash", "-lc"]

COPY requirements.txt requirements_docs.txt /tmp/

RUN conda init bash

RUN source /opt/conda/etc/profile.d/conda.sh && \
    conda activate base && \
    pip install --no-cache-dir -r /tmp/requirements.txt
