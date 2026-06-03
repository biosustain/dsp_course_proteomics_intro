#!/usr/bin/env bash

# Customise the terminal command prompt
printf "export PS1='\\[\\e[3;36m\\]\${PWD#/workspaces/} ->\\[\\e[0m\\] '\n" >> $HOME/.bashrc
export PS1='\[\e[3;36m\]${PWD#/workspaces/} ->\[\e[0m\] '

# Update Nextflow
nextflow self-update
nextflow -version

cat /usr/local/etc/vscode-dev-containers/first-run-notice.txt

conda init
conda activate base
# pip install gdown
# python 0_download_PXD040621_data.py
pip install -r requirements.txt

# sudo apt update
# sudo apt install -y azcopy

azcopy copy \
"https://dspteaching.blob.core.windows.net/msdatasets/PXD040621/mzML/*?sp=rwl&st=2026-05-29T10:56:55Z&se=2026-06-05T15:11:55Z&spr=https&sv=2026-02-06&sr=d&sig=gCopwuxco%2F8LQRohxJhxATxLszTnMSxpt%2B9gLA0KajU%3D&sdd=2" \
"/workspaces/dsp_course_proteomics_intro/data/PXD040621/mzML"
