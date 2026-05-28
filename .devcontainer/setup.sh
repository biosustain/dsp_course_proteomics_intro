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
