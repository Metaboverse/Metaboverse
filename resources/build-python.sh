#!/bin/bash

cd ${CLI_PATH}

rm -rf ${CONDA}/envs/pyinstaller

# Check if building on an M1 Mac 
if [[ $(uname -m) == "arm64" ]]; then
    CONDA_SUBDIR=osx-64 conda create -n pyinstaller python=3.9 -y
else
    conda create -n pyinstaller python=3.9 pyinstaller -y
fi

source "${CONDA_PATH}"
conda activate pyinstaller 

pip install pyinstaller
pip install -r ${CLI_PATH}/requirements.txt

pyinstaller ${CLI_PATH}/metaboverse-cli.spec