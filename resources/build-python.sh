#!/bin/bash

# Replace version number in $CLI_PATH/__version__.txt with the value in $VERSION 
cd ${CLI_PATH}
echo ${VERSION} > ${CLI_PATH}/__version__.txt


# Check if building on an M1 Mac 
if [[ $(uname -m) == "arm64" ]]; then
    CONDA_SUBDIR=osx-64 conda create -n pyinstaller python=3.9 -y
else
    conda create -n pyinstaller python=3.9 pyinstaller -y
fi

source "${CONDA_PATH}"
conda activate pyinstaller 

pip install pyinstaller
pip install -r requirements.txt

pyinstaller metaboverse-cli.spec