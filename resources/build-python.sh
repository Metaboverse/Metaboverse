#!/bin/bash

cd ${CLI_PATH}
rm -r dist/ build/


echo -e "\nInstalling pyinstaller environment..."
rm -rf ${CONDA}/envs/pyinstaller

# Check if building on an M1 Mac 
if [[ $(uname -m) == "arm64" ]]; then
    CONDA_SUBDIR=osx-64 conda create -n pyinstaller python=3.9 -y
else
    conda create -n pyinstaller python=3.9 -y
fi

source "${CONDA_PATH}"
conda activate pyinstaller 

echo -e "\nInstalling pyinstaller dependencies..."
pip install pyinstaller
pip install -r ${CLI_PATH}/requirements.txt

echo -e "\nBuilding the CLI..."
pyinstaller ${CLI_PATH}/metaboverse-cli.spec

conda deactivate