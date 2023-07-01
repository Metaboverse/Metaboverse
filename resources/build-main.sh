#!/bin/bash
###
### Requirements:
### - Conda and NPM need to be installed and accessible in your user path 
### - `parallel` command line tool needs to be installed on your system 
###       (`brew install parallel` or `sudo apt-get install parallel`)
### 
### Execute as sudo: `$ sudo bash build.sh`
### You may be prompted to login to Sourceforge for DB uploads
###


# Change this for each release 
export VERSION=0.10.1


# Check that these paths are correct 
export DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
export CONDA_PATH=~/miniconda3/etc/profile.d/conda.sh
export APP_PATH=$DIR/../app
export NODE_MODULES=$DIR/../app/node_modules
export MAIN_PATH=$DIR/../app/main.js
export CLI_PATH=$DIR/../cli
export CLI_DEST=$DIR/../app/python
export BUILD_PATH=$DIR/../build
# scp file.zip jsmith@frs.sourceforge.net:/home/frs/project/fooproject/release1
#https://sourceforge.net/projects/metaboverse/files/v0.10.1/mvrs/
export BD_DEST=j-berg@frs.sourceforge.net:/home/frs/project/metaboverse/v${VERSION}


# Parse optional argument "--build-db"
BUILD_DB=0
while (( "$#" )); do
  case "$1" in
    --build-db)
      BUILD_DB=1
      shift
      ;;
    --) # end argument parsing
      shift
      break
      ;;
    -*|--*=) # unsupported flags
      echo "Error: Unsupported flag $1" >&2
      exit 1
      ;;
  esac
done


# Build cli 
echo "Building the CLI..."
./build-python.sh
cp $CLI_PATH/dist/metaboverse-cli* $CLI_DEST


# Build electron app 
echo "Building the electron app..."
./build-electron.sh


# Code execution based on BUILD_DB flag
if [ ${BUILD_DB} -eq 1 ]; then
    echo "Building the database(s)..."
    mkdir -p ${BUILD_PATH}
    if [[ ${OS} ==*"MINGW"* ]]; then
        cp ${CLI_PATH}/dist/metaboverse-cli*.exe ${BUILD_PATH}/metaboverse-cli.exe
        export BUILD_EXE=${BUILD_PATH}/metaboverse-cli.exe
    else 
        cp ${CLI_PATH}/dist/metaboverse-cli* ${BUILD_PATH}/metaboverse-cli-nix
        export BUILD_EXE=${BUILD_PATH}/metaboverse-cli-nix
    fi
    ./build-db.sh
else
    echo "Not building database..."
fi


# Clean up 
rm -rf ${BUILD_PATH}


