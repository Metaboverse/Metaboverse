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
export VERSION=0.10.1b1


# Check that these paths are correct 
export DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# If the last subfolder of DIR is "resources", then we need to go up one more level
if [[ ${DIR} == *"resources"* ]]; then
    export DIR="$(dirname ${DIR})"
fi
echo -e "\nDIR: ${DIR}\n"

export CONDA=~/miniconda3
export CONDA_PATH=~/miniconda3/etc/profile.d/conda.sh
export APP_PATH=${DIR}/app
export NODE_MODULES=${DIR}/app/node_modules
export MAIN_PATH=${DIR}/app/main.js
export CLI_PATH=${DIR}/cli
export CLI_DEST=${DIR}/app/python
export BUILD_PATH=${DIR}/build
# scp file.zip jsmith@frs.sourceforge.net:/home/frs/project/fooproject/release1
#https://sourceforge.net/projects/metaboverse/files/v0.10.1/mvrs/
export BD_DEST=j-berg@frs.sourceforge.net:/home/frs/project/metaboverse/v${VERSION}


# Update version in app/package.json, cli/metaboverse_cli/__init__.py, CITATION.cff, and docs/conf.py
echo "v${VERSION}" > ${APP_PATH}/__version__.txt

# Extract major and minor version (e.g. if VERSION is "0.10.1", this gets "0.10")
MAJOR_MINOR_VERSION=$(echo "$VERSION" | cut -d'.' -f1,2)

# Get current date in YYYY-MM-DD format
CURRENT_DATE=$(date +'%Y-%m-%d')

# Modify the Python script
sed -i '' "s/^version = .*/version = '$MAJOR_MINOR_VERSION'/" ${DIR}/docs/conf.py
sed -i '' "s/^release = .*/release = '$VERSION'/" ${DIR}/docs/conf.py

# Modify the .cff file
sed -i '' "s/^version: .*/version: $VERSION/" ${DIR}/CITATION.cff
sed -i '' "s/^date-released: .*/date-released: $CURRENT_DATE/" ${DIR}/CITATION.cff

# Check if jq is installed
if ! command -v jq &> /dev/null
then
    echo "jq is not installed. Install it using 'brew install jq' and rerun the script."
    exit
fi

# Modify the main.js file
jq --arg VERSION "$VERSION" '.version = $VERSION' ${APP_PATH}/package.json > ${APP_PATH}/temp.json && mv ${APP_PATH}/temp.json ${APP_PATH}/package.json




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
echo -e "\nBuilding the CLI..."
chmod +x ${DIR}/resources/build-python.sh
${DIR}/resources/build-python.sh
chmod +wrx ${CLI_PATH}/metaboverse-cli*
cp ${CLI_PATH}/dist/metaboverse-cli* ${CLI_DEST}


# Build electron app 
echo -e "\nBuilding the electron app..."
chmod +x ${DIR}/resources/build-electron.sh
${DIR}/resources/build-electron.sh


# Code execution based on BUILD_DB flag
if [ ${BUILD_DB} -eq 1 ]; then
    echo -e "\nBuilding the database(s)..."
    mkdir -p ${BUILD_PATH}
    if [[ ${OS} == *"MINGW"* ]]; then
        cp ${CLI_PATH}/dist/metaboverse-cli*.exe ${BUILD_PATH}/metaboverse-cli.exe
        export BUILD_EXE=${BUILD_PATH}/metaboverse-cli.exe
    else 
        cp ${CLI_PATH}/dist/metaboverse-cli* ${BUILD_PATH}/metaboverse-cli-nix
        export BUILD_EXE=${BUILD_PATH}/metaboverse-cli-nix
    fi
    chmod +x ${DIR}/resources/build-db.sh
    ${DIR}/resources/build-db.sh
else
    echo -e "\nNot building database..."
fi


# Clean up 
echo -e "\nCleaning up..."
if [ ${BUILD_DB} -eq 1 ]; then
  rm -rf ${BUILD_PATH}
fi 
rm ${APP_PATH}/python/metaboverse-cli-darwin

echo -e "\nDone.\n"