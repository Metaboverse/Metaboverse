#!/bin/bash
###
### Requirements:
### - Conda and NPM need to be installed and accessible in your user path 
### - `parallel` command line tool needs to be installed on your system 
###       (`brew install parallel` or `sudo apt-get install parallel`)
### 
### Execute as sudo: `$ bash build.sh`
### You may be prompted to login to Sourceforge for DB uploads
### If you have permissions issues, try removing NPM cache: `$ sudo npm cache clean --force`
###

# Check if version argument is supplied
if [ $# -eq 0 ]; then
    echo "No version argument supplied"
    exit 1
fi

# Initialize build flag
BUILD_DB=false

# Process arguments
for arg in "$@"
do
    case $arg in
        --build)
        BUILD_DB=true
        shift # Remove --build_db from processing
        ;;
        *)
        export VERSION=${arg#v}
        ;;
    esac
done
echo "Compiling version: $VERSION"

# Install depedencies 
# If on macOS, install GNU parallel, jq using brew 
# If on Linux, install parallel, jq using apt-get
# If on Windows, install parallel, jq using choco
OS=$(uname -s)
if [[ $OS == *"Darwin"* ]]; then
    echo "Installing dependencies using brew..."
    brew install parallel jq -y
elif [[ $OS == *"Linux"* ]]; then
    echo "Installing dependencies using apt-get..."
    sudo apt-get install parallel jq -y
elif [[ $OS == *"MINGW"* ]]; then
    echo "Installing dependencies using choco..."
    choco install parallel jq -y
else
    echo "Unsupported OS: $OS"
    exit 1
fi

# Check if conda is installed
if ! command -v conda &> /dev/null
then
    echo "conda is not installed. Install it and rerun the script."
    exit
fi

# Check if NPM is installed 
if ! command -v npm &> /dev/null
then
    echo "npm is not installed. Install it and rerun the script."
    exit
fi


# Check that these paths are correct 
export DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# If the last subfolder of DIR is "resources", then we need to go up one more level
if [[ ${DIR} == *"resources"* ]]; then
    export DIR="$(dirname ${DIR})"
fi
echo -e "\nDIR: ${DIR}\n"

# Make sure user permissions are set for ease of downstream building 
#sudo chown -R $USER ${DIR}

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

OS=$(uname -s)
FILE1="${DIR}/docs/conf.py"
FILE2="${DIR}/CITATION.cff"

if [[ $OS == *"Darwin"* ]]; then
    sed -i '' "s/^version = .*/version = '$MAJOR_MINOR_VERSION'/" $FILE1
    sed -i '' "s/^release = .*/release = '$VERSION'/" $FILE1
    sed -i '' "s/^version: .*/version: $VERSION/" $FILE2
    sed -i '' "s/^date-released: .*/date-released: $CURRENT_DATE/" $FILE2
elif [[ $OS == *"Linux"* ]]; then
    sed -i "s/^version = .*/version = '$MAJOR_MINOR_VERSION'/" $FILE1
    sed -i "s/^release = .*/release = '$VERSION'/" $FILE1
    sed -i "s/^version: .*/version: $VERSION/" $FILE2
    sed -i "s/^date-released: .*/date-released: $CURRENT_DATE/" $FILE2
elif [[ $OS == *"MINGW"* ]]; then # Windows
    echo "Please consider using a Linux subsystem or Cygwin to use sed on Windows (or run this script using WSL first)."
    # Windows has a more complex environment for bash-like operations and may require a third-party software like Cygwin, WSL, or Git BASH.
else
    echo "Unsupported OS: $OS"
fi

# Check if jq is installed
if ! command -v jq &> /dev/null
then
    echo "jq is not installed. Install it using 'brew install jq' and rerun the script."
    exit
fi

# Modify the main.js file
jq --arg VERSION "$VERSION" '.version = $VERSION' ${APP_PATH}/package.json > ${APP_PATH}/temp.json && mv ${APP_PATH}/temp.json ${APP_PATH}/package.json


# Build cli 
echo -e "\nBuilding the CLI..."
#chmod 755 ${DIR}/resources/build-python.sh
${DIR}/resources/build-python.sh
#chmod 755 ${CLI_PATH}/metaboverse-cli*
cp ${CLI_PATH}/dist/metaboverse-cli* ${CLI_DEST}


# Build electron app 
echo -e "\nBuilding the electron app..."
#chmod 755 ${DIR}/resources/build-electron.sh
${DIR}/resources/build-electron.sh


# Code execution based on BUILD_DB flag
if [ "$BUILD_DB" = true ]; then
    echo -e "\nBuilding the database(s)..."
    mkdir -p ${BUILD_PATH}
    if [[ ${OS} == *"MINGW"* ]]; then
        cp ${CLI_PATH}/dist/metaboverse-cli*.exe ${BUILD_PATH}/metaboverse-cli.exe
        export BUILD_EXE=${BUILD_PATH}/metaboverse-cli.exe
    else 
        cp ${CLI_PATH}/dist/metaboverse-cli* ${BUILD_PATH}/metaboverse-cli-nix
        export BUILD_EXE=${BUILD_PATH}/metaboverse-cli-nix
    fi
    #chmod 755 ${DIR}/resources/build-db.sh
    ${DIR}/resources/build-db.sh
else
    echo -e "\nNot building database..."
fi


# Clean up 
echo -e "\nCleaning up..."
if [ "$BUILD_DB" = true ]; then
  rm -rf ${BUILD_PATH}
fi 
rm ${APP_PATH}/python/metaboverse-cli-darwin

echo -e "\nDone.\n"