#!/bin/bash
###
# Run as: $ bash ./update_version.sh 0.11.1
# Assuming you are the in /Metaboverse/resources/ folder
###
###
###

# Check if version argument is supplied
if [ $# -eq 0 ]; then
    echo "No version argument supplied"
    exit 1
fi

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



# Update version in app/package.json, cli/metaboverse_cli/__init__.py, CITATION.cff, and docs/conf.py
echo "v${VERSION}" > ${APP_PATH}/__version__.txt
echo "__version__='${VERSION}'" > ${CLI_PATH}/metaboverse_cli/__init__.py

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

echo "Version update to ${VERSION} complete."