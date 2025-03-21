#!/bin/bash

cd ${APP_PATH}
pwd


# Install node dependencies
echo -e "\nInstalling node dependencies..."
rm -rf ${NODE_MODULES}
npm install electron --save-dev -y
npm install electron-packager -g -y
npm install -y
npm audit fix -y
echo "========================================================================"
cat package.json
echo "========================================================================"
npm test


# Prep supplemental files
echo -e "\nPrepping supplemental files..."
cd ${APP_PATH}/data/
if [ -f test_data.zip ]; then
    echo "test_data.zip already exists. Removing the existing zip file."
    rm test_data.zip
fi
echo "Creating test_data.zip..."
zip -q -r test_data.zip test_data
#chmod 755 test_data.zip


# Build electron app
cd ${APP_PATH}
pwd


# Get OS name for darwin, win32, or linux
OS=$(uname -s)
if [[ $OS == *"Darwin"* ]]; then
    OS="darwin"
    ARCH="x64"
    LOGO="data/icon/nix/metaboverse_logo.icns"
    OS_LABEL="MacOS"  # Label for filename
elif [[ $OS == *"Linux"* ]]; then
    OS="linux"  # Keep original OS value for electron-packager
    ARCH="x64"
    SYS_ARCH=$(uname -m)  # Store the system architecture
    LOGO="data/icon/png/icon_1024x1024.png"
    OS_LABEL="Linux"  # Label for filename
elif [[ $OS == *"MINGW"* ]]; then
    OS="win32"  # Keep original OS value for electron-packager
    ARCH="x64"
    LOGO="data/icon/win32/metaboverse_logo.ico"
    OS_LABEL="Windows"  # Label for filename
else
    echo "Unsupported OS: $OS"
    exit 1
fi


# Build electron package
echo -e "\nBuilding the electron app..."
electron-packager ./ Metaboverse --platform=${OS} --arch=${ARCH} --icon=${LOGO} --build-version=${VERSION} --overwrite
cd ..


# Build release packages

#####
echo -e "\nBuilding the database(s)..."
#chmod -R 755 ${APP_PATH}/Metaboverse-${OS}-${ARCH}
# If ${DIR}/Metaboverse-${OS}-${ARCH}-${VERSION} already exists, delete it
if [ -d "${DIR}/Metaboverse-${OS_LABEL}-${ARCH}-${VERSION}" ]; then
    rm -rf ${DIR}/Metaboverse-${OS_LABEL}-${ARCH}-${VERSION}
fi
mv ${APP_PATH}/Metaboverse-${OS}-${ARCH} ${DIR}/Metaboverse-${OS_LABEL}-${ARCH}-${VERSION}
cp ${APP_PATH}/data/test_data.zip ${DIR}/Metaboverse-${OS_LABEL}-${ARCH}-${VERSION}


# Make OS-specific modifications to package
#if [[ ${OS} == "linux" ]]; then
    #chmod -R 755 ${DIR}/Metaboverse-${OS}-${ARCH}-${VERSION}/Metaboverse
    #chmod -R 755 ${DIR}/Metaboverse-${OS}-${ARCH}-${VERSION}/resources/app/python/metaboverse-cli-linux
#fi
#if [[ ${OS} == "darwin" ]]; then
    #chmod 755 ${DIR}/Metaboverse-${OS}-${ARCH}-${VERSION}/Metaboverse.app/Contents/Resources/app/python/metaboverse-cli-darwin
#fi


# Zip for distribution 
echo -e "\nZipping the release package..."
cd ${DIR}

# For Linux, use the system architecture for the filename
if [[ $OS == "linux" ]]; then
    ZIP_ARCH=$SYS_ARCH
else
    # Currently, the x64 Electron packager works fine for x64 and arm64, while arm64 does not seem to work at all. This is a temporary workaround
    ZIP_ARCH=$ARCH
    if [[ $(uname -m) == *"arm64" ]]; then 
        ZIP_ARCH="arm64"
    fi
fi

zip -q -r Metaboverse-${OS_LABEL}-${ZIP_ARCH}-${VERSION}.zip Metaboverse-${OS_LABEL}-${ARCH}-${VERSION}

echo -e "\nZipped release package:"
#chmod 755 ${DIR}/Metaboverse-${OS_LABEL}-${ZIP_ARCH}-${VERSION}.zip
echo -e "\nSHA256 checksum:"
shasum -a 256 ${DIR}/Metaboverse-${OS_LABEL}-${ZIP_ARCH}-${VERSION}.zip
echo -e "\nMD5 checksum:"
# If on linux, use md5sum, if on mac, use md5
if [[ $OS == "linux" ]]; then
    md5sum ${DIR}/Metaboverse-${OS_LABEL}-${ZIP_ARCH}-${VERSION}.zip
else
    md5 ${DIR}/Metaboverse-${OS_LABEL}-${ZIP_ARCH}-${VERSION}.zip
fi
echo -e "\n"

# Upload to Github
