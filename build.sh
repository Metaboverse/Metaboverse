#!/bin/bash

VERSION=0.10.1

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cd $DIR

# Install node dependencies
rm -rf app/node_modules/

cd app
pwd
npm install
npm audit fix
echo "========================================================================"
cat package.json
echo "========================================================================"

npm test

# Build python backend docker image 
cd python
cp ../../LICENSE .
cp ../../README.md .
cp ../__version__.txt . 
docker build -t metaboverse-cli .

# Prep supplemental files
cd $DIR
cd app/data/
rm test_data.zip
zip -r test_data.zip test_data

# Build electron app
cd $DIR
mv app/python/metaboverse-cli-linux .
mv app/python/metaboverse-cli-darwin .
cd app
pwd
electron-packager ./ Metaboverse --platform=win32 --icon=data/icon/win/metaboverse_logo.ico --electronVersion=13.6.9 --overwrite
cd ..

mv app/python/metaboverse-cli-windows.exe .
mv metaboverse-cli-darwin app/python/
cd app
pwd
electron-packager ./ Metaboverse --platform=darwin --icon=data/icon/mac/metaboverse_logo.icns --electronVersion=13.6.9 --overwrite
cd ..

mv app/python/metaboverse-cli-darwin .
mv metaboverse-cli-linux app/python/
cd app
pwd
electron-packager ./ Metaboverse --platform=linux --icon=data/icon/png/icon_1024x1024.png --electronVersion=13.6.9 --overwrite
cd ..

pwd
mv metaboverse-cli-darwin app/python/
mv metaboverse-cli-windows.exe app/python/

# Build release packages
mv app/Metaboverse-darwin-x64 ./Metaboverse-darwin-x64-${VERSION}
mv app/Metaboverse-linux-x64 ./Metaboverse-linux-x64-${VERSION}
mv app/Metaboverse-win32-x64 ./Metaboverse-win32-x64-${VERSION}

cp app/data/test_data.zip ./Metaboverse-darwin-x64-${VERSION}
cp app/data/test_data.zip ./Metaboverse-linux-x64-${VERSION}
cp app/data/test_data.zip ./Metaboverse-win32-x64-${VERSION}

chmod +x ./Metaboverse-linux-x64-${VERSION}/Metaboverse
chmod +x ./Metaboverse-linux-x64-${VERSION}/resources/app/python/metaboverse-cli-linux

chmod +x ./Metaboverse-darwin-x64-${VERSION}/Metaboverse.app/Contents/Resources/app/python/metaboverse-cli-darwin

zip -r ./Metaboverse-darwin-x64-${VERSION}.zip ./Metaboverse-darwin-x64-${VERSION}
zip -r ./Metaboverse-linux-x64-${VERSION}.zip ./Metaboverse-linux-x64-${VERSION}
zip -r ./Metaboverse-win32-x64-${VERSION}.zip ./Metaboverse-win32-x64-${VERSION}

shasum -a 256 ./Metaboverse-darwin-x64-${VERSION}.zip
shasum -a 256 ./Metaboverse-linux-x64-${VERSION}.zip
shasum -a 256 ./Metaboverse-win32-x64-${VERSION}.zip
