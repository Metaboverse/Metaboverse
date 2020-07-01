#!/bin/bash

VERSION=0.1.3b

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cd $DIR

rm app/package-lock.json
rm -rf app/node_modules/

cd app
pwd
npm install electron electron-packager mocha wine minimist zerorpc --save-dev
npm audit fix
npm install --save-dev --save-exact jsdom jsdom-global
npm audit fix
npm install fs path d3 jquery string-pixel-width d3-fisheye d3-save-svg jsonpickle plotly.js-dist filereader jpickle pickle save-svg-as-png unpickle update-electron-app xvfb jqueryui --save
npm audit fix
echo "========================================================================"
cat package.json
echo "========================================================================"

npm test
cd ..

mv app/python/metaboverse-cli-linux .
mv app/python/metaboverse-cli-mac .
cd app
pwd
electron-packager ./ Metaboverse --platform=win32 --icon=data/icon/metaboverse_logo.icns --overwrite
cd ..

mv app/python/metaboverse-cli-win.exe .
mv metaboverse-cli-mac app/python/
cd app
pwd
electron-packager ./ Metaboverse --platform=darwin --icon=data/icon/metaboverse_logo.icns --overwrite
cd ..

mv app/python/metaboverse-cli-mac .
mv metaboverse-cli-linux app/python/
cd app
pwd
electron-packager ./ Metaboverse --platform=linux --icon=data/icon/metaboverse_logo.icns --overwrite
cd ..

pwd
mv metaboverse-cli-mac app/python/
mv metaboverse-cli-win.exe app/python/

mv app/Metaboverse-darwin-x64 ./Metaboverse-darwin-x64-${VERSION}
mv app/Metaboverse-linux-x64 ./Metaboverse-linux-x64-${VERSION}
mv app/Metaboverse-win32-x64 ./Metaboverse-win32-x64-${VERSION}

zip -r ./Metaboverse-darwin-x64-${VERSION}.zip ./Metaboverse-darwin-x64-${VERSION}
zip -r ./Metaboverse-linux-x64-${VERSION}.zip ./Metaboverse-linux-x64-${VERSION}
zip -r ./Metaboverse-win32-x64-${VERSION}.zip ./Metaboverse-win32-x64-${VERSION}

shasum -a 256 ./Metaboverse-darwin-x64-${VERSION}.zip
shasum -a 256 ./Metaboverse-linux-x64-${VERSION}.zip
shasum -a 256 ./Metaboverse-win32-x64-${VERSION}.zip
