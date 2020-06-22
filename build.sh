#!/bin/bash

VERSION=0.1.0b

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cd $DIR

mv app/python/metaboverse-cli-linux .
mv app/python/metaboverse-cli-mac .
cd app
electron-packager ./ Metaboverse --platform=win32 --icon=data/icon/metaboverse_logo.icns --overwrite
cd ..

mv app/python/metaboverse-cli-win.exe .
mv metaboverse-cli-mac app/python/
cd app
electron-packager ./ Metaboverse --platform=darwin --icon=data/icon/metaboverse_logo.icns --overwrite
cd ..

mv app/python/metaboverse-cli-mac .
mv metaboverse-cli-linux app/python/
cd app
electron-packager ./ Metaboverse --platform=linux --icon=data/icon/metaboverse_logo.icns --overwrite
cd ..

mv metaboverse-cli-mac app/python/
mv metaboverse-cli-win.exe app/python/

mv app/Metaboverse-darwin-x64 ./Metaboverse-darwin-x64-${VERSION}
mv app/Metaboverse-linux-x64 ./Metaboverse-linux-x64-${VERSION}
mv app/Metaboverse-win32-x64 ./Metaboverse-win32-x64-${VERSION}

zip -r ./Metaboverse-darwin-x64-${VERSION}.zip ./Metaboverse-darwin-x64-${VERSION}
zip -r ./Metaboverse-linux-x64-${VERSION}.zip ./Metaboverse-linux-x64-${VERSION}
zip -r ./Metaboverse-win32-x64-${VERSION}.zip ./Metaboverse-win32-x64-${VERSION}
