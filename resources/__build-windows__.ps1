# Run as: .\build-windows.ps1 0.11.0 .. ..\app\ ..\cli\
# Assuming you are the in /Metaboverse/resources/ folder
param (
    [string]$VERSION,
    [string]$DIR,
    [string]$APP_PATH,
    [string]$CLI_PATH
)

$originalPath = Get-Location

# Convert to absolute paths
$DIR = Resolve-Path $DIR
$APP_PATH = Resolve-Path $APP_PATH
$CLI_PATH = Resolve-Path $CLI_PATH

# Build pyinstaller executable
Set-Location -Path $CLI_PATH

$FILE3 = "metaboverse_cli/__init__.py"

# Break down version to get major and minor part
$MAJOR_MINOR_VERSION = $VERSION.Split('.')[0] + '.' + $VERSION.Split('.')[1]

# Get the current date
$CURRENT_DATE = Get-Date -Format "yyyy-MM-dd"

# Substitute version in FILE3
(Get-Content -Path $FILE3) | Foreach-Object {
    $_ -replace "__version__='.*'", "__version__='$VERSION'"
} | Set-Content -Path $FILE3

conda create -n pyinstaller python=3.9 pyinstaller -y
conda activate pyinstaller 

pip install pyinstaller
pip install -r "requirements.txt"

pyinstaller "metaboverse-cli.spec"

conda deactivate

Set-Location $originalPath

# Prepare other files
Set-Location -Path $DIR
Get-Location

$FILE1 = "docs/conf.py"
$FILE2 = "CITATION.cff"

# Break down version to get major and minor part
$MAJOR_MINOR_VERSION = $VERSION.Split('.')[0] + '.' + $VERSION.Split('.')[1]

# Get the current date
$CURRENT_DATE = Get-Date -Format "yyyy-MM-dd"

# Substitute version and date in FILE1 and FILE2
(Get-Content -Path $FILE1) | Foreach-Object {
    $_ -replace "version = .*", "version = '$MAJOR_MINOR_VERSION'" -replace "release = .*", "release = '$VERSION'"
} | Set-Content -Path $FILE1

(Get-Content -Path $FILE2) | Foreach-Object {
    $_ -replace "version: .*", "version: $VERSION" -replace "date-released: .*", "date-released: $CURRENT_DATE"
} | Set-Content -Path $FILE2

# Move backend
Copy-Item -Path "$CLI_PATH\dist\metaboverse-cli-windows.exe" -Destination "$APP_PATH\python\metaboverse-cli-windows.exe"


Set-Location -Path $APP_PATH
Get-Location

# Install node dependencies
if (Test-Path "node_modules") {
    Remove-Item -Path "node_modules" -Force -Recurse
}
npm install electron --save-dev
npm install electron-packager -g
npm install
npm audit fix
Write-Host "========================================================================"
Get-Content package.json
Write-Host "========================================================================"
npm test

# Prep supplemental files
Set-Location -Path "${APP_PATH}/data/"
if (Test-Path "test_data.zip") {
    Remove-Item -Path "test_data.zip" -Force
}
Compress-Archive -Path "test_data" -DestinationPath "test_data.zip"

# Build electron app
Set-Location -Path $APP_PATH
Get-Location

# Get OS name for darwin, win32, or linux
$OS = "win32"
$ARCH = "x64"
$LOGO = "data/icon/win32/metaboverse_logo.ico"

# Build electron package
electron-packager ./ Metaboverse --platform=$OS --arch=$ARCH --icon=$LOGO --app-version=$VERSION --overwrite
Set-Location -Path ".."

# Build release packages
if (Test-Path "${DIR}/Metaboverse-${OS}-${ARCH}-${VERSION}") {
    Remove-Item -Path "${DIR}/Metaboverse-${OS}-${ARCH}-${VERSION}" -Force -Recurse
}
Move-Item -Path "${APP_PATH}/Metaboverse-${OS}-${ARCH}" -Destination "${DIR}/Metaboverse-${OS}-${ARCH}-${VERSION}"
Copy-Item -Path "${APP_PATH}/data/test_data.zip" -Destination "${DIR}/Metaboverse-${OS}-${ARCH}-${VERSION}"

# Zip for distribution 
Compress-Archive -Path "${DIR}/Metaboverse-${OS}-${ARCH}-${VERSION}" -DestinationPath "${DIR}/Metaboverse-${OS}-${ARCH}-${VERSION}.zip"

Write-Host "`nSHA256 checksum:"
Get-FileHash -Path "${DIR}/Metaboverse-${OS}-${ARCH}-${VERSION}.zip" -Algorithm SHA256
Write-Host "`nMD5 checksum:"
Get-FileHash -Path "${DIR}/Metaboverse-${OS}-${ARCH}-${VERSION}.zip" -Algorithm MD5
Write-Host "`n"

# Remove pyinstaller file 
Remove-Item -Path "$APP_PATH\python\metaboverse-cli-windows.exe" -ErrorAction Ignore

# Upload to Github
Set-Location $originalPath