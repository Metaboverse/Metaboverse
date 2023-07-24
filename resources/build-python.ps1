# Run as: .\build-python.ps1 ..\cli\
param (
    [string]$VERSION,
    [string]$CLI_PATH
)
$originalPath = Get-Location

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