# Run as: .\__build-windows__.ps1 0.11.0 .. ..\app\ ..\cli\
# Assuming you are the in /Metaboverse/resources/ folder
param (
    [string]$VERSION,
    [string]$DIR,
    [string]$APP_PATH,
    [string]$CLI_PATH
)

# Convert to absolute paths
$DIR = Resolve-Path $DIR
$APP_PATH = Resolve-Path $APP_PATH
$CLI_PATH = Resolve-Path $CLI_PATH

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