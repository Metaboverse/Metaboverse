# Run as: .\__build-windows__.ps1 0.11.0 .. ..\app\ ..\cli\
# Assuming you are the in /Metaboverse/resources/ folder
param (
    [string]$VERSION,
    [string]$DIR,
    [string]$APP_PATH,
    [string]$CLI_PATH
)

# Store the original location but don't change it
$originalPath = Get-Location

# Function to safely remove items with optimized performance
function Remove-ItemSafely {
    param (
        [string]$Path,
        [int]$MaxRetries = 2,
        [int]$RetryDelaySeconds = 1
    )
    
    if (-not (Test-Path $Path)) {
        # Path doesn't exist, so consider it a success
        return $true
    }
    
    # First attempt - just try to remove it directly
    try {
        Remove-Item -Path $Path -Force -Recurse -ErrorAction Stop
        return $true
    } catch {
        Write-Host "Initial removal attempt failed for $Path, trying with process termination..."
    }
    
    # If we're here, the first attempt failed - try more aggressive approach
    $retryCount = 0
    
    while ($retryCount -lt $MaxRetries) {
        try {
            # Only check for processes if it's a directory
            if (Test-Path $Path -PathType Container) {
                # Get only the most likely processes that could be locking files
                $processNames = @("electron", "Metaboverse", "node", "npm", "python")
                foreach ($procName in $processNames) {
                    Get-Process -Name $procName* -ErrorAction SilentlyContinue | 
                        Stop-Process -Force -ErrorAction SilentlyContinue
                }
            }
            
            # Quick garbage collection
            [System.GC]::Collect()
            [System.GC]::WaitForPendingFinalizers()
            
            # Try removal again
            Remove-Item -Path $Path -Force -Recurse -ErrorAction Stop
            Write-Host "Successfully removed: $Path"
            return $true
        } catch {
            $retryCount++
            if ($retryCount -lt $MaxRetries) {
                Write-Host "Retrying removal of $Path (Attempt $($retryCount+1) of $MaxRetries)..."
                Start-Sleep -Seconds $RetryDelaySeconds
            } else {
                Write-Host "Warning: Could not remove $Path after $MaxRetries attempts. Continuing anyway."
                # Return true anyway to not block the script
                return $true
            }
        }
    }
    
    return $false
}

# Function to terminate processes quickly
function Stop-RelevantProcesses {
    $processNames = @("electron", "Metaboverse", "node")
    
    foreach ($procName in $processNames) {
        Get-Process -Name $procName* -ErrorAction SilentlyContinue | 
            Stop-Process -Force -ErrorAction SilentlyContinue
    }
    
    # A shorter wait is usually sufficient
    Start-Sleep -Seconds 2
}

# Convert to absolute paths
$DIR = Resolve-Path $DIR
$APP_PATH = Resolve-Path $APP_PATH
$CLI_PATH = Resolve-Path $CLI_PATH

# Break down version to get major and minor part
$MAJOR_MINOR_VERSION = $VERSION.Split('.')[0] + '.' + $VERSION.Split('.')[1]

# Get the current date
$CURRENT_DATE = Get-Date -Format "yyyy-MM-dd"

# Update version numbers in all required files
Write-Host "Updating version numbers in files..."

try {
    # Update __version__.txt in app directory (missing from original script)
    Write-Host "Updating app version file..."
    $versionFilePath = Join-Path $APP_PATH "__version__.txt"
    "v$VERSION" | Set-Content -Path $versionFilePath
    Write-Host "Updated $versionFilePath with version v$VERSION"

    # Update CLI version
    Write-Host "Updating CLI version..."
    $FILE3 = Join-Path $CLI_PATH "metaboverse_cli/__init__.py"

    if (-not (Test-Path $FILE3)) {
        Write-Host "Error: $FILE3 not found" -ForegroundColor Red
    } else {
        # Substitute version in FILE3
        (Get-Content -Path $FILE3) | Foreach-Object {
            $_ -replace "__version__='.*'", "__version__='$VERSION'"
        } | Set-Content -Path $FILE3
        
        Write-Host "Updated $FILE3 with version $VERSION"
    }

    # Update documentation and citation files
    Write-Host "Updating documentation files..."
    
    $FILE1 = Join-Path $DIR "docs/conf.py"
    $FILE2 = Join-Path $DIR "CITATION.cff"

    # Update docs/conf.py
    if (-not (Test-Path $FILE1)) {
        Write-Host "Warning: $FILE1 not found, skipping update" -ForegroundColor Yellow
    } else {
        # Substitute version in FILE1
        (Get-Content -Path $FILE1) | Foreach-Object {
            $_ -replace "version = .*", "version = '$MAJOR_MINOR_VERSION'" -replace "release = .*", "release = '$VERSION'"
        } | Set-Content -Path $FILE1
        
        Write-Host "Updated $FILE1 with version $MAJOR_MINOR_VERSION and release $VERSION"
    }

    # Update CITATION.cff
    if (-not (Test-Path $FILE2)) {
        Write-Host "Warning: $FILE2 not found, skipping update" -ForegroundColor Yellow
    } else {
        # Substitute version and date in FILE2
        (Get-Content -Path $FILE2) | Foreach-Object {
            $_ -replace "version: .*", "version: $VERSION" -replace "date-released: .*", "date-released: $CURRENT_DATE"
        } | Set-Content -Path $FILE2
        
        Write-Host "Updated $FILE2 with version $VERSION and release date $CURRENT_DATE"
    }

    # Update package.json in the app directory if it exists
    $packageJsonPath = Join-Path $APP_PATH "package.json"
    if (Test-Path $packageJsonPath) {
        Write-Host "Updating package.json version..."
        
        # Load the JSON file
        $packageJson = Get-Content -Path $packageJsonPath -Raw | ConvertFrom-Json
        
        # Update the version
        $packageJson.version = $VERSION
        
        # Save the updated JSON
        $packageJson | ConvertTo-Json -Depth 10 | Set-Content -Path $packageJsonPath
        
        Write-Host "Updated package.json with version $VERSION"
    }

    Write-Host "Version update completed successfully" -ForegroundColor Green
} catch {
    Write-Host "Error updating version numbers: $_" -ForegroundColor Red
}

# Build pyinstaller executable
Push-Location -Path $CLI_PATH
Write-Host "Working in: $CLI_PATH"

# Ensure we're starting with a clean environment
Write-Host "Cleaning up existing environment..."
Remove-ItemSafely -Path "build"
Remove-ItemSafely -Path "dist"
Remove-ItemSafely -Path "$APP_PATH\python\metaboverse-cli-windows.exe"

# Create and activate conda environment
Write-Host "Setting up conda environment..."
conda create -n pyinstaller python=3.9 pyinstaller -y
conda activate pyinstaller 

# Install dependencies
Write-Host "Installing dependencies..."
pip install --upgrade pip
pip install pyinstaller
pip install -r "requirements.txt"

# Build the executable
Write-Host "Building executable..."
pyinstaller "metaboverse-cli.spec" --clean

# Verify the executable was created
$exePath = Join-Path $CLI_PATH "dist\metaboverse-cli-windows.exe"
if (-not (Test-Path $exePath)) {
    Write-Host "Error: Executable was not created at $exePath" -ForegroundColor Red
    exit 1
}

# Test the executable
Write-Host "Testing executable..."
try {
    $result = & $exePath --version
    Write-Host "Executable test successful: $result" -ForegroundColor Green
} catch {
    Write-Host "Error testing executable: $_" -ForegroundColor Red
    exit 1
}

conda deactivate

Pop-Location
Write-Host "Returned to: $(Get-Location)"

# Ensure the python directory exists
$pythonDir = Join-Path $APP_PATH "python"
if (-not (Test-Path $pythonDir)) {
    New-Item -ItemType Directory -Path $pythonDir -Force | Out-Null
}

# Move backend with verification
Write-Host "Moving executable to app directory..."
$sourcePath = Join-Path $CLI_PATH "dist\metaboverse-cli-windows.exe"
$destPath = Join-Path $APP_PATH "python\metaboverse-cli-windows.exe"

# Stop any running instances of the executable
Get-Process -Name "metaboverse-cli-windows" -ErrorAction SilentlyContinue | Stop-Process -Force

# Copy with retry logic
$maxRetries = 3
$retryCount = 0
$copySuccess = $false

while (-not $copySuccess -and $retryCount -lt $maxRetries) {
    try {
        Copy-Item -Path $sourcePath -Destination $destPath -Force
        if (Test-Path $destPath) {
            $copySuccess = $true
            Write-Host "Successfully copied executable to $destPath" -ForegroundColor Green
        } else {
            throw "File not found after copy"
        }
    } catch {
        $retryCount++
        Write-Host "Copy attempt $retryCount failed: $_" -ForegroundColor Yellow
        if ($retryCount -lt $maxRetries) {
            Start-Sleep -Seconds 2
        }
    }
}

if (-not $copySuccess) {
    Write-Host "Failed to copy executable after $maxRetries attempts" -ForegroundColor Red
    exit 1
}

Push-Location -Path $APP_PATH
Write-Host "Working in: $APP_PATH"

# Install node dependencies
Remove-ItemSafely -Path "node_modules"

npm install electron --save-dev
npm install electron-packager -g
npm install
npm audit fix
Write-Host "========================================================================"
Get-Content package.json
Write-Host "========================================================================"
npm test

# Prep supplemental files
Push-Location -Path "data"
Write-Host "Working in: $(Get-Location)"

if (Test-Path "test_data.zip") {
    Remove-Item -Path "test_data.zip" -Force
}
Compress-Archive -Force -Path "test_data" -DestinationPath "test_data.zip"

Pop-Location
Write-Host "Returned to: $(Get-Location)"

# Build electron app
Write-Host "Building electron app in: $(Get-Location)"

# Get OS name for darwin, win32, or linux
$OS = "win32"
$ARCH = "x64"
$LOGO = "data/icon/win32/metaboverse_logo.ico"

# Build electron package
electron-packager ./ Metaboverse --platform=$OS --arch=$ARCH --icon=$LOGO --app-version=$VERSION --overwrite

# Terminate processes more efficiently
Stop-RelevantProcesses

Pop-Location
Write-Host "Returned to: $(Get-Location)"

# Build release packages
$targetDir = "${DIR}/Metaboverse-${OS}-${ARCH}-${VERSION}"
$targetZip = "${DIR}/Metaboverse-${OS}-${ARCH}-${VERSION}.zip"

# Remove existing target directory if it exists
Remove-ItemSafely -Path $targetDir

# Try to move the packaged app to the target directory
$moveSuccess = $false
$retryCount = 0
$maxRetries = 2  # Reduced from 3 to 2 for speed

while (-not $moveSuccess -and $retryCount -lt $maxRetries) {
    try {
        Move-Item -Path "${APP_PATH}/Metaboverse-${OS}-${ARCH}" -Destination $targetDir -ErrorAction Stop
        $moveSuccess = $true
    } catch {
        $retryCount++
        Write-Host "Failed to move packaged app (Attempt $retryCount of $maxRetries): $_"
        
        # Terminate processes more efficiently
        Stop-RelevantProcesses
        
        if ($retryCount -lt $maxRetries) {
            Write-Host "Retrying in 2 seconds..."
            Start-Sleep -Seconds 2  # Reduced from 5 to 2 seconds
        } else {
            Write-Host "Trying to copy instead of move..."
            
            # As a fallback, try to copy the files instead of moving them
            try {
                New-Item -Path $targetDir -ItemType Directory -Force | Out-Null
                # Use robocopy for faster copying with better handling of locked files
                $robocopyResult = (Start-Process -FilePath "robocopy.exe" -ArgumentList "`"${APP_PATH}/Metaboverse-${OS}-${ARCH}`" `"$targetDir`" /E /NFL /NDL /NJH /NJS /nc /ns /np" -Wait -PassThru).ExitCode
                # Robocopy exit codes 0-7 indicate success with varying levels of copying activity
                if ($robocopyResult -ge 0 -and $robocopyResult -le 7) {
                    $moveSuccess = $true
                    Write-Host "Files copied successfully using robocopy."
                } else {
                    Write-Host "Robocopy failed with exit code $robocopyResult, trying PowerShell copy..."
                    Copy-Item -Path "${APP_PATH}/Metaboverse-${OS}-${ARCH}/*" -Destination $targetDir -Recurse -Force
                    $moveSuccess = $true
                }
            } catch {
                Write-Host "Failed to copy files: $_"
            }
        }
    }
}

if ($moveSuccess) {
    # Copy test data
    Copy-Item -Path "${APP_PATH}/data/test_data.zip" -Destination $targetDir -ErrorAction SilentlyContinue
    
    # Zip for distribution with error handling
    try {
        # Remove existing zip if it exists
        if (Test-Path $targetZip) {
            Remove-Item -Path $targetZip -Force
        }
        
        # Use 7-Zip if available for better handling of locked files
        $7zipPath = "C:\Program Files\7-Zip\7z.exe"
        if (Test-Path $7zipPath) {
            Write-Host "Using 7-Zip for compression..."
            & $7zipPath a -tzip $targetZip $targetDir
        } else {
            # Try to use the faster System.IO.Compression method first
            try {
                Write-Host "Using .NET compression..."
                Add-Type -AssemblyName System.IO.Compression.FileSystem
                [System.IO.Compression.ZipFile]::CreateFromDirectory($targetDir, $targetZip)
            } catch {
                Write-Host ".NET compression failed, falling back to PowerShell Compress-Archive..."
                # Force garbage collection before compression
                [System.GC]::Collect()
                [System.GC]::WaitForPendingFinalizers()
                Compress-Archive -Path "$targetDir\*" -DestinationPath $targetZip -Force
            }
        }
        
        # Calculate and display checksums
        if (Test-Path $targetZip) {
            Write-Host "`nSHA256 checksum:"
            Get-FileHash -Path $targetZip -Algorithm SHA256
            Write-Host "`nMD5 checksum:"
            Get-FileHash -Path $targetZip -Algorithm MD5
        } else {
            Write-Host "Warning: Zip file was not created successfully."
        }
    } catch {
        Write-Host "Error during compression or checksum calculation: $_"
    }
} else {
    Write-Host "Failed to prepare the package directory. Skipping compression and checksum steps."
}

# Script is complete
Write-Host "Build process completed."