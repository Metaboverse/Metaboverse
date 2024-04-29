# Release a new version of Metaboverse

## On Windows:
``` 
.\__build-windows__.ps1 0.11.0 .. ..\app\ ..\cli\
```

## On MacOS:
```
bash ./__build-nix__.sh 0.11.1
```

## On Linux:
```
bash ./__build-nix__.sh 0.11.1 --build
```
This will compile the Linux release, as well as the curated databases for each model organism to be uploaded to SourceForge.    

If running on Windows via WSL, you may encounter the following:
```
./__build-nix__.sh: line 20: $'\r': command not found
./__build-nix__.sh: line 21: $'\r': command not found
./__build-nix__.sh: line 33: syntax error near unexpected token `$'do\r''
'/__build-nix__.sh: line 33: `do
```

If so, this indicates that the build file has been edited or saved in a Windows environment, which uses carriage return and line feed ("\r\n") as line endings.     
To fix:
```
sudo apt update
sudo apt install dos2unix
dos2unix __build-nix__.sh
dos2unix build-python.sh
dos2unix build-electron.sh
dos2unix build-db.sh
```
or 
```
sed -i 's/\r$//' __build-nix__.sh
sed -i 's/\r$//' build-python.sh
sed -i 's/\r$//' build-electron.sh
sed -i 's/\r$//' build-db.sh
```

## Instructions
After each OS is compiled, you should see a corresponding zip archive, similar to `Metaboverse-{OS}-{release-version}.zip`. These need to be uploaded to the GitHub release.