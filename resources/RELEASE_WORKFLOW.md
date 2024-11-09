
## Test a release

### Tagged release
```
git tag v0.11.5-test
git push origin v0.11.5-test
```

### Undo failed tag
```
git tag -d v0.11.5-test
git push origin :refs/tags/v0.11.5-test
```

## Submit a release
1. Finish any changes to code 
2. Push changes to main branch
3. Update versions on dev branch (do no merge into main) - run `bash ./update_version.sh 0.11.5` or `.\update_version.ps1 0.11.5 .. ..\app\ ..\cli\` from `resources/` dir
4. Tag and push release
5. Edit release details
6. Build local version to build templates and upload to SourceForge - run `bash ./__build-nix__.sh 0.11.5 --build` or `.\__build-windows__.ps1 0.11.5 .. ..\app\ ..\cli\` from the `resources/` directory
7. At this point, you may delete the dev branch

### Tagged release
```
git tag v0.11.5
git push origin v0.11.5
```

### Undo failed tag
```
git tag -d v0.11.5
git push origin :refs/tags/v0.11.5
```