
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
3. Update versions on dev branch (do no merge into main)
4. Tag and push release
5. Edit release details

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