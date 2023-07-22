
var fs = require('fs');
var { ipcRenderer } = require('electron');

// Get app and user paths
function refresh_session() {
    ipcRenderer.invoke('get-paths').then((paths) => {
        console.log("Refreshing session data...")
        fs.copyFile(
            paths.sessionFileTemplatePath,
            paths.sessionFilePath,
            err => {
                if (err) throw err;
                console.log("Session data file was copied for this session");
            }
        );
    })
}