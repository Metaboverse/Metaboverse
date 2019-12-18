// components adapted from https://stackoverflow.com/a/8871836/9571488
window.addEventListener("load", function() {

  document.getElementById("element-upload").onchange = function (event) {

    var f = event.srcElement.files[0];
    console.log(f)

    if (!f.type.match('application/json')) {
        alert('Invalid JSON file');
      } else {

        var reader = new FileReader();
        reader.onloadend = function(e) {
        var result = JSON.parse(this.result);
        console.log(result);
        sessionStorage.setItem("saved_session", JSON.stringify(result));
      };
      reader.readAsText(f);

      }
    }
});
