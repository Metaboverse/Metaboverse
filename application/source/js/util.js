window.addEventListener("load", function() {

  document.getElementById("element-upload").onchange = function(event) {

    var reader = new FileReader();
    reader.readAsDataURL(event.srcElement.files[0]);

    reader.onload = function () {

      var fileContent = reader.result;

	  console.log(fileContent);

    }

}});
