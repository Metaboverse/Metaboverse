// Initialize species list
var reactome_api = "https://reactome.org/ContentService/data/species/all";

$.getJSON(reactome_api, function(data) {

  // Get species name and ID from Reactome API
  var abbreviation_dict = {}
  data.forEach(function(datum) {

      abbreviation_dict[datum["displayName"]] = datum["abbreviation"]


  });

  // Get species names (keys) as list
  speciesList = Object.getOwnPropertyNames(
    abbreviation_dict
  ).map(function(k) {
    return k;
  });
  speciesList.unshift("Select an organism..."); // Add select prompt to menu bar

  // Generate drop-down menu for species select
  var menu = document.getElementById("speciesMenu");
  for (var i = 0; i < speciesList.length; i++) {
    var option = document.createElement("option");
    option.innerHTML = speciesList[i];
    option.value = speciesList[i];
    menu.appendChild(option);
  };

});

// Change user selection based on input
function selectOrganism() {

  var selection = document.getElementById("speciesMenu").value;
  console.log(selection)
  secondary.style.display = 'block';
  tertiary.style.display = 'block';

};


// ********************** Drag & drop starts here ********************

let dropArea = document.getElementById("drop-area")

// Prevent default drag behaviors
;['dragenter', 'dragover', 'dragleave', 'drop'].forEach(eventName => {
dropArea.addEventListener(eventName, preventDefaults, false)
document.body.addEventListener(eventName, preventDefaults, false)
})

// Highlight drop area when item is dragged over it
;['dragenter', 'dragover'].forEach(eventName => {
dropArea.addEventListener(eventName, highlight, false)
})

;['dragleave', 'drop'].forEach(eventName => {
dropArea.addEventListener(eventName, unhighlight, false)
})

// Handle dropped files
dropArea.addEventListener('drop', handleDrop, false)

function preventDefaults (e) {
e.preventDefault()
e.stopPropagation()
}

function highlight(e) {
dropArea.classList.add('highlight')
}

function unhighlight(e) {
dropArea.classList.remove('active')
}

function handleDrop(e) {
var dt = e.dataTransfer
var files = dt.files

handleFiles(files)
}

function handleFiles(files) {
files = [...files]
initializeProgress(files.length)
files.forEach(uploadFile)
quinary.style.display = 'block';
}

function uploadFile(file, i) {
var xhr = new XMLHttpRequest()
var formData = new FormData()
xhr.open('POST', url, true)
xhr.setRequestHeader('X-Requested-With', 'XMLHttpRequest')

}

// ******************** Hide/Show starts here ******************

function toggleVis1() {
  quaternary.style.display = 'block';
}

function toggleVis2() {
  quinary.style.display = 'block';
}

// ****************** Progress Bar starts here ****************
