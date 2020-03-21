// d3.json('../data/net_data_toy1.json').then(data=>{
//     let metaGraph = new MetaGraph(data);
// })
var fs = require('fs')
var app = require('electron').remote.app
const exec = require('child_process').exec;
var userDataPath = app.getPath('userData');
var session_file = userDataPath + "/session_data.json"
let database_url = JSON.parse(fs.readFileSync(session_file).toString())["database_url"]
// get_session_info("database_url");
console.log("Database path: " + database_url)

// var data = JSON.parse(fs.readFileSync(database_url).toString());
// console.log(data)

function execute(command, callback) {

    exec(command, (error, stdout, stderr) => {
        callback(stdout);
        callback(stderr);
    });

};

execute("python python/motif_process.py "+database_url, (output) => {
    console.log(output);
    d3.json('../data/processed_data.json').then(data=>{
        let metaGraph = new MetaGraph(data);
    })
});


// let metaGraph = new MetaGraph(data);

// let content = new FormData();
// $.ajax({
//     type: "POST",
//     enctype: 'multipart/form-data',
//     url: "/initialize",
//     data: content,
//     processData: false, //prevent jQuery from automatically transforming the data into a query string
//     contentType: false,
//     cache: false,
//     dataType:'json',
//     success: function (response) {
//         console.log(response)
//         // clean_svg();
//         // let metaGraph = new MetaGraph(response);
//         // metaGraph.motifSearch();
//     },
//     error: function (error) {
//         console.log("error",error);
//     }
// });

// $("#import").click(function(){
//     $("#files").click();
// });

// d3.select("#files")
//     .on("change",()=>{
//         let form = $('#upload')[0];
//         let content = new FormData(form);
//         $.ajax({
//             type: "POST",
//             enctype: 'multipart/form-data',
//             url: "/import",
//             data: content,
//             processData: false, //prevent jQuery from automatically transforming the data into a query string
//             contentType: false,
//             cache: false,
//             dataType:'json',
//             success: function (response) {
//                 clean_svg();
//                 let metaGraph = new MetaGraph(response);
//                 metaGraph.draw_graph();
//             },
//             error: function (error) {
//                 console.log("error",error);
//             }
//         });
//     });

function clean_svg(){
    $('#networkSVG').remove();
    $('#container').append('<svg id="networkSVG" width="800" height="800"></svg>');
} 

