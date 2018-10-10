// creates a node.js based server used to render the incoming particle movements from the simulation 

const express = require('express');
const app = express();
const WebSocket = require('ws');

// reading from STDIN: https://stackoverflow.com/questions/20086849/how-to-read-from-stdin-line-by-line-in-node
var readline = require('readline');
var r1 = readline.createInterface({
  input: process.stdin,
  output: process.stdout,
  terminal: false
}); 

var fs = require('fs');
var path = require('path');

// open our web socket
const wss = new WebSocket.Server({port: 8079});

// send a particle position update  to the rendered instance
//   the update is in ndjson format
function sendUpdate(wss, update_str) {
    // https://github.com/websockets/ws#simple-server
    wss.clients.forEach(function each(client) {
     if(client !== wss && client.readyState == WebSocket.OPEN) {
      client.send(update_str);
      console.log(update_str);
     }
    });
}

// updates the frames folder
function update_frames() {
    fs.readdir('frames', function(err, files) {
        if(err) {
          console.error("could not list the directory", err); 
        } 
        files.forEach(function(file) {
            //console.log(file);
            app.get('/frames/'+file, (req, res) => res.sendFile((path.join(__dirname+'/frames/'+file))));
        });
    });
}

// get stdin data which will be passed to the rendered graph
// this will be the output of the executive
r1.on('line', function(line) {
  sendUpdate(wss, line); 
  update_frames();
});


app.get('/', (req, res) => res.sendFile(path.join(__dirname+'/index.html')));
app.get('/state.json', (req, res) => res.sendFile(path.join(__dirname+'/state.json')));
app.get('/state_frame.json', (req, res) => res.sendFile(path.join(__dirname+'/state_frame.json')));
app.get('/d3-3d.js', (req, res) => res.sendFile(path.join(__dirname+'/d3-3d.js')));

app.listen(3000, () => console.log('listening on port 3000'))
  
