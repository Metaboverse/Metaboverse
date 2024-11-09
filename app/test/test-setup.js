// test-setup.js
const jsdom = require('jsdom');
const { JSDOM } = jsdom;

const dom = new JSDOM('<!DOCTYPE html><html><body></body></html>', {
  url: 'http://localhost/',
  referrer: 'http://localhost/',
  contentType: 'text/html',
  includeNodeLocations: true,
  storageQuota: 10000000,
  resources: 'usable',
  runScripts: 'dangerously',
  beforeParse(window) {
    // Mock XMLHttpRequest to prevent actual HTTP requests
    window.XMLHttpRequest = class MockXMLHttpRequest {
      open() {}
      send() {}
      setRequestHeader() {}
    };
  }
});

// Set up global variables that would normally be available in the browser
global.window = dom.window;
global.document = dom.window.document;
global.navigator = {
  userAgent: 'node.js'
};

// Clean up JSDOM
process.on('exit', () => {
  dom.window.close();
});

module.exports = {
  dom
};