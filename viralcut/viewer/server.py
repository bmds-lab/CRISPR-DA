import http.server
from http.server import HTTPServer, BaseHTTPRequestHandler
from os import curdir,sep
import socketserver
import os
from io import BytesIO
#from viralcut import analysis 

PORT = 8080

extensions_map = {
    '.manifest': 'text/cache-manifest',
    '.html': 'text/html',
    '.png': 'image/png',
    '.jpg': 'image/jpg',
    '.svg': 'image/svg+xml',
    '.css': 'text/css',
    '.js':  'application/x-javascript'
}

def get_newick():
    with open('Euk2.txt', 'r') as fp:
        content = fp.readline()
        return content

routes = {
    'newick' : get_newick
}

class RequestHandler(BaseHTTPRequestHandler):
    def do_GET(self):
        is_route_registered = False
    
        if self.path == '/':
            self.path  = 'index.html'
            
        if self.path.startswith('/'):
            self.path = self.path[1:]

        if self.path in routes:
            is_route_registered = True
            content = routes[self.path]()
        
        try:
            sendReply = False
            
            mimeType = 'text/plain'
            for ext in extensions_map:
                if self.path.endswith(ext):
                    mimeType = extensions_map[ext]
                    break

            if not is_route_registered:
                f = open(os.path.join(curdir, self.path), 'rb')
            else:
                f = BytesIO(str.encode(content))
            
            self.send_response(200)
            self.send_header('Content-type', mimeType)
            self.end_headers()
            self.wfile.write(f.read())
            f.close()
            
            return
            
        except IOError:
            self.send_error(404,'File not found!')

Handler = RequestHandler# http.server.SimpleHTTPRequestHandler

Handler.extensions_map = extensions_map

httpd = socketserver.TCPServer(("", PORT), Handler)

print("serving at port", PORT)
httpd.serve_forever()

### from flask import Flask, send_from_directory
### 
### app = Flask(__name__,
###     static_url_path='', 
###     static_folder='static',
###     template_folder='templates'
### )
### 
### 
### @app.route("/")
### def root():
###     return send_from_directory('static', 'index.html')
### 
### if __name__ == '__main__':
###     app.run(host='0.0.0.0', port='8080', debug=True)