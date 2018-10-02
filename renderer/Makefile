all: launch_webserver

./node_modules:
	npm install express
	npm install ws

launch_webserver: ./node_modules state.json
	nodejs main

clean:
	rm -rf ./node_modules 
	rm -rf package-lock.json
