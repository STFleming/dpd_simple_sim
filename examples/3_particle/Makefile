current_dir=$(shell pwd)
LEVEL=../..

OBJS = $(LEVEL)/bin/SimSystem.o $(LEVEL)/bin/SpatialUnit.o $(LEVEL)/bin/Particle.o $(LEVEL)/bin/utils.o

.PHONY: all clean webserver replay

all: solve $(OBJS)

replay: webserver
	nodejs main

launch_webserver: webserver solve
	npm install express
	npm install ws
	./solve | nodejs main

webserver: solve 
	cp $(LEVEL)/renderer/*.js ./
	cp $(LEVEL)/renderer/index.html ./index.html
	mkdir -p ./frames

solve: $(OBJS) solve.cpp
	g++ -std=c++11 -O2 -I $(LEVEL)/inc solve.cpp -o solve $(LEVEL)/bin/*.o -ljsoncpp

$(LEVEL)/bin/%.o: $(LEVEL)/src/%.cpp $(LEVEL)/inc/%.hpp
	mkdir -p $(LEVEL)/bin
	g++ -c -std=c++11 -I $(LEVEL)/inc $< -o $@	

clean:
	rm -rf ./frames
	rm -rf ./bin
	rm -rf *.json
	rm -rf node_modules
	rm -f solve
	rm -f generate
	rm -f *.js *.html
