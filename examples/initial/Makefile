current_dir=$(shell pwd)
LEVEL=../..

OBJS = $(LEVEL)/bin/SimSystem.o $(LEVEL)/bin/SpatialUnit.o $(LEVEL)/bin/Particle.o $(LEVEL)/bin/utils.o

.PHONY: all clean webserver

all: generate solve $(OBJS)

state.json: generate
	./generate

launch_webserver: state.json webserver
	npm install express
	npm install ws
	nodejs main

webserver: state.json
	cp $(LEVEL)/renderer/*.js ./
	cp $(LEVEL)/renderer/*.html ./

solve: $(OBJS) solve.cpp
	g++ -std=c++11 -O2 -I $(LEVEL)/inc solve.cpp -o solve $(LEVEL)/bin/*.o -ljsoncpp

generate: $(OBJS) generate.cpp
	g++ -std=c++11 -O2 -I $(LEVEL)/inc generate.cpp -o generate $(LEVEL)/bin/*.o -ljsoncpp

$(LEVEL)/bin/%.o: $(LEVEL)/src/%.cpp $(LEVEL)/inc/%.hpp
	mkdir -p $(LEVEL)/bin
	g++ -c -std=c++11 -I $(LEVEL)/inc $< -o $@	

clean:
	rm -rf ./bin
	rm -rf *.json
	rm -rf node_modules
	rm -f solve
	rm -f generate
	rm -f *.js *.html
