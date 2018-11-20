current_dir=$(shell pwd)
LEVEL=../..

OBJS = $(LEVEL)/bin/SimSystem.o $(LEVEL)/bin/SpatialUnit.o $(LEVEL)/bin/Particle.o $(LEVEL)/bin/utils.o $(LEVEL)/bin/Vector3D.o

.PHONY: all clean webserver replay

all: solve $(OBJS)

run: solve
	./solve | nodejs $(LEVEL)/dpd-vis/visualiser.js

replay: 
	nodejs $(LEVEL)/dpd-vis/visualiser.js 

solve: $(OBJS) solve.cpp
	g++ -std=c++11 -O2 -I $(LEVEL)/inc solve.cpp -o solve $(LEVEL)/bin/*.o -ljsoncpp

$(LEVEL)/bin/%.o: $(LEVEL)/src/%.cpp $(LEVEL)/inc/%.hpp
	mkdir -p $(LEVEL)/bin
	g++ -c -std=c++11 -I $(LEVEL)/inc $< -o $@	

clean:
	rm -rf ./_frames
	rm -rf ./bin
	rm -rf state.json
	rm -rf expr.json
	rm -rf _meta.json
	rm -rf node_modules
	rm -f solve
	rm -f generate
	rm -f *.js *.html
