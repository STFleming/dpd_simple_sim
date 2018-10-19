current_dir=$(shell pwd)

OBJS = bin/SimSystem.o bin/SpatialUnit.o bin/Particle.o bin/utils.o bin/Vector3D.o

.PHONY: all clean

all: $(OBJS) dpd-vis/node_modules

dpd-vis/node_modules:
	make -C dpd-vis

bin/%.o: ./src/%.cpp ./inc/%.hpp
	mkdir -p ./bin
	g++ -c -std=c++11 -I ./inc $< -o $@	

clean:
	rm -rf ./bin
	rm -rf *.json
	rm -rf dpd-vis/node_modules 
	rm -rf dpd-vis/package-lock.json

