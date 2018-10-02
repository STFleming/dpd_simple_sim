current_dir=$(shell pwd)

OBJS = bin/SimSystem.o bin/SpatialUnit.o bin/Particle.o bin/utils.o

.PHONY: all clean

all: $(OBJS)

bin/%.o: ./src/%.cpp ./inc/%.hpp
	mkdir -p ./bin
	g++ -c -std=c++11 -I ./inc $< -o $@	

clean:
	rm -rf ./bin
	rm -rf *.json
