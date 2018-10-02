current_dir=$(shell pwd)

OBJS = bin/SimSystem.o bin/SpatialUnit.o bin/Particle.o bin/utils.o

.PHONY: all clean

all: test $(OBJS)

test: $(OBJS) src/test.cpp
	g++ -std=c++11 -O2 -I ./inc src/test.cpp -o ./bin/test ./bin/*.o

bin/%.o: ./src/%.cpp ./inc/%.hpp
	mkdir -p ./bin
	g++ -c -std=c++11 -I ./inc $< -o $@	

clean:
	rm -rf ./bin
	rm -rf *.json
