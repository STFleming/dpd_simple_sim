current_dir=$(shell pwd)
LEVEL=../..

OBJS = $(LEVEL)/bin/SimSystem.o $(LEVEL)/bin/SpatialUnit.o $(LEVEL)/bin/Particle.o $(LEVEL)/bin/utils.o

.PHONY: all clean

all: generate solve $(OBJS)

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
	rm -f solve
	rm -f generate
