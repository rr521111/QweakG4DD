CXXFLAGS = $(shell root-config --cflags) -I./include/ -ggdb -O0 -Wall
LIBS = $(shell root-config --glibs)
SRCS := \
	averageModel.cc \
	src/interpolatePEs.cc
	src/interpolatePOLs.cc

OFILES := $(addsuffix .o,$(basename $(SRCS)))

all: avgModel

avgModel: $(OFILES)
		g++ $(CXXFLAGS) -o $@ $^ $(LIBS)

%.o: %.cpp
		g++ $(CXXFLAGS) -o $@ $^

clean:
	rm -f avgModel *.o src/*.o
