SOURCES=$(shell ls *_nocmssw.cu)
EXECS=$(SOURCES:.cu=)
CUFLAGS=-std=c++11 -v -arch=sm_60

.SUFFIXES:

.SUFFIXES: .cu

all: $(EXECS)

.cu :
	/usr/local/cuda/bin/nvcc $(CUFLAGS) $^ -o $@

clean:
	rm $(EXECS)
