
FLAGS=-std=c++11 -g -Wall -O0
default: all

ldlt_seq: ldlt_seq.cc datagen.cc
	$(CXX) $(FLAGS) ldlt_seq.cc datagen.cc -o ldlt_seq -lm

ldlt_mpi: ldlt_mpi.cc datagen.cc
	export OMPI_MPICXX=$(CXX);mpic++ $(FLAGS) ldlt_mpi.cc datagen.cc -o ldlt_mpi

all: ldlt_seq ldlt_mpi

clean:
	rm ldlt_seq
	rm ldlt_mpi 
