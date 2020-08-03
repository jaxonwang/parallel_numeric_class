
FLAGS=-std=c++11 -g -Wall -O0
# FLAGS=-std=c++11 -Wall -O3
default: all

ldlt_seq: ldlt_seq.cc datagen.cc datagen.h
	$(CXX) $(FLAGS) ldlt_seq.cc datagen.cc -o ldlt_seq -lm

test: test.cc datagen.cc datagen.h
	$(CXX) $(FLAGS) test.cc datagen.cc -o test -lm

ldlt_leftlooking: ldlt_mpi.cc datagen.cc datagen.h
	export OMPI_MPICXX=$(CXX);mpic++ $(FLAGS) -DLEFTLOOKING ldlt_mpi.cc datagen.cc -o ldlt_leftlooking

ldlt_rightlooking: ldlt_mpi.cc datagen.cc datagen.h
	export OMPI_MPICXX=$(CXX);mpic++ $(FLAGS) ldlt_mpi.cc datagen.cc -o ldlt_rightlooking

all: ldlt_seq ldlt_leftlooking ldlt_rightlooking test

clean:
	rm ldlt_seq
	rm ldlt_leftlooking
	rm ldlt_rightlooking
	rm test
