
FLAGS=-std=c++11 -g -Wall -O0
# FLAGS=-std=c++11 -Wall -O3
default: all

ldlt_seq: ldlt_seq.cc datagen.cc datagen.h
	$(CXX) $(FLAGS) ldlt_seq.cc datagen.cc -o ldlt_seq -lm

test: test.cc datagen.cc datagen.h
	$(CXX) $(FLAGS) test.cc datagen.cc -o test -lm

ldlt_leftlooking: ldlt_mpi.cc datagen.cc datagen.h
	export MPICH_CC=$(CXX);mpic++ $(FLAGS) -pthread -DLEFTLOOKING ldlt_mpi.cc datagen.cc -o ldlt_leftlooking

ldlt_leftlooking_log: ldlt_mpi.cc datagen.cc datagen.h
	export MPICH_CC=$(CXX);mpic++ $(FLAGS) -pthread -DLEFTLOOKING ldlt_mpi.cc datagen.cc -o ldlt_leftlooking -L/home/wang/mpiP-3.4.1 -L/usr/lib/x86_64-linux-gnu -lmpiP -lbfd -liberty -lm

ldlt_rightlooking: ldlt_mpi.cc datagen.cc datagen.h
	export MPICH_CC=$(CXX);mpic++ $(FLAGS) -pthread ldlt_mpi.cc datagen.cc -o ldlt_rightlooking

all: ldlt_seq ldlt_leftlooking ldlt_rightlooking test ldlt_leftlooking_log

clean:
	rm ldlt_seq
	rm ldlt_leftlooking
	rm ldlt_rightlooking
	rm test
