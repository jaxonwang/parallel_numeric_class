
# FLAGS=-std=c++11 -g -Wall -O0
FLAGS=-std=c++11 -Wall -O3
MPICC=mpiicpc
default: all

ldlt_seq: ldlt_seq.cc datagen.cc datagen.h
	$(CXX) $(FLAGS) ldlt_seq.cc datagen.cc -o ldlt_seq -lm

test: test.cc datagen.cc datagen.h
	$(CXX) $(FLAGS) test.cc datagen.cc -o test -lm

ldlt_leftlooking: ldlt_mpi.cc datagen.cc datagen.h
	export MPICH_CC=$(CXX);$(MPICC) $(FLAGS) -pthread -DLEFTLOOKING ldlt_mpi.cc datagen.cc -o ldlt_leftlooking

ldlt_leftlooking_log: ldlt_mpi.cc datagen.cc datagen.h
	export MPICH_CC=$(CXX);$(MPICC) $(FLAGS) -pthread -DLEFTLOOKING ldlt_mpi.cc datagen.cc -o ldlt_leftlooking_log -L/home/wang/mpiP-3.4.1 -L/usr/lib/x86_64-linux-gnu -lmpiP -lbfd -liberty -lm -lunwind

ldlt_rightlooking: ldlt_mpi.cc datagen.cc datagen.h
	export MPICH_CC=$(CXX);$(MPICC) $(FLAGS) -pthread ldlt_mpi.cc datagen.cc -o ldlt_rightlooking

ldlt_rightlooking_log: ldlt_mpi.cc datagen.cc datagen.h
	export MPICH_CC=$(CXX);$(MPICC) $(FLAGS) -pthread ldlt_mpi.cc datagen.cc -o ldlt_rightlooking_log -L/home/wang/mpiP-3.4.1 -L/usr/lib/x86_64-linux-gnu -lmpiP -lbfd -liberty -lm -lunwind

all: ldlt_seq ldlt_leftlooking ldlt_rightlooking test ldlt_leftlooking_log ldlt_rightlooking_log

normal: ldlt_seq ldlt_leftlooking ldlt_rightlooking

clean:
	rm ldlt_seq
	rm ldlt_leftlooking
	rm ldlt_rightlooking
	rm test
