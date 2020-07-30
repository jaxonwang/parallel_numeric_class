
datagen.o:
	$(CXX) -std=c++11 -c datagen.cc -o datagen.o

ldlt_seq.o:
	$(CXX) -std=c++11 -c ldlt_seq.cc -o ldlt_seq.o

ldlt_seq: datagen.o ldlt_seq.o
	$(CXX) -std=c++11 ldlt_seq.o datagen.o -o ldlt_seq -lm

all: ldlt_seq

default: all

clean:
	rm datagen.o
	rm ldlt_seq.o
	rm ldlt_seq
