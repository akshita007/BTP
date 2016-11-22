CC=g++
CFLAGS=-I.

tiling_file=blocking.cpp
sequential_file=sequential.cpp
centroid_file=centroid_centric.cpp

number_of_points = 100
dimension=7
number_of_clusters=4

all: tiling sequential centroid run
tiling:
	$(CC) $(tiling_file) -fopenmp -o block.o

sequential: 
	$(CC) $(sequential_file) -fopenmp -o seqn.o

centroid : 
	$(CC) $(centroid_file) -fopenmp -o cntroid.o

run:
	./block.o ${number_of_points} ${dimension} ${number_of_clusters} 
	./seqn.o  ${number_of_points} ${dimension} ${number_of_clusters} 
	./cntroid.o  ${number_of_points} ${dimension} ${number_of_clusters} 

clean:
	rm block.o seqn.o cntroid.o
