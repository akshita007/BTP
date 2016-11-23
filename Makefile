CC=g++
CFLAGS=-I.

tiling_file=blocking.cpp
sequential_file=sequential.cpp
centroid_file=centroid_centric.cpp
thread_file=multithreading.cpp
input_file=Large_Input

number_of_points = 2075100
dimension=7
number_of_clusters= 10

all: tiling sequential centroid threading run
tiling:
	$(CC) $(tiling_file) -fopenmp -o block.o

sequential: 
	$(CC) $(sequential_file) -fopenmp -o seqn.o

centroid : 
	$(CC) $(centroid_file) -fopenmp -o centroid.o

threading:
	$(CC) $(thread_file) -fopenmp -o thread.o

run:
	./block.o  -n ${number_of_points} -d ${dimension} -k ${number_of_clusters} -i ${input_file} -t time -c clusters
	./seqn.o  -n ${number_of_points} -d ${dimension} -k ${number_of_clusters} -i ${input_file} -t time -c clusters
	./centroid.o  -n ${number_of_points} -d ${dimension} -k ${number_of_clusters} -i ${input_file} -t time -c clusters
	./thread.o  -n ${number_of_points} -d ${dimension} -k ${number_of_clusters} -i ${input_file}  -t time -c clusters

clean:
	rm *.o time clusters
