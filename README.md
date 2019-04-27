# Map-Reduce-with-dynamic-thread-balancing
This is a project on speeding up Map Reduce on multi-core, multi-node systems implemented in OpenMP and MPI. 

Map Reduce is a word count problem, where multiple input files are read and multiple files containing word count is output. 
There are reader, mapper, reducer and writer threads which can be parallelized and hence we get a speedup by running on multi-core
multi-node systems.

The job of each of these is:
Reader: Parse input files and identify words (separated by punctuations).
Mapper: Group identified words. Each unique word is a group, thus a hash function helps.
Reducer: Since multiple mapper threads are working in parallel, a reducer combines the work of all mappers.
Writer: The writer dumps the output of reducer in desired format.

Commands:
#OpenMP
To compile OpenMP code.
g++ -fopenmp OMP_wordcount.cc -o OMP_wordcount.o 
run command:
./OMP_wordcount.o 18 1 1 
First argument - number of input files to read. Input files are kept in Files directory.
Second argument - number of threads to run on.
Third argument -  the number of times the files are read (to increase the work load). 

#MPI
To compile MPI code:
mpiicpc -fopenmp  MPI_wordcount.cc -o MPI_wordcount.o -mt_mpi
run command:
mpiexec -n 2 ./MPI_wordcount.o 18 5 1 2
First argument - number of input files to read. Input files are kept in Files directory.
Second argument - number of threads that each process spawns. These are OpenMP threads.
Third argument -  the number of times the files are read (to increase the work load).
Fourth argument - the number of processes(on each process OMP code is run)
NOTE1: For mpiexec, the last argument passed should match with the argument passed to -n (otherwise segfault happens)
NOTE2: -n 2 means two nodes are used. If you do not have a cluster access, give -n 1.

