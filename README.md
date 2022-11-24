### CS451 Assignment 2

there are four files

1. gauss_pthread.c and its executable gauss_pthread.
2. gauss_openmp.c and its executable gauss_openmp.


how to compile

to compile the implementation of pthread
use the following command

gcc  gauss_pthreads.c -o gauss_pthreads -lpthread


to compile the implementation of openmp
use the following command

gcc  gauss_openmp.c -o gauss_openmp -fopenmp


how to execute

for pthreads

./gauss_pthreads 1000 0 4
 
 
 where 1000 is the size of matrix, 0 is the seed, and 4 is the number of threads
 
 
 
for openmp

./gauss_openmp 1000 0 4
 
 
 where 1000 is the size of matrix, 0 is the seed, and 4 is the number of threads
 
 
