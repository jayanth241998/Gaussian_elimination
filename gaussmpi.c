#include<stdio.h>
#include<mpi.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/times.h>
#include <sys/time.h>
#include <time.h>


/* Program Parameters */
#define MAXN 2500  /* Max value of N */
int N;  /* Matrix size */
/* Matrices and vectors */
float A[MAXN][MAXN], B[MAXN], X[MAXN];
int processors;
int rank;
time_t now;

unsigned int time_seed() {
  struct timeval t;
  struct timezone tzdummy;

  gettimeofday(&t, &tzdummy);
  return (unsigned int)(t.tv_usec);
}


void parameters(int argc, char **argv) {
	int seed = 0;  /* Random seed */
	// char uid[32]; /*User name */

	/* Read command-line arguments */
	srand(time_seed());  /* Randomize */

	if (argc == 3) {
		seed = atoi(argv[2]);
		srand(seed);
		printf("Random seed = %i\n", seed);
	} 
	if (argc >= 2) {
		N = atoi(argv[1]);
		if (N < 1 || N > MAXN) {
			printf("N = %i is out of range.\n", N);
			exit(0);
		}
	}
	else {
		printf("Usage: %s <matrix_dimension> [random seed]\n",
				argv[0]);    
		exit(0);
	}
	

	/* Print parameters */
	printf("\nMatrix dimension N = %i.\n", N);
	
}

void initialize_inputs() {
  int row, col;

  printf("\nInitializing...\n");
  for (col = 0; col < N; col++) {
    for (row = 0; row < N; row++) {
      A[row][col] = (float)rand() / 32768.0;
    }
    B[col] = (float)rand() / 32768.0;
    X[col] = 0.0;
  }

}

void print_inputs() {
  int row, col;

  if (N < 10) {
    printf("\nA =\n\t");
    for (row = 0; row < N; row++) {
      for (col = 0; col < N; col++) {
	printf("%5.2f%s", A[row][col], (col < N-1) ? ", " : ";\n\t");
      }
    }
    printf("\nB = [");
    for (col = 0; col < N; col++) {
      printf("%5.2f%s", B[col], (col < N-1) ? "; " : "]\n");
    }
  }
}


void gauss() {
 int norm, row, col,startingRow,endingRow,processor,k;   /*Normalization row, and zeroing
			* element row and col */
  float multiplier;
  // declaring request and status objects to store requests and status of each processors sends
  MPI_Request requests[processors];    
  MPI_Status  status[processors];
  
  //broadcasting the size of the matrix
  MPI_Bcast(&N,1,MPI_INT,0,MPI_COMM_WORLD);
  
   /*Gaussian elimination */
  for (norm = 0; norm < N - 1; norm++)
  { 
  
    //broadcasting the normalization row of A and B in each iteration
    MPI_Bcast((void *)&A[norm][0],N,MPI_FLOAT,0,MPI_COMM_WORLD);
    MPI_Bcast((void *)&B[norm],1,MPI_FLOAT,0,MPI_COMM_WORLD);
    
  
    MPI_Bcast(&norm,1,MPI_INT,0,MPI_COMM_WORLD);
    
     if(rank == 0)
     {
        
       // dividing the rows for each processors in processor 0 using static interleaving
       for(processor=1; processor < processors ; processor++) 
       {
         for(k = norm + 1 + processor; k < N ; k+= processors){
          MPI_Isend((void *)&A[k][0],N,MPI_FLOAT,processor,0,MPI_COMM_WORLD,&requests[processor]);
          MPI_Wait(&requests[processor],&status[processor]);
          MPI_Isend((void *)&B[k],1,MPI_FLOAT,processor,0,MPI_COMM_WORLD,&requests[processor]);
          MPI_Wait(&requests[processor],&status[processor]);
         
         
         }
       }
       
       //computing normalisation for processor 0
        
        for (row = norm+1; row < N; row+= processors) {
        multiplier = A[row][norm] / A[norm][norm];
        for (col = norm; col < N; col++) {
	     A[row][col] -= A[norm][col] * multiplier;
         }
         B[row] -= B[norm] * multiplier;
         }
         
          
      //recieving back the computed matrices portions from each processor
       for(processor=1; processor < processors ; processor++) {
        for(k = norm + 1 + processor ; k < N ; k+= processors){
       
          MPI_Recv((void *)&A[k][0],N,MPI_FLOAT,processor,0,MPI_COMM_WORLD,&status[rank]); 
          MPI_Recv((void *)&B[k],N,MPI_FLOAT,processor,0,MPI_COMM_WORLD,&status[rank]); 
        } 
       } 
       
       
             
     }
     
     else
     {
          
        // recieving the portion of rows in the matrix based on static interleaving
        for(k = norm + 1 + rank; k < N ; k+= processors){
          
         MPI_Recv((void *)&A[k],N,MPI_FLOAT,0,0,MPI_COMM_WORLD,&status[rank]); 
         MPI_Recv((void *)&B[k],N,MPI_FLOAT,0,0,MPI_COMM_WORLD,&status[rank]); 
         multiplier = A[k][norm] / A[norm][norm];
         for (col = norm; col < N; col++) {
	     A[k][col] -= A[norm][col] * multiplier;
	    
          }
          B[k] -= B[norm] * multiplier;
          MPI_Isend((void *)&A[k][0],N,MPI_FLOAT,0,0,MPI_COMM_WORLD,&requests[rank]);
          MPI_Wait(&requests[rank],&status[rank]);
          MPI_Isend((void *)&B[k],1,MPI_FLOAT,0,0,MPI_COMM_WORLD,&requests[rank]);
          MPI_Wait(&requests[rank],&status[rank]);
        }
  
       
     }
    
    
    
    // barrier at the end of each loop to ensure synchronization
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

void backsub()
{
int nrow,ncol;
 
  for (nrow = N - 1; nrow >= 0; nrow--) {
    X[nrow] = B[nrow];
    for (ncol = N-1; ncol > nrow; ncol--) {
      X[nrow] -= A[nrow][ncol] * X[ncol];
    }
    X[nrow] /= A[nrow][nrow];
  }


}


void print_X() {
  int row;

  if (N < 100) {
    printf("\nX = [");
    for (row = 0; row < N; row++) {
      printf("%5.2f%s", X[row], (row < N-1) ? "; " : "]\n");
    }
  }
}

void main(int arg,char** args)
{
 
  /* Timing variables */
  struct timeval etstart, etstop;  /* Elapsed times using gettimeofday() */
  double starttime,endtime; /* variables needed to calculated MPI_Wtime()*/
  struct timezone tzdummy;
  struct tms cputstart, cputstop;  /* CPU times for my processes */

   MPI_Init(&arg,&args);
   MPI_Comm_size(MPI_COMM_WORLD,&processors);
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    if(rank == 0)
    {  
      parameters(arg, args);  // initialising the parameters in processor 0
      initialize_inputs();    // initialising the input matrices 
      print_inputs();        // printing input matrices 
      gettimeofday(&etstart, &tzdummy);  // initialising date and time 
      times(&cputstart);
      starttime = MPI_Wtime();           // setting the time variable from MPI_Wtime()
    }
   
    gauss();                 // calling the gauss method outside processor 0
   
   
    MPI_Barrier(MPI_COMM_WORLD);  // ensures synchronization before backsubstitution
   
    if(rank == 0)
   {
     gettimeofday(&etstop, &tzdummy);
     times(&cputstop);
     endtime = MPI_Wtime();
     backsub();                  //back substitution
     print_X();                  // prints X matrix
     

      printf("\nElapsed time %g s.\n",(endtime-starttime));
      printf("(CPU times are accurate to the nearest %g ms)\n",
	 1.0/(float)CLOCKS_PER_SEC * 1000.0);
      printf("My total CPU time for parent = %g ms.\n",
	 (float)( (cputstop.tms_utime + cputstop.tms_stime) -
		  (cputstart.tms_utime + cputstart.tms_stime) ) /
	 (float)CLOCKS_PER_SEC * 1000);
     printf("My system CPU time for parent = %g ms.\n",
	 (float)(cputstop.tms_stime - cputstart.tms_stime) /
	 (float)CLOCKS_PER_SEC * 1000);
     printf("My total CPU time for child processes = %g ms.\n",
	 (float)( (cputstop.tms_cutime + cputstop.tms_cstime) -
		  (cputstart.tms_cutime + cputstart.tms_cstime) ) /
	 (float)CLOCKS_PER_SEC * 1000);
      /* Contrary to the man pages, this appears not to include the parent */
     printf("--------------------------------------------\n");
   }
   
   
   MPI_Finalize();

}
