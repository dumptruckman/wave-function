#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#include <math.h>

//#define M_PI 3.14159265359
#define MASTER 0
#define Init(x, y) MPI_Init(x, y)
#define Abort(x) MPI_Abort(MPI_COMM_WORLD, x)
#define Rank(x) MPI_Comm_rank(MPI_COMM_WORLD, x)
#define Size(x) MPI_Comm_size(MPI_COMM_WORLD, x)
#define Barrier() MPI_Barrier(MPI_COMM_WORLD)
#define Finalize() MPI_Finalize()
#define Scatter(x, y, z) MPI_Scatter(x, z, MPI_DOUBLE, y, z, MPI_DOUBLE, MASTER, MPI_COMM_WORLD)
#define Broadcast(x) MPI_Bcast(x, 1, MPI_INT, MASTER, MPI_COMM_WORLD)
#define ReduceSum(x, y) MPI_Reduce(x, y, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD)
#define ReduceMax(x, y) MPI_Reduce(x, y, 1, MPI_DOUBLE, MPI_MAX, MASTER, MPI_COMM_WORLD)
#define SendRecv(send, dest, recv, source, tag) MPI_Sendrecv(send, 1, MPI_DOUBLE, dest, tag, recv, 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE)
#define Time() MPI_Wtime()

// functions defined after main
void printArray(double *array, int length);
void printlnArray(double *array, int length);
void copyArray(double *source, double *dest, int length);
double * createRange(int n, int commSize, int rank);
double ** allocateResultSpace(int n);
double fInit(double **result, int j, int i, int n);
double f(double **result, int j, int i, int n);

void main() {
  int commSize;
  int myRank;

  Init(NULL, NULL);
  Size(&commSize);
  Rank(&myRank);

  if(commSize % 2 != 0){
    Abort(1);
  }

  int n = 6, end = 4, length = n/commSize;
  double *range = createRange(n); // creates our starting range
  double **result = allocateResultSpace(length); // allocates space for our results
  double *swap; // serves as a temporary pointer swap

  // fill the first 2 result rows with the initial values
  result[0] = range;
  copyArray(range, result[1], length);

  printf("Initial Values of Rank %d: ", rank);
  printlnArray(range, n);

  int i, j;
  // Generate results for the first 2 steps
  for (j = 0; j < 2; j++) {
    for (i = 0; i < length; i++) {
      result[j][i] = isArrayEnd(i, rank, commSize, length) == 1 ? 0 : fInit(result, j, i);
    }
    printlnArray(result[j], n);
  }

  double leftNeighborValue
  double rightNeighborValue;
  int leftNeighborRank;
  int rightNeighborRank;
  leftNeighborRank = rank==0? MPI_PROC_NULL: rank-1;
  rightNeighborRank = rank==commSize-1? MPI_PROC_NULL: rank+1;
  // Generate results from step 2 to step end.
  for (j = 2; j <= end; j++) {
    for (i = 0; i < n; i++) {
      result[2][i] = f(result, 2, i, n);
    }
    printlnArray(result[2], n);

    // Swap the pointers around save past 2 results only.
    if (j > 1) {
      swap = result[0];
      result[0] = result[1];
      result[1] = result[2];
      result[2] = swap;
    }
  }

  Finalize();
}

void printArray(double *array, int length) {
  printf("[");
  int i;
  for (i = 0; i < length; i++) {
    if (i != 0) {
      printf("  ");
    }
    printf("%g", array[i]);
  }
  printf("]");
}

void printlnArray(double *array, int length) {
  printArray(array, length);
  printf("\n");
}

void copyArray(double *source, double *dest, int length) {
  int i;
  for (i = 0; i < length; i++) {
    dest[i] = source[i];
  }
}

double * createRange(int n, int commSize, int rank) {
  int length = n/commSize;
  double *range = malloc((n/commSize) * sizeof(double));
  int i;
  for (i = 0; i < length; i++) {
    range[i] = (double) (i + rank*length) / (n - 1);
  }
  return range;
}

double ** allocateResultSpace(int n) {
  double **result = malloc(3 * sizeof(double));
  int i;
  for (i = 1; i <= 2; i++) {
    result[i] = calloc(n, sizeof(double));
  }
  return result;
}

double isArrayEnd(int i, int rank, int commSize, int length) {
  if ((i == 0 && rank == 0) || (rank == commSize - 1 && i == length - 1)) {
    return 1;
  }
  return 0;
}

double fInit(double **result, int j, int i) {
  return sin(M_PI * result[j][i]);
}

double f(double **result, int j, int i) {
  return 0.01 * (result[j-1][i-1] - 2 * result[j-1][i] + result[j-1][i+1]) + 2 * result[j-1][i] - result[j-2][i];
}
