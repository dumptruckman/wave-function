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
#define Gather(send, count, recv) MPI_Gather(send, count, MPI_DOUBLE, recv, count, MPI_DOUBLE, MASTER, MPI_COMM_WORLD)
#define Time() MPI_Wtime()

// functions defined after main
void printArray(double *array, int length);
void printlnArray(double *array, int length);
void copyArray(double *source, double *dest, int length);
double * createRange(int n, int commSize, int rank);
double ** allocateResultSpace(int n);
double isArrayEnd(int i, int rank, int commSize, int length);
double fInit(double **result, int j, int i);
double f(double **result, int j, int i, double stepSize, int length, double leftNeighborValue, double rightNeighborValue);

void main() {
  int commSize;
  int rank;

  Init(NULL, NULL);
  Size(&commSize);
  Rank(&rank);

  if(commSize % 2 != 0) {
    printf("commSize must be divisible by 2!\n");
    Abort(1);
  }

  int n = 20, end = 8;
  double stepSize = 1;
  if (n % commSize != 0) {
    printf("n must be divisible by commSize!\n");
    Abort(1);
  }
  int length = n / commSize;
  double *range = createRange(n, commSize, rank); // creates our starting range
  double **result = allocateResultSpace(length); // allocates space for our results
  double *swap; // serves as a temporary pointer swap

  double *combined;
  if (rank == MASTER) {
    combined = malloc(n * sizeof(double));
  }

  // fill the first 2 result rows with the initial values
  result[0] = range;
  copyArray(range, result[1], length);

  Gather(result[0], length, combined);
  if (rank == MASTER) {
    //printf("Step %d Rank %d :", step, rank);
    printf("Initial Values: ");
    printlnArray(combined, n);
  }
  //printf("Initial Values of Rank %d: ", rank);
  //printlnArray(range, length);

  int i, j, step;
  // Generate results for the first 2 steps
  for (j = 0; j < 2; j++) {
    for (i = 0; i < length; i++) {
      result[j][i] = isArrayEnd(i, rank, commSize, length) == 1 ? 0 : fInit(result, j, i);
    }
    Gather(result[j], length, combined);
    if (rank == MASTER) {
      //printf("Step %d Rank %d :", step, rank);
      printlnArray(combined, n);
    }
    //printf("Step %d Rank %d :", j, rank);
    //printlnArray(result[j], length);
  }

  double leftNeighborValue;
  double rightNeighborValue;
  int leftNeighborRank;
  int rightNeighborRank;
  leftNeighborRank = rank == 0 ? MPI_PROC_NULL : rank - 1;
  rightNeighborRank = rank == commSize - 1 ? MPI_PROC_NULL : rank + 1;
  // Generate results from step 2 to step end
  for (j = 2, step = 2; j <= end; j++, step++) {
    SendRecv(&result[1][length - 1], rightNeighborRank, &leftNeighborValue, leftNeighborRank, step);
    SendRecv(&result[1][0], leftNeighborRank, &rightNeighborValue, rightNeighborRank, step);
    for (i = 0; i < length; i++) {
      result[2][i] = isArrayEnd(i, rank, commSize, length) == 1 ? 0 : f(result, 2, i, stepSize, length, leftNeighborValue, rightNeighborValue);
    }
    Gather(result[2], length, combined);
    if (rank == MASTER) {
      //printf("Step %d Rank %d :", step, rank);
      printlnArray(combined, n);
    }

    // Rotate the rows to keep the last 2 calculations
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

double f(double **result, int j, int i, double stepSize, int length, double leftNeighborValue, double rightNeighborValue) {
  return stepSize * (
    (i == 0 ? leftNeighborValue : result[j-1][i-1])
    - 2 * result[j-1][i]
    + (i == length - 1 ? rightNeighborValue : result[j-1][i+1])
  ) + 2 * result[j-1][i] - result[j-2][i];
}
