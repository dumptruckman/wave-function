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
#define Time() MPI_Wtime()

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

double * createRange(int n) {
  double *range = malloc(n * sizeof(double));
  int i;
  for (i = 0.0; i < n; i++) {
    range[i] = (double) i / (n - 1);
  }
  return range;
}

double f_0(double *range, int i, int n) {
  if (i == 0 || i == n - 1) {
    return 0;
  }
  return sin(M_PI * range[i]);
}

void main() {

  int commSize;
  int myRank;

  Init(NULL, NULL);
  Size(&commSize);
  Rank(&myRank);

  int n = 6;
  double *range = createRange(n);

  printlnArray(range, n);

  int i;
  for (i = 0; i < n; i++) {
    range[i] = f_0(range, i, n);
  }

  printlnArray(range, n);

  Finalize();
}
