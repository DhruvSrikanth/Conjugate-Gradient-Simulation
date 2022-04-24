#include "cg_header.h"

int main(int argc, char* argv[]){
  // n is the physical domain size;
  int n = 75;
  printf("Solving Poisson Equation on %d x %d domain...\n", n, n);
  // Dimension of operator matrix and vectors is n^2
  int N = n*n;
  // Allocate full A matrix and vectors
  double ** A = matrix( N );
  double *  x = (double*) calloc(N, sizeof(double));
  double *  b = (double*) calloc(N, sizeof(double));
  printf("Dense Memory  = %.2lf MB\n", (N*N+5*N)*sizeof(double)/1024.0/1024.0);
  
  // Compute elements of 'A' matrix (Poisson Operator)
  fill_A(A, N);
  
  // Compute elements of boundary condition vector 'b'
  fill_b(b, N);
  
  // Run Dense CG Solve
  double start = get_time();
  cg_dense(A, x, b, N);
  double stop = get_time();
  printf("Dense Runtime = %.2lf seconds\n", stop-start);
  
  // Save Solution Vector to File
  save_vector(x,N, "dense.out");
  
  // Free A matrix
  matrix_free(A);
  
  // Free vectors
  free(x);
  free(b);
}
