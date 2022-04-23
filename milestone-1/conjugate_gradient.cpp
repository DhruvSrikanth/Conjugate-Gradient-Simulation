#include <stdio.h>
using namespace std;

#include <iostream>
using namespace std;

#include <algorithm>

#include <chrono>
using namespace std::chrono;

int contiguous_memory_index(int const& i, int const& j, int const& num_cols) {
    // Calculate the index of the contiguous memory matrix
    return j + (num_cols * i);
}

double* contiguous_memory_alloc(int const& n_rows, int const& n_cols) {
    // Allocate contiguous memory for the matrix
    double* mat = new double[n_rows * n_cols];
    // Initialize the matrix
    for (int i = 0; i < n_rows; i++) {
        for (int j = 0; j < n_cols; j++) {
            mat[contiguous_memory_index(i, j, n_cols)] = 0.0;
        }
    }
    return mat;
}

double* axpy(double const& a, double const* x, int const& b, double* y, int const& N) {
    // Perform res = a * x + b * y
}

double dot(double const* x, double* y, int const& N) {
    // Perform res = x * y
}

double* mat_vec_mult(double const* A, double const* x, int const& N) {
    // Perform res = A * x
}

double* conjugate_gradient(double* const& A, double* const& b, int const& n) {
    // Initialize variables
    int N = n * n;

    // Allocated memory for the following arrays follow a column major order
    double* x = contiguous_memory_alloc(N, 1); 
    double* r = contiguous_memory_alloc(N, 1);
    double* p = contiguous_memory_alloc(N, 1);
    double* z = contiguous_memory_alloc(N, 1);

    double alpha;
    double rsnew;

    // Tolerance for the solution to stop after convergence
    double tol = 1e-10;

    // r = b - Ax
    double* Ax = mat_vec_mult(A, x, N);
    r = axpy(1, b, -1, Ax, N);

    // p = r
    p = copy(r, r+N, p);

    // rsold = rT * r
    double rsold = dot(r, r, N);

    // Run simulation loop
    for (int i = 0; i < N; i++) {
        // z = A*p
        z = mat_vec_mult(A, p, N);

        // alpha = rsold / (p*z)
        alpha = rsold / dot(p, z, N);

        // x = x + alpha*p
        x = axpy(1, x, alpha, p, N);

        // r = r - alpha*z
        r = axpy(1, r, -alpha, z, N);
        
        // rsnew = rT*r
        rsnew = dot(r, r, N);

        // If the residual is small enough, stop
        if (rsnew < tol) {
            break;
        }
        
        // p = r + rsnew / rsold * p
        p = axpy(1, r, rsnew / rsold, p, N);

        // rsold = rsnew
        rsold = rsnew;
    }

    return x;

}

int main(int argc, char** argv) {

    // Initialize variables
    int n = stoi(argv[1]);

    double* A = contiguous_memory_alloc(n*n, n*n);
    double* b = contiguous_memory_alloc(n*n, 1);

    cout << "Simulation Parameters:" << endl;
    cout << "n = " << n << "\n" << endl;

    // Esimate memory usage
    cout << "Estimated memeory usage = " << n*n*n*n*sizeof(double)/1e6 << " MB" << "\n" << endl;

    // Start timer
    auto t1 = high_resolution_clock::now();

    // Run simulation
    double* x = conjugate_gradient(A, b, n);
    cout << "x = " << x << "\n" << endl;

    // Stop timer
    auto t2 = high_resolution_clock::now();

    // Print time taken
    cout << "Time taken = " << 1e6*(duration_cast<microseconds>(t2 - t1)).count() << " seconds" << "\n" << endl;
}