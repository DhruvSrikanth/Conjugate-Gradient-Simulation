#include <stdio.h>
using namespace std;

#include <iostream>
using namespace std;

#include <algorithm>
#include <math.h>
#include<assert.h>
#include <fstream>

#include <chrono>
using namespace std::chrono;


void write_to_file(double* result, string const& filename, int const& n_iter, int const& N) { 
    // Allocate memory for the file
    ofstream file;
    file.open(filename);

    // Compute dimensions of the 1D array to represent as 2D array
    int n = sqrt(N);

    // Write the timestep to the file
    file << "[" << n_iter << "], ";

    // Write the data to the file
    file << "[";
    for (int i = 0; i < n; i++) {
        if (i == 0) {
            file << "[";
        } 
        else {
            file << ", [";
        }
        for (int j = 0; j < n; j++) {
            double res = result[i*n + j];
            if (j == N - 1) {
                file << res;
            } 
            else {
                file << res << ", ";
            }
        }
        file << "]";
    }
    file << "]";

    // Release the memory for the file
    file.close();
}

double find_b(int i, int j, int n) {
    double delta = 1.0 / double(n);

    double x = -0.5 + delta + delta * j;
    double y = -0.5 + delta + delta * i;

    // Check if within a circle
    double radius = 0.1;
    if ( x*x + y*y < radius*radius ) {
        return delta * delta / 1.075271758e-02;
    }
    else {
        return 0.0;
    }
}

void fill_b(double* b, int const& N) {
    int n = sqrt(N);
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++ ) {
            b[i*n + j] = find_b(i,j,n);
        }
    }
}

void axpy(double* res, double const& a, double* const& x, int const& b, double* const& y, int const& N) {
    // Perform res = a * x + b * y
    int n = sqrt(N);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int idx = i*n + j;
            res[idx] = a * x[idx] + b * y[idx];
        }
    }
}

double dot(double* const& x, double* const& y, int const& N) {
    // Perform res = x * y
    int n = sqrt(N);
    double res = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int idx = i*n + j;
            res += x[idx] * y[idx];
        }
    }
    return res;
}

double* poisson_on_the_fly(double* const& w, int const& N) {
    double t1 = 0.0;
    double t2 = 0.0;
    double t3 = 0.0;
    double t4 = 0.0;
    double t5 = 0.0;
    int n = sqrt(N);
    double* res = new double[N];

    for (int i = 0; i < N; i++) {
        t3 = 4*w[i];

        if (i < n) {
            t1 = 0.0;
        } else {
            t1 = w[i - n];
        }

        if (i < 1) {
            t2 = 0.0;
        } else {
            t2 = w[i - 1];
        }

        if (i >= N - 1) {
            t4 = 0.0;
        } else {
            t4 = w[i + 1];
        }

        if (i >= N - n) {
            t5 = 0.0;
        } else {
            t5 = w[i + n];
        }

        res[i] = t3 - t1 - t2 - t4 - t5;
    }

    return res;
}

void conjugate_gradient(double* const& b, double* x, int const& n) {
    // Initialize variables
    int N = n * n;
    // Tolerance for the solution to stop after convergence
    double tol = 1e-10;

    // Write RHS to file
    write_to_file(b, "./output/b.txt", 0, N);

    double* r = new double[N];
    // r = -Ax + b
    double* Ax = poisson_on_the_fly(x, N);
    axpy(r, -1.0, Ax, 1.0, b, N);


    double* p = new double[N];
    // p = r
    copy(r, r + N*sizeof(double), p);
    
    double* z = new double[N];    

    // rsold = rT * r
    double rsold = dot(r, r, N);

    // Run simulation loop
    for (int i = 0; i < N; i++) {

        // Start timer
        auto ts = high_resolution_clock::now();

        // z = A*p
        z = poisson_on_the_fly(p, N);

        // alpha = rsold / (p*z)
        double alpha = rsold / dot(p, z, N);

        // x = x + alpha*p
        axpy(x, 1.0, x, alpha, p, N);
        

        // r = r - alpha*z
        axpy(r, 1.0, r, -alpha, z, N);

        // rsnew = rT*r
        double rsnew = dot(r, r, N);

        // If the residual is small enough, stop
        if (sqrt(rsnew) <= tol) {
            cout << "Converged after " << i << " iterations" << endl;
            break;
        }

        // p = r + rsnew / rsold * p
        axpy(p, 1.0, r, rsnew/rsold, p, N);

        rsold = rsnew;

        // Stop timer
        auto te = high_resolution_clock::now();

        // Print time taken
        auto duration = duration_cast<microseconds>(te - ts);
        cout << "Iteration: " << i << " - Grind Rate: " << int(1/(1e-6*duration.count())) << " iter/sec" << endl;

        // if (i % 100 == 0) {
        //     // Write solution to file
        //     char filename[100];
        //     int dummy_var = sprintf(filename, "./output/output_x_%d.txt", i);
        //     write_to_file(x, filename, i, N);
        // }
    }

    // Write solution to file
    write_to_file(x, "./output/x.txt", N, N);
}

int main(int argc, char** argv) {

    // Initialize variables
    int n = stoi(argv[1]);
    int N = n * n;

    double* b = new double[N];
    fill_b(b, N);

    // Result vector
    double* x = new double[N];

    cout << "Simulation Parameters:" << endl;
    cout << "n = " << n << "\n" << endl;
    cout << "Estimated memeory usage = " << 5*N*sizeof(double)/1e6 << " MB" << "\n" << endl;

    // Start timer
    auto s = high_resolution_clock::now();

    // Run simulation
    conjugate_gradient(b, x, n);

    // Stop timer
    auto e = high_resolution_clock::now();

    // Print time taken
    auto duration = duration_cast<microseconds>(e - s);
    cout << "\nTime taken = " << 1e-6*duration.count() << " seconds" << endl;

}