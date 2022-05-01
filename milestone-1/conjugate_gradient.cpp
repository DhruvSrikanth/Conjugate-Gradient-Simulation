#include <iostream>
using namespace std;

#include <math.h>
#include <fstream>

#include <chrono>
using namespace std::chrono;

void write_to_file(double *result, string filename, int n_iter, int N) { 
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

void fill_b(double *b, int N) {
    int n = sqrt(N);
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++ ) {
            b[i*n + j] = find_b(i,j,n);
        }
    }
}

double dotp(double *x, double *y, int N) {
    // Perform res = x * y
    double res = 0.0;
    for (int i = 0; i < N; i++) {
        res += (x[i] * y[i]);
    }
    return res;
}

void poisson_on_the_fly(double *v, double *w, int N) {
    double t1 = 0.0;
    double t2 = 0.0;
    double t3 = 0.0;
    double t4 = 0.0;
    double t5 = 0.0;
    int n = sqrt(N);
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

        v[i] = t3 - t1 - t2 - t4 - t5;
    }
}

void axpy(double *output, double alpha, double *v1, double beta, double *v2, int N) {
    for (int i = 0; i < N; i++) {
        v1[i] *= alpha;
        v2[i] *= beta;
        output[i] = v1[i] + v2[i];
    }
}

void conjugate_gradient(double *b, double *x, int n) {
    // Initialize variables
    int N = n * n;

    // Tolerance for the solution to stop after convergence
    double tol = 1e-10;

    // Write RHS to file
    write_to_file(b, "./output/b.txt", 0, N);

    double r[N];
    double p[N];
    double z[N];

    // Temporary variables
    double Ax[N]; 

    // r = b - Ax
    poisson_on_the_fly(Ax, x, N);
    axpy(r, 1.0, b, -1.0, Ax, N);

    // p = r
    for (int i = 0; i < N; i++) {
        p[i] = r[i];
    }

    // rsold = rT * r
    double rsold = dotp(r, r, N);

    for (int i = 1; i < N + 1; i++) {
        // Start timer
        auto ts = high_resolution_clock::now();

        // z = A*p
        poisson_on_the_fly(z, p, N);

        // alpha = rsold / (p*z)
        double alpha = rsold / dotp(p, z, N);

        // x = x + alpha*p
        axpy(x, 1.0, x, alpha, p, N);

        // r = r - alpha*z
        axpy(r, 1.0, r, -alpha, z, N);
        
        // rsnew = rT*r
        double rsnew = dotp(r, r, N);

        // If the residual is small enough, stop
        if (sqrt(rsnew) <= tol) {
            cout << "Converged after " << i << " iterations" << endl;
            break;
        }

        // p = r + rsnew / rsold * p
        axpy(p, 1.0, r, rsnew / rsold, p, N);

        // Stop timer
        auto te = high_resolution_clock::now();

        // Print time taken
        auto duration = duration_cast<microseconds>(te - ts);
        cout << "Iteration: " << i << " - Grind Rate: " << int(1/(1e-6*duration.count())) << " iter/sec" << endl;

        // For generating a movie
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

    double b[N];
    fill_b(b, N);

    // Result vector
    double x[N];
    for (int i = 0; i < N; i++) {
        x[i] = 0.0;
    }

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