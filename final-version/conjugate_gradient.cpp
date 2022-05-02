#include <iostream>
using namespace std;

#include <math.h>
#include <fstream>
#include <assert.h>

#include <chrono>
using namespace std::chrono;

#include <mpi.h>

// Dimension of node arrangement in the MPI grid
#define DIMENSION 1

int global_start (int mype, int N_local) {
    return mype * N_local;
}

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

void collect_and_write_array(double *arr_to_collect, string filename, int n_iter, int mype, int N_global, int nprocs, MPI_Comm comm1d) {
    int N_local = N_global / nprocs;

    double local_buffer[N_local];

    // For processor 0
    if (mype == 0) {
        // Processor 0's contribution to the global output
        double global_out[N_global];
        for (int i = 0; i < N_local; i++) {
            global_out[i] = arr_to_collect[i];
        }

        // Every other processor's contribution to the global output
        for (int proc = 1; proc < nprocs; proc++) {
            MPI_Recv(local_buffer, N_local, MPI_DOUBLE, proc, 0, comm1d, MPI_STATUS_IGNORE);
            int global_start_index = global_start(proc, N_local);
            for (int i = 0; i < N_local; i++) {
                global_out[global_start_index + i] = local_buffer[i];
            }
        }
        // Write to file
        write_to_file(global_out, filename, n_iter, N_global); 
    }

    // For other processors
    else {
        // send to processor 0
        for (int i = 0; i < N_local; i++) {
            local_buffer[i] = arr_to_collect[i];
        }
        MPI_Send(local_buffer, N_local, MPI_DOUBLE, 0, 0, comm1d);
    }
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

void fill_b(double *b, int N, int mype, int nprocs) {
    // TODO: Fill the b array with the CORRECT values
    int global_start_index = global_start(mype, N);
    int global_i;
    int global_j;
    int n_global = sqrt(N * nprocs);
    int n = sqrt(N);
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++ ) {
            b[i*n + j] = find_b(global_i, global_j, n_global);
        }
    }
}

void initialize_x(double *x, int N) {
    for(int i = 0; i < N; i++) {
        x[i] = 0.0;
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

double parallel_dotp(double *x, double *y, int N, MPI_Comm comm1d) {
    double global_dotp = 0.0;
    double local_dotp = 0.0;
    local_dotp = dotp(x, y, N);
    MPI_Allreduce(&local_dotp, &global_dotp, 1, MPI_DOUBLE, MPI_SUM, comm1d);
    return global_dotp;
}

void axpy(double *output, double alpha, double *v1, double beta, double *v2, int N) {
    for (int i = 0; i < N; i++) {
        v1[i] *= alpha;
        v2[i] *= beta;
        output[i] = v1[i] + v2[i];
    }
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
    // MPI initialization
    int mype, nprocs;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &mype);

    // MPI Cartesian Grid Creation (1D)
    int dims[DIMENSION], periodic[DIMENSION], coords[DIMENSION];
    MPI_Comm comm1d;
    dims[0] = nprocs;
    periodic[0] = 1;

    // Create Cartesian Communicator
    MPI_Cart_create(MPI_COMM_WORLD, DIMENSION, dims, periodic, 1, &comm1d);

    // Extract this MPI rank's N-dimensional coordinates from its place in the MPI Cartesian grid
    MPI_Cart_coords(comm1d, mype, DIMENSION, coords);

    // Determine 1D neighbor ranks for this MPI rank
    int left, right;
    MPI_Cart_shift(comm1d, 0, 1, &left, &right);
    // cout << "Rank " << mype << " has left neighbor " << left << " and right neighbor " << right << endl;

    // Initialize variables
    int n_global = stoi(argv[1]);
    string solver_type = argv[2];

    // Check whether the physical domain can be divided evenly among the MPI ranks
    assert(n_global % nprocs == 0);

    int N_global = n_global * n_global;
    int N_local = N_global / nprocs;
    int n_local = sqrt(N_local);

    double b_local[N_local];
    fill_b(b_local, N_local, mype, nprocs);

    // Save each processor's data to a file
    // char filename[100];
    // int dummy_var = sprintf(filename, "./output/mype_%d.txt", mype);
    // write_to_file(b_local, filename, 0, N_local);


    // Result vector
    double x_local[N_local];
    initialize_x(x_local, N_local);
    
    if (mype == 0) {
        cout << "Simulation Parameters:" << endl;
        cout << "Solver type = " << solver_type << endl;
        cout << "n = " << n_global << "\n" << endl;
        cout << "Estimated memeory usage = " << 5*N_global*sizeof(double)/1e6 << " MB" << "\n" << endl;
    }
    
    // Start timer
    auto s = high_resolution_clock::now();

    // Run simulation
    // conjugate_gradient(b_local, x_local, n_local);

    // Stop timer
    auto e = high_resolution_clock::now();

    // Print time taken
    auto duration = duration_cast<microseconds>(e - s);
    cout << "\nTime taken = " << 1e-6*duration.count() << " seconds" << endl;

    // MPI finalization
    MPI_Finalize();

}