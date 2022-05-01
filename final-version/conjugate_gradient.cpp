#include <iostream>
using namespace std;

#include <math.h>
#include <fstream>

#include <chrono>
using namespace std::chrono;

#include <mpi.h>

#define DIMENSION 1

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

void collect_and_write(int const& mype, double* arr_to_send, int const& N_global, int const& nprocs, MPI_Comm const& comm1d, string const& filename, int const& n_iter) {
    // write b out
    int N_local = N_global / nprocs;
    double* local_buffer = new double[N_local];
    

    if (mype == 0) {
        write_to_file(arr_to_send, "./output/b_0.txt", 0, N_local);
    }
    for (int i = 0; i < N_local; i++) {
        if (arr_to_send[i] != 0.0) {
            cout << "arr_to_send[" << i << "] = " << arr_to_send[i] << endl;
            cout << "mype = " << mype << endl;
        }
    }
        

    // For processor 0
    if (mype == 0) {
        // Processor 0's contribution to the global output
        double* global_out = new double[N_global];
        for (int i = 0; i < N_global; i++) {
            global_out[i] = 0.0;
        }
        
        int global_start = N_local * mype;
        int local_i = 0;
        for (int i = global_start; i < global_start + N_local; i++) {
            local_i = i - global_start;
            global_out[i] = arr_to_send[local_i];
        }

        // Every other processor's contribution to the global output
        for (int x = 1; x < nprocs; x++) {
            MPI_Recv(local_buffer, N_local, MPI_DOUBLE, x, 0, comm1d, MPI_STATUS_IGNORE);
            for (int i = global_start; i < global_start + N_local; i++) {
                local_i = i - global_start;
                global_out[i] = local_buffer[local_i];
            }
        }
        // Write RHS to file
        write_to_file(global_out, filename, n_iter, N_global);
        
    }

    // For other processors
    else {
        // Initialize local output as a contiguous array to send
        double* send_buffer = new double[N_local];
        int local_i = 0;
        int global_start = N_local * mype;
        for (int i = global_start; i < global_start + N_local; i++) {
                local_i = i - global_start;
                send_buffer[local_i] = arr_to_send[local_i];
        }
        MPI_Send(send_buffer, N_local, MPI_DOUBLE, 0, 0, comm1d);
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

void fill_b(double* b, int const& N_local, int const& mype, int const& nprocs) {
    int n = sqrt(N_local);
    int N_global = N_local * nprocs;
    // start - N/nprocs * mype
    // end - N/nprocs * (mype + 1)
    int global_start = N_local * mype;
    int global_i = 0;
    int global_j = 0;
    for(int i = 0; i < n; i++) {
        global_i = global_start + i;
        for(int j = 0; j < n; j++ ) {
            global_j = j;
            b[i*n + j] = find_b(global_i,global_j,n);
        }
    }
}

double dotp(double* const& x, double* const& y, int const& N) {
    // Perform res = x * y
    double res = 0.0;
    for (int i = 0; i < N; i++) {
        res += x[i] * y[i];
    }
    return res;
}

double* poisson_on_the_fly(double* const& w, int const& N) {
    double* v = new double[N];
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
    return v;
}

void vec_scale(double* v, double const& alpha, int const& N) {
    for (int i = 0; i < N; i++) {
        v[i] *= alpha;
    }
}

void vec_add(double* output, double* const& v1, double* const& v2, double const& op, int const& N) {
    for (int i = 0; i < N; i++) {
        output[i] = v1[i] + op * v2[i];
    }
}

void axpy(double* output, double* const& v1, double const& a, double* const& v2, double const& b, int const& N) {
    for (int i = 0; i < N; i++) {
        output[i] = (a * v1[i]) + (b * v2[i]);
    }
}

void conjugate_gradient(double* const& b, double* & x, int const& n, int const& mype, int const& nprocs, MPI_Comm const& comm1d) {
    // Initialize variables
    int N = n * n;
    int N_local = N / nprocs;

    // Tolerance for the solution to stop after convergence
    double tol = 1e-10;
    collect_and_write(mype, b, N, nprocs, comm1d, "./output/b.txt", 0);
    

    return ;

    
    double* r = new double[N];
    double* p = new double[N];
    double* z = new double[N];

    // Temporary variables
    double* Ax = poisson_on_the_fly(x, N);

    // r = b - Ax
    // vec_add(r, b, Ax, -1.0, N);
    axpy(r, b, 1.0, Ax, -1.0, N);

    // Free Ax since we wont use it after this
    free(Ax);

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
        z = poisson_on_the_fly(p, N);

        // alpha = rsold / (p*z)
        double alpha = rsold / dotp(p, z, N);

        // x = x + alpha*p
        vec_scale(p, alpha, N);
        vec_add(x, x, p, 1.0, N);
        // axpy(x, x, 1.0, p, alpha, N);

        // r = r - alpha*z
        // vec_scale(z, alpha, N);
        // vec_add(r, r, z, -1.0, N);
        axpy(r, r, 1.0, z, -1.0*alpha, N);
        
        // rsnew = rT*r
        double rsnew = dotp(r, r, N);

        // If the residual is small enough, stop
        if (sqrt(rsnew) <= tol) {
            cout << "Converged after " << i << " iterations" << endl;
            break;
        }

        // p = r + rsnew / rsold * p
        // vec_scale(p, rsnew / rsold, N);
        // vec_add(p, r, p, 1.0, N);
        axpy(p, r, 1.0, p, (double) rsnew/rsold, N);

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

    // MPI Cartesian Grid Creation
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
    int n = stoi(argv[1]);
    int N = n * n;
    int N_local = N / nprocs;

    double* b = new double[N_local];
    fill_b(b, N_local, mype, nprocs);

    for (int i = 0; i < N_local; i++) {
        if (b[i] != 0.0) {
            cout << "arr_to_send[" << i << "] = " << b[i] << endl;
            cout << "mype = " << mype << endl;
        }
    }


    // Result vector
    double* x = new double[N_local];
    
    for (int i = 0; i < N_local; i++) {
        x[i] = 0.0;
    }

    cout << "Simulation Parameters:" << endl;
    cout << "n = " << n << "\n" << endl;
    cout << "Estimated memeory usage = " << 5*N*sizeof(double)/1e6 << " MB" << "\n" << endl;

    // Start timer
    auto s = high_resolution_clock::now();

    // Run simulation
    conjugate_gradient(b, x, n, mype, nprocs, comm1d);

    // Stop timer
    auto e = high_resolution_clock::now();

    // MPI finalization
    MPI_Finalize();

    // Print time taken
    auto duration = duration_cast<microseconds>(e - s);
    cout << "\nTime taken = " << 1e-6*duration.count() << " seconds" << endl;
}