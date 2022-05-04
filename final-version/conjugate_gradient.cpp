#include <iostream>
using namespace std;

#include <math.h>
#include <cmath>
#include <fstream>
#include <assert.h>

#include <chrono>
using namespace std::chrono;

#include <mpi.h>

// Dimension of node arrangement in the MPI grid
#define DIMENSION 1

bool check_array_elements(double *array, int N) {
    for (int i = 0; i < N; i++) {
        // char element[100];
        // int dummy_var = sprintf(element, "%s", array[i]);
        // cout << element << " ";
        // if (element ==  "(null)") {
        //     cout << "NaN detected at index " << i << endl;
        //     return false;
        // }

        if (!array[i]) {
            cout << "NaN detected at index " << i << endl;
            return false;
        } else {
            cout << "not nan" << endl;
            cout << "arr[i]: " << array[i] << endl;
        }

        // cout << "arr[i]: " << array[i] << endl;
    }
    return true;
}

int global_start (int proc, int N_local) {
    return proc * N_local;
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
            if (j == n - 1) {
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
    // For processor 0
    if (mype == 0) {
        // Processor 0's contribution to the global output
        double global_out[N_global];
        for (int i = 0; i < N_local; i++) {
            global_out[i] = arr_to_collect[i];
        }

        // Every other processor's contribution to the global output
        double local_buffer[N_local];
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
        MPI_Send(arr_to_collect, N_local, MPI_DOUBLE, 0, 0, comm1d);
    }
}

void distribute_ghost_cells(double *arr_to_distribute, double *left_ghost_cells, double *right_ghost_cells, int N, int left, int right, MPI_Comm comm1d) {
    // Send and receive ghost cell
    MPI_Status status;
    MPI_Sendrecv(arr_to_distribute, N, MPI_DOUBLE, left, 99, &right_ghost_cells, N, MPI_DOUBLE, right, MPI_ANY_TAG, comm1d, &status);
    MPI_Sendrecv(arr_to_distribute, N, MPI_DOUBLE, right, 99, &left_ghost_cells, N, MPI_DOUBLE, left, MPI_ANY_TAG, comm1d, &status);
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
    int N_global = N * nprocs;
    int global_start_index = global_start(mype, N);
    int n = sqrt(N);
    int n_global = sqrt(N_global);
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++ ) {
            int local_1d = i*n + j;
            int global_1d = global_start_index + local_1d;
            int global_i = floor(global_1d / n_global);
            int global_j = global_1d % n_global;
            b[local_1d] = find_b(global_i, global_j, n_global);
            // b[i*n + j] = find_b(i, j, n);
        }
    }
}

void initialize_array(double *array, int N) {
    for(int i = 0; i < N; i++) {
        array[i] = 0.0;
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

void poisson_on_the_fly(double *v, double *w, int N, int mype, int nprocs, int left, int right, MPI_Comm comm1d, bool p_flag) {
    int N_global = N * nprocs;
    int n_global = sqrt(N_global);
    int n = sqrt(N);
    int global_start_index = global_start(mype, N);
    int global_i;

    // ghost cells
    double left_ghost_cells[N];
    double right_ghost_cells[N];
    if (p_flag) {
        initialize_array(left_ghost_cells, N);
        initialize_array(right_ghost_cells, N);
        distribute_ghost_cells(w, left_ghost_cells, right_ghost_cells, N, left, right, comm1d);
    }

    for (int i = 0; i < N; i++) {
        // Terms to compute vector values "on the fly"
        double t1 = 0.0;
        double t2 = 0.0;
        double t3 = 0.0;
        double t4 = 0.0;
        double t5 = 0.0;

        global_i = global_start_index + i;

        t3 = 4*w[i];
        if (i - n >= 0) {
            t1 = w[i - n];
        }
        else if (global_i - n >= 0 && p_flag) {
            // get ghost cell w[global_i - n_global]
            t1 = left_ghost_cells[N + i - n];
        }
        if (i - 1 >= 0) {
            t2 = w[i - 1];
        }
        else if (global_i - 1 >= 0 && p_flag) {
            // get ghost cell w[global_i - 1]
            t2 = left_ghost_cells[N + i - 1];
        }
        if (i + 1 < N) {
            t4 = w[i + 1];
        }
        else if (global_i + 1 < N_global && p_flag) {
            // get ghost cell w[global_i + 1]
            t4 = right_ghost_cells[i + 1 - N];
        }
        if (i + n < N) {
            t5 = w[i + n];
        }
        else if (global_i + n < N_global && p_flag) {
            // get ghost cell w[global_i + n_global]
            t5 = right_ghost_cells[i + n - N];
        }
        v[i] = t3 - t1 - t2 - t4 - t5;
    }
}

void conjugate_gradient(double *b, double *x, int n, int mype, int nprocs, int left, int right, MPI_Comm comm1d, bool p_flag) {
    // Initialize variables
    int N = n * n;
    int N_global = N * nprocs;
    
    // Tolerance for the solution to stop after convergence
    double tol = 1e-10;

    // Write RHS to file
    // write_to_file(b, "./output/b.txt", 0, N); // serial
    collect_and_write_array(b, "./output/b.txt", 0, mype, N_global, nprocs, comm1d); // parallel
    // return ;

    double r[N];
    double p[N];
    double z[N];

    // Temporary variables
    double Ax[N]; 

    // r = b - Ax
    poisson_on_the_fly(Ax, x, N, mype, nprocs, left, right, comm1d, p_flag);
    axpy(r, 1.0, b, -1.0, Ax, N);

    // p = r
    for (int i = 0; i < N; i++) {
        p[i] = r[i];
    }
    
    // rsold = rT * r
    double rsold = parallel_dotp(r, r, N, comm1d);

    for (int i = 1; i < N + 1; i++) {
        // Start timer
        auto ts = high_resolution_clock::now();

        // z = A*p
        poisson_on_the_fly(z, p, N, mype, nprocs, left, right, comm1d, p_flag);      

        // alpha = rsold / (p*z)
        double alpha = rsold / parallel_dotp(p, z, N, comm1d);
        
        // x = x + alpha*p
        axpy(x, 1.0, x, alpha, p, N);

        // r = r - alpha*z
        axpy(r, 1.0, r, -alpha, z, N);
        
        // rsnew = rT.r
        double rsnew = parallel_dotp(r, r, N, comm1d);

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
    }
    // Write solution to file
    collect_and_write_array(x, "./output/x.txt", N-1, mype, N_global, nprocs, comm1d);
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
    periodic[0] = 0;

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
    bool p_flag = false;
    if (solver_type == "parallel"){
        p_flag = true;
    }

    // Check whether the physical domain can be divided evenly among the MPI ranks
    assert(n_global % nprocs == 0);

    int N_global = n_global * n_global;
    int N_local = N_global / nprocs;
    int n_local = sqrt(N_local);

    
    // Source vector
    double b_local[N_local];
    fill_b(b_local, N_local, mype, nprocs);

    // Result vector
    double x_local[N_local];
    initialize_array(x_local, N_local);
    
    if (mype == 0) {
        cout << "Simulation Parameters:" << endl;
        cout << "Solver type = " << solver_type << endl;
        cout << "n = " << n_global << "\n" << endl;
        cout << "Estimated memeory usage = " << 5*N_global*sizeof(double)/1e6 << " MB" << "\n" << endl;
    }
    
    // Start timer
    auto s = high_resolution_clock::now();

    // Run simulation
    conjugate_gradient(b_local, x_local, n_local, mype, nprocs, left, right, comm1d, p_flag);

    // Stop timer
    auto e = high_resolution_clock::now();

    // Print time taken
    auto duration = duration_cast<microseconds>(e - s);
    // cout << "\nTime taken = " << 1e-6*duration.count() << " seconds" << endl;
    // cout << "\nAverage Grind Rate = " << int(n_global*n_global/(1e-6*duration.count())) << " iter/sec" << endl;

    // MPI finalization
    MPI_Finalize();

}