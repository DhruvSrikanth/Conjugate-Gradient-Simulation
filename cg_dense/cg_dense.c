#include "cg_header.h"

// Serial Conjugate Gradient Solver Function for Ax = b
// A must be symmetric and positive definite
// A = 2D operator matrix
// x = 1D solution vector
// b = 1D vector
// N = dimension
void cg_dense(double ** A, double * x, double * b, long N)
{
    // r = -A*x + b
    double * r = (double *) malloc( N*sizeof(double));
    matvec(r, A, x, N);
    axpy(-1.0, r, 1.0, b, N);

    //p = r;
    double * p = (double *) malloc( N*sizeof(double));
    memcpy(p, r, N*sizeof(double));

    //rsold = r' * r;
    double rsold = dotp(r, r, N);

    // z
    double * z = (double *) malloc( N*sizeof(double));

    long iter = 0;
    for( iter = 0; iter < N; iter++ )
    {
        //z = A * p;
        matvec(z, A, p, N);

        //alpha = rsold / (p' * z);
        double alpha = rsold / dotp(p, z, N);

        //x = x + alpha * p;
        axpy(1.0, x, alpha, p, N);

        //r = r - alpha * z;
        axpy(1.0, r, -alpha, z, N);

        double rsnew = dotp(r,r,N);

        if( sqrt(rsnew) < 1.0e-10 )
            break;

        //p = (rsnew / rsold) * p + r;
        axpy(rsnew/rsold, p, 1.0, r, N);

        rsold = rsnew;
    }
    printf("CG converged in %ld iterations.\n", iter);

    free(r);
    free(p);
    free(z);
}

// General Matrix Vector Product for v = M * w
// v, w = 1D vectors
// M = 2D matrix
// N = dimension
void matvec( double * v, double ** M, double * w, long N )
{
    // Set solution vector to 0
    memset( v, 0, N*sizeof(double));

    for (long i = 0; i < N; i++)
        for (long j = 0; j < N; j++)
            v[i] += (M[i][j] * w[j]);
}

// Specific "On the fly" Matrix Vector Product for v = A * w
// where 'A' is the 2D Poisson operator matrix
// that is applied without explicit storage
// v, w = 1D vectors
// N = dimension

// Dot product of c = a * b
// c = result scalar that's returned
// a, b = 1D Vectors
// N = dimension
double dotp( double * a, double * b, long N)
{
    double c = 0.0;
    for( long i = 0; i < N; i++ )
        c += a[i]*b[i];

    return c;
}

// Scale and add of two vectors (axpy)
// Solves w = alpha * w + beta * v (overwrites what's in w)
// alpha, beta = scalars
// w, v = 1D vectors
// N = dimension
void axpy( double alpha, double * w, double beta, double * v, long N)
{
    for( long i = 0; i < N; i++ )
        w[i] = alpha * w[i] + beta * v[i];
}


// Fills an Explicit Fully Stored 2D Poisson Operator Matrix 'A'
// A = 2D Matrix
// N = dimension
void fill_A(double ** A, long N)
{
    long n = sqrt(N);

    for( long i = 0; i < N; i++ )
    {
        for( long j = 0; j < N; j++ )
        {
            // Compute x-coordinate of index
            long x = j % n;

            // main diagonal
            if(i == j)
                A[i][j] = 4.0;
            // left diagonal
            else if (i == j+1 )
            {
                A[i][j] = -1.0;
                if(x == n-1)
                    A[i][j] = 0.0;
            }
            // right diagonal
            else if (i == j-1 )
            {
                A[i][j] = -1.0;
                if(x == 0)
                    A[i][j] = 0.0;
            }
            // far left diagonal
            else if(j+n == i)
                A[i][j] = -1.0;
            // far right diagonal
            else if(j-n == i)
                A[i][j] = -1.0;
            // Otherwise, A is 0
            else
                A[i][j] = 0.0;
        }
    }
}

// Sets a circular source term for the right hand side 'b' vector
// Takes as input 2D spatial domain indices
// i, j = indices
// n = physical dimension
double find_b(long i, long j, long n)
{
    double delta = 1.0 / (double) n;

    double x = -.5 + delta + delta * j;
    double y = -.5 + delta + delta * i;

    // Check if within a circle
    double radius = 0.1;
    if( x*x + y*y < radius*radius )
        return delta * delta / 1.075271758e-02;
    else
        return 0.0;
}

// Fills a 1-D RHS vector specifying boundary conditions
// by calling the get_b method
// b = 1-D RHS vector
// N = dimension (length of b)
void fill_b(double * b, long N)
{
    long n = sqrt(N);
    for( long i = 0; i < n; i++ )
        for( long j = 0; j < n; j++ )
        {
            long idx = i*n + j;
            b[idx] = find_b(i,j,n);
        }
}
