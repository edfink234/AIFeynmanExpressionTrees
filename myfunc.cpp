// myfunc.cpp
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <Eigen/Dense>
#include <emscripten/emscripten.h>

using namespace Eigen;
using std::complex;

double lastResidual = 0;

// Example helper: hyperbolic secant
extern "C"
{
    EMSCRIPTEN_KEEPALIVE
    double sech(double x)
    {
        return 1.0 / cosh(x);
    }
}
extern "C" {
    EMSCRIPTEN_KEEPALIVE
    double getLastResidual() {
        return lastResidual;
    }
}
// Export a function to compute eigenvalues from a complex matrix.
// We assume that the input matrix is passed as a pointer to double
// representing 2*N*N doubles, where each complex number is interleaved:
// [Re_0, Im_0, Re_1, Im_1, ..., Re_{N*N-1}, Im_{N*N-1}]
// N is provided as an additional argument.
//
// This function computes the eigenvalues using Eigen’s ComplexEigenSolver
// and returns a pointer to a newly allocated array of 2*N doubles (real and imaginary parts)
// representing the eigenvalues. (Memory management on the JS side is your responsibility.)

extern "C"
{
    EMSCRIPTEN_KEEPALIVE
    double* computeEigenspectrum(double* inMat, int N)
    {
        // Create an Eigen complex matrix of size N x N.
        MatrixXcd M(N, N);
        // Fill the matrix by reading from the interleaved input array.
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                int idx = 2 * (i * N + j);
                double re = inMat[idx];
                double im = inMat[idx + 1];
                M(i, j) = complex<double>(re, im);
            }
        }
//        std::cout << "done" << std::endl;
        // Compute the eigenvalues.
        ComplexEigenSolver<MatrixXcd> ces;
//        std::cout << "Matrix M size: " << M.rows() << " x " << M.cols() << std::endl;
//        std::cout << "M(0,0) = " << M(0,0) << std::endl;  // check that data looks sane
        ces.compute(M);
        
        if (ces.info() != Eigen::Success)
        {
            // Eigenvalue computation failed
            return nullptr;
        }
        
        
        auto lambdas = ces.eigenvalues();
        auto vectors = ces.eigenvectors();

        // === STEP 1: Create index list [0, 1, ..., N-1]
        std::vector<int> indices(N);
        std::iota(indices.begin(), indices.end(), 0);

        // === STEP 2: Sort indices by descending real part of eigenvalues
        std::sort(indices.begin(), indices.end(), [&](int i, int j) {
            return lambdas[i].real() > lambdas[j].real();
        });

        // === STEP 3: Allocate output buffer
        int totalDoubles = 2 * N + 2 * N * N;
        double* output = (double*) malloc(totalDoubles * sizeof(double));
        if (!output) return nullptr;

        // === STEP 4: Write sorted eigenvalues
        for (int i = 0; i < N; i++) {
            int idx = indices[i];
            std::complex<double> lambda = lambdas[idx];
            output[2 * i]     = lambda.real();
            output[2 * i + 1] = lambda.imag();
        }

        // === STEP 5: Write sorted eigenvectors (column-major: each col is a vec)
        for (int j = 0; j < N; j++) {
            int idx = indices[j];  // get sorted index
            VectorXcd v = vectors.col(idx);
            for (int i = 0; i < N; i++) {
                int outIdx = 2 * N + 2 * (j * N + i);  // after eigenvalues
                output[outIdx]     = v(i).real();
                output[outIdx + 1] = v(i).imag();
            }
        }
        
        // Compute accuracy: ||Mv - lambda*v|| for each eigenpair
        double max_residual = 0.0;
        for (int j = 0; j < N; j++) {
            VectorXcd v = ces.eigenvectors().col(j);
            complex<double> lambda = ces.eigenvalues()(j);
            VectorXcd residual = M * v - lambda * v;
            double norm = residual.norm();  // 2-norm
            if (norm > max_residual) max_residual = norm;
        }
        lastResidual = max_residual;
        std::cout << "max residual =" << lastResidual << '\n';
        // Store the max residual in a static variable so we can retrieve it in JS

        return output;
    }
}


//emcc myfunc.cpp -O2 -s WASM=1 -s MODULARIZE=1 -s EXPORT_NAME="createModule" \
//  -s EXPORTED_FUNCTIONS='["_sech", "_computeEigenspectrum", "_getLastResidual", "_malloc", "_free"]' -s EXPORT_ES6=1 \
//  -s EXPORTED_RUNTIME_METHODS="['ccall', 'cwrap']" \
//  -s TOTAL_MEMORY=1073741824 \
//  -s INITIAL_MEMORY=268435456 \
//  -s STACK_SIZE=10485760 \
//  -I/opt/homebrew/opt/eigen/include/eigen3 \
//  -o myfunc.js
