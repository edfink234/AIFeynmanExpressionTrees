// myfunc.cpp
#include <cmath>
#include <complex>
#include <cstdlib>
#include <Eigen/Dense>
#include <emscripten/emscripten.h>

using namespace Eigen;
using std::complex;

// Example helper: hyperbolic secant
extern "C" {
EMSCRIPTEN_KEEPALIVE
double sech(double x) {
    return 1.0 / cosh(x);
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

extern "C" {
EMSCRIPTEN_KEEPALIVE
double* computeEigenspectrum(double* inMat, int N) {
    // Create an Eigen complex matrix of size N x N.
    MatrixXcd M(N, N);
    // Fill the matrix by reading from the interleaved input array.
    for (int i = 0; i < N; i++){
      for (int j = 0; j < N; j++){
         int idx = 2 * (i * N + j);
         double re = inMat[idx];
         double im = inMat[idx + 1];
         M(i, j) = complex<double>(re, im);
      }
    }
    
    // Compute the eigenvalues.
    ComplexEigenSolver<MatrixXcd> ces;
    ces.compute(M);
    
    // Allocate an output array of 2*N doubles.
    double* output = (double*) malloc(2 * N * sizeof(double));
    if (!output) return nullptr; // allocation failure
    
    // Write eigenvalues into the output array.
    // (Here, we simply output them in the same order as Eigen returns them.
    // You could sort them if required.)
    for (int i = 0; i < N; i++){
       complex<double> lambda = ces.eigenvalues()(i);
       output[2 * i]     = lambda.real();
       output[2 * i + 1] = lambda.imag();
    }
    
    return output;
}
}

//emcc myfunc.c -Os -s WASM=1 -s EXPORT_ES6=1 -s MODULARIZE=1 -s EXPORT_NAME="createModule" -s EXPORTED_FUNCTIONS='["_sech"]' -s EXPORTED_RUNTIME_METHODS="['ccall', 'cwrap']" -o myfunc.js
