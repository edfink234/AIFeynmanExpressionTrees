#include <math.h>

// This is your performance-critical function that you wish to call from JS.
// For example, a function that computes f(x) = sin(x) * cos(x).

double sech(double x) {
    return 1.0/cosh(x);
}
   
//emcc myfunc.cpp -O2 -s WASM=1 -s MODULARIZE=1 -s EXPORT_NAME="createModule" \
//  -s EXPORTED_FUNCTIONS='["_sech", "_computeEigenspectrum", "_malloc", "_free"]' -s EXPORT_ES6=1 \
//  -s EXPORTED_RUNTIME_METHODS="['ccall', 'cwrap']" \
//  -I/opt/homebrew/opt/eigen/include/eigen3 \
//  -o myfunc.js
