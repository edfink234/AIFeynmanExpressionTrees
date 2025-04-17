#include <math.h>

// This is your performance-critical function that you wish to call from JS.
// For example, a function that computes f(x) = sin(x) * cos(x).

double sech(double x) {
    return 1.0/cosh(x);
}
   
