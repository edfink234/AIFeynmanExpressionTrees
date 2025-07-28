// myfunc.cpp
#include <cmath>
#include <complex>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <Eigen/Dense>
#include <emscripten/emscripten.h>
#include <boost/numeric/odeint.hpp>
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>

using namespace Eigen;
using std::complex;

/* --- contiguous storage helpers --------------------------------- */
#define RE(i)  (psi_buffer[(i)])          /* 0 … N-1           */
#define IM(i)  (psi_buffer[(Npsi + (i))]) /* N … 2N-1          */
double lastResidual = 0;
double lastNLSResidual = 0.0;
/* ──────────────────── global simulation state ─────────────────── */
static double current_time = 0.0;
static double total_time   = 10.0;
static double dt           = 1e-3;
static double x_particle   = 0.0;
static double v_particle   = 0.0;
/* ───────────────────────── ODE parameters ───────────────────────── */
static double m_param     = 1.0;
static double Omega_param = 0.2;
static double A_param     = 1.0;
static double b_param     = 1.0;
/* ───────────────────────── PDE parameters ───────────────────────── */
static double dx_pde      = 0.05;
static double Omega_pde   = 0.2;
static double A_pde       = 1.0;
static double b_pde       = 1.0;
static double x_min_pde   = 0.0;

/* simple ψ-buffer so JS can peek at the field                 */
static int    Npsi        = 1024; // grid size
static double *psi_buffer = nullptr; // real‐imag interleaved

/* ---------------------------------------------------------------- */
// tiny helpers the JS GUI expects
extern "C"
{
    EMSCRIPTEN_KEEPALIVE double sech(double x)    { return 1.0 / std::cosh(x);}
    EMSCRIPTEN_KEEPALIVE double getCurrentTime()  { return current_time; }
    EMSCRIPTEN_KEEPALIVE double getXParticle()    { return x_particle; }
    EMSCRIPTEN_KEEPALIVE double getVParticle()    { return v_particle; }
    EMSCRIPTEN_KEEPALIVE double *getPsiPointer()  { return psi_buffer; }
}

/* acceleration a(x) = (-Ω² x + 2 A b sech²(bx) tanh(bx)) / m */
static inline double accel(double x)
{
    const double s  = sech(b_param * x);
    const double t  = std::tanh(b_param * x);
    return (-Omega_param * Omega_param * x +
            2.0 * A_param * b_param * s * s * t) / m_param;
}

/* ───────────── dummy integrators – replace with real code ───────── */
static inline void doODEstep()
{
    double k1x = v_particle;
    double k1v = accel(x_particle);

    double k2x = v_particle + 0.5 * dt * k1v;
    double k2v = accel(x_particle + 0.5 * dt * k1x);

    double k3x = v_particle + 0.5 * dt * k2v;
    double k3v = accel(x_particle + 0.5 * dt * k2x);

    double k4x = v_particle + dt * k3v;
    double k4v = accel(x_particle + dt * k3x);

    x_particle += (dt / 6.0) * (k1x + 2*k2x + 2*k3x + k4x);
    v_particle += (dt / 6.0) * (k1v + 2*k2v + 2*k3v + k4v);
}

/* central-difference Laplacian with periodic BC */
void laplacian(int N, double dx, const double *Re, const double *Im, double *outRe, double *outIm)
{
    const double dx2 = dx*dx;
    for (int i = 0; i < N; ++i)
    {
        int ip = (i + 1) % N;
        int im = (i - 1 + N) % N;
        outRe[i] = ( Re[ip] - 2.0*Re[i] + Re[im] ) / dx2;
        outIm[i] = ( Im[ip] - 2.0*Im[i] + Im[im] ) / dx2;
    }
}

/* one RK-4 step of  i u_t = -½ u_xx + V u - |u|² u  */
static inline void doPDEstep()
{
    const int N = Npsi;
    const double h = dt;

    auto V_at = [](double x)
    {
        return 0.5*Omega_pde*Omega_pde*x*x
             + A_pde*std::pow(sech(b_pde*x), 2.0);
    };

    auto pde_rhs = [&](const Eigen::VectorXd &ReIn, const Eigen::VectorXd &ImIn, Eigen::VectorXd &ReOut, Eigen::VectorXd &ImOut)
    {
        static Eigen::VectorXd LapR(Npsi), LapI(Npsi);
        laplacian(Npsi, dx_pde, ReIn.data(), ImIn.data(), LapR.data(), LapI.data());

        for (int i = 0; i < Npsi; ++i)
        {
            double x   = x_min_pde + i*dx_pde;
            double Vx  = V_at(x);
            double mod2 = ReIn[i]*ReIn[i] + ImIn[i]*ImIn[i];

            double inRe = -0.5*LapR[i] - mod2*ReIn[i] + Vx*ReIn[i];
            double inIm = -0.5*LapI[i] - mod2*ImIn[i] + Vx*ImIn[i];

            ReOut[i] = -inIm;              /* multiply by +i */
            ImOut[i] =  inRe;
        }
    };

    /* copy ψ ⇒ Eigen vectors */
    Eigen::VectorXd uR(Npsi), uI(Npsi);
    for (int i = 0; i < Npsi; ++i) { uR[i] = RE(i);  uI[i] = IM(i); }
    
    /* buffers */
    Eigen::VectorXd k1R(Npsi), k1I(Npsi),
                        k2R(Npsi), k2I(Npsi),
                        k3R(Npsi), k3I(Npsi),
                        k4R(Npsi), k4I(Npsi),
                        tmpR(Npsi), tmpI(Npsi);

    /* k1 */
    pde_rhs(uR, uI, k1R, k1I);

    /* u2 = u + dt/2 * k1 */
    tmpR = uR + 0.5*dt * k1R;
    tmpI = uI + 0.5*dt * k1I;
    pde_rhs(tmpR, tmpI, k2R, k2I);

    /* u3 = u + dt/2 * k2 */
    tmpR = uR + 0.5*dt * k2R;
    tmpI = uI + 0.5*dt * k2I;
    pde_rhs(tmpR, tmpI, k3R, k3I);

    /* u4 = u + dt * k3 */
    tmpR = uR + dt * k3R;
    tmpI = uI + dt * k3I;
    pde_rhs(tmpR, tmpI, k4R, k4I);

    /* new u = u + dt/6*(k1 + 2k2 + 2k3 + k4) */
    uR += dt/6.0 * (k1R + 2.0*k2R + 2.0*k3R + k4R);
    uI += dt/6.0 * (k1I + 2.0*k2I + 2.0*k3I + k4I);

    /* write back to ψ-buffer */
    for (int i = 0; i < Npsi; ++i) { RE(i) = uR[i];  IM(i) = uI[i]; }
}



extern "C"
{
    EMSCRIPTEN_KEEPALIVE

    void stepForPlotInterval(double maxMillis)
    {
        using clock = std::chrono::high_resolution_clock;
        auto start = clock::now();

        while (current_time < total_time)
        {
            doODEstep(); // RK4 particle update
            printf("u[%d] = %lf + %lfi\n", Npsi/3, RE(Npsi/3), IM(Npsi/3));
            doPDEstep(); // Split-step or RK4 PDE update
            current_time += dt;

            auto now = clock::now();
            double elapsed = std::chrono::duration<double, std::milli>(now - start).count();
            if (elapsed >= maxMillis)
            {
                break;
            }
        }
    }

    /* allow JS to reset / configure the solver on each run */
    EMSCRIPTEN_KEEPALIVE
    void setSimParameters(double new_dt, double new_T, int new_Npsi, double m, double Omega, double Atrap, double btrap, double x0, double v0, double dx, double x_min, double OmegaPDE, double Apde, double bpde, double A_sol)
    {
        dt = new_dt;
        total_time = new_T;
        current_time = 0.0; //restarting simulation here
        m_param     = m;
        Omega_param = Omega;
        A_param     = Atrap;
        b_param     = btrap;
        x_particle  = x0;
        v_particle  = v0;
        dx_pde    = dx;
        x_min_pde = x_min;
        Omega_pde = OmegaPDE;
        A_pde     = Apde;
        b_pde     = bpde;
        
//        printf("Simulation Parameters:\n");
//        printf("  dt = %lf\n", dt);
//        printf("  total_time = %lf\n", total_time);
//        printf("  current_time = %lf\n", current_time);
//        printf("  m_param = %lf\n", m_param);
//        printf("  Omega_param = %lf\n", Omega_param);
//        printf("  A_param = %lf\n", A_param);
//        printf("  b_param = %lf\n", b_param);
//        printf("  x_particle = %lf\n", x_particle);
//        printf("  v_particle = %lf\n", v_particle);
//        printf("  dx_pde = %lf\n", dx_pde);
//        printf("  x_min_pde = %lf\n", x_min_pde);
//        printf("  Omega_pde = %lf\n", Omega_pde);
//        printf("  A_pde = %lf\n", A_pde);
//        printf("  b_pde = %lf\n", b_pde);

//        if (new_Npsi != Npsi) // re-allocate ψ if grid changes
//        {
            if (psi_buffer)
            {
                free(psi_buffer);
            }
            Npsi = new_Npsi;
            psi_buffer = (double*) malloc(sizeof(double) * 2 * Npsi);
//        }
        /* initialise ψ to a sech – real part only, imag = 0 */
        for (int i = 0; i < Npsi; ++i)
        {
            double xpos  = x_min + i * dx_pde;              // ← use x_min
            double arg   = A_sol * (xpos - x0);
            double amp   = A_sol * sech(arg);
            double phase = v0  * xpos;
            RE(i) = amp * std::cos(phase);
            IM(i) = amp * std::sin(phase);
        }
        printf("u[%d] = %lf + %lfi\n", Npsi/3, psi_buffer[(2*(Npsi/3))], psi_buffer[(2*(Npsi/3))+1]);
    }
}

extern "C"
{
    EMSCRIPTEN_KEEPALIVE
    double getLastResidual()
    {
        return lastResidual;
    }
    EMSCRIPTEN_KEEPALIVE
    double getNLSResidual()
    {
        return lastNLSResidual;
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
        std::sort(indices.begin(), indices.end(), [&](int i, int j)
        {
            return lambdas[i].real() > lambdas[j].real();
        });

        // === STEP 3: Allocate output buffer
        int totalDoubles = 2 * N + 2 * N * N;
        double* output = (double*) malloc(totalDoubles * sizeof(double));
        if (!output)
        {
            return nullptr;
        }

        // === STEP 4: Write sorted eigenvalues
        for (int i = 0; i < N; i++)
        {
            int idx = indices[i];
            std::complex<double> lambda = lambdas[idx];
            output[2 * i]     = lambda.real();
            output[2 * i + 1] = lambda.imag();
        }

        // === STEP 5: Write sorted eigenvectors (column-major: each col is a vec)
        for (int j = 0; j < N; j++)
        {
            int idx = indices[j];  // get sorted index
            VectorXcd v = vectors.col(idx);
            for (int i = 0; i < N; i++)
            {
                int outIdx = 2 * N + 2 * (j * N + i);  // after eigenvalues
                output[outIdx]     = v(i).real();
                output[outIdx + 1] = v(i).imag();
            }
        }
        
        // Compute accuracy: ||Mv - lambda*v|| for each eigenpair
        double max_residual = 0.0;
        for (int j = 0; j < N; j++)
        {
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

// Define a functor for Levenberg-Marquardt
template<typename T>
struct NLSResidualFunctor
{
    using Scalar = T;
    using InputType = Eigen::Matrix<T, Eigen::Dynamic, 1>;
    using ValueType = Eigen::Matrix<T, Eigen::Dynamic, 1>;
    using JacobianType = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

    enum
    {
        InputsAtCompileTime = Eigen::Dynamic,
        ValuesAtCompileTime = Eigen::Dynamic
    };
    
    int N;
    double dx;
    double* V;
    double g;
    double omega;

    NLSResidualFunctor(int N_, double dx_, double* V_, double g_, double omega_)
        : N(N_), dx(dx_), V(V_), g(g_), omega(omega_) {}

    int inputs() const { return 2 * N; }       // Real and Imag parts
    int values() const { return 2 * N; }       // Residual vector size

    // Computes F(x) = residual vector
    int operator()(const Eigen::VectorXd& U, Eigen::VectorXd& F) const
    {
        Eigen::VectorXd Ur = U.head(N);
        Eigen::VectorXd Ui = U.tail(N);
        Eigen::VectorXd D2Ur(N), D2Ui(N);

        for (int i = 0; i < N; ++i)
        {
            int ip = (i + 1) % N;
            int im = (i - 1 + N) % N;
            D2Ur[i] = (Ur[ip] - 2 * Ur[i] + Ur[im]) / (dx * dx);
            D2Ui[i] = (Ui[ip] - 2 * Ui[i] + Ui[im]) / (dx * dx);
        }

        for (int i = 0; i < N; ++i)
        {
            double r = Ur[i], im = Ui[i];
            double U2 = r*r + im*im;
            double common = g * U2 + V[i] + omega;
            F[i] = -0.5 * D2Ur[i] + common * r;
            F[N + i] = -0.5 * D2Ui[i] + common * im;
        }
        return 0;
    }
    
    int df(const InputType& U, JacobianType& J) const
    {
        Eigen::VectorXd Ur = U.head(N);
        Eigen::VectorXd Ui = U.tail(N);

        Eigen::MatrixXd D2 = Eigen::MatrixXd::Zero(N, N);
        for (int i = 0; i < N; ++i)
        {
            D2(i, i) = -2;
            D2(i, (i + 1) % N) = 1;
            D2(i, (i - 1 + N) % N) = 1;
        }
        D2 /= (dx * dx);

        Eigen::VectorXd diagJ11(N), diagJ22(N), diagJ12(N);
        for (int i = 0; i < N; ++i)
        {
            double r = Ur[i], im = Ui[i];
            diagJ11[i] = g * (3 * r * r + im * im) + V[i] + omega;
            diagJ22[i] = g * (r * r + 3 * im * im) + V[i] + omega;
            diagJ12[i] = 2 * g * r * im;
        }

        Eigen::MatrixXd J11 = (-0.5 * D2).eval() + diagJ11.asDiagonal().toDenseMatrix();
        Eigen::MatrixXd J22 = (-0.5 * D2).eval() + diagJ22.asDiagonal().toDenseMatrix();
        Eigen::MatrixXd J12 = diagJ12.asDiagonal();

        // Assemble full 2N x 2N Jacobian
        J.topLeftCorner(N, N)     = J11;
        J.topRightCorner(N, N)    = J12;
        J.bottomLeftCorner(N, N)  = J12;
        J.bottomRightCorner(N, N) = J22;

        return 0;
    }
};

extern "C"
{
    EMSCRIPTEN_KEEPALIVE
    double* refineCpp(double* U_init, double* V, int N, double dx, double g, double omega, int useLM)
    {
        Eigen::VectorXd U(2 * N);
        for (int i = 0; i < 2 * N; ++i)
        {
            U[i] = U_init[i];
        }
        
        Eigen::VectorXd F(2 * N);
        NLSResidualFunctor<double> functor(N, dx, V, g, omega);
        
        if (useLM)
        {
            //Eigen::NumericalDiff<NLSResidualFunctor<double>> numDiff(functor);
            //Eigen::LevenbergMarquardt<Eigen::NumericalDiff<NLSResidualFunctor<double>>> lm(numDiff);
            Eigen::LevenbergMarquardt<NLSResidualFunctor<double>> lm(functor);
            lm.parameters.maxfev = 10000;
            lm.parameters.xtol = 1e-10;
            lm.minimize(U);
            functor(U, F);
        }
        else
        {
            // Simple Newton with fixed iterations and learning rate
            Eigen::VectorXd DU(2 * N);
            for (int iter = 0; iter < 10; ++iter)
            {
                functor(U, F);
                lastNLSResidual = F.norm();
                std::cout << "iter " << iter << ", " << "||F|| = " << lastNLSResidual
                << std::endl;
                if (lastNLSResidual < 1e-10) break;
                Eigen::MatrixXd J(2 * N, 2 * N);
                functor.df(U, J);
                DU = J.ldlt().solve(-F);  // Or .ldlt() if symmetric
                U += DU;
            }
        }

        lastNLSResidual = F.norm();  // L2 norm of residual
        // Allocate result
        double* result = (double*) malloc(sizeof(double) * 2 * N);
        if (!result) return nullptr;
        for (int i = 0; i < 2 * N; ++i) result[i] = U[i];
        return result;
    }
}

//emcc myfunc.cpp -O2 -s WASM=1 -s MODULARIZE=1 -s EXPORT_NAME="createModule" \
//  -s EXPORTED_FUNCTIONS='["_sech", "_setSimParameters", "_computeEigenspectrum", "_getLastResidual", "_refineCpp", "_stepForPlotInterval",  "_getCurrentTime", "_getXParticle", "_getVParticle", "_getPsiPointer", "_malloc", "_free"]' \
//  -s EXPORT_ES6=1 \
//  -s EXPORTED_RUNTIME_METHODS="['ccall', 'cwrap']" \
//  -s TOTAL_MEMORY=1073741824 \
//  -s INITIAL_MEMORY=268435456 \
//  -s STACK_SIZE=10485760 \
//  -I/opt/homebrew/opt/eigen/include/eigen3 \
//  -L/opt/homebrew/Cellar/boost/1.84.0 \
//  -I/opt/homebrew/Cellar/boost/1.84.0/include \
//  -o myfunc.js
