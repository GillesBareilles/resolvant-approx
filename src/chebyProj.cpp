#include <mex.h>
#include <complex>
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include "utilsMex.hpp"

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

// mex GCC='/usr/bin/g++-4.7' CXXFLAGS="\$CXXFLAGS -std=c++11" -I"./include/" -largeArrayDims src/chebyEval.cpp

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /* Computes
  $$ (A - \eta I)^{-1} \sum_{j=0}^{cheb_order_i} alphas(i, j) T_j(-2\eta(A-\eta I)^{-1}-I)b
  $$
  where $T_n$ are first order chebyshev polynomials
  With suitable coefficients alphas, can be used to approximate
  $$(inv(M)*K-tau_i*I)^{-1} * b$$ for a familly of complex shifts $\tau_i$
  */
    // Expected input :
        // - K : sparse double matrix
        // - M : sparse double matrix (optionally)
        // - b : double vector
        // - eta : real negative value
        // - alphas : *complex* matrix - coefficients of the chebyshev projections
        // - cheb_orders : integer vector - order of chebyshev approximation
    // Output :
        // - cheb_approx : complex array containing partial fraction expansion

    Eigen::SparseMatrix<double> K, M;
    Eigen::MatrixXd b;
    double eta;
    Eigen::MatrixXcd alphas;
    Eigen::VectorXd cheb_orders;

    if (nrhs == 5) {
      readSparse(prhs[0], K);
      setSparseIdentity(M, K.rows());
      readRealMatrix(prhs[1], b);
      eta = mxGetScalar(prhs[2]);
      readComplexMatrix(prhs[3], alphas);
      readRealVector(prhs[4], cheb_orders);
    } else if (nrhs == 6) {
      readSparse(prhs[0], K);
      readSparse(prhs[1], M);
      readRealMatrix(prhs[2], b);
      eta = mxGetScalar(prhs[3]);
      readComplexMatrix(prhs[4], alphas);
      readRealVector(prhs[5], cheb_orders);
    } else {
      mexErrMsgIdAndTxt( "MATLAB:chebyProj", "5 or 6 arguments expected, e.g. [ out ] = chebyProj(K, M, b, eta, alphas, cheb_degrees);");
    }

    // TODO : expand to dense matrices

    /* Check for proper input and output arguments */
    if (K.rows() != K.cols()) {
      mexErrMsgIdAndTxt( "MATLAB:chebyProj",
      "K must be a square matrix");
    }

    if (M.rows() != M.cols()) {
      mexErrMsgIdAndTxt( "MATLAB:chebyProj",
      "M must be a square matrix");
    }

    if ((K.rows() != M.rows()) || (K.cols() != b.rows())) {
      mexErrMsgIdAndTxt( "MATLAB:chebyProj",
      "K, M and b must be compatible");
    }

    if (eta > 0) {
      mexErrMsgIdAndTxt( "MATLAB:chebyProj",
      "eta must be negative");
    }

    if (alphas.rows() != cheb_orders.size()) {
      mexErrMsgIdAndTxt( "MATLAB:chebyProj",
      "alphas must have as many rows as cheb_orders elements (number of poles)");
    }

    for (int i=0;i<cheb_orders.size();i++) {
      if (std::floor(cheb_orders(i)) != cheb_orders(i)) {
        std::string message = "non integer chebyshev order for pole " + std::to_string(i);
        mexErrMsgIdAndTxt( "MATLAB:chebyProj", &(message.front()));
      }
      if (cheb_orders(i) > alphas.cols()) {
        std::string message = "pole " + std::to_string(i) + " : Asking for chebyshev approx of order " + std::to_string((int) std::round(cheb_orders(i))) + " but only " + std::to_string(alphas.cols()) + " coeffs given (at best).\nWill proceed to degree " + std::to_string(alphas.cols()) + " approximation for this pole.";
        mexWarnMsgTxt(&(message.front()));
        cheb_orders(i) = alphas.cols();
      }
      if (cheb_orders(i) < 0) {
        std::string message = "Pole " + std::to_string(i) + " : Asking for chebyshev approx of order " + std::to_string((int) std::round(cheb_orders(i))) + " < 0.\nWill proceed to degree " + std::to_string(alphas.cols()) + " approximation for this pole.";
        mexWarnMsgTxt(&(message.front()));
        cheb_orders(i) = alphas.cols();
      }
    }

    // Setting the coefficient matrix to approximate each pole to the required chebyshev degree

    for (int i=0; i<alphas.rows(); i++) {
      for (int j=cheb_orders[i]+1; j<alphas.cols(); j++) {
        alphas(i,j) = 0;
      }
    }

    // Declare solvers, variables
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solverRealShiftedSyst;
    solverRealShiftedSyst.compute(K-eta*M);

    Eigen::MatrixXcd temp;
    temp.resize(b.rows(), b.cols());
    temp.setZero();

    Eigen::MatrixXcd cheb_approx;
    cheb_approx.resize(b.rows(), alphas.rows() * b.cols());
    cheb_approx.setZero();

    Eigen::MatrixXd temp_real = M*b;
    b = solverRealShiftedSyst.solve(temp_real);

    Eigen::MatrixXcd T0 = b;

    temp_real = M*b;
    Eigen::MatrixXcd current_sol = solverRealShiftedSyst.solve(temp_real);
    Eigen::MatrixXcd T1 = -2*eta*current_sol -T0;

    for (int i=0;i<alphas.rows();i++) {
      cheb_approx.middleCols(i*b.cols(), b.cols()) = alphas(i,0)*T0 + alphas(i,1)*T1; // P(:, i*b.cols()+1:(i+1)*b.cols())
    }

    for (int j=2;j<alphas.cols();j++) {
      temp = T1;
      current_sol = solverRealShiftedSyst.solve(M*T1);
      T1 = -4*eta*current_sol -2*T1 - T0;
      T0 = temp;

      for (int i=0;i<alphas.rows();i++) {
        cheb_approx.middleCols(i*b.cols(), b.cols()) += alphas(i,j)*T1;
      }
    }


    plhs[0] = writeComplexMatrix(cheb_approx);

    return;
}
