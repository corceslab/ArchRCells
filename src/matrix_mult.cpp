#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
SEXP mat_mult_dens_dens_cpp(
    const Eigen::Map<Eigen::MatrixXd> &A,
    const Eigen::Map<Eigen::MatrixXd> &B 
    ){
  
  Eigen::MatrixXd C = A * B;
  return Rcpp::wrap(C);

}

// [[Rcpp::export]]
SEXP mat_mult_dens_sp_cpp(
    const Eigen::Map<Eigen::MatrixXd> &A,
    const Eigen::SparseMatrix<double> &B
    ){
  
  Eigen::MatrixXd C = A * B;
  return Rcpp::wrap(C);

}

// [[Rcpp::export]]
SEXP mat_mult_sp_dens_cpp(
    const Eigen::SparseMatrix<double> & A,
    const Eigen::Map<Eigen::MatrixXd> & B 
    ){
  
    Eigen::MatrixXd C = A * B;
    return Rcpp::wrap(C);

}

// [[Rcpp::export]]
SEXP mat_mult_sp_sp_cpp(
    const Eigen::SparseMatrix<double> & A,
    const Eigen::SparseMatrix<double> & B
    ){
  
    const Eigen::SparseMatrix<double> C = A * B;
    return Rcpp::wrap(C);

}

