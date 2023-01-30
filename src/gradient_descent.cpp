#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::depends('RcppArmadillo')]]
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

using namespace Rcpp;
using namespace std;
using namespace arma;

// [[Rcpp::export]]
arma::uvec col_mins_id_cpp(arma::mat &X) {
   
   arma::urowvec mins = arma::index_min(X, 0);
   arma::uvec mins2 = arma::conv_to<uvec>::from(mins);
   
   return(mins2);

}

// [[Rcpp::export]]
arma::uvec sp_col_mins_id_cpp(arma::sp_mat &X) {
   
  arma::urowvec mins = arma::index_min(X, 0);
  arma::uvec mins2 = arma::conv_to<uvec>::from(mins);

  for (int i = 0; i < X.n_cols; ++i){
    if(mins2[i] > 0){
      if(X(0, i) == 0){
        mins2[i] = 0;
      }
    }
  }   

  return(mins2);

}

// [[Rcpp::export]]
arma::sp_mat update_B_cpp(
	arma::sp_mat &A, 
	arma::sp_mat B, 
	const arma::sp_mat &K, 
	int &it,
	int t = 0,
	bool verbose = true
	){

  //precompute some gradient terms
  arma::sp_mat t1s = A * A.t();
  arma::sp_mat t2s = K * A.t();
  arma::mat t1 = arma::conv_to<mat>::from(t1s);
  arma::mat t2 = arma::conv_to<mat>::from(t2s);

  //Progress
  Progress p(it, verbose);

  //update rows of A for itven number of iterations
  while(t < it){

  	R_CheckUserInterrupt();
  	p.increment(); 

		arma::mat G = 2 * (K * B * t1 - t2);

		//get colmin indices
		arma::uvec col_mins = col_mins_id_cpp(G);
		arma::uvec col_ids = arma::regspace<uvec>(0, G.n_cols - 1);

		//locations
		arma::umat locations = join_rows(col_mins, col_ids);

		// vector of 1s
		arma::vec values;
		values = values.ones(G.n_cols);

		//initialize matrix
		arma::sp_mat e(locations.t(), values, G.n_rows, G.n_cols, true, true);

		//update
		B += (2. / (t + 2)) * (e - B);
		t += 1;
	    
	}

  return(B);

}

// [[Rcpp::export]]
arma::sp_mat update_B_cpp_sp(
	arma::sp_mat &A, 
	arma::sp_mat B, 
	const arma::sp_mat &K, 
	int &it,
	int t = 0,
	bool verbose = true
	){

  //precompute some gradient terms
  arma::sp_mat t1 = A * A.t();
  arma::sp_mat t2 = K * A.t();

  //Progress
  Progress p(it, verbose);

  //update rows of A for itven number of iterations
  while(t < it){

  	R_CheckUserInterrupt();
  	p.increment(); 

		arma::sp_mat G = 2 * (K * B * t1 - t2);

		//get colmin indices
		arma::uvec col_mins = sp_col_mins_id_cpp(G);
		arma::uvec col_ids = arma::regspace<uvec>(0, G.n_cols - 1);

		//locations
		arma::umat locations = join_rows(col_mins, col_ids);

		// vector of 1s
		arma::vec values;
		values = values.ones(G.n_cols);

		//initialize matrix
		arma::sp_mat e(locations.t(), values, G.n_rows, G.n_cols, true, true);

		//update
		B += (2. / (t + 2)) * (e - B);
		t += 1;
	    
	}

  return(B);

}

// [[Rcpp::export]]
arma::mat initialize_A_cpp(const int &nr, const int &nc){

  // Initialize
  arma::mat A(nr, nc);

  // Fill Matrix
  for (int i = 0; i < nc; ++i){
  
  	arma::vec rv = randu(nr);

  	A.col(i) = rv / sum(rv);

  }

  return(A);

}

// [[Rcpp::export]]
arma::uvec initialize_A_idx(const arma::sp_mat &t1, const arma::sp_mat &t2, bool verbose = true){

  // Initialize
  arma::uvec min_idx(t2.n_cols);

  //Progress
  Progress p(t2.n_cols, verbose);

  // Determine A
  for (int i = 0; i < t2.n_cols; ++i){

  	p.increment(); 

    //Random Vector
    arma::vec rv = randu(t2.n_rows);
    rv = rv / sum(rv);

    //Multiply
    arma::vec G = 2 * (t1 * rv - t2.col(i));

    //Min
    min_idx(i) = arma::index_min(G);
  
  }

  return(min_idx);

}

// [[Rcpp::export]]
arma::sp_mat update_A_cpp(
	arma::sp_mat A, 
	arma::sp_mat &B, 
	const arma::sp_mat &K, 
	int &it,
	int t = 0,
	bool verbose = true
	){

  //precompute some gradient terms
  arma::sp_mat t2s = (K * B).t();
  arma::sp_mat t1s = t2s * B;
  arma::mat t1 = arma::conv_to<mat>::from(t1s);
  arma::mat t2 = arma::conv_to<mat>::from(t2s);

  //Progress
  Progress p(it, verbose);

  //update rows of A for itven number of iterations
  while(t < it){

  	R_CheckUserInterrupt();
  	p.increment(); 

		arma::mat G = 2 * (t1 * A - t2);

		//get colmin indices
		arma::uvec col_mins = col_mins_id_cpp(G);
		arma::uvec col_ids = arma::regspace<uvec>(0, G.n_cols - 1);

		//locations
		arma::umat locations = join_rows(col_mins, col_ids);

		// vector of 1s
		arma::vec values;
		values = values.ones(G.n_cols);

		//initialize matrix
		arma::sp_mat e(locations.t(), values, G.n_rows, G.n_cols, true, true);

		//update
		A += (2. / (t + 2)) * (e - A);
		t += 1;

	}

	return(A);
	
}

// [[Rcpp::export]]
arma::sp_mat update_A_cpp_sp(
	arma::sp_mat A, 
	arma::sp_mat &B, 
	const arma::sp_mat &K, 
	int &it,
	int t = 0,
	bool verbose = true
	){

  //precompute some gradient terms
  arma::sp_mat t2 = (K * B).t();
  arma::sp_mat t1 = t2 * B;

  //Progress
  Progress p(it, verbose);

  //update rows of A for itven number of iterations
  while(t < it){

  	R_CheckUserInterrupt();
  	p.increment(); 

		arma::sp_mat G = 2 * (t1 * A - t2);

		//get colmin indices
		arma::uvec col_mins = sp_col_mins_id_cpp(G);
		arma::uvec col_ids = arma::regspace<uvec>(0, G.n_cols - 1);

		//locations
		arma::umat locations = join_rows(col_mins, col_ids);

		// vector of 1s
		arma::vec values;
		values = values.ones(G.n_cols);

		//initialize matrix
		arma::sp_mat e(locations.t(), values, G.n_rows, G.n_cols, true, true);

		//update
		A += (2. / (t + 2)) * (e - A);
		t += 1;

	}

	return(A);
	
}


