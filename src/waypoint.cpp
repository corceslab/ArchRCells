#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
using namespace arma;

// [[Rcpp::depends('RcppArmadillo')]]

// [[Rcpp::export]]
arma::uvec determine_waypoints_cpp(arma::colvec &v, int &idx, int &N, int &it, int &nc, bool verbose){
   
	//zero
	//arma::uvec zero;
	//zero.zeros(0);

 	//Initialize
 	arma::uvec iter_set(it);
 	iter_set(0) = idx - 1;

 	//Dist Matrix
 	arma::mat dists(N, it);
 	dists.col(0) = abs(v - v[idx - 1]);

 	for(int k = 0; k < it - 1; k++){
 		
 		if((k > 0) & verbose){
	 		if(k % 50 == 0){
	 			Rcout << "Component " << nc << " : Finding Waypoints " << k << " of " << it << std::endl;
	 		} 			
 		}

 		//subset by
 		arma::uvec k_idx = arma::conv_to<uvec>::from(regspace(0, k));

 		//Minimum distances across the current set
 		arma::vec min_dists = arma::min(dists.cols(k_idx), 1);

 		//Max
 		int idx_max = arma::index_max(min_dists);

 		//Add
 		dists.col(k + 1) = abs(v - v[idx_max]);
 		iter_set(k + 1) = idx_max;

 	}

 	return(iter_set);

}
