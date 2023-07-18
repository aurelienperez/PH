#include <RcppArmadillo.h>
#include <vector>
#include <cmath>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat Van_loan(double x, arma::mat S, arma::colvec s, arma::rowvec alpha) {
  int dimension = S.n_rows;
  arma::mat s_alpha = s * alpha;
  
  arma::mat A = arma::join_rows(S, s_alpha);
  arma::mat B = arma::join_rows(arma::zeros<arma::mat>(dimension, dimension), S);
  arma::mat C = arma::join_cols(A, B);
  
  arma::mat result = arma::expmat(C * x);
  
  return result.submat(0, dimension, dimension - 1, 2 * dimension - 1);
}



// [[Rcpp::export]]
arma::mat expm(arma::mat inputMatrix) {
  arma::mat result = arma::expmat(inputMatrix);
  return result;
}

// [[Rcpp::export]]
arma::rowvec denominateur(const Rcpp::List& expS, const arma::colvec& s, const arma::rowvec& alpha) {
  int n=expS.size();
  arma::rowvec result(n);
  for (int i = 0; i < n; ++i) {
      arma::mat matrix = Rcpp::as<arma::mat>(expS[i]);
	  result[i] = arma::accu(alpha * matrix * s);
  }

  return result;
}

// Fonction A_EM
// [[Rcpp::export]]
arma::rowvec A_EM(const arma::vec& alpha, const arma::vec& weight, const arma::colvec& s, const Rcpp::List& expS, const arma::vec& denom)
{
  int dimension = alpha.size();
  arma::rowvec result(dimension);

  for (int k = 0; k < dimension; ++k) {
    double sum = 0.0;
    for (int i = 0; i < expS.size(); ++i) {
      arma::mat matrix = Rcpp::as<arma::mat>(expS[i]);
	  arma::rowvec extractedRow = matrix.row(k);
	  double dotprod = arma::dot(extractedRow , s);
      sum += weight[i] * dotprod / denom[i];
    }
    result[k] = alpha[k] * sum ;
  }

  return result;
}

// [[Rcpp::export]]
arma::rowvec B_EM(const double& dimension, const Rcpp::List& Gy,const arma::vec& weight, const arma::vec& denom)
{
  arma::rowvec result(dimension);

  for (int k = 0; k < dimension; ++k) {
    double sum = 0.0;
    for (int i = 0; i < Gy.size(); ++i) {
      arma::mat matrix = Rcpp::as<arma::mat>(Gy[i]);
	  double extractedDouble = matrix(k,k);
      sum += weight[i] * extractedDouble / denom[i];
    }
    result[k] =  sum ;
  }

  return result;
}

// [[Rcpp::export]]
arma::mat C_EM(const double& dimension,const arma::mat& S, const Rcpp::List& Gy,const arma::vec& weight, const arma::vec& denom,const arma::vec& s,const Rcpp::List& expS,const arma::rowvec& alpha)
{
  arma::mat result(dimension,dimension);

  for (int k = 0; k < dimension; ++k) {
	  for (int l = 0; l < dimension; ++l) {
		  if (k!=l){
			double sum = 0.0;
			for (int i = 0; i < Gy.size(); ++i) {
			  arma::mat matrix = Rcpp::as<arma::mat>(Gy[i]);
			  double extractedDouble = matrix(l,k);
			  sum += weight[i] * extractedDouble / denom[i];
			}
			result(l,k) = S(k,l) * sum ;}
		  else{
			double sum = 0.0;
			for (int i = 0; i < Gy.size(); ++i) {
			  arma::mat matrix = Rcpp::as<arma::mat>(expS[i]);
			  arma::colvec extractedCol = matrix.col(k);
			  double dotprod = arma::dot(alpha,extractedCol);
			  sum += weight[i] * dotprod / denom[i];
			}
			result(l,k) = s(k) * sum ;}
		  }
  }
  

  return result;
}


// [[Rcpp::export]]
arma::vec Exit(const arma::mat& S){
	int n=S.n_rows;
	arma::vec result(n);
	for (int i=0;i<n;++i){
		result(i)=-arma::accu(S.row(i));
	}
	return result;
}



// [[Rcpp::export]]
double logvraisemblance(const arma::vec& weight,const arma::vec& obs, const arma::mat& S, const arma::rowvec& alpha, const arma::colvec& s){

	int n=obs.n_elem;
	double sum=0;
	for (int k=0; k<n; ++k){
		sum += weight(k) * log(arma::accu(alpha * arma::expmat(S*obs(k)) * s));
	}
	return sum;
}

// [[Rcpp::export]]
arma::vec Densite(const arma::rowvec& alpha,const arma::mat& S,const arma::colvec& s,const arma::vec& y){
	int n=y.n_elem;
	arma::vec result(n);
	for (int i; i<n;++i){
		result(i) = arma::accu(alpha * arma::expmat(S*y(i)) * s);
	}
	return result;
}

double lambdaGompertz(double t, double beta){
	return exp(t*beta);
}

double invgGompertz(double t, double beta){
	return (lambdaGompertz(t,beta)-1)/beta;
}

// [[Rcpp::export]]
arma::vec DensiteGompertz(const arma::rowvec& alpha,const arma::mat& S,const arma::colvec& s,const arma::vec& y,const double& beta){
	int n=y.n_elem;
	arma::vec result(n);
	for (int i; i<n;++i){
		result(i) = lambdaGompertz(y(i),beta) * arma::accu(alpha * arma::expmat(S*invgGompertz(y(i),beta)) * s);
	}
	return result;
}