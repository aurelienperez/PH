#include <RcppArmadillo.h>
#include <vector>
#include <cmath>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

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
arma::mat matrix_vanloan(arma::mat A, arma::mat B, arma::mat C) {
  int dimension = A.n_rows;

  arma::mat M1 = arma::join_rows(A, C);
  arma::mat M2 = arma::join_rows(arma::zeros<arma::mat>(dimension, dimension), B);
  arma::mat result = arma::join_cols(M1, M2);
  
  
  return result;
}

//' Calcule les elements S^n / n! jusqu'a une taille donnee

// [[Rcpp::export]]
void vector_of_matrices(std::vector<arma::mat> & vect, const arma::mat & S, double a, int vect_size) {
  arma::mat I;
  I.eye(size(S));
  
  arma::mat P = I + S * (1 / a);
  
  vect.push_back(I);
  
  for (int k{1}; k <= vect_size; ++k) {
    vect.push_back((P * (1.0 / k) ) * vect[k - 1]);
  }
}


//' Calcule exp(Sx) 

// [[Rcpp::export]]
arma::mat m_exp_sum(double x, int n, const std::vector<arma::mat> & pow_vector, double a) {
  arma::mat res_mat = pow_vector[0];
  
  for (int i{1}; i <= n; ++i) {
    res_mat = res_mat + pow_vector[i] * exp(i * std::log(a * x));
  }
  res_mat = res_mat * exp(-a * x);
  
  return res_mat;
}


//' Calcule A^(2^n)

// [[Rcpp::export]]
void pow2_matrix(int n , arma::mat & A) {
  arma::mat aux_mat(size(A));
  
  for (int i{1}; i <= n; ++i) {
    aux_mat = A * A;
    A = aux_mat;
  }
}

//' Maximum de la diagonale de la matrice
// [[Rcpp::export]]
double max_diagonal(const arma::mat & A) {
  double maximum{A(0,0)};
  for (int i{0}; i < A.n_rows; ++i) {
    if (A(i,i) > maximum) {
      maximum = A(i,i);
    }
  }
  return maximum;
}
//' n tel que P(N > n) = h ou N suit une loi de Poisson
// [[Rcpp::export]] 
int find_n(double h, double lambda) {
  int n{0};
  double cum_prob{0.0};
  
  do {
    cum_prob += R::dpois(n, lambda, false);
    ++n;
  } while (cum_prob < 1.0 - h);
  
  return (n - 1);
}


// [[Rcpp::export]]
arma::mat expSxUNI(const double x, const double m, const std::vector<arma::mat> aux_vect, const double a){
	arma::mat J;
	if (x * a <= 1.0) {
      J = m_exp_sum(x, m, aux_vect, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      J = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
      
      pow2_matrix(n, J);
    }
	return J;
}

//' Calcule les statistiques A, B et C de l'etape E de l'algo EM
// [[Rcpp::export]]
Rcpp::List ABC(const arma::vec & alpha,const arma::mat & S, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight) {
  unsigned p{S.n_rows}; 


  arma::mat exit_vect = Exit(S);
  
  arma::vec A = arma::zeros(p);
  arma::vec B = arma::zeros(p);
  arma::mat C = arma::zeros(p,p + 1);
  
  arma::mat avector(1,p);
  arma::mat bvector(p,1);
  arma::mat cmatrix(p,p);
  arma::mat aux_exp(p,p);
  
  arma::mat aux_mat(1,1);
  
  arma::mat J(2 * p,2 * p);
  arma::mat s_prod_alpha(p,p);
  s_prod_alpha = exit_vect * alpha.t();

  J = matrix_vanloan(S, S, s_prod_alpha);
  
  double a = max_diagonal(J * (-1));
  
  // Ca fait 8
  int m{find_n(0.00001, 1)};

  
  std::vector<arma::mat> aux_vect;
  
  vector_of_matrices(aux_vect, J, a, m);
  

  double density{0.0};
  
  // E-step
  // Uncensored data
  for (int k{0}; k < obs.size(); ++k) {

    
    double x{obs[k]};
    
    J = expSxUNI(x, m, aux_vect, a);
    
    for (int i{0}; i < p; ++i) {
      for (int j{0}; j < p; ++j) {
        aux_exp(i,j) = J(i,j);
        cmatrix(i,j) = J(i,j + p);
      }
    }
	avector = alpha.t() * aux_exp;
    bvector = aux_exp * exit_vect;
    aux_mat = alpha.t() * bvector;
    density = aux_mat(0,0);
    
    // E-step
    for (int i{0}; i < p; ++i) {
      A(i) += alpha[i] * bvector(i,0) * weight[k] / density;
      C(i,p) += avector(0,i) * exit_vect(i,0) * weight[k] / density;
      B(i) += cmatrix(i,i) * weight[k] / density;
      for (int j{0}; j < p; ++j) {
        C(i,j) += S(i,j) * cmatrix(j,i) * weight[k] / density;
      }
    }
  }
   Rcpp::List result(3);
   result[0] = Rcpp::wrap(A);
  result[1] = Rcpp::wrap(B);
  result[2] = Rcpp::wrap(C);
  return result;
}







//' Densite loi PH
// [[Rcpp::export]]
arma::vec Densite(const arma::rowvec& alpha,const arma::mat& S,const arma::colvec& s,const arma::vec& y){
	int n=y.n_elem;
	arma::vec result(n);
	
	double a = max_diagonal(S * (-1));
    int m{8};
    std::vector<arma::mat> aux_vect;
    vector_of_matrices(aux_vect, S, a, m);

	for (int i; i<n;++i){
		result(i) = arma::accu(alpha *  expSxUNI(y(i), m, aux_vect, a) * s);
	}
	return result;
}

//' Fonction de survie loi PH
// [[Rcpp::export]]
arma::vec Survie(const arma::rowvec& alpha,const arma::mat& S,const arma::vec& y){
	int n=y.n_elem;
	arma::vec result(n);

	double a = max_diagonal(S * (-1));
    int m{8};
    std::vector<arma::mat> aux_vect;
    vector_of_matrices(aux_vect, S, a, m);

	for (int i; i<n;++i){
		result(i) = arma::accu(alpha * expSxUNI(y(i), m, aux_vect, a) );
	}
	return result;
}

//' intensite Gompertz
double lambdaGompertz(double t, double beta){
	return exp(t*beta);
}

//' g-1 Gompertz
double invgGompertz(double t, double beta){
	return (lambdaGompertz(t,beta)-1)/beta;
}

//' Densite loi IPH-Gompertz
// [[Rcpp::export]]
arma::vec DensiteGompertz(const arma::rowvec& alpha,const arma::mat& S,const arma::colvec& s,const arma::vec& y,const double& beta){
	int n=y.n_elem;
	arma::vec result(n);

	double a = max_diagonal(S * (-1));
    int m{8};
    std::vector<arma::mat> aux_vect;
    vector_of_matrices(aux_vect, S, a, m);

	for (int i; i<n;++i){
		result(i) = lambdaGompertz(y(i),beta) * arma::accu(alpha * expSxUNI(invgGompertz(y(i),beta), m, aux_vect, a) * s);
	}
	return result;
}

//' Log-Densite loi IPH-Gompertz
// [[Rcpp::export]]
arma::vec logDensiteGompertz(const arma::rowvec& alpha,const arma::mat& S,const arma::colvec& s,const arma::vec& y,const double& beta){
	int n=y.n_elem;
	arma::vec result(n);

	double a = max_diagonal(S * (-1));
    int m{8};
    std::vector<arma::mat> aux_vect;
    vector_of_matrices(aux_vect, S, a, m);

	for (int i; i<n;++i){
		result(i) = beta * y(i) + std::log(arma::accu(alpha * expSxUNI(invgGompertz(y(i),beta), m, aux_vect, a) * s));
	}
	return result;
}


//' Fonction survie loi IPH-Gompertz
// [[Rcpp::export]]
arma::vec SurvieGompertz(const arma::rowvec& alpha,const arma::mat& S,const arma::vec& y,const double& beta){
	int n=y.n_elem;
	arma::vec result(n);

	double a = max_diagonal(S * (-1));
    int m{8};
    std::vector<arma::mat> aux_vect;
    vector_of_matrices(aux_vect, S, a, m);

	for (int i; i<n;++i){
		result(i) =  arma::accu(alpha * expSxUNI(invgGompertz(y(i),beta), m, aux_vect, a) );
	}
	return result;
}


//' intensite logLogistique
double lambdaLlogis(double t, double gamma,double theta){
	return theta*pow(t,(theta-1))/(pow(t,theta)+pow(gamma,theta));
}

//' g-1 logLogistique
double invgLlogis(double t, double gamma,double theta){
	return log(pow((t/gamma),theta)+1);
}

//' Densite loi IPH-logLogistique
// [[Rcpp::export]]
arma::vec DensiteLlogis(const arma::rowvec& alpha,const arma::mat& S,const arma::colvec& s,const arma::vec& y,const double& gamma,const double& theta){
	int n=y.n_elem;
	arma::vec result(n);

	double a = max_diagonal(S * (-1));
    int m{8};
    std::vector<arma::mat> aux_vect;
    vector_of_matrices(aux_vect, S, a, m);

	for (int i; i<n;++i){
		result(i) = lambdaLlogis(y(i),gamma,theta) * arma::accu(alpha * expSxUNI(invgLlogis(y(i),gamma,theta), m, aux_vect, a) * s);
	}
	return result;
}

//' Log-Densite loi IPH-logLogistique
// [[Rcpp::export]]
arma::vec logDensiteLlogis(const arma::rowvec& alpha,const arma::mat& S,const arma::colvec& s,const arma::vec& y,const double& gamma,const double& theta){
	int n=y.n_elem;
	arma::vec result(n);

	double a = max_diagonal(S * (-1));
    int m{8};
    std::vector<arma::mat> aux_vect;
    vector_of_matrices(aux_vect, S, a, m);

	for (int i; i<n;++i){
		result(i) = std::log(theta)+(theta-1)*std::log(y(i)) - std::log(pow(y(i),theta)+pow(gamma,theta)) + std::log(arma::accu(alpha * expSxUNI(invgLlogis(y(i),gamma,theta), m, aux_vect, a) * s));
	}
	return result;
}


//' Fonction survie loi IPH-logLogistique
// [[Rcpp::export]]
arma::vec SurvieLlogis(const arma::rowvec& alpha,const arma::mat& S,const arma::vec& y,const double& gamma,const double& theta){
	int n=y.n_elem;
	arma::vec result(n);

	double a = max_diagonal(S * (-1));
    int m{8};
    std::vector<arma::mat> aux_vect;
    vector_of_matrices(aux_vect, S, a, m);

	for (int i; i<n;++i){
		result(i) =  arma::accu(alpha * expSxUNI(invgLlogis(y(i),gamma,theta), m, aux_vect, a) );
	}
	return result;
}

//' Densite loi IPH-Gompertz conditionnellement covariables
// [[Rcpp::export]]
arma::vec DensiteGompertzs(const arma::rowvec& alpha,const arma::mat& S,const arma::vec& y,const double& beta, const arma::vec& mx){
	int n=y.n_elem;
	arma::vec s=Exit(S);
	arma::vec result(n);

	double a = max_diagonal(S * (-1));
    int m{8};
    std::vector<arma::mat> aux_vect;
    vector_of_matrices(aux_vect, S, a, m);

	for (int i; i<n;++i){
		result(i) = mx(i) * lambdaGompertz(y(i),beta) * arma::accu(alpha * expSxUNI(mx(i) * invgGompertz(y(i),beta), m, aux_vect, a) * s);
	}
	return result;
}

//' Log-Densite loi IPH-Gompertz conditionnellement covariables
// [[Rcpp::export]]
arma::vec logDensiteGompertzs(const arma::rowvec& alpha,const arma::mat& S,const arma::vec& y,const double& beta, const arma::vec& l_mx){
	int n=y.n_elem;
	arma::vec s=Exit(S);
	arma::vec result(n);

	double a = max_diagonal(S * (-1));
    int m{8};
    std::vector<arma::mat> aux_vect;
    vector_of_matrices(aux_vect, S, a, m);

	for (int i; i<n;++i){
		result(i) = l_mx(i) + y(i)*beta + std::log(arma::accu(alpha * expSxUNI(exp(l_mx(i)) * invgGompertz(y(i),beta), m, aux_vect, a) * s));
	}
	return result;
}

// [[Rcpp::export]]
arma::vec logDensiteGompertzs2(const arma::rowvec& alpha,const arma::mat& S,const arma::vec & y,const arma::mat& X,const double gamma0,const arma::vec &  gamma, arma::vec & beta){
	int n=y.n_elem;
	arma::vec s=Exit(S);
	arma::vec result(n);

	double a = max_diagonal(S * (-1));
    int m{8};
    std::vector<arma::mat> aux_vect;
    vector_of_matrices(aux_vect, S, a, m);

	for (int i; i<n;++i){
		double exp_x_beta = exp(arma::accu(X.row(i) * beta));
		double gamma_sum = arma::accu(X.row(i) * gamma);
		double exp_gamma0_sum = exp(gamma0 + gamma_sum);


		result(i) = exp_x_beta  + exp_gamma0_sum * y(i) + std::log(arma::accu(alpha * expSxUNI(exp_x_beta * invgGompertz(y(i),exp_gamma0_sum), m, aux_vect, a) * s));
	}
	return result;
}

// [[Rcpp::export]]
arma::vec SurvieGompertzs(const arma::rowvec& alpha,const arma::mat& S,const arma::vec& y,const double& beta, const arma::vec& mx){
	int n=y.n_elem;
	arma::vec s=Exit(S);
	arma::vec result(n);

	double a = max_diagonal(S * (-1));
    int m{8};
    std::vector<arma::mat> aux_vect;
    vector_of_matrices(aux_vect, S, a, m);

	for (int i; i<n;++i){
		result(i) =  arma::accu(alpha * expSxUNI(mx(i) * invgGompertz(y(i),beta), m, aux_vect, a));
	}
	return result;
}

// [[Rcpp::export]]
arma::vec derivBetaGompertz(const arma::rowvec& alpha,const arma::mat& S,const arma::colvec& s,const arma::vec& y,const double& beta){
	int n=y.n_elem;
	arma::vec result(n);

	double a = max_diagonal(S * (-1));
    int m{8};
    std::vector<arma::mat> aux_vect;
    vector_of_matrices(aux_vect, S, a, m);

	for (int i; i<n;++i){
		result(i) = (1+lambdaGompertz(y(i),beta)*(beta * y(i) -1))/(beta*beta) * arma::accu(alpha * S * expSxUNI(invgGompertz(y(i),beta), m, aux_vect, a) * s);
	}
	return result;
}

// [[Rcpp::export]]
arma::mat derivExpMat(const arma::mat& S,const double& t, double i, double j){
	unsigned p{S.n_rows}; 
	i=i-1;
	j=j-1;
	arma::mat I;
	I = arma::eye(p, p);

	arma::mat Eij =  I.col(i) * I.row(j);
	arma::mat J = matrix_vanloan(S, S, Eij);

	

	double a = max_diagonal(S * (-1));
    int m{8};
    std::vector<arma::mat> aux_vect;
    vector_of_matrices(aux_vect, J, a, m);
	arma::mat cmatrix(p,p);
	J = expSxUNI(t, m, aux_vect, a);
	for (int l{0}; l < p; ++l) {
      for (int k{0}; k < p; ++k) {
        cmatrix(l,k) = J(l,k + p);
      }
    }


	return cmatrix;
}

// [[Rcpp::export]]
arma::vec trans_y(const arma::vec & y,const arma::mat& X,const double gamma0,const arma::vec &  gamma, arma::vec & beta){
	int n=y.n_elem;
	arma::vec result(n);

	for (int i{0}; i < n; ++i) {


		double exp_x_beta = exp(accu(X.row(i) * beta));
		double gamma_sum = arma::accu(X.row(i) * gamma);
		double exp_gamma0_sum = exp(gamma0 + gamma_sum);
		result(i) = exp_x_beta * invgGompertz(y(i),exp_gamma0_sum);
    }


	return result;
}

// [[Rcpp::export]]
arma::mat PIx(const arma::rowvec& alpha,const arma::mat& S,const arma::vec& y,const double& beta){
	int n=y.n_elem;
	int p=alpha.n_elem;
	arma::mat result(n,p);

	double a = max_diagonal(S * (-1));
    int m{8};
    std::vector<arma::mat> aux_vect;
    vector_of_matrices(aux_vect, S, a, m);

	for (int i; i<n;++i){
			result.row(i) =  alpha * expSxUNI(invgGompertz(y(i),beta), m, aux_vect, a)/arma::accu(alpha * expSxUNI(invgGompertz(y(i),beta), m, aux_vect, a) );
	}
	return result;
}
