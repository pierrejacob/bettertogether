#include <RcppEigen.h>
using namespace Rcpp;
using namespace std;


// made for export to R
// [[Rcpp::export]]
NumericVector plummer_full_posterior_(NumericMatrix theta1s, NumericMatrix theta2s, NumericVector nhpv,
                                      NumericVector Npart, NumericVector ncases, NumericVector Npop_normalized,
                                      double theta2_mean_prior, double theta2_sd_var){
  int n = theta1s.nrow();
  NumericVector evals(n);
  std::fill(evals.begin(), evals.end(), 0);
  double mu, logmu;
  for (int i = 0; i < n; i++){
    for (int j = 0; j < 13; j++){
      // prior module 1
      if (theta1s(i,j) < 0 || theta1s(i,j) > 1){
        evals(i) = -INFINITY;
      }
    }
    if (evals(i) != -INFINITY){
      // prior module 2
      evals(i) += (-0.5/theta2_sd_var) * (theta2s(i,0) - theta2_mean_prior) * (theta2s(i,0) - theta2_mean_prior);
      evals(i) += (-0.5/theta2_sd_var) * (theta2s(i,1) - theta2_mean_prior) * (theta2s(i,1) - theta2_mean_prior);
      for (int j = 0; j < 13; j++){
        // likelihood module 1
        // k log p + (n-k) log(1-p)
        evals(i) += nhpv(j) * log(theta1s(i,j)) + (Npart(j) - nhpv(j)) * log(1 - theta1s(i,j));
        // likelihood module 2
        logmu = theta2s(i,0) + theta1s(i,j) * theta2s(i,1) + Npop_normalized(j);
        mu = exp(logmu);
        // ncases ~ Poisson(mu)
        // ncases * log(mu) - mu
        evals(i) += ncases(j) * logmu - mu;
      }
    }
  }
  return evals;
}

// [[Rcpp::export]]
NumericVector plummer_module2_conditional_(NumericVector theta1, NumericMatrix theta2s,
                                          NumericVector ncases, NumericVector Npop_normalized){
  int n = theta2s.nrow();
  NumericVector evals(n);
  std::fill(evals.begin(), evals.end(), 0);
  double mu, logmu;
  // double theta2_var = (theta2_sd_prior*theta2_sd_prior);
  for (int i = 0; i < n; i++){
    if (evals(i) != -INFINITY){
      // prior module 2
      // evals(i) += (-0.5/theta2_var) * (theta2s(i,0) - theta2_mean_prior) * (theta2s(i,0) - theta2_mean_prior);
      // evals(i) += (-0.5/theta2_var) * (theta2s(i,1) - theta2_mean_prior) * (theta2s(i,1) - theta2_mean_prior);
      for (int j = 0; j < 13; j++){
        // likelihood module 2
        logmu = theta2s(i,0) + theta1(j) * theta2s(i,1) + Npop_normalized(j);
        mu = exp(logmu);
        // ncases ~ Poisson(mu)
        // ncases * log(mu) - mu
        evals(i) += ncases(j) * logmu - mu;
      }
    }
  }
  return evals;
}

// [[Rcpp::export]]
NumericVector plummer_module2_loglikelihood_(NumericMatrix thetas1, NumericMatrix theta2s,
                                           NumericVector ncases, NumericVector Npop_normalized){
  int n = theta2s.nrow();
  NumericVector evals(n);
  std::fill(evals.begin(), evals.end(), 0);
  double mu, logmu;
  for (int i = 0; i < n; i++){
    if (evals(i) != -INFINITY){
      // prior module 2
      // evals(i) += (-0.5/theta2_sd_var) * (theta2s(i,0) - theta2_mean_prior) * (theta2s(i,0) - theta2_mean_prior);
      // evals(i) += (-0.5/theta2_sd_var) * (theta2s(i,1) - theta2_mean_prior) * (theta2s(i,1) - theta2_mean_prior);
      for (int j = 0; j < 13; j++){
        // likelihood module 2
        logmu = theta2s(i,0) + thetas1(i,j) * theta2s(i,1) + Npop_normalized(j);
        mu = exp(logmu);
        // ncases ~ Poisson(mu)
        // ncases * log(mu) - mu
        evals(i) += ncases(j) * logmu - mu;
      }
    }
  }
  return evals;
}
