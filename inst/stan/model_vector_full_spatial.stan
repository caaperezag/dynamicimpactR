//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.


functions {
  
  vector[] cumsum_vector(vector[] input_vector) {
    
    int N_SIZE = dims(input_vector)[1];

    vector[dims(input_vector)[2]] result_vector[N_SIZE];
    
    result_vector[1] = input_vector[1];

    for(t in 2:(N_SIZE) ) {
       result_vector[t] = result_vector[t - 1] + input_vector[t];

    }
    return(result_vector);
  }
  
}

data {
  
  
  
  int<lower=0> N; // total of observations
  int<lower=0> N_before;
  
  int<lower=0> K; // number of stations
  int<lower=0> P; // number of regresors
  int<lower=1> J; // number of spatial kernels
  
  vector[K] Y[N];
  vector[P] X[N];
  
  vector[2] COORDINATES[K];

  real coordinates_lower;
  real coordinates_upper;
  
  int<lower=0, upper=1> use_predefined_stations_var; // use a predefined (pass by the user) as the observation variance of stations
  matrix[K , K]  predefined_stations_var; // Rstan does not allow to leave this parameter empty, even when is not in use.
  
  
  
}

transformed data {

  int N_after = N - N_before; // number of times after the intervention
  
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  cov_matrix[use_predefined_stations_var ? 0 : K] sigma_entry_obs_stations;
  // vector[J*P*K] theta_vec[N_before]; //theta es (J*P) X K
  vector[P*K] theta_vec[J, N_before]; //theta es (J*P) X K
  // cov_matrix[J*P*K] theta_vec_cov_matrix;
  // cov_matrix[P*K] theta_vec_cov_matrix[J];
  cov_matrix[P*K] theta_vec_cov_matrix;
  
  vector<lower=coordinates_lower, upper=coordinates_upper>[2] kernels[J];
  
  cov_matrix[2] sigma_kernels;

  vector[P*K] theta_vec_no_espatial[N_before];
  cov_matrix[P*K] theta_vec_cov_matrix_no_spatial;

  real<lower=0, upper=1> scale_spatial_param;

}

transformed parameters {
  
  vector[K] pi_mixture[J];


    // pi_mixture =   rep_vector(0, J) ;

    for(j in 1:J) {

      for(k in 1:K) {

        pi_mixture[j][k] =  exp(multi_normal_lpdf( COORDINATES[k] | kernels[j], sigma_kernels));
        // pi_mixture[j] =  pi_mixture[j] + multi_normal_lpdf( COORDINATES[k] | kernels[j], sigma_kernels);

      }
      
    }
      
  
}

model {
  
  vector[K] mu[N_before];
  
  // scale_spatial_param ~ gamma(0.01, 0.01);
  // scale_spatial_param ~ normal(0, 100);
  scale_spatial_param ~ beta(0.01, 0.01);
  

  if(use_predefined_stations_var == 0) {
    sigma_entry_obs_stations ~ inv_wishart(1.0*K, diag_matrix(rep_vector(100, K)) );
  }
  
  theta_vec_cov_matrix ~ inv_wishart(1.0*P*K, diag_matrix(rep_vector(100, P*K)) );
  
  for(j in 1:J) {
       // theta_vec_cov_matrix[J] ~ inv_wishart(1.0*P*K, diag_matrix(rep_vector(100, P*K)) );
       
       theta_vec[j, 1]          ~  multi_normal(rep_vector(0, P*K) , theta_vec_cov_matrix);
       theta_vec[j, 2:N_before] ~  multi_normal(theta_vec[j, 1:(N_before-1)] , theta_vec_cov_matrix);
  }
  
  // theta_vec_cov_matrix ~ inv_wishart(1.0*P*K, diag_matrix(rep_vector(100, P*K)) );
  
  
  theta_vec_cov_matrix_no_spatial ~ inv_wishart(1.0*P*K, diag_matrix(rep_vector(100, P*K)) );
  
  theta_vec_no_espatial[1] ~  multi_normal(rep_vector(0, K*P) , theta_vec_cov_matrix_no_spatial);
  theta_vec_no_espatial[2:N_before] ~  multi_normal(theta_vec_no_espatial[1:(N_before-1)] , theta_vec_cov_matrix_no_spatial);
  
  //  https://mc-stan.org/docs/2_22/stan-users-guide/multivariate-outcomes.html
  for (t in 1:N_before) {
    // mu[t] = X[t] * theta[t];
    // mu[t] = (X[t] .* theta[t]) * rep_vector(1, P);
    // print(X[t] .* theta[t]);
    // mu[t] = (to_matrix(theta_vec[t], P, K)') * X_FIT[t]  ;
    // mu[t] = rep_vector(0, K);
    for(j in 1:J) {

       mu[t] =  pi_mixture[j] .* ((to_matrix(theta_vec[j, t], P, K)')*X[t]);
      // mu[t] =  pi_mixture[j] * ((to_matrix(theta_vec[j, t], P, K)')*X[t]);
    }
    
    // mu[t] =  (scale_spatial_param*mu[t]) + ((to_matrix(theta_vec_no_espatial[t], P, K)')*X[t]);
    mu[t] =  (scale_spatial_param*mu[t]) + ( (1-scale_spatial_param) * (to_matrix(theta_vec_no_espatial[t], P, K)')*X[t]);
  }
  
  
  sigma_kernels ~ inv_wishart(1.0*2, diag_matrix(rep_vector(100, 2)) );
  
  if(use_predefined_stations_var) {
    
     Y[1:N_before] ~ multi_normal(mu[1:N_before] , predefined_stations_var);
     
   } else {
     
     Y[1:N_before] ~ multi_normal(mu[1:N_before] , sigma_entry_obs_stations);
     
   }
}

