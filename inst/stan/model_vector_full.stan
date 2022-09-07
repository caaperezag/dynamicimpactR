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

  vector[K] Y[N];
  vector[P] X[N];


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
  vector[P*K] theta_vec[N_before];
  cov_matrix[P*K] theta_vec_cov_matrix;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {

  vector[K] mu[N_before];

  //  https://mc-stan.org/docs/2_22/stan-users-guide/multivariate-outcomes.html
  for (t in 1:N_before) {
    mu[t] = (to_matrix(theta_vec[t], P, K)') * X[t]  ;
  }

  if(use_predefined_stations_var == 0) {
    sigma_entry_obs_stations ~ inv_wishart(1.0*K, diag_matrix(rep_vector(100, K)) );
  }

  theta_vec_cov_matrix ~ inv_wishart(1.0*K*P, diag_matrix(rep_vector(100, K*P)) );

  theta_vec[1]          ~  multi_normal(rep_vector(0, K*P) , theta_vec_cov_matrix);
  theta_vec[2:N_before] ~  multi_normal(theta_vec[1:(N_before-1)] , theta_vec_cov_matrix);

  if(use_predefined_stations_var) {

     Y[1:N_before] ~ multi_normal(mu[1:N_before] , predefined_stations_var);

   } else {

     Y[1:N_before] ~ multi_normal(mu[1:N_before] , sigma_entry_obs_stations);

   }
}


