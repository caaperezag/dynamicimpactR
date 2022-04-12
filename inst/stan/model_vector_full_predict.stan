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
  
}

generated quantities {
  
  vector[K] Y_pred[N]; // matricial version of y
  vector[K*P] theta_vec_pred[N]; // matricial version of theta
  vector[K] mu[N];
  
  vector[K] difference[N];
  vector[K] cumsum_difference[N]; // array of cumsum matrices, it starts of the initial time.
  vector[K] cumsum_only_after[N]; // array of cumsum matrices, it starts after the intervention.
  vector[K] arco_only_after[N]; //
  vector[N] arco_only_after_aggregated;
  
  
  Y_pred[1:N_before] = Y[1:N_before];
  theta_vec_pred[1:N_before] = theta_vec[1:N_before];
  
  for (t in (N_before+1):N) {
    print("n_before: ", N_before);
    print("t: ", t);
    theta_vec_pred[t] = theta_vec[N_before];
  }

  for (t in 1:N) {
    mu[t] = (to_matrix(theta_vec_pred[t], P, K)') * X[t]  ;
  }
  
  if(use_predefined_stations_var) {
    
    Y_pred = multi_normal_rng(mu, predefined_stations_var);

  } else {
    Y_pred = multi_normal_rng(mu, sigma_entry_obs_stations);
  }
  
  for (t in 1:N) {
    difference[t] = Y[t] - Y_pred[t];
  }
  
  cumsum_only_after[1:N_before] = rep_array(rep_vector(0, K), N_before);
  cumsum_only_after[(N_before+1):N] = cumsum_vector(difference[(N_before+1):N]);
  
  cumsum_difference = cumsum_vector(difference);
  
  arco_only_after[1:N_before] =  rep_array(rep_vector(0, K), N_before);
  arco_only_after_aggregated[1:N_before] = rep_vector(0, N_before);
  
  for(t in (N_before+1):N) {
    arco_only_after[t] = (1.0/(t-N_before+1.0)) * difference[t];
    arco_only_after_aggregated[t] = mean(arco_only_after[t]);
  }
  
  
}
