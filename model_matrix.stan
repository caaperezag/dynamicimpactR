//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.

functions {

  // taken from https://jrnold.github.io/ssmodels-in-stan/stan-functions.html
  // caclulates the kronecker product of two matrices as defined in https://en.wikipedia.org/wiki/Kronecker_product
  matrix kronecker_prod(matrix A, matrix B) {

    matrix[rows(A) * rows(B), cols(A) * cols(B)] C = rep_matrix(0.0, rows(A) * rows(B), rows(A) * rows(B) );
    int m;
    int n;
    int p;
    int q;

    m = rows(A);
    n = cols(A);
    p = rows(B);
    q = cols(B);

    if( (rows(A) == 1) && (cols(A) == 1) ) {
        return (A[1][1] * B);
    }

    if( (rows(B) == 1) && (cols(B) == 1) ) {
        return (B[1][1] * A);
    }

    for (i in 1:m) {

      for (j in 1:n) {
        int row_start;
        int row_end;
        int col_start;
        int col_end;
        row_start = (i - 1) * p + 1;
        row_end = (i - 1) * p + p;
        col_start = (j - 1) * q + 1;
        col_end = (j - 1) * q + q;
        C[row_start:row_end, col_start:col_end] = A[i, j] * B;

    }
  }

  return C;
}

  // log density of matrix normal distribution for diference y_real - y_est
  real log_matrix_normal_density(matrix y_real, matrix y_est, matrix V, matrix U) {

    real det_V = determinant(V);
    real det_U = determinant(U);


    // use inverse_spd for speed as the V and U matrices must be semi-positive definite by definition
    matrix[dims(V)[1], dims(V)[2]] inv_V  = inverse_spd(V);
    matrix[dims(U)[1], dims(U)[2]] inv_U  = inverse_spd(U);

    int n = dims(y_real)[1]; // R
    int p = dims(y_real)[2]; // K


    matrix[p, p] temp_matrix = (-1/2.0) * (inv_V * ((y_real - y_est)') * inv_U * (y_real - y_est));


    real result_value = trace(temp_matrix) - log( pow(2.0*pi(), 1.0*n*p/2.0) ) - log(pow(det_V, 1.0*n/2.0)) - log(pow(det_U, 1.0*p/2.0));

    return(result_value);
  }

   // takes as input an array of matrices and returns the an array of "cumsum matrices"
   matrix[] cumsum_matrix(matrix[] input_matrix) {

    matrix[dims(input_matrix)[2], dims(input_matrix)[3]] result_matrix[dims(input_matrix)[1]];

    for(t in 1:(dims(input_matrix)[1]) ) {

      if(t == 1) {
        result_matrix[t] = input_matrix[t];
      } else {
        result_matrix[t] = result_matrix[t-1] + input_matrix[t];
      }

    }
    return(result_matrix);
  }

}

data {

  int N; // total of observations
  int N_before;


  int K; // number of stations
  int P; // number of regresors
  int R; // number of sensors



  matrix[R, K] y[N]; // array of response matrices
  matrix[R, P] X[N]; // array of input matrices



  // model parameters
  int<lower=0, upper=1> use_discount_factor; // currently no in use

  int<lower=0, upper=1> share_stations_var; // use the same covariance matrix for the response and evoluation

  int<lower=0, upper=1> keep_theta_static_for_prediction; // stop the evoluation of theta after the intervention

  // TODO fix-later it gives an error in R
  // matrix[use_predefined_sensors_var ? R : 0, use_predefined_sensors_var ? R : 0] predefined_sensors_var;
  int<lower=0, upper=1> use_predefined_sensors_var; // use a predefined (pass by the user) as the observation variance of sensors.
  matrix[R , R]  predefined_sensors_var; // Rstan does not allow to leave this parameter empty, even when is not in use.


  int<lower=0, upper=1> use_predefined_stations_var; // use a predefined (pass by the user) as the observation variance of stations
  matrix[K , K]  predefined_stations_var; // Rstan does not allow to leave this parameter empty, even when is not in use.

}

transformed data {

  int N_after = N - N_before; // number of times after the intervention

  vector[R*K] y_vector[N]; // requiered to use the kronecker product, rstan does not support the normal matrix variate distribution

  for(t in 1:N) {
    y_vector[t] = to_vector(y[t]); // vector version of the observation to use the multivariate normal distribution
  }

}

parameters {




  cov_matrix[use_predefined_stations_var ? 0 : K] sigma_entry_obs_stations; // Sigma
  cov_matrix[use_predefined_sensors_var ? 0 : R]  sigma_entry_obs_sensores; // v_t

  cov_matrix[P] level_sigma_variables; // w_t
  cov_matrix[share_stations_var ? 0 : K] level_sigma_stations;


  vector[P*K] theta[N]; // theta is defined as vector, to use the multivariate normal distribution
}


transformed parameters {
  vector[R*K]  mu[N_before];


  // initilize mu.

  // mu is also a vector to. Because X is matrix  theta is converted to matrix form
  // before the multiplication, the result is then converted to a vector
  mu[1] = to_vector( X[1] * to_matrix(theta[1], P, K) );

   // repeat for all times before the intervention
   for (t in 2:N_before) {

      mu[t] =   to_vector( ( X[t] * to_matrix(theta[t], P, K)) );

   }

}

model {

  // only use a distribution for sigma_entry_obs_stations if the user dones not pass a matrix
  if(use_predefined_stations_var == 0) {
    sigma_entry_obs_stations ~ inv_wishart(1.0*K, diag_matrix(rep_vector(100, K)) );
  }

  // only use a distribution for sigma_entry_obs_sensores if the user dones not pass a matrix
  if(use_predefined_sensors_var == 0) {
    sigma_entry_obs_sensores   ~ inv_wishart(1.0*R, diag_matrix(rep_vector(100, R)) );
  }


  level_sigma_variables     ~ inv_wishart(1.0*P, diag_matrix(rep_vector(100, P)) );

  // only use a distribution for sigma_entry_obs_sensores if share_stations_var is 0 (False)
  if(share_stations_var == 0) {
    level_sigma_stations ~ inv_wishart(1.0*K, diag_matrix(rep_vector(100, K)) );
  }


  // repeat for all times before the intervention
  for(t in 1:N_before) {

    if(t <= 1) {

    // initilization
    theta[t] ~ multi_normal(rep_vector(0.0, P*K), kronecker_prod(diag_matrix(rep_vector(1, K)), diag_matrix(rep_vector(1, P))));

    } else {

       // use the same stations variance for the evolution if this option is true
      if(share_stations_var) {

        // use the stations variance pass by the user if this option is true
        if(use_predefined_stations_var) {

          // kronecker_prod is use to be able to fit the  multivariate normal distribution
          theta[t] ~ multi_normal(theta[t-1], kronecker_prod(predefined_stations_var, level_sigma_variables));
        } else {
          // kronecker_prod is use to be able to fit the  multivariate normal distribution
          theta[t] ~ multi_normal(theta[t-1], kronecker_prod(sigma_entry_obs_stations, level_sigma_variables));
        }


      } else {
        // kronecker_prod is use to be able to fit the  multivariate normal distribution
        theta[t] ~ multi_normal(theta[t-1], kronecker_prod(level_sigma_stations,     level_sigma_variables));
      }

    }
  }

  // repeat for all times before the intervention
  for (t in 2:N_before) {


    // use the sensors variance pass by the user if this option is true
    if(use_predefined_sensors_var) {

      // use the stations variance pass by the user if this option is true
      if(use_predefined_stations_var) {
        // kronecker_prod is use to be able to fit the  multivariate normal distribution
        y_vector[t] ~  multi_normal(mu[t], kronecker_prod(predefined_stations_var,  predefined_sensors_var));
      } else {
        // kronecker_prod is use to be able to fit the  multivariate normal distribution
        y_vector[t] ~  multi_normal(mu[t], kronecker_prod(sigma_entry_obs_stations, predefined_sensors_var));
      }


    } else {

      // use the stations variance pass by the user if this option is true
      if(use_predefined_stations_var) {

        // kronecker_prod is use to be able to fit the  multivariate normal distribution
        y_vector[t] ~  multi_normal(mu[t], kronecker_prod(predefined_stations_var, sigma_entry_obs_sensores));
      } else {
        // kronecker_prod is use to be able to fit the  multivariate normal distribution
        y_vector[t] ~  multi_normal(mu[t], kronecker_prod(sigma_entry_obs_stations, sigma_entry_obs_sensores));
      }

    }

  }

}



