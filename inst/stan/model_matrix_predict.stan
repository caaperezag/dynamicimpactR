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


// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {

  cov_matrix[use_predefined_stations_var ? 0 : K] sigma_entry_obs_stations; // Sigma
  cov_matrix[use_predefined_sensors_var ? 0 : R]  sigma_entry_obs_sensores; // v_t

  cov_matrix[P] level_sigma_variables; // w_t

  // TODO mejorar esto para la prediccion
  // cov_matrix[share_stations_var ? 0 : K] level_sigma_stations;


  vector[P*K] theta[N]; // theta is defined as vector, to use the multivariate normal distribution


}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {

}


generated quantities {

  //moved from transformed parameters
  vector[R*K]  mu[N_before];
  //

  // these variables are computed after MCMC
  matrix[R, K] y_pred[N]; // matricial version of y
  matrix[P, K] theta_pred[N]; // matricial version of theta
  matrix[R, K] mu_pred[N]; // matricial version of mu
  matrix[R, K] difference[N]; // array of diference matrices
  matrix[R, K] cumsum_difference[N]; // array of cumsum matrices, it starts of the initial time.
  matrix[R, K] cumsum_only_after[N]; // array of cumsum matrices, it starts after the intervention.


  matrix[R, K] arco_only_after[N]; // Artificial counter factual
  //vector[N] arco_only_after_aggregated;


  // this was moved from transformed parameters

  // mu is also a vector to. Because X is matrix  theta is converted to matrix form
  // before the multiplication, the result is then converted to a vector
  mu[1] = to_vector( X[1] * to_matrix(theta[1], P, K) );

   // repeat for all times before the intervention
   for (t in 2:N_before) {

      mu[t] =   to_vector( ( X[t] * to_matrix(theta[t], P, K)) );

   }


  // get all paramters in their correct dimension
  for(t in 1:N_before) {

    theta_pred[t] = to_matrix(theta[t], P, K);

    mu_pred[t] = to_matrix(mu[t], R, K);

    // repeat the processs in the models block.
    // here we need also need to use multi_normal_rng to reproduce the variation
    if(use_predefined_sensors_var) {

      if(use_predefined_stations_var) {

        y_pred[t] = to_matrix(multi_normal_rng(to_vector( mu_pred[t] ), kronecker_prod(predefined_stations_var, predefined_sensors_var)), R, K);

      } else {

        y_pred[t] = to_matrix(multi_normal_rng(to_vector( mu_pred[t] ), kronecker_prod(sigma_entry_obs_stations, predefined_sensors_var)), R, K);

      }


    } else {

      if(use_predefined_stations_var) {

        y_pred[t] = to_matrix(multi_normal_rng(to_vector( mu_pred[t] ), kronecker_prod(predefined_stations_var, sigma_entry_obs_sensores)), R, K);

      } else {

        y_pred[t] = to_matrix(multi_normal_rng(to_vector( mu_pred[t] ), kronecker_prod(sigma_entry_obs_stations, sigma_entry_obs_sensores)), R, K);

      }

    }

    // diference variable
    difference[t] = y[t] - y_pred[t];
    // the cumsum_only_after is the cero matrix for all times before the intervention
    cumsum_only_after[t] = to_matrix(rep_vector(0, R*K), R, K);
  }

  // predictions
  for(t in (N_before+1):N) {

    // the variable keep_theta_static_for_prediction controls if theta is evolves after the intervention or not.
    if(keep_theta_static_for_prediction) {

      // if theta does not evolves after the prediction, make it equal to theta in the time before the intervention.
      theta_pred[t] = theta_pred[N_before];

    } else {
      // repeat the process defined model block
      if(share_stations_var) {

        if(use_predefined_stations_var) {

          theta_pred[t] = to_matrix(multi_normal_rng(to_vector(theta_pred[t-1]), kronecker_prod(predefined_stations_var, level_sigma_variables)), P, K);

        } else {

          theta_pred[t] = to_matrix(multi_normal_rng(to_vector(theta_pred[t-1]), kronecker_prod(sigma_entry_obs_stations, level_sigma_variables)), P, K);

        }


      } else {
        theta_pred[t] = to_matrix(multi_normal_rng(to_vector(theta_pred[t-1]), kronecker_prod(level_sigma_stations,     level_sigma_variables)), P, K);
      }
    }

    mu_pred[t] = (X[t] * theta_pred[t]);

     // repeat the process defined in the model block
     if(use_predefined_sensors_var) {

      if(use_predefined_stations_var) {

        y_pred[t] = to_matrix(multi_normal_rng(to_vector( mu_pred[t] ), kronecker_prod(predefined_stations_var,  predefined_sensors_var)), R, K);

      } else {

        y_pred[t] = to_matrix(multi_normal_rng(to_vector( mu_pred[t] ), kronecker_prod(sigma_entry_obs_stations, predefined_sensors_var)), R, K);

      }

    } else {

      if(use_predefined_stations_var) {

        y_pred[t] = to_matrix(multi_normal_rng(to_vector( mu_pred[t] ), kronecker_prod(predefined_stations_var,  sigma_entry_obs_sensores)), R, K);

      } else {

        y_pred[t] = to_matrix(multi_normal_rng(to_vector( mu_pred[t] ), kronecker_prod(sigma_entry_obs_stations, sigma_entry_obs_sensores)), R, K);

      }

    }

    // diference variable
    difference[t] = y[t] - y_pred[t];

  }

  // compute cumsums
  cumsum_difference                 = cumsum_matrix(difference);
  cumsum_only_after[(N_before+1):N] = cumsum_matrix(difference[(N_before+1):N]);


  for(t in (N_before+1):N) {
    //arco_only_after[t] = (1.0/(t-N_before+1.0)) * cumsum_only_after[t];
    arco_only_after[t] = (1.0/(t-N_before)) * cumsum_only_after[t];
    //arco_only_after_aggregated[t] = mean(arco_only_after[t]);
  }
}


