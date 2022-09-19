
MODULES_IC_SIMULATION <- modules::module({
  
    import("keep", "karray")
    sum_of_normals_var <- function(coefficients, cov_matrix) {
      
      N <- length(coefficients)
      
      
      if(length(dim(cov_matrix)) != 2) {
        stop("cov_matrix must be a 2D array")
      } else {
        
        if( any(dim(cov_matrix) != N ) ) {
          
          stop(paste0("cov_matrix must be a ", N,"X",N, " array"))
          
        }
        
      }
      
      part_1 <- sum(coefficients*coefficients*diag(cov_matrix))
      
      part_2 <- 0
      
      for(j in 1:N) {
        
        for(i in 1:j) {
          if(i == j) {next}
          coefficients[i]*coefficients[j]*cov_matrix[i,j]
        }
        
      }
      
      part_2 <- 2*part_2
      
      
      final_result <- part_1 + part_2
      
      return(final_result)
      
    }
    
    
    
    make_theta_based <- function(model_result, n_simul, y_scaled_data=NULL, use_percent=TRUE, m_weights=NULL, MODULES_SCALE=NULL) {
      
      # browser()

      n_time <- dim(model_result$y_t_after)[1]
      n_est  <-  dim(model_result$y_t_after)[3]
      n_cont <-  dim(model_result$y_t_after)[2]
      n_pred_vars <- dim(model_result$M_t)[2]


      
      N_before <- dim(model_result$y_t_before)[1]
      N_full <- dim(model_result$y_t_full)[1]
      

      result_array <- karray(NA_real_, dim=c(n_simul, n_time, n_cont, n_est))
      percent_array <- karray(NA_real_, dim=c(n_simul, n_time, n_cont, n_est))
      error_array <- karray(NA_real_, dim=c(n_simul, n_time, n_cont, n_est))
      theta_array <- karray(NA_real_, dim=c(n_simul, n_time, dim(model_result$M_t)[2], dim(model_result$M_t)[3]))
      
      
      if(is.null(m_weights)) {
        m_weights = rep(1, n_est)
      }
      
      if(is.null(dim(m_weights))) {
        m_weights <- rep(m_weights, n_time) |> 
          matrix(nrow=n_time, ncol=n_est, byrow = T)
      }
      
      if(use_percent) {
        m_y <- model_result$E_UP[N_before:N_full,,]/model_result$y_t_full[N_before:N_full,,]
      } else {
        m_y <- model_result$E_UP[N_before:N_full,,]
      }
      

      
      beta <- diag(model_result$discount**(-1/2), dim(model_result$C_t)[2])
      
      
      for(t in 1:n_time) {
        

        
        
        #W_t <- (beta %*% model_result$C_t[N_before+t-2,,] %*% beta) - model_result$C_t[N_before+t-2,,]
        W_t <- (beta %*% model_result$C_t[N_before-1,,] %*% beta) - model_result$C_t[N_before-1,,]
        theta_old <-       MixMatrix::rmatrixnorm(n = n_simul, 
                                                  mean =  model_result$M_t[N_before+t-1,,] ,
                                               
                                                   # U =  model_result$C_t[N_before+t-1,,], #rows
                                                   U =  W_t, #rows
                                                   V = model_result$S_t[N_before-1,,] # cols
                                                 
          )

        print(paste0("t:", t))
        for(i in 1:n_simul) {

          # browser()

          # W_t <- (beta %*% model_result$C_t[N_before+t-1,,] %*% beta) - model_result$C_t[N_before+t-1,,]   
          W_t <- (beta %*% model_result$C_t[N_before,,] %*% beta) - model_result$C_t[N_before,,]   

          #print(paste0("i:", i))

          theta_array[i,t,,] <-  MixMatrix::rmatrixnorm(n = 1, 
                                                        mean = theta_old[,,i] , 
                                                        # U =  model_result$C_t[N_before+t,,], #rows
                                                        U =  W_t, #rows
                                                        #V = model_result$S_t[N_before+t,,] # cols
                                                        V = model_result$S_t[N_before,,]
                                          
                                                        
          )
          
          # temp_matrix <-  matrix(model_result$X_t[N_before+t,], ncol=1)
          
          temp_matrix <-  model_result$X_t[N_before+t,,]

          result_array[i,t,,] <- MixMatrix::rmatrixnorm(n = 1, 
                                                       mean =  ( (temp_matrix) %*% theta_array[i,t,,] ),
                                                       #U = model_result$V_t[N_before+t] |> diag(1), #rows
                                                       U =  model_result$V_t[N_before+t,,], #rows
                                                       #V = model_result$S_t[N_before+t,,] # cols,
                                                       V = model_result$S_t[N_before,,]
                                                       )
          

          
          temp_real_y <- model_result$y_t_after[t,,]
          

          
          percent_array[i, t,, ] <- (result_array[i,t,,] - temp_real_y)/temp_real_y
          error_array[i, t,, ] <-  (result_array[i,t,,] - temp_real_y)


          
        }
        

      }
      
      # browser()
      
      if(!is.null(y_scaled_data)) {
        
        for(i in 1:n_simul) {

          error_array[i,,, ] <- error_array[i,,, ] |> 
                               # MODULES_SCALE$unscale_matrix(result_list=y_scaled_data)
                               MODULES_SCALE$unscale_array_3d(result_list=y_scaled_data)
          
        }
        
        
      }
      
      
      return(list(
        y_Array = result_array,
        error_array= error_array,
        theta_array = theta_array,
        percent_array = percent_array
      ))
      
    }
    
    apply_across_time <-  function(simul_model_raw, m_fum) {
      
      array_result <- karray(NA_real_, dim = dim(simul_model_raw))
      
      n_simul <-  dim(simul_model_raw)[1]
      n_time  <-  dim(simul_model_raw)[2]
      n_est   <-  dim(simul_model_raw)[4]
      n_cont  <-  dim(simul_model_raw)[3]
      
      for(idx in 1:n_simul) {
        
        for(i in 1:n_est) {

          for(j in 1:n_cont) {

            array_result[idx,,j,i] <- simul_model_raw[idx,,j,i] |> m_fum()

          }
          
          
        }
        
        
      }
      
      return(array_result)
      
    }
    
    apply_across_time_single <-  function(simul_model_raw, m_fum) {
      
      
      
      n_simul <-  dim(simul_model_raw)[1]
      n_time  <-  dim(simul_model_raw)[2]
     
      
      array_result <- matrix(NA_real_, nrow=n_simul, ncol=n_time)
      
      for(idx in 1:n_simul) {

        array_result[idx,] <- simul_model_raw[idx,] |>  m_fum()

      }
      
      return(array_result)
      
    }
    
    extract_ic <- function(simul_model_cumsum, m_confidence_level) {
      
      
      
      n_time <- dim(simul_model_cumsum)[2]
      n_est  <- dim(simul_model_cumsum)[3]
      
      ic_result_lower <- matrix(NA_real_, nrow=n_time, ncol=n_est)
      ic_result_upper <- matrix(NA_real_, nrow=n_time, ncol=n_est)
      
      
      
      for(idx in 1:n_est) {
        
        temp_result_ic <- simul_model_cumsum[,,idx] |> apply(2, bayestestR::hdi, ci=m_confidence_level)
        
        ic_result_lower[,idx] <- temp_result_ic |>  sapply(function(x){x$CI_low }, USE.NAMES = F, simplify = T)
        ic_result_upper[,idx] <- temp_result_ic |>  sapply(function(x){x$CI_high}, USE.NAMES = F, simplify = T)
        
      }
      
      
      return(list(
        lower=ic_result_lower,
        upper=ic_result_upper
      ))
      
    }
     
    
    extract_ic_single <- function(mean_simul_model_cumsum_raw,  m_confidence_level) {
      
      
      temp_result_ic <- mean_simul_model_cumsum_raw |> apply(2, bayestestR::hdi, ci=m_confidence_level)
      
      ic_result_lower <- temp_result_ic |>  sapply(function(x){x$CI_low }, USE.NAMES = F, simplify = T)
      ic_result_upper <- temp_result_ic |>  sapply(function(x){x$CI_high}, USE.NAMES = F, simplify = T)
      
      
      return(list(
        lower=ic_result_lower,
        upper=ic_result_upper
      ))
      
    }
    
    make_sumation <- function(simul_model_raw, alpha, m_weights=NULL) {
      
      browser()
      
      n_simul <-  dim(simul_model_raw)[1]
      n_time  <-  dim(simul_model_raw)[2]
      n_est   <-  dim(simul_model_raw)[4]
      n_cont  <-  dim(simul_model_raw)[3]
      
    
      
      if(is.null(m_weights)) {
        m_weights = rep(1, n_est)
      }
      
      if(is.null(dim(m_weights))) {
        m_weights <- rep(m_weights, n_time) |> 
          matrix(nrow=n_time, ncol=n_est, byrow = T)
      }
      

      simul_model_cumsum <- simul_model_raw |>  apply_across_time(cumsum)
      

      m_confidence_level <- 1-alpha


      hdi_ic <- simul_model_cumsum |> extract_ic(m_confidence_level)

      mean_inter  <- simul_model_cumsum |> apply(c(2,3,4), mean)
      median_inter  <- simul_model_cumsum |> apply(c(2,3,4), stats::median)
      quantile_inter  <- simul_model_cumsum |> apply(c(2,3,4), function(x) {
        x  |> stats::quantile(probs=m_confidence_level)  |> as.numeric()
      } )
      lower_inter <- hdi_ic$lower
      upper_inter <- hdi_ic$upper
      
      
      mean_simul_model_cumsum_raw <- simul_model_raw |>
                                     apply(c(1, 2), mean) |>  
                                     apply_across_time_single(cumsum)
      
      
      aggrate_inter <- mean_simul_model_cumsum_raw |>  extract_ic_single(m_confidence_level)
      
      lower_aggregate <-   aggrate_inter$lower
      upper_aggregate <-   aggrate_inter$upper
      mean_aggregate  <-   mean_simul_model_cumsum_raw |> apply(2, mean)
      median_aggregate  <-   mean_simul_model_cumsum_raw |> apply(2, stats::median)
      quantile_aggregate  <- mean_simul_model_cumsum_raw |> apply(2, function(x) {
        x  |> stats::quantile(probs=m_confidence_level)  |> as.numeric()
      })

      
      
      return(
        list(
          mean_inter = mean_inter,
          median_inter = median_inter,
          quantile_inter = quantile_inter,
          lower_limit = lower_inter,
          upper_limit = upper_inter,

          cumsum_result = simul_model_cumsum,
          cumsum_aggregate = mean_simul_model_cumsum_raw,
          t = 1:n_time,
          
          aggregate = data.frame(
            lower_limit = lower_aggregate,
            upper_limit = upper_aggregate,
            mean_aggregate = mean_aggregate,
            median_aggregate = median_aggregate,
            quantile_aggregate = quantile_aggregate,
            t = 1:n_time
          )
          
        )
      )
      
      
    }
    
    
    
    
   
    
    
  }
  
)