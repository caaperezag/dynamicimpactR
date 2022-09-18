
MODULE_IMPACT <-  modules::module(
  
  
  {
  
    import("keep", "karray")
    import("abind", "abind")
    import('stats')
    
    MATRIX_LIST <- c("ml", "identity")

    # DEFAULT_DISCOUNTS <- seq(from = 0.7, to = 0.97, length.out=100)
    DEFAULT_DISCOUNTS <- seq(from = 0, to = 1, length.out=1000)

    get_variable_matrix <- function(matrix_type, Y_data, event_initial= NULL) {
    
    
      if( !(matrix_type %in% MATRIX_LIST) ) {
        stop("unknow matrix type")
      }
      
      if(is.null(event_initial)) {
        event_initial <-  dim(Y_data)[1]
      } else {
        event_initial <- min(dim(Y_data)[1], event_initial)
      }
      
      
      n_dim <- Y_data |> dim()  |> length()
        
      if( (n_dim %in% c(2, 3)) ) {
        
        if(matrix_type == "ml") {
          
          if(n_dim == 3) {
            
            result <-  estamate_ml_from_array( Y_data [1:event_initial,,] )$U
            
          } else if (n_dim == 2 ) {
            Y_data <- Y_data |> 
              array(dim=c(dim(Y_data)[1], 1, dim(Y_data)[2]))
            
            result <-  estamate_ml_from_array( Y_data [1:event_initial,,] )$v
          } 
          
          
        } else  {

          if(n_dim == 2 ) {
            result <-  diag(dim(Y_data)[2])
          } else if(n_dim == 3) {

              if(dim(Y_data)[2] == 1) {
                result <-  diag(dim(Y_data)[3])
              } else {
                result <-  diag(dim(Y_data)[2])
              }

          }

        }
        
        
        
      } else {
        
        stop("Y must be a 2D or 3D array")
        
      }
      

      return(result)
    }
  
    
    matrix_trace <- function(m_matrix) {
      
      return(sum(diag(m_matrix)))
      
    }
    
    get_model_error_abs <- function(error_array) {
      
      acc <- sum(abs(error_array))
      
      # acc <- 0
      # 
      # for(t in (1:dim(error_array)[1]) ) {
      #   
      #   acc <- acc + sum(abs(error_array[t,,]))
      #   
      # }
      
      # browser()
      
      return(acc)
      
      
    }
    
    
    
    estamate_ml_from_array <- function(y_array) {
      
      temp_list <- list()

      array_size  <- y_array |> dim() |> length()
      
      for(t in 1:dim(y_array)[1]) {
        
        if(array_size == 2) {
          temp_list[[t]] <- array(y_array[t,], dim= c(dim(y_array)[2], 1)  )
        } else {
           temp_list[[t]] <- y_array[t,,]
        }

      }
      
      result <- MixMatrix::MLmatrixnorm(temp_list, tol = 1e-8, max.iter = 1000)
      
      return(result)

    }
    
    matrix_to_array_rep <- function(m_matrix, m_size) {
      
      result_array <-  array(NA_real_, dim = c(m_size, dim(m_matrix)[1], dim(m_matrix)[2]))
      
      for(i in 1:m_size) {
        
        result_array[i,,] <- m_matrix
      }
      
      return(result_array)
      
    }
    
    
    
    build_data_list_from_data <- function(data_list_input) {
      
      # browser()
      data_list <- list()
      
      data_list$P <- dim(data_list_input$matrix_x)[3]
      data_list$Q <- dim(data_list_input$matrix_y)[2]
      
      
      Y <- data_list_input$matrix_y[,1,]
      X <- data_list_input$matrix_x[,1,]
      
      
      
      
      # browser()
      
      data_list$K <- dim(Y)[2]
      
      
      
      data_list$N        <- dim(Y)[1]
      data_list$N_before <- data_list_input$N_before
      data_list$N_after  <- data_list$N - data_list$N_before
      
      
      data_list$y <- karray(Y, dim = dim(Y))
      data_list$X <- karray(X, dim = dim(X))
      
      # browser()
      data_list$X_before <- data_list$X[1:data_list$N_before,]
      data_list$X_after  <- data_list$X[(data_list$N_before+1):data_list$N,]
      
      data_list$y_before <- data_list$y[1:data_list$N_before,]
      data_list$y_after  <- data_list$y[(data_list$N_before+1):data_list$N,]
      
      data_list$G_t  <- diag(data_list$P)
      
      data_list$G_t  <- matrix_to_array_rep(data_list$G_t, m_size=data_list$N )
      
      # browser()
      
      # V_observation <- estamate_ml_from_array( data_list_input$matrix_y[1:data_list$N_before,,] )$U
      # data_list$V_t <- matrix_to_array_rep(V_observation, m_size=data_list$N )
      # data_list$V_t <- karray(data_list$V_t, dim=dim(data_list$V_t))
      
      # if(dim(data_list_input$matrix_y)[2] > 1) {
      #   V_observation <- estamate_ml_from_array( data_list_input$matrix_y[1:data_list$N_before,,] )$U
      #   data_list$V_t <- matrix_to_array_rep(V_observation, m_size=data_list$N )
      #   data_list$V_t <- karray(data_list$V_t, dim=dim(data_list$V_t))
      # } else {
      #   # V_observation <- estamate_ml_from_array( data_list_input$matrix_y[1:data_list$N_before,,] )$var
      #   V_observation  <- 1
      #   data_list$V_t <-  rep(V_observation, data_list$N)
      #   
      # }
      
      # V_observation <-  estamate_ml_from_array( data_list_input$matrix_y[1:data_list$N_before,,] )$var
      V_observation  <- 1
      data_list$V_t <-  rep(V_observation, data_list$N)
      
      # browser()
      

      # browser()
      
      # V_observation <- estamate_ml_from_array( data_list_input$matrix_y[1:data_list$N_before,,] )$U
      
      
      
      W_t_level <- karray(rWishart(1,data_list$P ,
                                   diag(data_list$P )), 
                          dim=c(data_list$P,data_list$P,1))[,,1]
      
      # browser()
      
      data_list$W_t <- matrix_to_array_rep(W_t_level, m_size=data_list$N )
      data_list$W_t <- karray(data_list$W_t, dim=dim(data_list$W_t))
      
      # data_list$initial_ct <- diag(1000, data_list$P)
      data_list$initial_ct <- diag(20, data_list$P)
      data_list$initial_mt <- matrix(0 , nrow=data_list$P, ncol=data_list$K)
      
      return(data_list)
      
    }
    
    
    update_C_t_star <- function(G_t, C_t, W_t=NULL, discount=NULL, inv_1_2_discount=NULL, use_discount_factor=FALSE) {
      
      # browser()
      
      if(use_discount_factor) {
        
        if(is.null(inv_1_2_discount)) {
          
          beta <- diag(discount, dim(G_t)[1])
          
          # inv_1_2_discount <- matlib::mpower(beta, -1/2)
          
          if(dim(G_t)[2] > 1) {
          
          inv_1_2_discount <- matlib::mpower(diag(discount, dim(G_t)[1]), -1/2)
          
          } else {
            inv_1_2_discount <- discount**(-1/2)
          }
        


          
        } 
        
        C_star_t <- inv_1_2_discount %*% G_t %*% C_t %*% t(G_t) %*% inv_1_2_discount
        
      } else {
        C_star_t <-  W_t + (G_t %*% C_t %*% t(G_t) )
      }
      
      return(C_star_t)
    }
    
    
    run_model <- function(data_list, discount=NULL ) {
      
      if(!is.null(discount)) {
        
        return(
          run_model_single(data_list, discount)
        )
        
      }
      
      model_list <- list()
      
      N_before <- data_list$N_before
      
      for(m_discount in DEFAULT_DISCOUNTS) {
        
        # tryCatch(
          # {
            
            temp_list <- list()

            
            temp_list$model    <- run_model_single(data_list, m_discount)
            temp_list$error    <- temp_list$model$E_UP[1:N_before,,] |>   get_model_error_abs()
            temp_list$discount <- m_discount
            
            model_list[[length(model_list) + 1]] <- temp_list
          # },  error = function(e) {}
        # )
        
      }

      # para agregar el factor de bayes
      model_list  <- model_list  |> compute_bayes_factor_for_model_list(y_before=data_list$y_before)

      
      # browser()
      
      best_model <- get_best_models(model_list, score_variable='bayes')
      # print(best_model$discount)
      best_model <- best_model$model
      
      return(best_model)
      
    }
    
    run_model_single <- function(data_list, discount=NULL) {
      
      # browser()
      
      # print(update_C_t_star)
      
      G_t <- data_list$G_t
      V_t <- data_list$V_t
      W_t <- data_list$W_t
      
      if(is.null(V_t)) {
        V_t <- 1
      }
      
      
      
      N <- data_list$N
      N_before <- data_list$N_before
      N_after <- data_list$N_after
      
      K <- data_list$K
      P <- data_list$P
      Q <- data_list$Q
      R  <- Q
      
      inv_1_2_discount = NULL
      use_discount_factor = NULL
      
      # if(is.null(dim(V_t))) {
      #   
      #   if( length(V_t) != N_before ) {
      #     V_t <- rep(V_t, N_before)
      #   }
      #   
      # }
      
      # browser()

      if(!is.null(discount)) {
        
        if(dim(G_t)[2] > 1) {
          
          inv_1_2_discount <- matlib::mpower(diag(discount, dim(G_t)[2]), -1/2)
          
        } else {
          inv_1_2_discount <- discount**(-1/2)
        }
        
        use_discount_factor <- TRUE
        
      }
      
      # browser()
      
      F_t_before <- abind(NA_real_ |> array(dim = c(1, 
                                                 dim(data_list$X_before)[2], 
                                                 dim(data_list$X_before)[3])), data_list$X_before, along = 1)
      F_t_before  <- karray(F_t_before, dim=dim(F_t_before))
  
      y_t_before <- abind(NA_real_ |> array(dim = c(1, 
                                                    dim(data_list$y_before)[2], 
                                                    dim(data_list$y_before)[3])), data_list$y_before, along = 1)
      y_t_before  <- karray(y_t_before, dim=dim(y_t_before))
      
      G_t <- abind(array(NA_real_, dim = c(1, dim(G_t)[2], dim(G_t)[3])), G_t, along = 1)
      G_t <- karray(G_t, dim=dim(G_t))
      
      V_t <- abind(NA_real_ |> array(dim = c(1, dim(V_t)[2], dim(V_t)[3])), V_t, along = 1)
      V_t <- karray(V_t, dim=dim(V_t))
      
      W_t <- abind(array(NA_real_, dim = c(1, dim(W_t)[2], dim(W_t)[3])), W_t, along = 1)
      W_t <- karray(W_t, dim=dim(W_t))
      
      C_star_t <- karray(NA_real_, dim = c(N_before+N_after+1, P, P))
      C_t      <- karray(NA_real_, dim = c(N_before+N_after+1, P, P))
      
      M_star_t <- karray(NA_real_, dim = c(N_before+N_after+1, P, K))
      M_t      <- karray(NA_real_, dim = c(N_before+N_after+1, P, K))
      
      Y_UP   <- karray(NA_real_, dim = c(N_before+N_after+1, R, K))
      Y_DOWN <- karray(NA_real_, dim = c(N_before+N_after+1, R, R))
     
      A_t  <- karray(NA_real_, dim = c(N_before+N_after+1, P, R))
      
      E_UP <- karray(NA_real_, dim = c(N_before+N_after+1, R, K))
      
      S_t  <- karray(NA_real_, dim = c(N_before+N_after+1, K, K))
      
      
      C_t[1,,] <- data_list$initial_ct # matriz con la diagonal con valores grandes, con cov cero
      M_t[1,,] <- data_list$initial_mt #|> t() # para que coicida conn quintana

      S_t[1,,] <- diag(10, K) # por ahora con esto basta
      

      n_t <- 1
      
      B_t <- diag(discount**(-1/2), P)

      for(t in 2:(N_before+1)) {
        
        
        # browser()
        
        
        #F_t <-  matrix(F_t_before[t,], ncol=1)

        temp_matrix <-  F_t_before[t,,]
        

        C_star_t[t,,] <- update_C_t_star(G_t=G_t[t,,], 
                                         C_t=C_t[t-1,,],
                                         W_t= W_t[t,,], 
                                         discount=discount, 
                                         inv_1_2_discount=inv_1_2_discount, 
                                         use_discount_factor=use_discount_factor)

        # W_t[t,,] <-  (B_t %*%  C_t[t-1,,]  %*%  B_t) - C_t[t-1,,]
        
      
  
       
        M_star_t[t,,] <- G_t[t,,] %*% M_t[t-1,,]
        
        Y_DOWN[t,,] <- V_t[t,,] +   (temp_matrix %*% C_star_t[t,,] %*% t(temp_matrix) )

        Y_UP[t,,] <- temp_matrix %*% M_star_t[t,,] 

        Y_DOWN_inverse <- matlib::inv(Y_DOWN[t,,])

        A_t[t,,] <-  C_star_t[t,,]  %*% t(temp_matrix) %*% Y_DOWN_inverse

        C_t[t,,] <- C_star_t[t,,] -  (A_t[t,,] %*% Y_DOWN[t,,] %*% t(A_t[t,,]))

        E_UP[t,,] <- y_t_before[t,,] -  Y_UP[t,,]

        W_t[t,,] <-  (B_t %*%  C_t[t-1,,]  %*%  B_t) - C_t[t-1,,]

        M_t[t,,] <- M_star_t[t,,] + (A_t[t,,] %*% E_UP[t,,])
        
        
        n_t_old <- n_t
        
        n_t <- n_t + 1
        
        # S_t[t,,] <- ((1/n_t)* (discount*n_t_old*S_t[t-1,,])) +  ( t(E_UP[t,,]) %*% Y_DOWN_inverse %*% E_UP[t,,] ) 
        S_t[t,,] <- (S_t[t-1,,]) +  ( t(E_UP[t,,]) %*% Y_DOWN_inverse %*% E_UP[t,,] ) 
        #S_t[t,,] <- (1/n_t)* ( (n_t_old*S_t[t-1,,]) +  (( E_UP[t,,] %*%  t(E_UP[t,,]) ) *  (1/Y_DOWN[t]) ))
          
   
          
        
        
        
      }
      
      # browser()
      
      F_t_before <- F_t_before[2:(N_before+1),, ]
      y_t_before <- y_t_before[2:(N_before+1),, ]
      
      G_t <- G_t[2:(N_before+N_after+1),,]
      V_t <- V_t[2:(N_before+N_after+1),,]
      W_t <- W_t[2:(N_before+N_after+1),,]
      
      C_star_t <- C_star_t[2:(N_before+N_after+1), ,]
      C_t      <- C_t[2:(N_before+N_after+1), ,]
      
      
      M_star_t <- M_star_t[2:(N_before+N_after+1),,]
      M_t      <- M_t[2:(N_before+N_after+1),,]
      
      
      Y_UP   <- Y_UP[2:(N_before+N_after+1),, ]
      Y_DOWN <- Y_DOWN[2:(N_before+N_after+1),,]
      
      A_t <- A_t[2:(N_before+N_after+1),, ]
      
      
      E_UP <- E_UP[2:(N_before+N_after+1),, ]
      
      S_t <- S_t[2:(N_before+N_after+1),, ]
      
      
      F_t_after <- data_list$X_after
      y_t_after <- data_list$y_after

      browser()
      
      for(t in 1:N_after) {
        
        # browser()
        
        # W_t[t+N_before,,] <-  (B_t %*%  C_t[t+N_before-1,,]  %*%  B_t) - C_t[t+N_before-1,,]
        
        Y_DOWN[N_before+t,,] <- V_t[N_before+t,,] +   F_t_after[t,,] %*% C_star_t[N_before,,] %*% t(F_t_after[t,,]) 
    
        Y_UP[N_before+t,,] <- F_t_after[t,,] %*% M_star_t[N_before,,] 
        
        Y_DOWN_inverse <- matlib::inv(Y_DOWN[t,,])

        A_t[t+N_before,,] <-  C_star_t[N_before,,]  %*% t(F_t_after[t,,]) %*% Y_DOWN_inverse

        C_t[t+N_before,,] <- C_star_t[N_before,,] -  A_t[t,,] %*% Y_DOWN[t,,] %*% t(A_t[t,,])

        W_t[t+N_before,,] <-  (B_t %*%  C_t[t+N_before-1,,]  %*%  B_t) - C_t[t+N_before-1,,]

        E_UP[t+N_before,,] <- y_t_after[t,,] -  Y_UP[t+N_before,,]

        M_t[t+N_before,,] <- M_star_t[N_before,,] + A_t[t+N_before,,] %*% E_UP[t+N_before,,]
        
        n_t_old <- n_t
        n_t <- n_t + 1
        #S_t[t+N_before,,] <- (1/n_t)* ( (n_t_old*S_t[t+N_before-1,,]) + ((  E_UP[t+N_before,,] %*%  t(E_UP[t+N_before,,]))*(1/Y_DOWN[t+N_before]) ))

        S_t[t+N_before,,] <- (S_t[t-1,,]) +  ( t(E_UP[t+N_before,,]) %*% Y_DOWN_inverse %*% E_UP[t+N_before,,] ) 
      
        
        
        
        
        
          
          
        
      }
      
      
      # browser()
      
      result <- list(

        W_t = W_t,

        G_t = G_t,
        
        Y_DOWN=Y_DOWN,
        Y_UP=Y_UP,
        
        A_t=A_t,
        C_t=C_t,
        
        E_UP=E_UP,
        M_t=M_t,
        
        M_star_t=M_t,
        #C_star_t=C_star_t,
        
        y_t_before = y_t_before,
        y_t_after  = y_t_after,
        y_t_full   = abind(y_t_before, y_t_after, along=1),
        S_t = S_t,
        discount = discount,
        n_t = n_t,
        
        X_t = data_list$X,
        
        # df_plot=df_plot
        # 
        
        V_t = V_t
        
      )
      
      # browser()
      
      return(result)
      
      
    }
    
    
    get_log_likehood_from_model <- function(model) {
      
      # browser()

      N_BEOFRE <- dim(model$y_t_before)[1]
      log_likehood <- rep(NA_real_, N_BEOFRE)
      
      
      # if(!is.null(model$sigma)) {
      #   sigma <- model$sigma
      # }
      
      # if(!is.null(sigma)) {
      #   
      #   if(length(sigma) == 1) {
      #     sigma <- rep(N_BEOFRE, sigma)
      #   }
      # }
      
      
      for(t in 1:N_BEOFRE) {
        # log_likehood[t] <- log_matrix_variable_normal_density(y_est = model$Y_UP[t,,], 
        #                                                       y_real= model$y_t_before[t,] |> matrix(ncol=1), 
        #                                                       U = model$Y_DOWN[t,,], 
        #                                                       V = sigma)

        # browser()
        #print(t)
        # browser()
        log_likehood[t] <- MixMatrix::dmatrixnorm(x = model$Y_UP[t,,] |>  matrix(ncol=1), 
                                                  mean= model$y_t_before[t,] |>  matrix(ncol=1), 
                                                  U = model$S_t[t,,], # cols
                                                  V = model$V_t[t] |> diag(1), #rows
                                                  log = TRUE) 
        
        
      }
      
      return(log_likehood)
      
    }
    
    
    
    log_matrix_variable_normal_density <- function(y_est, y_real, U, V) {
      
      # browser()
      
      n = dim(y_real)[1] # R
      p = dim(y_real)[2] # K
      
      
      inv_U <- matlib::inv(U)
      inv_V <- matlib::inv(V)
      
      det_U <- matlib::Det(U)
      det_V <- matlib::Det(V)
      
      
      temp_matrix <- -(1/2) * inv_V %*% ( t(y_est - y_real) ) %*% inv_U %*% (y_est - y_real);
      
      result <- matrix_trace(temp_matrix) - log( (2*pi)**(n*p/2) ) - log( det_V**(n/2) ) - log( det_U**(p/2) )
      
      
      return(result)
      
    }
    
    
    
    get_best_models <-  function(model_list, score_variable) {

        # browser()

        if(score_variable == 'bayes') {
          direction="max"

        } else {
          direction="min"
        }

        score_vectors <- sapply(model_list, function(x) {
    
          return(x[[score_variable]])
          
        }, simplify = TRUE, USE.NAMES = FALSE)

        if(direction == "max") {
          m_index <- which.max(score_vectors)
        } else {
          m_index <- which.min(score_vectors)
        }  
        
       # browser()

        return(model_list[[m_index]])
      
    }
    
    get_null_model_array <- function(ml_result, y_t_before) {
      
      
      N_BEOFRE <- dim(y_t_before)[1]
      
      log_likehood <- rep(NA_real_, N_BEOFRE)
      
      for(t in 1:N_BEOFRE) {
        
        # log_likehood[t] <- log_matrix_variable_normal_density(y_est  = ml_result$mean, 
        #                                                       y_real = y_t_before[t,,], 
        #                                                       U      = ml_result$U, 
        #                                                       V      = ml_result$V)
        #                                                       
        # browser()
        log_likehood[t] <- MixMatrix::dmatrixnorm(x = ml_result$mean |>  matrix(ncol=1), 
                                                  mean= y_t_before[t,] |>  matrix(ncol=1), 
                                                  U = ml_result$U, #rows
                                                  V = ml_result$V,
                                                  log = TRUE) # cols
                                                              
      }
      
      
      return(log_likehood)
      
    }
    
    get_bayes_factor <- function(log_array_model, log_array_null_model) {
      
      # browser()
      
      
      result <- sum(log_array_model) - sum(log_array_null_model)
      
      
      #result <- exp(result)
      
      if(is.nan(result) | is.infinite(result)) {
        result <- NA
      }
      
      return(result)
      
    }
    
    
    compute_bayes_factor_for_model_list <- function(model_result_list, y_before) {
      
      # browser()
      
      ml_result <- estamate_ml_from_array(y_before) 
      
      # sigma <- ml_result$V
      
      log_likehood_array_null_model <- get_null_model_array(ml_result, y_before)
      
      for(t in 1:length(model_result_list)) {
        
          # print(t)
        
        #if(t == 2) {
          # browser()
        #}
        
        if( !all ( is.na(model_result_list[[t]]) ) ) {
          
          # browser()
          
          model_result_list[[t]]$log_likehood_array <- get_log_likehood_from_model(model=model_result_list[[t]]$model)

          # browser()
          
          model_result_list[[t]]$bayes       <- get_bayes_factor(log_array_model= model_result_list[[t]]$log_likehood_array, 
                                                                 log_array_null_model = log_likehood_array_null_model)
          
        } else {
          
          model_result_list[[t]]$log_likehood_array <- NA
          model_result_list[[t]]$bayes <- NA
          
        }
        
        
      }
      
      # browser()
      
      return(model_result_list)
      
    }
  
  
  }
  
  

)