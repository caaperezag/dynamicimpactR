MODULES_SCALE <- modules::module({
  
    
    import('stats')
    
    
    scale_matrix <- function(input_matrix) {

      is_3d  <-  length(dim(input_matrix)) == 3

      if(is_3d) {

        input_matrix  <- input_matrix[,1,]

      }
      
      result_matrix <- input_matrix |> scale() 
      
      temp_matrix <- result_matrix
      attr(temp_matrix,"scaled:center") <- NULL
      attr(temp_matrix,"scaled:scale") <- NULL

      if(is_3d) {
        temp_matrix  <- temp_matrix  |> 
                        array( dim=c(dim(input_matrix)[1], 
                                     1 ,
                                     dim(input_matrix)[2]) )
      }
      
      result_list <- list(
        'scaled_matrix' = temp_matrix,
        'original_means' = attr(result_matrix, "scaled:center"),
        'original_scale' = attr(result_matrix, "scaled:scale"),
        'is_3d' = is_3d
      )
      
      
      
      return(result_list)
      
    }
    
    unscale_matrix <- function(input_matrix, result_list) {
      
      
      # n_row <- result_list$scaled_matrix |> nrow()
      # n_col <- result_list$scaled_matrix |> ncol()
      # 
      
      # browser()
      
      is_3d  <-  length(dim(input_matrix)) == 3

      if(is_3d) {
        # browser()
        input_matrix  <- input_matrix[,1,]

      }

      n_row <- input_matrix |> nrow()
      n_col <- input_matrix |> ncol()
      
      result_matrix <- matrix(NA_real_, nrow=n_row, ncol = n_col)
      
      for(idx in 1:n_col) {
        
        # browser()
        
        temp_mean <- result_list$original_means[idx] |> as.numeric()
        temp_sd   <- result_list$original_scale[idx] |> as.numeric()
        
        result_matrix[,idx] <- (input_matrix[,idx]*temp_sd) + temp_mean
        
      }

      
      if(is_3d) {
        result_matrix  <- result_matrix  |> 
                          array(dim=c(dim(input_matrix)[1], 
                                1,
                                dim(input_matrix)[2]))
      }
      
      return(result_matrix)
      
    }
    
    
    unscale_array_3d <- function(input_array, result_list, m_diff_array=NULL) {
      
      # browser()
      
      # browser()
      m_result <- array(NA_real_, dim=dim(input_array))
      
      m_simul <- dim(m_result)[1]
      
      
      for(i in 1:m_simul) {
        
        m_result[i,,] <- input_array[i,,] |> unscale_matrix(result_list)

        if(!is.null(m_diff_array)) {
           # browser()

          m_result[i,,]  <- m_diff_array - m_result[i,,]
        }
        
      }
      
      
      return(m_result)
      
    }
    
    # input_array is in the original scale, 
    unscale_cumsum  <-  function(input_array, result_list, start_event= 0, m_diff_array=NULL) {
      
      # browser()
      
      m_result <- array(NA_real_, dim=dim(input_array))
      
      m_simul <- dim(m_result)[1]
      m_elem <- dim(m_result)[3]
      
      m_unscaled_array <- unscale_array_3d(input_array=input_array, 
                                           result_list=result_list,
                                           m_diff_array = m_diff_array)
      
      for(i in 1:m_simul) {
        
        
        for(j in 1:m_elem) {
          
          temp_array <- m_unscaled_array[i,,j]
          
          if(start_event > 0) {
            # browser()
            temp_array[1:start_event] <- 0
          }
          
          
          m_result[i,,j] <- temp_array |> cumsum()
        }
        
        
      }
      
      return(m_result)
      
    }

    
    
    
}
  

  
  
  
)