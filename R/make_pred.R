MODULE_PRED <- modules::module(
  
  {
    
    #import("doParallel")
    #import("foreach")
    #
    
    import("parallel")
   
    vector_array_to_matrix_matrix <- function(data_array, row_dim=NULL, col_dim=NULL) {
      
      # 
      array_dim <- dim(data_array)[length(data_array)]
      
      if(is.null(row_dim)) {
        row_dim <- sqrt(array_dim) 
        
        temp_result <- all.equal(round(row_dim), row_dim)
        
        if(!is.logical(temp_result)) {
          stop("The matrix is not square, you need spesify the at least the row dimension")
        }
        
        row_dim <- round(row_dim)
        
      }
      
      if(is.null(col_dim)) {
        col_dim = row_dim
      }
      
      result_aray <- array(NA_real_, dim = c(dim(data_array)[1: (length(data_array)-1) ], row_dim, col_dim))
      
      temp_df <- tibble::tibble(.rows = length(data_array)-1)
      
      for( m_idx in 1:(length(data_array)-1) ) {
        
        temp_df[[m_idx]] <- 1:(dim(data_array)[m_idx])
        
      }
      
      temp_df <-  temp_df |> tidyr::complete()
      
      for(m_idx in 1:nrow(temp_df)) {
        
        temp_index <- temp_df[m_idx,]
        
        result_aray[temp_index, ] <- data_array[temp_index, ] |> matrix(nrow = row_dim, ncol = col_dim)
        
      }
      
      return(result_aray)
    }
    
    get_single_simul_data <- function(extracted_data, index) {
      
      
      max_simul_index <- dim(extracted_data[[1]])[1]
      
      if(index > max_simul_index) {
        stop("The index is greater than the total number of simulations")
      }
      
      
      for(m_name in names(extracted_data)) {
        
      }
      
    }
     
    make_pred <- function(model, initial_event_time, ncores=NULL, end_event_time=NULL) {
      
      end_time  <- model$get_end_time()
      
      if(is.null(ncores))  {
        ncores <- round(parallel::detectCores()/2)
        ncores <- getOption("cl.cores", ncores)
      }
      
      if(is.null(end_event_time)) {
        end_event_time <- end_time
      }
      
      
      cl <- makeCluster(ncores)
      
      clusterExport(cl, "xx")
      
      stopCluster(cl)
      
      

      
      
      
    }
    
  }
)