UTILS <- modules::module({
  
  apply_across_simulation <- function(m_3d_array) {
    
    
    data_result <- m_array |> apply(c(2, 3), mean)
    
    return(data_result)

  }
  
  
  matrix_to_df <- function(input_matrix, row_names, colum_names) {
    
    new_df <- data.frame(
      value = as.vector(input_matrix),
      row = rep(row_names, length(input_matrix)/length(row_names)),
      col = rep(colum_names, each = length(row_names))
    )
    
    return(new_df)
    
  }
  
  adjust_3d_array <- function(x) {
    
    if(length(dim(x)) == 3) {
      return(x)
    }
    
    if(length(dim(x)) == 4) {
      
      if(dim(x)[1] == 1) {
        
        temp_result <- x[1,,,]
        
        if(length(dim(temp_result)) == 2) {
          temp_result <- array(temp_result, dim = c(dim(temp_result)[1], 1, dim(temp_result)[2]))
        }
        
        return(temp_result)
      }
    }
    
    stop("The 3d array dimensions are inconsistent")
    
  }
  
  array_3d_to_df <- function(array_3d, row_names, colum_names, time_index) {
    
    array_3d <- adjust_3d_array(array_3d)
    
    
    result_df <- list()
    
    N <- dim(array_3d)[1]
    
    if(N != length(time_index)) {
      stop("The 3d array first dimension is diferent form the length of time_index")
    }
    
    for(i in 1:N) {
      
      index <- length(result_df)+1
      
      result_df[[index]] <- matrix_to_df(input_matrix=array_3d[i,,],
                                         row_names=row_names,
                                         colum_names=colum_names)
      
      result_df[[index]]$t <- time_index[i]
      
    }
    
    
    result_df <- do.call(rbind, result_df)
    
    return(result_df)
    
  }

  # taken from :
  # https://stackoverflow.com/questions/26553602/release-memory-by-gc-in-silence
  gc_quiet <- function(quiet = TRUE, ...) {
    if(quiet) invisible(gc()) else gc(...)
  }
  
  
})