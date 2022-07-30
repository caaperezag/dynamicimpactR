MODULES_SCALE <- modules::module({


    import('stats')
    import("keep")


    scale_3d_array  <- function(input_array, use_log) {

       result_array  <- array(NA_real_, dim=dim(input_array)) |> keep::as.karray()

       indexer  <- dim(input_array)[2]

       n_cols  <- dim(input_array)[3]

       scaled_array  <- matrix(NA_real_, indexer, n_cols) |> keep::as.karray()
       center_array  <- matrix(NA_real_, indexer, n_cols) |> keep::as.karray()



       for(idx in 1:indexer ) {

        result_scale  <- scale_matrix(input_array[,idx,], use_log)
        result_array[,idx,]  <- result_scale$scaled_matrix

        scaled_array[idx,] <- result_scale$original_scale
        center_array[idx,] <- result_scale$original_means

       }


       result_list <- list(
        'scaled_matrix' = result_array,
        'original_means' = center_array,
        'original_scale' = scaled_array,
        'is_3d' = TRUE,
        'is_array' = TRUE,
        "is_log" = use_log
      )



    }

    scale_matrix <- function(input_matrix, use_log) {

      is_3d  <-  length(dim(input_matrix)) == 3

      if(is_3d) {

        input_matrix  <- input_matrix[,1,]

      }

      if(use_log) {
        if(any(input_matrix <= 0) ) {
          stop("when using log elements must be positive")
        }

        input_matrix  <- log(input_matrix)
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
        # 'is_3d' = is_3d,
        'is_3d' = FALSE, # this function is only when the original array is 2d an then reshaped to 3d.
        'is_array' = FALSE,
        "is_log" = use_log
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

      if(result_list$use_log) {

        result_matrix  <- result_matrix  |> exp()
      }

      return(result_matrix)

    }

    unscale_array_3d  <- function(input_array, result_list, m_diff_array=NULL) {

      if(result_list$is_3d) {


        m_result <- array(NA_real_, dim=dim(input_array)) |> keep::as.karray()

        indexer <- dim(input_array)[3]

        for(idx in 1:indexer) {


            temp_result_list  <- list(
                'original_means' = result_list$original_means[idx,],
                'original_scale' = result_list$original_scale[idx,],
                'is_3d' = FALSE,
                "is_log" = result_list$use_log

            )

            m_result[,,idx,]  <- unscale_array_3d_from_2d_array(input_array = input_array[,,idx,] ,
                                                                result_list = temp_result_list,
                                                                m_diff_array=m_diff_array[,idx,])

        }


        return(m_result)

      } else {

        result  <- unscale_array_3d_from_2d_array(input_array=input_array,
                                       result_list=result_list,
                                       m_diff_array=m_diff_array)

        return(result)

      }


    }

    unscale_array_3d_from_2d_array <- function(input_array, result_list, m_diff_array=NULL) {

      # browser()

      # browser()
      m_result <- array(NA_real_, dim=dim(input_array))

      m_simul <- dim(m_result)[1]


      for(i in 1:m_simul) {

        m_result[i,,] <- input_array[i,,] |> unscale_matrix(result_list)

        if(!is.null(m_diff_array)) {
           # browser()

          # if(length(m_diff_array) == 2) {}

          m_result[i,,]  <- m_diff_array - m_result[i,,]
        }

      }


      return(m_result)

    }

    unscale_cumsum  <- function(input_array, result_list, start_event= 0, m_diff_array=NULL) {

      if(result_list$is_3d) {


        m_result <- array(NA_real_, dim=dim(input_array)) |> keep::as.karray()

        indexer <- dim(input_array)[3]

        for(idx in 1:indexer) {


            temp_result_list  <- list(
                'original_means' = result_list$original_means[idx,],
                'original_scale' = result_list$original_scale[idx,],
                'is_3d' = FALSE,
                "is_log" = result_list$use_log

            )

            m_result[,,idx,]  <- unscale_cumsum_from_2d_array(input_array = input_array[,,idx,],
                                                              result_list = temp_result_list,
                                                              start_event= start_event,
                                                              m_diff_array=m_diff_array[,idx,])

        }


        return(m_result)

      } else {

        result  <- unscale_cumsum_from_2d_array(input_array = input_array,
                                                result_list = result_list,
                                                start_event= start_event,
                                                m_diff_array=m_diff_array)

        return(result)

      }



    }

    # input_array is in the original scale,
    unscale_cumsum_from_2d_array   <-  function(input_array, result_list, start_event= 0, m_diff_array=NULL) {

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
