MODULE_SUMMARY <- modules::module({
  
  import("dplyr")
  import("tidyr")
  import("stats")
  
  # PREDIFINED_CONFIDENCE  <- 0.5
  
  ics_to_data_frame <- function(lower_matrix, upper_matrix, event_min, 
                            event_max, variables_names, 
                            lower_limit_name="lower",
                            upper_limit_name="upper") {
    
    
    
    var_names <- variables_names
    current_time_index <- event_min:event_max
    
    m_df_lower <- lower_matrix
    m_df_upper <- upper_matrix
    
    
    colnames(m_df_lower) <- var_names
    colnames(m_df_upper) <- var_names

    # browser()
    
    m_df_lower <- m_df_lower |>
      tibble::as_tibble()  |>
      mutate(time_index=current_time_index) |> 
      pivot_longer(-time_index, names_to="variable", values_to = lower_limit_name)
    
    
    m_df_upper <- m_df_upper |>
      tibble::as_tibble()  |>
      mutate(time_index=current_time_index) |> 
      pivot_longer(-time_index, names_to="variable", values_to = upper_limit_name)
    
    
    result_df <- m_df_lower |> left_join(m_df_upper,  by=c("time_index", "variable"))
    
    
    return(result_df)
    
    
  }

  

  merge_multiuples_df  <- function(df_list, dates_df=NULL) {

    result_df  <- df_list[[1]]

    if(length(df_list) > 1) {

      for(i in 2:length(df_list)) {

        result_df  <- result_df |> left_join(df_list[[i]], by=c("time_index", "variable"))

      }

    }

    if(!is.null(dates_df)) {
      
      result_df <- result_df |> left_join(dates_df, by="time_index")
      
    }

    return(result_df)

  }

  get_averange_df <- function(variable_array, event_min, event_max, variables_names, prefix) {
    
    
    matrix_mean <- variable_array[,event_min:event_max,] |> 
                   apply(c(2,3), function(x){  mean(x, na.rm=T) })
    
    
    median_matrix <- variable_array[,event_min:event_max,] |> 
                     apply(c(2,3), function(x){  median(x, na.rm=T)  })
    
    averange_df <- ics_to_data_frame(lower_matrix = matrix_mean, 
                                     upper_matrix = median_matrix, 
                                     event_min = event_min, 
                                     event_max = event_max, 
                                     variables_names = variables_names, 
                                     lower_limit_name = paste0(prefix,"mean"),
                                     upper_limit_name = paste0(prefix,"median"))
    
    return(averange_df)
    
  }

  make_df_from_simul_resut  <- function(simul_result, event_min, event_max, dates_df, variable_name="global", m_quantile=0.9) {

    names(simul_result)  <- names(simul_result)  |> 
                            stringr::str_remove_all("_inter")  |> 
                            stringr::str_remove_all("_aggregate")  


    event_min_date  <- event_min |> get_date_from_index(dates_df)
    event_max_date  <- event_max |> get_date_from_index(dates_df)
    
    
    quantile_cumsum_name  <- paste0("quantile_",m_quantile,"_cumsum")
    quantile_arco_name  <- paste0("quantile_",m_quantile,"_arco")
      
    event_max  <- min(event_max, dim(simul_result[[1]])[1])

    time_index <- event_min:event_max

    result_df  <- data.frame(time_index = time_index)

    result_df$event_min  <- event_min
    result_df$event_max  <- event_max

    result_df$event_min_date  <- event_min_date
    result_df$event_max_date  <- event_max_date

    result_df$cumsum_lower  <- simul_result$lower
    result_df$cumsum_upper  <- simul_result$upper

    result_df$cumsum_mean  <- simul_result$mean
    result_df$cumsum_median  <- simul_result$median 

    result_df$arco_mean  <- (1/time_index)*simul_result$mean
    result_df$arco_median  <- (1/time_index)*simul_result$median 

    result_df$lower_arco  <- (1/time_index)*simul_result$lower
    result_df$upper_arco  <- (1/time_index)*simul_result$upper


    result_df$lower_arco  <- (1/time_index)*simul_result$lower
    result_df$upper_arco  <- (1/time_index)*simul_result$upper


    result_df[[quantile_cumsum_name]]  <- simul_result$quantile
    result_df[[quantile_arco_name]]  <- (1/time_index)*simul_result$quantile


    return(result_df)


  }
  
  extract_single_from_simul  <- function(simul_result, event_min, event_max, idx1=NULL, idx2=NULL) {

    
    

    if(is.null(idx1)) {

      if(is.null(idx2)) {
          simul_result  <- simul_result$aggregate
      }

      simul_result$cumsum_result  <-  NULL
      simul_result$cumsum_aggregate  <-  NULL
      simul_result$t  <-  NULL
        
    } else if(is.null(idx2)) {
        simul_result$aggregate  <- NULL

        simul_result$cumsum_result  <-  NULL
        simul_result$cumsum_aggregate  <-  NULL
        simul_result$t  <-  NULL
        
        for(m_elem in names(simul_result)) {
          simul_result[[simul_result]]  <- simul_result[[simul_result]][,idx1]
        }

    } else {

      simul_result$aggregate  <- NULL

      simul_result$cumsum_result  <-  NULL
      simul_result$cumsum_aggregate  <-  NULL
      simul_result$t  <-  NULL

      for(m_elem in names(simul_result)) {
          simul_result[[simul_result]]  <- simul_result[[simul_result]][,idx1, idx2]
      }
      
    }



    return(simul_result)

  }

  get_impact_manual  <- function(m_model, event_min, event_max, 
                              variables_names, 
                              variable_index = NULL,
                              dates_df = NULL, ci=0.9, discount=NULL) {

    N <- dim(m_model$X_data)[1]
    
    if( (event_min < 1) |  (event_max < 1) ) {
      stop("All the events times must be greater than zero")
    }
    
    if(event_max < 0) {
      stop("All the initial evnets times must be greater than zero")
    }
    
    if(event_min >= event_max) {
      stop("the initial time must be greater than the end time")
    }


    m_model$.fit(event_initial=event_min, discount=discount)

    
    simul_result  <- m_model$.__enclos_env__$private$.simul_sumation




  }
 
  # TODO mejorar esto, mucho codigo repetido
  get_impact_stan_matrix  <- function(m_model, event_min, event_max, 
                              variables_names, dates_df = NULL, ci=0.9) {


    # browser()

    N <- dim(m_model$X_data)[1]
    
    if( (event_min < 1) |  (event_max < 1) ) {
      stop("All the events times must be greater than zero")
    }
    
    if(event_max < 0) {
      stop("All the initial events times must be greater than zero")
    }
    
    if(event_min >= event_max) {
      stop("the initial time must be greater than the end time")
    }
    
    m_model$predict(event_initial=event_min)

    
    
   i_lower <- m_model$.__enclos_env__$private$.extracted_data$cumsum_only_after[,event_min:event_max,,] |> 
              apply(c(2,3,4), 
                     function(x){  bayestestR::hdi(x, c=ci)$CI_low  }
              )
    UTILS$gc_quiet()
   
   i_upper <- m_model$.__enclos_env__$private$.extracted_data$cumsum_only_after[,event_min:event_max,,] |> 
              apply(c(2,3,4), 
                    function(x){  
                      bayestestR::hdi(x, c=ci)$CI_high  
              })
    UTILS$gc_quiet()
   
   
   i_lower_arco <- m_model$.__enclos_env__$private$.extracted_data$arco_only_after[,event_min:event_max,,] |> 
                   apply(c(2,3,4), 
                         function(x){  bayestestR::hdi(x, c=ci)$CI_low  }
                   )
    UTILS$gc_quiet()
   
   i_upper_arco <- m_model$.__enclos_env__$private$.extracted_data$arco_only_after[,event_min:event_max,,] |> 
                   apply(c(2,3,4), 
                         function(x){  
                           bayestestR::hdi(x, c=ci)$CI_high  
                         })
    UTILS$gc_quiet()


    i_arco_quantile <- m_model$.__enclos_env__$private$.extracted_data$arco_only_after[,event_min:event_max,,] |> 
                             apply(c(2,3,4), 
                                  function(x){  
                                   x |> quantile(ci, na.rm=T) |> as.numeric() 
                                  })
    UTILS$gc_quiet()

    i_cumsum_quantile <- m_model$.__enclos_env__$private$.extracted_data$cumsum_only_after[,event_min:event_max,,] |> 
                               apply(c(2,3,4), 
                                    function(x){  
                                     x |> quantile(ci, na.rm=T) |> as.numeric() 
                                    })
    UTILS$gc_quiet()

   

   m_df_result_list  <- list()

  #indexer = dim(m_model$.__enclos_env__$private$.extracted_data$cumsum_only_after)[3]
  indexer = length(m_model$vector_name)

   for(idx in 1:indexer) {

     i_quantile_arco_cumsum <- ics_to_data_frame(lower_matrix = i_arco_quantile[,idx,], 
                                                upper_matrix = i_cumsum_quantile[,idx,], 
                                                event_min = event_min, 
                                                event_max = event_max, 
                                                variables_names = variables_names, 
                                                lower_limit_name=paste0("quantile_", ci, "_arco"),
                                                upper_limit_name=paste0("quantile_", ci, "_cumsum")
                                                )
     
     

     data_frame_result_cumsum <- ics_to_data_frame(lower_matrix = i_lower[,idx,], 
                                                 upper_matrix = i_upper[,idx,], 
                                                 event_min = event_min, 
                                                 event_max = event_max, 
                                                 variables_names = variables_names, 
                                                 lower_limit_name="cumsum_lower",
                                                 upper_limit_name="cumsum_upper")
   
   
    data_frame_result_arco <- ics_to_data_frame(lower_matrix = i_lower_arco[,idx,], 
                                                upper_matrix = i_upper_arco[,idx,], 
                                                event_min = event_min, 
                                                event_max = event_max, 
                                                variables_names = variables_names, 
                                                lower_limit_name="lower_arco",
                                                upper_limit_name="upper_arco")
 
   
    averange_df <-  get_averange_df(variable_array=m_model$.__enclos_env__$private$.extracted_data$cumsum_only_after[,,idx,], 
                                    event_min=event_min, 
                                    event_max=event_max, 
                                    variables_names=variables_names,
                                    prefix ="cumsum_")
   
   
    averange_arco_df <-  get_averange_df(variable_array=m_model$.__enclos_env__$private$.extracted_data$arco_only_after[,,idx,], 
                                          event_min=event_min, 
                                          event_max=event_max, 
                                          variables_names=variables_names,
                                          prefix = "arco_")

    temp_result_df  <- merge_multiuples_df(list(
      data_frame_result_cumsum,
      data_frame_result_arco,
      averange_df,
      averange_arco_df,
      i_quantile_arco_cumsum
    ), dates_df=dates_df)

    temp_result_df$vector  <- m_model$vector_name[idx]


    m_df_result_list[[length(m_df_result_list) + 1]]  <- temp_result_df


   }

    result_df  <- m_df_result_list  |> bind_rows()

    return(result_df)

  }


  
  get_impact_stan <- function(m_model, event_min, event_max, 
                              variables_names, dates_df = NULL, ci=0.9) {
    
    N <- dim(m_model$X_data)[1]

    variables_names  <- c(variables_names, "global")
    
    if( (event_min < 1) |  (event_max < 1) ) {
      stop("All the events times must be greater than zero")
    }
    
    if(event_max < 0) {
      stop("All the initial events times must be greater than zero")
    }
    
    if(event_min >= event_max) {
      stop("the initial time must be greater than the end time")
    }
    
    m_model$predict(event_initial=event_min)

    
    
   i_lower <- m_model$.__enclos_env__$private$.extracted_data$cumsum_only_after[,event_min:event_max,] |> 
              apply(c(2,3), 
                     function(x){  bayestestR::hdi(x, c=ci)$CI_low  }
              )
    UTILS$gc_quiet()

    i_lower_global <- m_model$.__enclos_env__$private$.extracted_data$cumsum_only_after[,event_min:event_max,] |> 
              apply(c(1, 2), function(x) {  mean(x, na.rm=T) } )  |> 
              apply(c(2), 
                     function(x){  bayestestR::hdi(x, c=ci)$CI_low  }
              )  |> matrix(ncol=1)

    UTILS$gc_quiet()
   
   i_upper <- m_model$.__enclos_env__$private$.extracted_data$cumsum_only_after[,event_min:event_max,] |> 
              apply(c(2,3), 
                    function(x){  
                      bayestestR::hdi(x, c=ci)$CI_high  
              })

    UTILS$gc_quiet()

    i_upper_global <- m_model$.__enclos_env__$private$.extracted_data$cumsum_only_after[,event_min:event_max,] |> 
              apply(c(1, 2), function(x) {  mean(x, na.rm=T) } )  |> 
              apply(c(2), 
                    function(x){  
                      bayestestR::hdi(x, c=ci)$CI_high  
              })  |>  matrix(ncol=1)
   
   
   i_lower_arco <- m_model$.__enclos_env__$private$.extracted_data$arco_only_after[,event_min:event_max,] |> 
                   apply(c(2,3), 
                         function(x){  bayestestR::hdi(x, c=ci)$CI_low  }
                   )
    UTILS$gc_quiet()

    i_lower_arco_global <- m_model$.__enclos_env__$private$.extracted_data$arco_only_after[,event_min:event_max,] |> 
                   apply(c(1, 2), function(x) {  mean(x, na.rm=T) } )  |> 
                   apply(c(2), 
                         function(x){  bayestestR::hdi(x, c=ci)$CI_low  }
                   )  |>  matrix(ncol=1)
    UTILS$gc_quiet()
   
   i_upper_arco <- m_model$.__enclos_env__$private$.extracted_data$arco_only_after[,event_min:event_max,] |> 
                   apply(c(2,3), 
                         function(x){  
                           bayestestR::hdi(x, c=ci)$CI_high  
                         })
    UTILS$gc_quiet()

    i_upper_arco_global <- m_model$.__enclos_env__$private$.extracted_data$arco_only_after[,event_min:event_max,] |> 
                   apply(c(1, 2), function(x) {  mean(x, na.rm=T) } )  |> 
                   apply(c(2), 
                         function(x){  
                           bayestestR::hdi(x, c=ci)$CI_high  
                         }) |>  matrix(ncol=1)
    UTILS$gc_quiet()


    i_arco_quantile <- m_model$.__enclos_env__$private$.extracted_data$arco_only_after[,event_min:event_max,] |> 
                             apply(c(2,3), 
                                  function(x){  
                                   x |> quantile(ci, na.rm=T) |> as.numeric() 
                                  })
    UTILS$gc_quiet()

    i_arco_quantile_global <- m_model$.__enclos_env__$private$.extracted_data$arco_only_after[,event_min:event_max,] |> 
                             apply(c(1, 2), function(x) {  mean(x, na.rm=T) } )  |> 
                             apply(c(2), 
                                  function(x){  
                                   x |> quantile(ci, na.rm=T) |> as.numeric() 
                                  }) |>  matrix(ncol=1)
    UTILS$gc_quiet()

    i_cumsum_quantile <- m_model$.__enclos_env__$private$.extracted_data$cumsum_only_after[,event_min:event_max,] |> 
                               apply(c(2,3), 
                                    function(x){  
                                     x |> quantile(ci, na.rm=T) |> as.numeric() 
                                    })
    UTILS$gc_quiet()

    i_cumsum_quantile_global <- m_model$.__enclos_env__$private$.extracted_data$cumsum_only_after[,event_min:event_max,] |> 
                               apply(c(1, 2), function(x) {  mean(x, na.rm=T) } )  |> 
                               apply(c(2), 
                                    function(x){  
                                     x |> quantile(ci, na.rm=T) |> as.numeric() 
                                    }) |>  matrix(ncol=1)
    UTILS$gc_quiet()


    i_quantile_arco_cumsum <- ics_to_data_frame(lower_matrix = cbind(i_arco_quantile, i_arco_quantile_global), 
                                                upper_matrix = cbind(i_cumsum_quantile, i_cumsum_quantile_global), 
                                                event_min = event_min, 
                                                event_max = event_max, 
                                                variables_names = variables_names, 
                                                lower_limit_name=paste0("quantile_", ci, "_arco"),
                                                upper_limit_name=paste0("quantile_", ci, "_cumsum")
                                                )
    
    UTILS$gc_quiet()
   
    # browser()
   data_frame_result_cumsum <- ics_to_data_frame(lower_matrix = cbind(i_lower, i_lower_global), 
                                                 upper_matrix = cbind(i_upper, i_upper_global), 
                                                 event_min = event_min, 
                                                 event_max = event_max, 
                                                 variables_names = variables_names, 
                                                 lower_limit_name="cumsum_lower",
                                                 upper_limit_name="cumsum_upper")
   
   
   data_frame_result_arco <- ics_to_data_frame(lower_matrix = cbind(i_lower_arco, i_lower_arco_global), 
                                               upper_matrix = cbind(i_upper_arco, i_upper_arco_global), 
                                               event_min = event_min, 
                                               event_max = event_max, 
                                               variables_names = variables_names, 
                                               lower_limit_name="lower_arco",
                                               upper_limit_name="upper_arco")
 
   temp_dim  <- m_model$.__enclos_env__$private$.extracted_data$cumsum_only_after  |>  dim()
   temp_dim[3]  <- 1

   averange_df <-  get_averange_df(variable_array=
                                   abind::abind(
                                    m_model$.__enclos_env__$private$.extracted_data$cumsum_only_after,
                                    m_model$.__enclos_env__$private$.extracted_data$cumsum_only_after  |> 
                                    apply(c(1, 2), function(x) {  mean(x, na.rm=T) } )  |> 
                                    array(dim = temp_dim),
                                    along=length(temp_dim)
                                   ), 
                                   event_min=event_min, 
                                   event_max=event_max, 
                                   variables_names=variables_names,
                                   prefix ="cumsum_")
   
   
   averange_arco_df <-  get_averange_df(variable_array=
                                        abind::abind(
                                          m_model$.__enclos_env__$private$.extracted_data$arco_only_after, 
                                          m_model$.__enclos_env__$private$.extracted_data$arco_only_after  |> 
                                          apply(c(1, 2), function(x) {  mean(x, na.rm=T) } )  |> 
                                          array(dim = temp_dim),
                                          along=length(temp_dim)
                                        ),
                                        event_min=event_min, 
                                        event_max=event_max, 
                                        variables_names=variables_names,
                                        prefix = "arco_")

    result_df  <- merge_multiuples_df(list(
      data_frame_result_cumsum,
      data_frame_result_arco,
      averange_df,
      averange_arco_df,
      i_quantile_arco_cumsum
    ), dates_df=dates_df)

    return(result_df)
    
    
  }
  
  get_date_from_index <- function(m_index, dates_df) {
    
    # browser()

    if(is.null(dates_df)) {
      return(m_index)
    }

    result  <- dates_df |> dplyr::filter(time_index == m_index) |> pull(Date) |> first()  |> lubridate::as_date()

    return(result)
    
  }
  
  # get_date_from_index <- Vectorize(get_date_from_index,  vectorize.args = c("m_index"), SIMPLIFY = TRUE,  USE.NAMES = FALSE)
  
  
  get_index_from_date <- function(m_date, dates_df) {

    if(is.null(dates_df)) {
      return(m_date)
    }
    
    result  <- dates_df |> dplyr::filter(Date == m_date) |> pull(time_index) |> first()

    return(result)
    
  }
  
  # get_index_from_date <- Vectorize(get_index_from_date, vectorize.args = c("m_date"), SIMPLIFY = TRUE,USE.NAMES = FALSE)
  

  print_single_summary  <- function(result) {

    cat("\n\n")

    if(lubridate::is.Date(result$event_min_date) &
       lubridate::is.Date(result$event_max_date)) {

      cat(
        paste0(
          "Date: From ", result$event_min_date, " To ", result$event_max_date, "\n"
        )
      )

    }

    cat(
      paste0(
        "Index: From ", result$event_min, " To ", result$event_max, "\n\n"
      )
    )



    #arco_result  <- 
    temp_result  <- 
    result$impact_df |> 
                    dplyr::filter(time_index == result$event_max)  |> 
                    select(
                      # time_index, 
                      variable, 
                      lower_arco, upper_arco, arco_median, #arco_mean,
                      cumsum_lower, cumsum_upper, cumsum_median #, cumsum_mean
                    
                    )  |> 
                    mutate_if(is.numeric, function(x){round(x,2)})  
    
    stargazer::stargazer(temp_result, type = 'text', summary = FALSE)

    #stargazer::stargazer(arco_result,  type = 'text', summary = FALSE)
    #print(arco_result)

  }

  result_summary  <- function(result_list) {


    for(elem_name in names(result_list)) {

      current_data  <- result_list[[elem_name]]

      print_single_summary(current_data)

      #print(current_data)

    }

    #return(result_list)


  }
  
  get_get_multiple_impacts_stan <- function(model, impact_list, dates_df=NULL, ci=0.9, is_matrix_model=FALSE) {


    UTILS$gc_quiet()

    if(is.null(dates_df)) {
       dates_df  <- model$get_dates_df()
    }
    

    result_list  <- list()

    for(index in 1:length(impact_list)) {

      m_text  <- paste0(
        "Computing impact: ",
        index,"/", length(impact_list)
      )

      print(m_text)

      UTILS$gc_quiet()

      #current_index  <- length(result_list) + 1

      # browser()
      if(lubridate::is.Date(impact_list[[index]][1])) {
        
        current_event_min_date <- impact_list[[index]][1] |> lubridate::as_date()
        current_event_min <- current_event_min_date |> get_index_from_date(dates_df)
          
      } else {
        
        current_event_min <- impact_list[[index]][1]
        current_event_min_date <- current_event_min |> get_date_from_index(dates_df)
  
      }
      
      if(lubridate::is.Date(impact_list[[index]][2])) {
        
        current_event_max_date <- impact_list[[index]][2] |> lubridate::as_date()
        current_event_max <- current_event_max_date |> get_index_from_date(dates_df)
        
      } else {
        
        current_event_max <- impact_list[[index]][2]
        current_event_max_date <- current_event_max |> get_date_from_index(dates_df)
        
      }
      
      # browser()
      elem_name  <- paste0("[",
                           current_event_min_date |> as.character(),
                           ",",
                           current_event_max_date |> as.character(),
                           "]")


      if(is_matrix_model) {

        impact_df  <- get_impact_stan_matrix(
            m_model=model, 
            event_min=current_event_min, 
            event_max=current_event_max, 
            variables_names=model$variables_names, 
            dates_df = dates_df, 
            ci=ci
        )



      } else {

        impact_df  <- get_impact_stan(
            m_model=model, 
            event_min=current_event_min, 
            event_max=current_event_max, 
            variables_names=model$variables_names, 
            dates_df = dates_df, 
            ci=ci
        )

      }

      

      UTILS$gc_quiet()

      result_list[[elem_name]]  <- list(
        "impact_df"  = impact_df,

        "event_min_date" = current_event_min_date,
        "event_max_date" = current_event_max_date,

        "event_max" = current_event_max,
        "event_min" = current_event_min

      )



    }

    rm(model)
    UTILS$gc_quiet()

    return(result_list)


    
    
  }
  
  
  
})
  
  
  

