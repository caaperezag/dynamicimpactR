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
      tibble::as.tibble()  |>
      mutate(time_index=current_time_index) |> 
      pivot_longer(-time_index, names_to="variable", values_to = lower_limit_name)
    
    
    m_df_upper <- m_df_upper |>
      tibble::as.tibble()  |>
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

    if(is.null(dates_df)) {
      
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
  
 
  
  get_impact_stan <- function(m_model, event_min, event_max, 
                              variables_names, dates_df = NULL, ci=0.95) {
    
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
    
    m_model$predict(event_initial=event_min)

    
    
   i_lower <- m_model$.__enclos_env__$private$.extracted_data$cumsum_only_after[,event_min:event_max,] |> 
              apply(c(2,3), 
                     function(x){  bayestestR::hdi(x, c=ci)$CI_low  }
              )
    UTILS$gc_quiet()
   
   i_upper <- m_model$.__enclos_env__$private$.extracted_data$cumsum_only_after[,event_min:event_max,] |> 
              apply(c(2,3), 
                    function(x){  
                      bayestestR::hdi(x, c=ci)$CI_high  
              })
    UTILS$gc_quiet()
   
   
   i_lower_arco <- m_model$.__enclos_env__$private$.extracted_data$arco_only_after[,event_min:event_max,] |> 
                   apply(c(2,3), 
                         function(x){  bayestestR::hdi(x, c=ci)$CI_low  }
                   )
    UTILS$gc_quiet()
   
   i_upper_arco <- m_model$.__enclos_env__$private$.extracted_data$arco_only_after[,event_min:event_max,] |> 
                   apply(c(2,3), 
                         function(x){  
                           bayestestR::hdi(x, c=ci)$CI_high  
                         })
    UTILS$gc_quiet()


    i_arco_quantile <- m_model$.__enclos_env__$private$.extracted_data$arco_only_after[,event_min:event_max,] |> 
                             apply(c(2,3), 
                                  function(x){  
                                   x |> quantile(ci) |> as.numeric() 
                                  })
    UTILS$gc_quiet()

    i_cumsum_quantile <- m_model$.__enclos_env__$private$.extracted_data$cumsum_only_after[,event_min:event_max,] |> 
                               apply(c(2,3), 
                                    function(x){  
                                     x |> quantile(ci) |> as.numeric() 
                                    })
    UTILS$gc_quiet()


    i_quantile_arco_cumsum <- ics_to_data_frame(lower_matrix = i_arco_quantile, 
                                                upper_matrix = i_cumsum_quantile, 
                                                event_min = event_min, 
                                                event_max = event_max, 
                                                variables_names = variables_names, 
                                                lower_limit_name=paste0("quantile_", ci, "_arco"),
                                                upper_limit_name=paste0("quantile_", ci, "_cumsum")
                                                )
    UTILS$gc_quiet()
   
   
    # browser()
   data_frame_result_cumsum <- ics_to_data_frame(lower_matrix = i_lower, 
                                                 upper_matrix = i_upper, 
                                                 event_min = event_min, 
                                                 event_max = event_max, 
                                                 variables_names = variables_names, 
                                                 lower_limit_name="cumsum_lower",
                                                 upper_limit_name="cumsum_upper")
   
   
   data_frame_result_arco <- ics_to_data_frame(lower_matrix = i_lower_arco, 
                                               upper_matrix = i_upper_arco, 
                                               event_min = event_min, 
                                               event_max = event_max, 
                                               variables_names = variables_names, 
                                               lower_limit_name="lower_arco",
                                               upper_limit_name="upper_arco")
 
   
   averange_df <-  get_averange_df(variable_array=m_model$.__enclos_env__$private$.extracted_data$cumsum_only_after, 
                                   event_min=event_min, 
                                   event_max=event_max, 
                                   variables_names=variables_names,
                                   prefix ="cumsum_")
   
   
   averange_arco_df <-  get_averange_df(variable_array=m_model$.__enclos_env__$private$.extracted_data$arco_only_after, 
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
    
    result$impact_df |> 
                    dplyr::filter(time_index == result$event_max)  |> 
                    select(
                      # time_index, 
                      variable, 
                      lower_arco, upper_arco, arco_median, #arco_mean,
                      cumsum_lower, cumsum_upper, cumsum_median #, cumsum_mean
                    
                    )  |> 
                    mutate_if(is.numeric, function(x){round(x,2)})  |> 
                    stargazer::stargazer( type = 'text', summary = FALSE)

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
  
  get_get_multiple_impacts_stan <- function(model, impact_list, dates_df=NULL, ci=0.9) {


    UTILS$gc_quiet()

    if(is.null(dates_df)) {
       dates_df  <- model$get_dates_df()
    }
    

    result_list  <- list()

    for(index in 1:length(impact_list)) {

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

      impact_df  <- get_impact_stan(
            m_model=model, 
            event_min=current_event_min, 
            event_max=current_event_max, 
            variables_names=model$variables_names, 
            dates_df = dates_df, 
            ci=ci
      )

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
  
  
  

