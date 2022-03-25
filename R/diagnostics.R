
DIAGNOSTICS <- modules::module({

    import("ggplot2")
    import("stats")
  
    filter_df  <- function(df, current_parameter_name) {

        df_result  <- df  |> 
                      dplyr::filter(parameter  |>    
                                    stringr::str_detect(current_parameter_name) )  |>  
                      dplyr::mutate(param_name = parameter  |>  stringr::str_extract("[a-zA-Z_0-9]*"))  |> 
                      dplyr::filter(  param_name == current_parameter_name   )  |> 
                      dplyr::select(-param_name)  

        return(df_result)

    }
    
    get_chain <- function(stan_fit, chain) {
      
      summary_list  <- stan_fit  |> build_summary_df() 
      
      # browser()
      
      if(chain > (length(summary_list) + 1) ) {
        stop("Chain of out index")
      }
      
      if(chain < 1) {
        m_data_df  <- summary_list$global
      } else {
        m_data_df  <-  summary_list[[paste0("C", chain)]]
      }
      
      if(nrow(m_data_df) < 1) {
        stop("The parameter dataframe is empty.")
      }
      
      return(m_data_df)
      
    }

  
    box_plot <- function(stan_fit, 
                         variable_name = c("Rhat", "n_eff", "mean",
                                            "sd", "2.5%", "25%", "50%",
                                            "75%", "97.5%"),
                         parameter_name=NA_character_, chain=0) {

      variable_name = match.arg(variable_name)

      m_valid_parameters  <- c(NA_character_,stan_fit  |>  names() |> stringr::str_extract("[a-zA-Z_0-9]*") |> unique())

      if( !(parameter_name %in% m_valid_parameters) ) {
          stop(
            paste0(
              "The only valid parameters names are: ",
              paste0(m_valid_parameters, collapse=","),
              "."
            )
          )
        }
      
      m_data_df <- get_chain(stan_fit, chain) |> 
                   dplyr::mutate( parameter = parameter |>  stringr::str_extract("[a-zA-Z_0-9]*") ) 
      
      
      if(!is.na(parameter_name)) {
        m_data_df <-  m_data_df |>  
                      dplyr::filter(parameter == parameter_name)
      }
      
                   
      parameter_name  <- ifelse(is.na(parameter_name), "", parameter_name)
      m_title  <- paste0(parameter_name, variable_name, collapse=" - ")
      
      
      # browser()
      
      m_plot  <- ggplot(data=m_data_df, mapping=aes(x=parameter, y=(!!sym(variable_name)) )) + 
                 stat_boxplot(geom = "errorbar", width = 0.2) + 
                 geom_boxplot() +
                 theme_bw() +
                 ggtitle(m_title)
      
      return(m_plot)
      
    }
    
    plot_ics  <- function(stan_fit, variable_name = c("Rhat", "n_eff", "mean",
                                                       "sd", "2.5%", "25%", "50%",
                                                       "75%", "97.5%"), parameter_name, chain=0) {

        # browser()
        # 
        
        m_valid_parameters  <- stan_fit  |>  names() |> stringr::str_extract("[a-zA-Z_0-9]*") |> unique()
        
        variable_name = match.arg(variable_name)

        if( !(parameter_name %in% m_valid_parameters) ) {
          stop(
            paste0(
              "The only valid parameter names are: ",
              paste0(m_valid_parameters, collapse=","),
              "."
            )
          )
        }

        m_data_df <- get_chain(stan_fit, chain)

        # browser()
        
        suppressWarnings({
          
          m_data_df_grouped  <- m_data_df |> 
            dplyr::mutate( parameter = parameter |>  stringr::str_extract("[a-zA-Z_0-9]*") ) |> 
            dplyr::filter(parameter == parameter_name) |> 
            dplyr::group_by(time_index, parameter) |> 
            dplyr::summarise_all( list(~ mean(., na.rm=TRUE), 
                                       ~ median(., na.rm = TRUE), 
                                       ~ min(., na.rm = TRUE),
                                       ~ max(., na.rm = TRUE)))
          
        })
        
        # browser()
        
        if( !(variable_name %in% names(m_data_df)) ) {
            stop(paste("The variable", 
                       variable_name, 
                       "is not in the dataframe, the avalible variaes are", 
                       stan_fit, chainnames(m_data_df_grouped)))
        }

        min_name  <- paste0(variable_name, "_min")
        max_name  <- paste0(variable_name, "_max")
        mean_name <- paste0(variable_name, "_mean")
        median_name <- paste0(variable_name, "_median")
        m_title  <- paste0(parameter_name, " - ", variable_name)

        pd <- position_dodge(0.1)

        
        # browser()
        
        m_plot  <-  ggplot(m_data_df_grouped, aes(x=time_index, y=!!sym(mean_name))) +
                    geom_errorbar(aes(ymin=!!sym(min_name), 
                                      ymax=!!sym(max_name)), colour="black", width=.1, position=pd) +
                    geom_line(position=pd, colour="blue") +
                    # geom_point(position=pd, size=3) +
                    theme_bw() +
                    ggtitle(m_title)
                   

        return(m_plot)

    }

    
    add_time_index  <- function(m_df) {

        # browser()
        # 
        
        m_df$comma_count <- m_df$parameter |> stringr::str_count(pattern = ",")
        
        count_0 <-  m_df |> 
                    dplyr::filter( (comma_count == 0) || (comma_count > 2)) |> 
                    dplyr::mutate(time_index = NA_integer_)
                    
        
        count_2 <- m_df |> 
                   dplyr::filter(comma_count == 2)
        count_2$time_index  <-  count_2$parameter |> 
                                stringr::str_extract( "(?<=\\[)[0-9]+,[0-9]+,[0-9]+[^\\]\\[]*(?=])")  |>                                 stringr::str_extract(",[0-9]+,")  |> 
                                stringr::str_remove_all(",")  |> 
                                as.integer()
        
        # browser()
        
        count_1 <- m_df |> 
                   dplyr::filter(comma_count == 1)
        count_1$time_index  <-  count_1$parameter |> 
                                stringr::str_extract( "(?<=\\[)[0-9]+,[0-9]+[^\\]\\[]*(?=])")  |> 
                                stringr::str_extract("[0-9]+,") |> 
                                stringr::str_remove_all(",")  |> 
                                as.integer()
                
        # m_df <- m_df |> dplyr::select(-comma_count)
        
        
        result_df = dplyr::bind_rows(
          count_0,
          count_1,
          count_2,
        ) |> dplyr::select(-comma_count)
        
        return(result_df)

    }
  
    prepare_df  <- function(m_data) {

     m_data  <- m_data  |> 
                data.frame()  |> 
                dplyr::add_rownames(var = "parameter")  |> 
                add_time_index()


     return(m_data)

    }

    build_summary_df  <- function(stan_fit) {

      # browser()

      m_summary  <- rstan::summary(stan_fit)

      result_list  <- list()

      result_list[["global"]]  <- m_summary$summary |> prepare_df()

      for(i in 1:(dim(m_summary$c_summar)[2]) ) {

        result_list[[paste0("C", i)]]  <- m_summary$summary_c[,,i] |> prepare_df()

      }

      return(result_list)

    }
  
})