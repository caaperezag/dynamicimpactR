MODULE_PLOT_EXTRA <- modules::module({
  
 
  import("stats")
  import("dplyr")
  import("tidyr")
  import("ggplot2")
  import("patchwork")

  MIN_DATE_DIST <- 5
  
  get_break_dates  <- function(max_date, event_initial, min_date=1, length.out=5) {
    
    m_break_dates  <- seq(min_date,max_date, length.out=length.out) |> floor()
    
    
    
    remove_index <- abs(m_break_dates - event_initial) > MIN_DATE_DIST
    remove_index <-  (remove_index)
    
    m_break_dates <- m_break_dates[remove_index]
    
    m_break_dates  <- c(event_initial, m_break_dates)
    
    return(m_break_dates)
    
    
  }
  
  plot_model <- function(m_model, event_initial, variable_name) {
    
    m_model$predict(event_initial)
    
    if(variable_name == "global") {
      base_df <- m_model$.__enclos_env__$private$.plot_df_aggregate
    } else {
      base_df <- m_model$.__enclos_env__$private$.plot_df |> filter(variable == variable_name)
    }
    
    # browser()
    arco_df <- base_df |> 
      filter(type == "cumsum" ) |> 
      mutate(type = "ArCo", class = "ArCo",
             across( c(value, lower_limit , upper_limit), ~ .x/(time_index-event_initial) )  )
    
    
    full_df <- bind_rows(base_df, arco_df)
    
    
    current_plot <-  full_df |> plot_df(event_initial = event_initial,
                                        dates_df=NULL)
    
    
    return(current_plot)
  }
  
  
  plot_df  <- function(plot_df, event_initial, dates_df=NULL, break_dates = NULL) {
    
    # browser()
    # 
    
    
    if(is.null(dates_df)) {
      dates_df <- data.frame(time_index = plot_df$time_index |> unique() |> sort())
      dates_df$Date <- dates_df$time_index
    }
    
    temp_df <-  plot_df |> 
      mutate(
        class   = if_else(type == "cumsum", "cumsum", class)
      ) |> 
      rename(lower = lower_limit, upper = upper_limit) |> 
      mutate(lower = if_else( (time_index <= event_initial) , NA_real_, lower)) |> 
      mutate(upper = if_else(  (time_index <= event_initial) , NA_real_, upper)) |> 
      mutate(value = if_else(  (time_index <= event_initial) &
                                 (!(class %in% c("input", "prediction", "real","error"))) , NA_real_, value)) |> 
      mutate(lower = if_else(  class == "prediction" , NA_real_, lower)) |> 
      mutate(upper = if_else(  class == "prediction" , NA_real_, upper)) |> 
      left_join(dates_df)
    # browser()
    
    temp_df <-  temp_df |> bind_rows(
      
      temp_df |> dplyr::filter(class == "real") |> mutate(class = "Y ")
    )
    
    
    # browser()
    # tibble::view(temp_df)
    
    plot_df <- temp_df
    
    max_date  <- plot_df$time_index  |> max()
    
    
    if(is.null(break_dates)) {
      # break_dates = waiver()
      
      break_dates  <- get_break_dates(max_date=max_date, 
                                      event_initial=event_initial) 
      
      # break_dates_index <- break_dates
      
      break_dates <- dates_df |> dplyr::filter(time_index %in% break_dates) |> pull(Date)
      
    }
    
    data_event_initial <- dates_df |> dplyr::filter(time_index == event_initial) |> pull(Date) |> first()
    
    # browser()
    
    plot_df  <- plot_df  |> 
      rename(param = class) |> 
      mutate(
        param = param   |> recode(real = "Y", 
                                  error = "Impact",
                                  # prediction = "Contra factual",
                                  prediction = "Predicción",
                                  cumsum = "cumulative impact",
                                  ArCo = "ArCo",
                                  input = "X"
        )
      )   |>  
      mutate( face_break = 
                case_when(
                  param %in% c("Y","Predicción") ~ "Real series and prediction",
                  param %in% c("Impact") ~ "Impact",
                  param %in% c("cumulative impact") ~ "cumulative impact",
                  param %in% c("ArCo") ~ "ArCo",
                  #param %in% c("X", "Y") ~ "Input"
                  param %in% c("X", "Y ") ~ "Input"
                  
                ) |> factor(levels = c("Real series and prediction", "Input", "Impact", "cumulative impact", "ArCo"))
              
              
      ) |> 
      mutate(
        
        face_break2 = case_when(
          face_break %in% c("Real series and prediction", "Input")  ~ "Modelo",
          face_break %in% c("Impact", "cumulative impact", "ArCo") ~ "Efecto"
        ) |> factor(levels=c("Modelo", "Efecto"))
        
      )
    
    m_scale_break <- scale_x_continuous(breaks = break_dates )
    
    if(lubridate::is.Date(break_dates)) {
      
      m_scale_break <- scale_x_date(labels = scales::date_format("%m/%d"), breaks = break_dates ) 
      
    }
    
    plot_result1  <- plot_df  |> 
      dplyr::filter(face_break2 == "Modelo")  |> 
      ggplot(aes(x=Date, y=value, col=param)) +
      geom_ribbon( aes(ymin = lower, ymax = upper), fill = "grey90" ) +
      geom_line() +
      geom_vline(xintercept = data_event_initial, linetype = "dashed") +
      geom_hline(yintercept = 0, linetype = "dashed") +
      facet_wrap(~face_break, ncol=2, scales='fixed') +
      #scale_colour_manual(values=c("#006666", "#ff7f00", "#ff7f00", "#e41a1c")) +
      #scale_color_discrete(breaks = c("X", "Y", "Predicción"), values=c("#006666", "#ff7f00", "#ff7f00",  "#e41a1c") )+ 
      # "#3b80b9"
      # scale_colour_manual(values=c("#006666", "#ff7f00",  "#3b80b9",  "#ff7f00") )+ 
      scale_colour_manual(breaks = c("X", "Y", "Predicción"), values=c("#006666", "#ff7f00",  "#3b80b9",  "#ff7f00"), na.value="#ff7f00" )+ 
      #scale_colour_manual(values=c("#006666", "#ff7f00", "#ff7f00",  "#e41a1c")) +
      # scale_x_continuous(breaks = break_dates ) +
      m_scale_break +
      theme_bw()  +
      theme(
        legend.title = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text = element_text(size=14),
        legend.position="bottom",
        strip.text = element_text(size=11)#,
        # strip.text.x = element_blank()
      ) + ylab("") + xlab("")
    
    
    plot_result2  <- plot_df  |> 
      dplyr::filter(face_break2 == "Efecto",
                    Date >= data_event_initial
      )  |> 
      ggplot(aes(x=Date, y=value, col=param)) +
      geom_ribbon( aes(ymin = lower, ymax = upper), fill = "grey90" ) +
      geom_line() +
      geom_vline(xintercept = data_event_initial, linetype = "dashed") +
      geom_hline(yintercept = 0, linetype = "dashed") +
      facet_wrap(~face_break, ncol=3, scales='free_y') +
      scale_colour_manual(values=c("#e41a1c", "#4daf4a", "#9d57a7", "#3b80b9")) +
      # scale_x_date(labels = scales::date_format("%m/%d"), breaks = break_dates ) +
      m_scale_break +
      theme_bw()  +
      theme(
        legend.title = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text = element_text(size=14),
        legend.position="bottom",
        strip.text = element_text(size=11)#,
        # strip.text.x = element_blank()
      ) + ylab("")  + xlab("Fecha")
    
    # plot_result  <- (plot_result1/plot_result2) + plot_layout(guides = "collect") & theme(legend.position='bottom')
    plot_result  <- patchwork:::`/.ggplot`(plot_result1,plot_result2) + plot_layout(guides = "collect") & theme(legend.position='bottom')
    
    
    return(plot_result)
    
  }
  

  
  
  
})