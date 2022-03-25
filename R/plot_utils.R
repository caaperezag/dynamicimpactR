PLOT_UTILS <- modules::module({
  
  import('ggplot2')

  date_1 <- lubridate::ymd(paste("2020", 1, 1, sep= ' '))
  date_2 <- lubridate::ymd(paste("2020", 5, 1, sep= ' '))
  
  dates_df <- data.frame(
    time_index = 1:122,
    Fecha = (date_1:date_2) |> lubridate::as_date()
  ) 
  
  BREAK_DATES  <- c(
    lubridate::ymd(paste("2020", 1, 1, sep= ' ')),
    lubridate::ymd(paste("2020", 2, 1, sep= ' ')),
    lubridate::ymd(paste("2020", 3, 1, sep= ' ')),
    
    lubridate::ymd(paste("2020", 3, 23, sep= ' ')),

    #lubridate::ymd(paste("2020", 4, 1, sep= ' ')),
    lubridate::ymd(paste("2020", 5, 1, sep= ' '))
    )

  plot_individual <- function(plot_df, plot_variable, event_initial, vector_name, exclude_input=F) {
    
    temp_unique_vars  <- plot_df$variable  |> unique()

    if( !(plot_variable %in% temp_unique_vars) ) {

        stop(paste0(
          "The plot variable must be one of: ",
          paste0(temp_unique_vars, collapse,  ",")
        ))
    }

    m_mean <- plot_df |> 
              dplyr::filter(type == "prediction", 
                            class == 'real',
                            !is_impact) |> 
              dplyr::pull(value) |> 
              mean(na.rm=T)
    
    plot_df_line_df <- data.frame(
      type = c('Impacto Acumulado', 'Ajuste', 'Impacto') |> 
             factor(levels=c("Ajuste", 'Impacto', 'Impacto Acumulado')) ,
      value = c(0, m_mean, 0)
    )

    if(!exclude_input) {

       # browser()

        m_df_add_1  <- plot_df  |> 
                       dplyr::filter(class == 'input') |> 
                       dplyr::mutate(type="Entrada")

        
        m_df_add_2  <- plot_df  |> 
                       dplyr::filter(class == 'real')   |> 
                       dplyr::mutate(type="Entrada")

        # browser()

        plot_df  <- dplyr::bind_rows(plot_df, m_df_add_1, m_df_add_2)

    } else {

      plot_df <- plot_df   |> 
                           dplyr::filter(class != 'input') 

    }

    
    # browser()
    plot_df_individual <- plot_df |> 
       dplyr::left_join(dates_df, by='time_index')  |> 
      dplyr::mutate(lower_limit = dplyr::if_else(class != 'real', lower_limit, NA_real_),
                    upper_limit = dplyr::if_else(class != 'real', upper_limit, NA_real_)) |> 
      dplyr::mutate(lower_limit = dplyr::if_else( (class == 'error') & (type == 'error'), NA_real_, lower_limit),
                    upper_limit = dplyr::if_else( (class == 'error') & (type == 'error'), NA_real_, upper_limit)) |> 
      dplyr::mutate(lower_limit = dplyr::if_else(type != 'prediction', lower_limit, NA_real_),
                    upper_limit = dplyr::if_else(type != 'prediction', upper_limit, NA_real_))  |> 
      dplyr::mutate(value = dplyr::if_else( (!is_impact) & (type == 'cumsum'), NA_real_, value))  |> 
      dplyr::mutate(value = dplyr::if_else(  (type == 'prediction') & (class == 'input') , NA_real_, value))  |>
      dplyr::mutate(type = type |> 
                      as.character() |> 
                      dplyr::recode('prediction'='Ajuste',
                                    'cumsum' =  'Impacto Acumulado',
                                    'error' = 'Impacto',
                                    "input"= "Entrada")) |>    
      dplyr::mutate(class = class  |> factor(levels = c('error','real', "input", 'prediction') , ordered = F) )  |>   

      dplyr::filter(variable == plot_variable) |> 
      ggplot(aes(x=Fecha, y=value, col=class)) + 
      geom_ribbon( aes(ymin = lower_limit, ymax = upper_limit), fill = "grey90" ) +
      geom_line() + 
      # geom_vline(xintercept = event_initial, linetype = "dashed") +
      geom_vline(xintercept = dates_df$Fecha[event_initial], linetype = "dashed") +
      # geom_hline(yintercept = 0, linetype = "dashed") +
      geom_hline(data = plot_df_line_df, aes(yintercept=value), linetype = "dashed") +
      facet_wrap(~type, ncol=1, scales='free_y') +
      ggtitle(paste0(vector_name, ' - ' ,plot_variable)) +
      theme_bw() +
      theme(
        legend.title = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0),
        # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position="bottom",
        strip.text = element_text(size=11)
      )  + ylab('')
    
    
    return(plot_df_individual)

  }

  plot_aggregate <- function(plot_df_aggregate, event_initial, vector_name, exclude_input=F) {
    
    # browser()
    # 
    plot_df_aggregate <- 
    plot_df_aggregate |> 
      # dplyr::filter() |> 
      # dplyr::filter(type != 'error') |> 
      # dplyr::filter(class != 'input') |> 
      dplyr::mutate(class = class |> as.character()) |> 
      dplyr::left_join(dates_df, by='time_index') 

    if(!exclude_input) {

       # browser()

        m_df_add_1  <- plot_df_aggregate  |> 
                       dplyr::filter(class == 'input') |> 
                       dplyr::mutate(type="Entrada")

        
        m_df_add_2  <- plot_df_aggregate  |> 
                       dplyr::filter(class == 'real')   |> 
                       dplyr::mutate(type="Entrada")

        # browser()

        plot_df_aggregate  <- dplyr::bind_rows(plot_df_aggregate, m_df_add_1, m_df_add_2)

    } else {

      plot_df_aggregate <- plot_df_aggregate   |> 
                           dplyr::filter(class != 'input') 

    }

    m_plot <- plot_df_aggregate |> 
              dplyr::mutate(lower_limit = dplyr::if_else(class != 'real', lower_limit, NA_real_),
                            upper_limit = dplyr::if_else(class != 'real', upper_limit, NA_real_)) |> 
              dplyr::mutate(lower_limit = dplyr::if_else( (class == 'error') & (type == 'error'), NA_real_, lower_limit),
                            upper_limit = dplyr::if_else( (class == 'error') & (type == 'error'), NA_real_, upper_limit)) |> 
              dplyr::mutate(lower_limit = dplyr::if_else(type != 'prediction', lower_limit, NA_real_),
                            upper_limit = dplyr::if_else(type != 'prediction', upper_limit, NA_real_))  |> 
              dplyr::mutate(value = dplyr::if_else( (!is_impact) & (type == 'cumsum'), NA_real_, value))  |> 
              dplyr::mutate(value = dplyr::if_else(  (type == 'prediction') & (class == 'input') , NA_real_, value))  |>
              dplyr::mutate(type = type |> 
                              as.character() |> 
                              dplyr::recode('prediction'='Ajuste',
                                            'cumsum' =  'Impacto Acumulado',
                                            'error' = 'Impacto',
                                            "input"= "Entrada")) |> 
              # dplyr::mutate(class = class  |> factor(levels = c('error','real', 'prediction')  |> rev(), ordered = T) )  |>  
              #dplyr::mutate(class = class  |> factor(levels = c('error','real', 'prediction')  |> rev(), ordered = F) )  |>     
              dplyr::mutate(class = class  |> factor(levels = c('error','real', "input", 'prediction') , ordered = F) )  |>                         
              ggplot(aes(x=Fecha, y=value, col=class)) + 
              # geom_ribbon( aes(ymin = lower_limit, ymax = upper_limit), fill = "grey70" ) +
              geom_ribbon( aes(ymin = lower_limit, ymax = upper_limit), fill = "grey90" ) +
              geom_line(aes(linetype=class)) + 
              geom_vline(xintercept = dates_df$Fecha[event_initial], linetype = "dashed") +
              geom_hline(yintercept = 0, linetype = "dashed") +
              facet_wrap(~type, ncol=1, scales='free_y') +
              ggtitle(vector_name) +
              theme_bw() +
              theme(
                legend.title = element_blank(),
                axis.text=element_text(size=12),
                axis.title=element_text(size=14),
                legend.text = element_text(size=14),
                # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0),
                # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                legend.position="bottom",
                strip.text = element_text(size=11)
              ) + ylab('') +
              scale_linetype_manual(values=c("solid", "solid", 'solid', 'solid')) +
              #scale_linetype_manual(values=c("dashed", "solid", 'solid')) +
              scale_colour_brewer(palette = "Set1") +
              #scale_x_date(labels = scales::date_format("%m/%d"))
              scale_x_date(labels = scales::date_format("%m/%d"), breaks = BREAK_DATES )
              #scale_linetype_manual(values=c("solid", "solid", 'solid'))
    
    
    return(m_plot)
    
  }
  
  
  
  
})