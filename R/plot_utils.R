PLOT_UTILS <- modules::module({

  import('ggplot2')

  plot_individual <- function(plot_df, plot_variable, event_initial, vector_name, dates_df=NULL, break_dates=NULL, exclude_input=F) {

    temp_unique_vars  <- plot_df$variable  |> unique()

    dates_df_is_null  <- is.null(dates_df)

    if(is.null(break_dates)) {
      break_dates = waiver()
    }

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


        m_df_add_1  <- plot_df  |>
                       dplyr::filter(class == 'input') |>
                       dplyr::mutate(type="Entrada")


        m_df_add_2  <- plot_df  |>
                       dplyr::filter(class == 'real')   |>
                       dplyr::mutate(type="Entrada")


        plot_df  <- dplyr::bind_rows(plot_df, m_df_add_1, m_df_add_2)

    } else {

      plot_df <- plot_df   |>
                           dplyr::filter(class != 'input')

    }

    if(is.null(dates_df)) {
      dates_df  <- data.frame(
        time_index =  plot_df$time_index  |> unique(),
        Date =  plot_df$time_index  |> unique()

      )
    }

    m_xintercep   <- dates_df  |> 
                     filter(time_index==event_initial)  |> 
                     pull(Date)  |> 
                     first()

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
      ggplot(aes(x=Date, y=value, col=class)) +
      geom_ribbon( aes(ymin = lower_limit, ymax = upper_limit), fill = "grey90" ) +
      geom_line() +
      geom_vline(xintercept = m_xintercep, linetype = "dashed") +
      geom_hline(data = plot_df_line_df, aes(yintercept=value), linetype = "dashed") +
      facet_wrap(~type, ncol=1, scales='free_y') +
      ggtitle(paste0(vector_name, ' - ' ,plot_variable)) +
      theme_bw() +
      theme(
        legend.title = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.position="bottom",
        strip.text = element_text(size=11)
      )  + ylab('')


    if(!dates_df_is_null) {
        plot_df_individual  <- plot_df_individual + scale_x_date(labels = scales::date_format("%m/%d"), breaks = break_dates )
    }


    return(plot_df_individual)

  }

   

  plot_aggregate <- function(plot_df_aggregate, event_initial, vector_name, dates_df=NULL, break_dates=NULL, exclude_input=F) {

    # browser()
    #

    dates_df_is_null  <- is.null(dates_df)

    if(dates_df_is_null) {
      dates_df  <- data.frame(
        time_index =  plot_df_aggregate$time_index  |> unique(),
        Date =  plot_df_aggregate$time_index  |> unique()

      )
    }

    plot_df_aggregate <-
    plot_df_aggregate |>
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

    

    if(is.null(break_dates)) {
      break_dates = waiver()
    }

    m_xintercep   <- dates_df  |> 
                     filter(time_index==event_initial)  |> 
                     pull(Date)  |> 
                     first()


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
              dplyr::mutate(class = class  |> factor(levels = c('error','real', "input", 'prediction') , ordered = F) )  |>
              ggplot(aes(x=Date, y=value, col=class)) +
              geom_ribbon( aes(ymin = lower_limit, ymax = upper_limit), fill = "grey90" ) +
              geom_line(aes(linetype=class)) +
              geom_vline(xintercept = m_xintercep, linetype = "dashed") +
              geom_hline(yintercept = 0, linetype = "dashed") +
              facet_wrap(~type, ncol=1, scales='free_y') +
              ggtitle(vector_name) +
              theme_bw() +
              theme(
                legend.title = element_blank(),
                axis.text=element_text(size=12),
                axis.title=element_text(size=14),
                legend.text = element_text(size=14),
                legend.position="bottom",
                strip.text = element_text(size=11)
              ) + ylab('') +
              scale_linetype_manual(values=c("solid", "solid", 'solid', 'solid')) +
              scale_colour_brewer(palette = "Set1") 

    if(!dates_df_is_null) {
        m_plot  <- m_plot + scale_x_date(labels = scales::date_format("%m/%d"), breaks = break_dates )
    }

    return(m_plot)

  }




})
