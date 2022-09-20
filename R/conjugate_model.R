#' @export
ConjugateModel <- R6::R6Class('ConjugateModel',
                           inherit = BaseImpactModel,
                           
                           public = list(
                             discount_factor = NULL,
                             n_simul = 1000,
                             initialize  = function(name='model impact', event_initial, X_data, Y_data, 
                                                    vector_name, variables_names, confidence_level, 
                                                    n_simul, discount_factor=NULL, dates=NULL, log_x=FALSE, log_y=FALSE) {

                               # super$initialize(name, event_initial, X_data, Y_data, vector_name, variables_names, confidence_level, dates)
                               
                               super$initialize(name=name, event_initial=event_initial, 
                                                X_data=X_data, Y_data=Y_data,
                                                vector_name=vector_name, variables_names=variables_names,
                                                confidence_level=confidence_level, 
                                                log_x=log_x, log_y=log_y, dates=dates)

                               self$n_simul <- n_simul

                              if(length(discount_factor) > 1) {
                                stop("The discount factor must an scalar")
                              }

                              if(!is.null(discount_factor)) {

                                if( (discount_factor < 0) | (discount_factor > 1)  ) {
                                  stop("The discount factor must in (0, 1)") 
                                }

                               }

                               self$discount_factor  <- discount_factor

                               

                               
                             },


                             fit = function() {

                               # browser()

                               private$.fit()

                             },
                             predic = function(event_initial=NULL) {

                               private$.fit(event_initial=event_initial)

                             },
                             plot_individual = function(plot_variable) {

                               # browser()

                               # if(is.na(private$.simul_result) | is.na(private$.simul_sumation)) {
                               #   # TODO agregar despues
                               #   private$.build_simulation()
                               # }

                               plot_df_individual  <- PLOT_UTILS$plot_individual(
                                 plot_variable = plot_variable,
                                 plot_df = private$.plot_df,
                                 event_initial = self$event_initial,
                                 vector_name = self$vector_name)

                               return(plot_df_individual)

                             },

                             plot_aggregate = function() {

                               plot_df_aggregate  <- PLOT_UTILS$plot_aggregate(
                                 plot_df_aggregate = private$.plot_df_aggregate,
                                 event_initial = self$event_initial,
                                 vector_name = self$vector_name)

                               return(plot_df_aggregate)

                             }
                           ),

                           private = list(

                             .residuals = NA_real_,
                             .fitted_model = NA,
                             .plot_df = NA,
                             .plot_df_scaled = NA,
                             .plot_df_aggregate = NA,
                             .simul_result = NA,
                             .simul_resul_scaled = NA,
                             .simul_sumation = NA,
                             .simul_sumation_scaled = NA,
                             .fitted_model_orignal = NA,
                             .original_variance = NA_real_,


                             .build_data_list = function(event_initial=NULL) {

                               # browser()

                               event_initial  <- private$.get_event_initial(event_initial)

                               m_list <- list(

                                 N_array = dim(self$X_data)[1],
                                 N_before = event_initial,

                                 parameter_names = self$vector_name,
                                 estaciones_names = self$variables_names,

                                 N_contaminates = 1,
                                 N_estacioens = dim(self$X_data)[3],

                                 matrix_x = self$X_data,
                                 matrix_y = self$Y_data
                               )

                               m_list <-  m_list |>  MODULE_IMPACT$build_data_list_from_data()


                               return(m_list)

                             },
                             .build_individual_ic = function(plot_df, variable_name, variable_type) {

                               alpha <- 1-self$confidence_level

                               # result_df <- private$.plot_df |>
                               #              # dplyr::filter(type=variable_type, variable=variable_name)
                               #              dplyr::filter(variable=variable_name)


                               result_df <-  plot_df |>
                                 dplyr::filter(variable==variable_name)

                               N <- plot_df$time_index |> unique() |> length()


                               temp_join <- data.frame(

                                 time_index = plot_df$time_index |> unique() |> sort(),
                                 var_estimation = rep(NA_real_, N)

                               )

                               m_coefficients = (variable_name == self$variables_names)  |> as.numeric()

                               # browser()
                               for(index in 1:N) {
                                 # browser()
                                 temp_join$var_estimation[index] =
                                   MODULES_IC_SIMULATION$sum_of_normals_var(coefficients=m_coefficients,
                                                                            cov_matrix = private$.fitted_model$S_t[index,,])
                                 # browser()
                                 # temp_join$var_estimation[index] =
                                 # MODULES_IC_SIMULATION$sum_of_normals_var(coefficients=m_coefficients,
                                 #                                          cov_matrix = kronecker( private$.fitted_model$S_t[index,,], private$.original_variance) )





                               }
                               # browser()

                               result_df <- result_df |>
                                 dplyr::left_join(temp_join, by = 'time_index') |>
                                 dplyr::mutate(var_estimation = dplyr::if_else(class == 'real', NA_real_, var_estimation) ) |>
                                 dplyr::mutate(var_estimation = dplyr::if_else( (class == 'input') & (type == "prediction"), NA_real_, var_estimation) ) |>
                                 dplyr::mutate(
                                   upper_limit = value + (qnorm(1-(alpha/2)) * sqrt(var_estimation) ),
                                   lower_limit = value + (qnorm((alpha/2))   * sqrt(var_estimation) )
                                 )

                               return(result_df)

                             },
                             .build_aggregate_ic = function(plot_df_aggregate) {

                               # browser()

                               result_df <- plot_df_aggregate


                               alpha <- 1-self$confidence_level

                               N <- plot_df_aggregate$time_index |> unique() |> length()

                               temp_join <- data.frame(

                                 time_index = plot_df_aggregate$time_index |> unique() |> sort(),
                                 var_estimation = rep(NA_real_, N)

                               )

                               m_coefficients = rep(1,length(self$variables_names))/length(self$variables_names)

                               # browser()
                               for(index in 1:N) {
                                 # browser()

                                 # temp_join$var_estimation[index] =
                                 #   MODULES_IC_SIMULATION$sum_of_normals_var(coefficients=m_coefficients,
                                 #                                            cov_matrix = kronecker( private$.fitted_model$S_t[index,,], private$.original_variance) )

                                 temp_join$var_estimation[index] =
                                   MODULES_IC_SIMULATION$sum_of_normals_var(coefficients=m_coefficients,
                                                                            cov_matrix = private$.fitted_model$S_t[index,,] )

                               }

                               result_df <- result_df |>
                                 dplyr::left_join(temp_join, by = 'time_index') |>
                                 dplyr::mutate(var_estimation = dplyr::if_else(class == 'real', NA_real_, var_estimation) ) |>
                                 dplyr::mutate(var_estimation = dplyr::if_else( (class == 'input') & (type == "prediction"), NA_real_, var_estimation) ) |>
                                 dplyr::mutate(
                                   upper_limit = value + (qnorm(1-(alpha/2)) * sqrt(var_estimation) ),
                                   lower_limit = value + (qnorm((alpha/2))   * sqrt(var_estimation) )
                                 )




                             },

                             .build_plot_scaled_df = function() {

                               df_list <- list()


                               for(idx in 1:length(self$variables_names)) {
                                 # browser()
                                 m_df <- data.frame(

                                   variable=self$variables_names[idx],
                                   time_index = 1:dim(self$X_data)[1],
                                   value = private$.fitted_model$Y_UP[, idx, 1],
                                   type = "prediction",
                                   class = "prediction"

                                 )

                                 m_df_real <- data.frame(

                                   variable=self$variables_names[idx],
                                   time_index = 1:dim(self$X_data)[1],
                                   # value = private$.fitted_model_orignal$y_t_full[, idx],
                                   value = private$.fitted_model$y_t_full[, idx],
                                   type = "prediction",
                                   class = 'real'

                                 )

                                 m_df_real$upper_limit <- NA
                                 m_df_real$lower_limit <- NA
                                 m_df_real$var_estimation <- NA

                                 m_df_real_input <- data.frame(

                                   variable=self$variables_names[idx],
                                   time_index = 1:(dim(self$X_data)[1]),
                                   value = self$X_data[, 1, idx],
                                   type = "prediction",
                                   class = 'input'

                                 )

                                 m_df_real_input$upper_limit <- NA
                                 m_df_real_input$lower_limit <- NA
                                 m_df_real_input$var_estimation <- NA

                                 m_df_error <- data.frame(

                                   variable=self$variables_names[idx],
                                   time_index = 1:dim(self$X_data)[1],
                                   # value = private$.fitted_model_orignal$y_t_full[, idx] - private$.fitted_model_orignal$Y_UP[, idx, 1],
                                   value = private$.fitted_model$y_t_full[, idx] - private$.fitted_model_orignal$Y_UP[, idx, 1],
                                   # value = private$.fitted_model_orignal$E_UP[, idx],
                                   type = "error",
                                   class = 'error'

                                 )



                                 m_df$is_impact       <- m_df$time_index       >= self$event_initial
                                 m_df_error$is_impact <- m_df_error$time_index >= self$event_initial
                                 m_df_real$is_impact  <- m_df_real$time_index  >= self$event_initial
                                 m_df_real_input$is_impact  <- m_df_real_input$time_index  >= self$event_initial

                                 df_list[[length(df_list) + 1]] <- m_df_real

                                 df_list[[length(df_list) + 1]] <- m_df |>
                                   private$.build_individual_ic(variable_name=self$variables_names[idx],
                                                                variable_type="prediction")

                                 df_list[[length(df_list) + 1]] <- m_df_error |>
                                   private$.build_individual_ic(variable_name=self$variables_names[idx],
                                                                variable_type="error")  |>
                                   private$.add_cumsum_to_df()


                                 df_list[[length(df_list) + 1]] <- m_df_real_input
                               }

                               private$.plot_df_scaled <- do.call(rbind, df_list)


                               private$.plot_df_scaled$type <- private$.plot_df_scaled$type |>
                                 factor(levels=c("prediction", 'error', 'cumsum'))

                               return(private$.plot_df_scaled)

                             },

                             .build_plot_df = function(event_initial=NULL) {

                              browser()

                               df_list <- list()
                               df_list_aggregate <- list()

                               event_initial = private$.get_event_initial(event_initial)

                               for(m_index_grops  in 1:length(self$vector_name)) {
                                 # browser()

                                 for(m_index in 1:(length(self$variables_names)+1)) {

                                  is_global <- m_index > length(self$variables_names)

                                  if(is_global) {
                                   m_variable_name <-  'global'
                                   idx = 1:length(self$variables_names)
                                  } else {
                                    idx = m_index
                                    m_variable_name <-  self$variables_names[idx]
                                  }

                                  m_df <- data.frame(

                                   variable=self$variables_names[idx],
                                   time_index = 1: (dim(self$X_data)[1]),
                                   value = private$.fitted_model$Y_UP[,m_index_grops,idx],
                                   type = "prediction",
                                   class = "prediction",
                                   vector = vector_name[m_index_grops]

                                 )

                                 m_df_real <- data.frame(

                                   variable=self$variables_names[idx],
                                   time_index = 1: (dim(self$X_data)[1]),
                                   value = private$.fitted_model$y_t_full[, m_index_grops, idx],
                                   type = "prediction",
                                   class = 'real',
                                   vector = vector_name[m_index_grops]

                                 )

                                 m_df_real_input <- data.frame(

                                   variable=self$variables_names[idx],
                                   time_index = 1: (dim(self$X_data)[1]),
                                   value = private$.original_x[,m_index_grops, idx],
                                   type = "prediction",
                                   class = 'input',
                                   vector = vector_name[m_index_grops]

                                 )


                                 m_df_real$upper_limit <- NA
                                 m_df_real$lower_limit <- NA
                                 m_df_real$var_estimation <- NA

                                 m_df_real_input$upper_limit <- NA
                                 m_df_real_input$lower_limit <- NA
                                 m_df_real_input$var_estimation <- NA

                                 # browser()
                                 m_df_error <- data.frame(

                                   variable=self$variables_names[idx],
                                   time_index = 1: (dim(self$X_data)[1]),
                                   value =  private$.fitted_model$y_t_full[, m_index_grops, idx] - private$.fitted_model$Y_UP[, m_index_grops, idx],
                                   type = "error",
                                   class = 'error',
                                   vector = vector_name[m_index_grops]

                                 )

                                 m_df$is_impact             <- m_df$time_index       >= self$event_initial
                                 m_df_error$is_impact       <- m_df_error$time_index >= self$event_initial
                                 m_df_real$is_impact        <- m_df_real$time_index  >= self$event_initial
                                 m_df_real_input$is_impact  <- m_df_real_input$time_index  >= self$event_initial


                                 if(is_global) {

                                    df_list_aggregate[[length(df_list_aggregate) + 1]] <- m_df_real
                                    df_list_aggregate[[length(df_list_aggregate) + 1]] <- m_df_real_input

                                    df_list_aggregate[[length(df_list_aggregate) + 1]] <- m_df |>
                                      private$.build_individual_ic(variable_name=self$variables_names[idx],
                                                                    variable_type="prediction")

                                    df_list_aggregate[[length(df_list_aggregate) + 1]] <- m_df_error |>
                                      private$.build_individual_ic(variable_name=self$variables_names[idx],
                                                                    variable_type="error")  |>
                                      private$.add_cumsum_to_df()

                                 } else {

                                  df_list[[length(df_list) + 1]] <- m_df_real
                                  df_list[[length(df_list) + 1]] <- m_df_real_input

                                  df_list[[length(df_list) + 1]] <- m_df |>
                                    private$.build_individual_ic(variable_name=self$variables_names[idx],
                                                                  variable_type="prediction")

                                  df_list[[length(df_list) + 1]] <- m_df_error |>
                                    private$.build_individual_ic(variable_name=self$variables_names[idx],
                                                                  variable_type="error")  |>
                                    private$.add_cumsum_to_df()
                                  
                                 }




                                 }

                                 
                               }

                               private$.plot_df <- do.call(rbind, df_list)


                               private$.plot_df$type <- private$.plot_df$type |>
                                 factor(levels=c("prediction", 'error', 'cumsum'))

                               return(private$.plot_df)

                             },

                             .fit  =  function(event_initial=NULL, discount=NULL) {


                                if(is.null(discount)) {
                                  discount = self$discount_factor
                                }


                                m_list <- private$.build_data_list(event_initial)

                                
                                private$.fitted_model <- m_list |>  MODULE_IMPACT$run_model(discount = discount)

                               

                                private$.build_simulation()

                                private$.build_plot_df()
                                # private$.build_plot_scaled_df()
                                private$.build_plot_df_aggregate()


                                return(NULL)



                             },

                             .build_plot_df_aggregate = function() {

                               # if(is.na(private$.plot_df)) {
                               #   private$.build_plot_df()
                               # }

                               # browser()

                               # private$.plot_df_aggregate <- private$.plot_df |>
                               # private$.plot_df_aggregate <- private$.plot_df_scaled |>
                               private$.plot_df_aggregate <- private$.plot_df |>
                                 dplyr::filter(type != 'cumsum') |>
                                 dplyr::group_by(time_index, type, class) |>
                                 dplyr::summarise(
                                   value   = mean(value) #,
                                   # median_value = median(value)
                                 ) |>
                                 dplyr::ungroup() |>
                                 private$.build_aggregate_ic()


                               private$.plot_df_aggregate$is_impact <- private$.plot_df_aggregate$time_index >= self$event_initial

                               # browser()
                               private$.plot_df_aggregate <-  private$.plot_df_aggregate |>
                                 private$.add_cumsum_to_df_aggregate()

                               private$.plot_df_aggregate$type <- private$.plot_df_aggregate$type |>
                                 factor(levels=c("prediction", 'error', 'cumsum'))


                               return(private$.plot_df_aggregate)

                             },
                             .add_cumsum_to_df = function(plot_df) {

                               result_df <- plot_df |>
                                 dplyr::filter(type == 'error') |>
                                 dplyr::filter(is_impact)

                               variable_name <- result_df$variable |> dplyr::first()



                               m_variable <- (variable_name == self$variables_names)  |> as.numeric()  |> which.max()


                               print('*******************************')
                               print(variable_name)
                               print(m_variable)
                               print('*******************************')

                               # browser()

                               cumsum_df <- data.frame(
                                 variable = variable_name,
                                 time_index = (self$event_initial+1):(dim(private$.simul_sumation$mean_inter)[1] + self$event_initial ),
                                 value = -private$.simul_sumation$mean_inter[, m_variable],
                                 type = 'cumsum',
                                 class = 'error',
                                 is_impact = TRUE,
                                 var_estimation=NA_real_,
                                 # upper_limit = -private$.simul_sumation$upper_limit[, m_variable],
                                 # lower_limit = -private$.simul_sumation$lower_limit[, m_variable],
                                 upper_limit = -private$.simul_sumation$lower_limit[, m_variable],
                                 lower_limit = -private$.simul_sumation$upper_limit[, m_variable]
                               )

                               # result_df <- rbind(result_df, plot_df)
                               result_df <- rbind(cumsum_df, plot_df)

                               return(result_df)

                             },
                             .add_cumsum_to_df_aggregate = function(plot_df) {
                               # browser()

                               result_df <- plot_df |>
                                 dplyr::filter(type == 'error') |>
                                 dplyr::filter(is_impact)


                               cumsum_df <- data.frame(

                                 time_index = (self$event_initial+1):(dim(private$.simul_sumation$mean_inter)[1] + self$event_initial ),
                                 value = -private$.simul_sumation$aggregate$mean_aggregate,
                                 type = 'cumsum',
                                 class = 'error',
                                 is_impact = TRUE,
                                 var_estimation=NA_real_,
                                 upper_limit = -private$.simul_sumation$aggregate$upper_limit,
                                 lower_limit = -private$.simul_sumation$aggregate$lower_limit
                               )

                               # result_df <- rbind(result_df, plot_df)
                               result_df <- rbind(cumsum_df, plot_df)

                               return(result_df)

                             },
                             .unscaled_fitted_model=function(model) {

                               temp_fitted <- model
                               # browser()
                               temp_fitted$Y_UP <- temp_fitted$Y_UP |>
                                 MODULES_SCALE$unscale_array_3d(result_list=private$.scaled_data_y)

                               temp_fitted$E_UP <- temp_fitted$E_UP |>
                                 MODULES_SCALE$unscale_array_3d(result_list=private$.scaled_data_y)

                               temp_fitted$y_t_before <- temp_fitted$y_t_before |>
                                 MODULES_SCALE$unscale_array_3d(result_list=private$.scaled_data_y)

                               temp_fitted$y_t_after <- temp_fitted$y_t_after |>
                                 MODULES_SCALE$unscale_array_3d(result_list=private$.scaled_data_y)

                               # temp_fitted$y_t_full <- temp_fitted$y_t_full |>
                               #   MODULES_SCALE$unscale_array_3d(result_list=private$.scaled_data_y)

                               temp_fitted$y_t_after  <- abind::abind(temp_fitted$y_t_before, temp_fitted$y_t_after, along=1)

                               temp_fitted$X_t <- temp_fitted$X_t |>
                                 MODULES_SCALE$unscale_array_3d(result_list=private$.scaled_data_x)

                               return(temp_fitted)

                             },
                             
                             .build_simulation = function() {

                               alpha <- 1-self$confidence_level

                               private$.simul_result <- MODULES_IC_SIMULATION$make_theta_based(
                                 model_result = private$.fitted_model,
                                 n_simul = self$n_simul,
                                 y_scaled_data=private$.scaled_data_y,
                                 use_percent=FALSE,
                                 m_weights=NULL,
                                 MODULES_SCALE=MODULES_SCALE)

            
                               private$.simul_sumation <- MODULES_IC_SIMULATION$make_sumation(
                                 simul_model_raw = private$.simul_result$error_array,
                                 alpha = alpha,
                                 m_weights=NULL)

                               private$.fitted_model_orignal <- private$.fitted_model
                               private$.fitted_model <- private$.fitted_model |> private$.unscaled_fitted_model()

                             },


                             .make_ic_simulation = function(cumsum_df, event_time, m_variable, alpha) {



                               result_df <- cumsum_df

                               if(m_variable=='AGGREGATE') {

                                 aggregate_df <- sumation_result$aggregate

                                 aggregate_df$t <- aggregate_df$t + self$event_initial

                                 result_df <- result_df |> dplyr::left_join(
                                   aggregate_df,
                                   by=c('time_index'='t')
                                 )

                               } else {

                                 m_var_index <- (m_variable == self$variables_names) |> which.max()

                                 individual_df <- data.frame(
                                   time_index = sumation_result$t + self$event_initial,
                                   # variable = m_variable,
                                   lower_limit = sumation_result$lower_limit[, m_var_index],
                                   upper_limit = sumation_result$upper_limit[, m_var_index],
                                   value = sumation_result$cumsum_result[, m_var_index]
                                 )

                                 result_df <- result_df |> dplyr::left_join(
                                   individual_df, by='time_index'
                                 )

                               }


                             }

                           )
)
