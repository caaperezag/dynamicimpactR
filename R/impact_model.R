#' R6 class for StanModelVector
#'
#'
#'
#' @details
#' This class is used to fit Bayesian model when the response is a vector.
#'
#' @section Public fields:
#' * `name`: name of the model.
#' * `event_initial`: The time when the event starts.
#' * `vector_name`: name of the response variable (only used in plots).
#' * `variables_names`: name of each of the components of vector (only used in plots).
#' * `confidence_level`: level of the credible interval.
#' * `n_simul`: number of the iteration for MCMC, half of the iteration will be burned, default 2000.
#' * `n_chains`: number of chains for MCMC, default 4.
#' * `n_cores`: number of cores used in MCMC, if NA(default) the number of CPU cores will be used.
#' * `thin`: thinning used in MCMC, default 1.
#' @export
StanModelVector <- R6::R6Class('StanModelVector',
                           inherit = BaseImpactModel,
                           public = list(
                             n_simul = NA_integer_,
                             n_chains = NA_integer_,
                             n_cores = NA_integer_,
                             thin = NA_integer_,

                             #' @details
                             #' Create a new `StanModelVector` object.
                             #'
                             #' @param  name: name of the model.
                             #' @param  Y_data: response matrix.
                             #' @param  X_data: covariates matrix.
                             #' @param  event_initial: The time when the event starts.
                             #' @param  vector_name: name of the response variable (only used in plots).
                             #' @param  variables_names: name of each of the components of vector (only used in plots).
                             #' @param  confidence_level: level of the credible interval.
                             #' @param  n_simul: number of the iteration for MCMC, half of the iteration will be burned, default 2000.
                             #' @param  n_chains: number of chains for MCMC, default 4.
                             #' @param  n_cores: number of cores used in MCMC, if NA(default) the number of CPU cores will be used.
                             #' @param  thin: thinning used in MCMC, default 1.
                             #' @param  predefined_cov_matrix: If you want used a predefined covariance matrix you can used this variable.
                             #' @param  predefined_cov_matrix_type: The covariance matrix type, if = mcmc  the covariance matrix is is estimated using MCMC, if ml the covariance matrix is estimated used maximun likehood and then "plugin" as constant in MCMC, if = identity the covariance matrix is "plugin"
                             initialize  = function(name='Model Vector', event_initial,
                                                    X_data, Y_data, vector_name="vector model",
                                                    variables_names, confidence_level,
                                                    n_simul=2000, n_chains=4, n_cores=NA,
                                                    predefined_cov_matrix_type = c("mcmc", "ml", "identity"),
                                                    predefined_cov_matrix = NULL,
                                                    stan_fit=NA_real_,
                                                    thin=1,
                                                    log_x=FALSE, log_y=FALSE,
                                                    dates=NULL) {

                               #browser()

                               self$thin  <- thin

                               super$initialize(name=name, event_initial=event_initial,
                                                X_data=X_data, Y_data=Y_data,
                                                vector_name=vector_name, variables_names=variables_names,
                                                confidence_level=confidence_level,
                                                log_x=log_x, log_y=log_y, dates=dates)

                               if(!is.na(stan_fit)) {
                                 private$.stan_result <- stan_fit
                               }


                               # browser()
                               self$n_simul <- n_simul
                               self$n_chains <- n_chains

                               private$.model_path <- stanmodels$model_vector_full
                               private$.predict_model_path  <- stanmodels$model_vector_full_predict


                               # if(length(dim(Y_data)) == 2) {
                               #   self$.original_y <-  Y_data |> array(dim = c(dim(Y_data)[1], 1, dim(Y_data)[2]))
                               # }
                               #
                               # if(length(dim(X_data)) == 2) {
                               #   self$.original_x  <-  X_data |> array(dim = c(dim(X_data)[1], 1, dim(X_data)[2]))
                               # }

                               # private$.scaled_data_x  <- X_data |> MODULES_SCALE$scale_matrix()
                               # private$.scaled_data_y  <- Y_data |> MODULES_SCALE$scale_matrix()
                               #
                               # self$X_data <-  private$.scaled_data_x$scaled_matrix
                               # self$Y_data <-  private$.scaled_data_y$scaled_matrix


                               # browser()

                               if(!is.na(dim(private$.original_y)[3])) {
                                 private$.N_elem <- dim(private$.original_y)[3]
                               } else {
                                 private$.N_elem <- dim(private$.original_y)[2]
                               }


                               if(!is.na(dim(X_data)[3])) {
                                 private$.N_pred_var <- dim(X_data)[3]
                               } else {

                                 private$.N_pred_var <- dim(X_data)[2]
                               }

                               # private$.N_pred_var <- dim(X_data)[2]
                               # private$.N_pred_var <- dim(X_data)[3]

                               # private$.N_elem <- dim(Y_data)[3]
                               private$.N_time <- dim(private$.original_y)[1]

                               if(!is.na(n_cores)) {
                                 n_cores = parallel::detectCores()
                               }

                               self$n_cores <- n_cores



                              private$.predefined_cov_matrix_type = match.arg(predefined_cov_matrix_type)

                              private$.predefined_cov_matrix <-  predefined_cov_matrix






                             },
                             #' @details
                             #' Fit the model using MCMC.
                             fit = function() {

                               # browser()


                               stan_data = private$.get_stan_data(
                                 # self$get_end_time(),
                                 private$.get_event_initial_or_end_time()
                               )


                               options(mc.cores = self$n_cores)


                               if(is.na(private$.stan_result)) {

                                 # browser()
                                 private$.stan_result <- rstan::sampling(private$.model_path,
                                                                         data = stan_data,
                                                                         chains = self$n_chains,
                                                                         iter=self$n_simul,
                                                                         thin=self$thin)

                               }



                               private$.extracted_data <- rstan::extract(private$.stan_result)


                             },
                             #' @details
                             #' make a prediction using MCMC.
                             #' @param event_initial time when event of interest start, if NULL the event_initial set in the initialization is used.
                             predict = function(event_initial=NULL) {

                               event_initial = private$.get_event_initial(event_initial)

                               stan_data = private$.get_stan_data(event_initial)


                               result_redraw <- rstan::gqs(private$.predict_model_path,
                                                           draws = as.matrix(private$.stan_result),
                                                           data=stan_data)

                               # rm(stan_data)

                               UTILS$gc_quiet()

                               extracted_data2 <- result_redraw |> rstan::extract()

                               rm(result_redraw)

                               UTILS$gc_quiet()

                               private$.extracted_data  <- NA_real_

                               UTILS$gc_quiet()

                               private$.extracted_data <- rstan::extract(private$.stan_result)

                               UTILS$gc_quiet()

                               private$.extracted_data <- private$.extracted_data  |> append(extracted_data2)

                               rm(extracted_data2)

                               UTILS$gc_quiet()

                               private$.build_plot_df(event_initial)

                               UTILS$gc_quiet()

                               return(NULL)

                             },

                             #' @details
                             #' plot a single variable.
                             #' @param plot_variable name of the variable to plot, the name must be in the vector variables_names.
                             #' @return a ggplot object.
                             plot_individual = function(plot_variable, event_initial=NULL) {

                               UTILS$gc_quiet()

                               private$.can_plot(event_initial)

                               event_initial = private$.get_event_initial(event_initial)

                               plot_df_individual  <- PLOT_UTILS$plot_individual(
                                 plot_variable = plot_variable,
                                 plot_df = private$.plot_df,
                                 event_initial = event_initial,
                                 vector_name = self$vector_name)

                               return(plot_df_individual)


                             },

                             #' @details
                             #' plot an aggregate of all variables.
                             #' @return a ggplot object.
                             plot_aggregate = function(event_initial=NULL) {

                               UTILS$gc_quiet()

                               private$.can_plot(event_initial)

                               event_initial = private$.get_event_initial(event_initial)

                               plot_df_aggregate  <- PLOT_UTILS$plot_aggregate(
                                 plot_df_aggregate = private$.plot_df_aggregate,
                                 event_initial = event_initial,
                                 vector_name = self$vector_name,
                                 dates_df = self$get_dates_df()
                                )

                               return(plot_df_aggregate)

                             },

                             #' @details
                             #' default plot method
                             #' @return a ggplot object.
                             plot = function() {

                              UTILS$gc_quiet()

                              event_initial  <- private$.get_event_initial_or_end_time(NULL)

                              private$.can_plot(event_initial)

                              current_plot  <-  MODULE_PLOT_EXTRA$plot_model(self, event_initial, "global")

                              return(current_plot)
                             },

                             #' @details
                             #' plot a diagnosis measure.
                             #' @param parameter_name name of the parameter of the model to plot, if parameter_name=NA_character_ and plot_type = boxplot all of the parameters are ploted.
                             #' @param measure_name name of the measure plot must one of: "Rhat", "n_eff", "mean".
                             #' @param plot_type The possible values are interval for a confidence interval or boxplot.
                             #' @return a ggplot object.
                             plot_diagnosis = function(parameter_name="theta_vec",
                                                       measure_name="Rhat",
                                                       plot_type=c("interval", "boxplot"),
                                                       event_initial=NULL
                                                       ) {


                              event_initial = private$.get_event_initial(event_initial)


                              plot_type = match.arg(plot_type)

                              if(plot_type == "interval") {

                                plot_result  <- DIAGNOSTICS$plot_ics(private$.stan_result,
                                                                     parameter_name=parameter_name,
                                                                     variable_name=measure_name) +
                                                geom_vline(xintercept = event_initial, linetype="dashed", color="red")

                              } else {

                                  plot_result  <- DIAGNOSTICS$box_plot(private$.stan_result,
                                                                       parameter_name=parameter_name,
                                                                       variable_name=measure_name) +
                                                  geom_vline(xintercept = event_initial, linetype="dashed", color="red")
                              }

                              return(plot_result)
                             },


                             summary = function(dates_list, confidence_level=NA) {

                                private$.build_summary(
                                  dates_list=dates_list,
                                  confidence_level=confidence_level
                                )

                                MODULE_SUMMARY$result_summary(private$.summary_result )


                             }



                           ),
                           private = list(
                             .N_pred_var = NA_integer_,
                             .summary_result = NA_real_,
                             .stan_result = NA_real_,
                             .predefined_cov_matrix = NA_real_,
                             .use_predefined_stations_var = 0,
                             .residuals = NA_real_,
                             .model_path = NA_character_,
                             .N_elem = NA_integer_,
                             .N_time = NA_integer_,
                             .extracted_data = NA_real_,
                             .plot_df = NA,
                             .plot_df_aggregate = NA,
                             .scaled_data_x = NA_real_,
                             .scaled_data_y = NA_real_,
                             .predict_model_path = NA_character_,
                             .predefined_cov_matrix_type = NA_real_,

                             .get_predefined_cov_matrix  = function(event_initial) {

                              event_initial = private$.get_event_initial(event_initial)
                              private$.use_predefined_stations_var = 1

                              predefined_cov_matrix  <- private$.predefined_cov_matrix

                              if(private$.predefined_cov_matrix_type != "mcmc") {

                                 if(!is.matrix(predefined_cov_matrix)) {

                                   predefined_cov_matrix  <- MODULE_IMPACT$get_variable_matrix(
                                     matrix_type=private$.predefined_cov_matrix_type,
                                     Y_data=private$.scaled_data_y$scaled_matrix,
                                     event_initial=event_initial
                                   )


                                 }




                              }

                              if(is.null(predefined_cov_matrix)) {

                                  # in this case the matrix is not use in stan
                                  predefined_cov_matrix = matrix(0,
                                                                 private$.N_elem,
                                                                 private$.N_elem)

                                  private$.use_predefined_stations_var = 0

                              }

                              return(predefined_cov_matrix)

                             },

                             .get_stan_data = function(event_initial) {

                              #browser()

                               event_initial = private$.get_event_initial(event_initial)

                               temp_predefined_var  <-  private$.get_predefined_cov_matrix(event_initial)

                               stan_data = list(
                                 N = private$.N_time,
                                 N_before = event_initial,
                                 K = private$.N_elem,
                                 P = private$.N_pred_var,
                                 Y = self$Y_data[,1,],
                                 X = self$X_data[,1,],
                                 use_predefined_stations_var = private$.use_predefined_stations_var,
                                 #predefined_stations_var = private$.predefined_cov_matrix,
                                 predefined_stations_var = temp_predefined_var,
                                 use_log_x = 1,
                                 use_log_y = 1
                               )

                               return(stan_data)

                             },
                             .build_df = function(m_matrix, station_names, event_initial) {

                               event_initial = private$.get_event_initial(event_initial)

                               df <- data.frame(m_matrix)
                               colnames(df) <- station_names
                               df$t <- 1:nrow(df)

                               df <- df |> dplyr::pivot_longer(-t, names_to=self$vector_name, values_to="value")

                               df$before = ifelse(df$t > event_initial, "after", "before")

                               return(df)
                             },



                             .get_ic_from_variable = function(m_array, global=FALSE) {

                               temp_array <- m_array

                               if(global) {

                                 temp_array <- m_array |>  apply(c(1,2), mean)

                               }

                               temp_ic <-
                                 apply(temp_array, c(2),
                                       function(x, ci) {
                                         result <- bayestestR::hdi(x, ci=ci)
                                         return ( c(result$CI_low, result$CI_high))
                                       }, ci=self$confidence_level)



                               result_df <- data.frame(value = temp_array |> apply(2, median),
                                                       time_index = 1:private$.N_time )


                               result_df$lower_limit <-  temp_ic[1,]
                               result_df$upper_limit <-  temp_ic[2,]




                               return(result_df)

                             },



                             .build_plot_df = function(event_initial=NULL) {

                              #browser()

                               event_initial = private$.get_event_initial(event_initial)

                               df_list <- list()
                               df_list_aggregate <- list()

                               y_pred_unscaled <- private$.extracted_data$Y_pred |>
                                 MODULES_SCALE$unscale_array_3d(result_list = private$.scaled_data_y)

                               difference_unscaled <- private$.extracted_data$Y_pred |>
                               #difference_unscaled <- -private$.extracted_data$difference  |>
                                 MODULES_SCALE$unscale_array_3d(
                                   result_list = private$.scaled_data_y,
                                   m_diff_array = private$.original_y[,1,]
                                   #m_diff_array = NULL
                                 )



                               cumsum_unscaled <- private$.extracted_data$Y_pred |>
                               # cumsum_unscaled <- private$.extracted_data$difference |>
                                 MODULES_SCALE$unscale_cumsum(result_list = private$.scaled_data_y,
                                                              start_event = event_initial,
                                                              m_diff_array = private$.original_y[,1,]
                                                              #m_diff_array = NULL
                                                              )





                               for(m_index in 1:(length(self$variables_names)+1)  ) {

                                 is_global <- m_index > length(self$variables_names)



                                 # idx = ifelse(is_global, 1:length(self$variables_names), m_index )

                                 if(is_global) {
                                   m_variable_name <-  'global'
                                   idx = 1:length(self$variables_names)
                                 } else {
                                   idx = m_index
                                   m_variable_name <-  self$variables_names[idx]
                                 }


                                 # m_df_stan_pred <- private$.extracted_data$Y_pred[,,idx]  |>
                                 m_df_stan_pred <- y_pred_unscaled[,,idx]  |>
                                   private$.get_ic_from_variable(global=is_global)
                                 m_df_stan_pred$variable <- m_variable_name
                                 m_df_stan_pred$type <- "prediction"
                                 m_df_stan_pred$class <- "prediction"




                                 if(is_global) {

                                   # browser()

                                   m_df_real <- data.frame(

                                     variable=m_variable_name,
                                     time_index = 1:(dim(self$X_data)[1]),
                                     # value = self$Y_data |> apply(c(1), mean),
                                     value = private$.original_y |> apply(c(1), mean),
                                     type = "prediction",
                                     class = 'real'

                                   )

                                   m_df_real_input <- data.frame(

                                     variable=m_variable_name,
                                     time_index = 1:(dim(self$X_data)[1]),
                                     # value = self$Y_data |> apply(c(1), mean),
                                     value = private$.original_x |> apply(c(1), mean),
                                     type = "prediction",
                                     class = 'input'

                                   )

                                 } else {
                                   # browser()

                                   m_df_real <- data.frame(

                                     variable=m_variable_name,
                                     time_index = 1:(dim(self$X_data)[1]),
                                     # value = self$Y_data[,1, idx],
                                     value = private$.original_y[,1, idx],
                                     type = "prediction",
                                     class = 'real'

                                   )

                                   m_df_real_input <- data.frame(

                                     variable=m_variable_name,
                                     time_index = 1:(dim(self$X_data)[1]),
                                     # value = self$Y_data[,1, idx],
                                     value = private$.original_x[,1, idx],
                                     type = "prediction",
                                     class = 'input'

                                   )
                                 }

                                 m_df_real_input$upper_limit <- NA
                                 m_df_real_input$lower_limit <- NA


                                 m_df_real$upper_limit <- NA
                                 m_df_real$lower_limit <- NA


                                 m_df_stan_error <- difference_unscaled[,,idx]  |>
                                   private$.get_ic_from_variable(global=is_global)
                                 m_df_stan_error$variable <- m_variable_name
                                 m_df_stan_error$type <- "error"
                                 m_df_stan_error$class <- "error"


                                 # m_df_stan_error <- private$.extracted_data$cumsum_only_after[,,idx] |>
                                 m_df_stan_error_cusum <- cumsum_unscaled[,,idx] |>
                                   private$.get_ic_from_variable(global=is_global)
                                 m_df_stan_error_cusum$variable <- m_variable_name
                                 m_df_stan_error_cusum$type <- "cumsum"
                                 m_df_stan_error_cusum$class <- "error"

                                 if(is_global) {

                                   # browser()

                                   df_list_aggregate[[length(df_list_aggregate) + 1]] <- m_df_stan_pred
                                   df_list_aggregate[[length(df_list_aggregate) + 1]] <- m_df_real
                                   df_list_aggregate[[length(df_list_aggregate) + 1]] <- m_df_stan_error
                                   df_list_aggregate[[length(df_list_aggregate) + 1]] <- m_df_stan_error_cusum
                                   df_list_aggregate[[length(df_list_aggregate) + 1]] <- m_df_real_input

                                 } else {

                                   df_list[[length(df_list) + 1]] <- m_df_stan_pred
                                   df_list[[length(df_list) + 1]] <- m_df_real
                                   df_list[[length(df_list) + 1]] <- m_df_stan_error
                                   df_list[[length(df_list) + 1]] <- m_df_stan_error_cusum
                                   df_list[[length(df_list) + 1]] <- m_df_real_input


                                 }

                               }

                               # browser()

                               private$.plot_df <- do.call(rbind, df_list)

                               private$.plot_df$type <- private$.plot_df$type |>
                                 factor(levels=c("prediction", 'error', 'cumsum'))

                               private$.plot_df$is_impact  <- private$.plot_df$time_index  >= event_initial


                               private$.plot_df_aggregate <- do.call(rbind, df_list_aggregate)
                               private$.plot_df_aggregate$type <- private$.plot_df_aggregate$type |>
                                 factor(levels=c("prediction", 'error', 'cumsum'))

                               private$.plot_df_aggregate$is_impact  <-
                                 private$.plot_df_aggregate$time_index  >= event_initial


                               return(private$.plot_df)

                             },

                             .build_summary = function(dates_list, confidence_level) {

                                if(is.na(confidence_level)) {
                                  confidence_level = self$confidence_level
                                }

                                if(is.na(confidence_level)) {
                                  stop("The confidence level can't be NA")
                                }

                                dates_df  <- self$get_dates_df()

                                private$.summary_result  <- MODULE_SUMMARY$get_get_multiple_impacts_stan(
                                      model=self,
                                      impact_list=dates_list,
                                      dates_df=dates_df,
                                      ci=confidence_level
                                )

                                return(private$.summary_result)

                             },

                             .can_plot = function(event_initial) {

                                can_plot = TRUE

                                if(!is.data.frame(private$.plot_df)) {
                                  can_plot = FALSE
                                }

                                if(!is.data.frame(private$.plot_df_aggregate)) {
                                  can_plot = FALSE
                                }

                                if(!can_plot) {

                                  message("The model will make the predictions before ploting")

                                  self$predict(event_initial=event_initial)
                                }



                                return(can_plot)

                             }


                           )

)
