#' R6 class for SpatialVectorModel
#'
#'
#'
#' @details
#' This class is used to fit Bayesian model when the response is a vector with associated coordiantes for each value.
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
#' * `n_spatial_kernels`: number of spatial kernels to use, default 2.
#' @export
SpatialModel <- R6::R6Class('SpatialVectorModel',
                        inherit = StanModelVector,

                        #' @details
                        #' Create a new `SpatialModel` object.
                        #'
                        #' @param  name: name of the model.
                        #' @param  y_data: response matrix.
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
                        #' @param  `coordinates`: matrix of coordinates, must be in a planar projection.
                        public = list(
                              n_spatial_kernels = NA_integer_,
                              initialize  = function(name='Model Vector spatial', event_initial,
                                                     X_data, Y_data, vector_name,
                                                     variables_names, confidence_level,
                                                     n_simul=20000, n_chains=4, n_cores=NA,
                                                     predefined_cov_matrix_type = c("mcmc", "ml", "identity"),
                                                     predefined_cov_matrix = NULL,
                                                     stan_fit=NA_real_,
                                                     n_spatial_kernels=2,
                                                     thin=1,
                                                     dates=NULL,
                                                     log_x=FALSE, log_y=FALSE,
                                                     coordinates) {


                                  #self$thin  <- thin

                                  super$initialize(name=name, event_initial=event_initial,
                                                   X_data=X_data, Y_data=Y_data, vector_name=vector_name,
                                                   variables_names=variables_names, confidence_level=confidence_level,
                                                   n_simul=n_simul, n_chains=n_chains, n_cores=n_cores,
                                                   predefined_cov_matrix=predefined_cov_matrix,
                                                   predefined_cov_matrix_type=predefined_cov_matrix_type,
                                                   stan_fit=stan_fit,
                                                   thin=thin,
                                                   dates=dates,
                                                   log_x=log_x, log_y=log_y
                                                   )



                                  private$.model_path <- stanmodels$model_vector_full_spatial
                                  private$.predict_model_path  <- stanmodels$model_vector_full_spatial_predict

                                  if(n_spatial_kernels < 1) {

                                    stop("The number of spatial kernels must be a atleast 2")

                                  }

                                  self$n_spatial_kernels <- n_spatial_kernels

                                  # browser()


                                  if(length(dim(coordinates)) != 2) {
                                    stop("coordinates must be a 2D array")
                                  }


                                  if( dim(coordinates)[1] !=  private$.N_elem) {

                                    stop("The size of the first dimension of coordiantes must be equal to the number of variables")

                                  }

                                  if( dim(coordinates)[2] !=  2) {

                                    stop("The size of the second dimension of coordinates must be 2, only 2D coordiantes are supported")

                                  }

                                  private$.coordiates = coordinates
                                  private$.scaled_coordinates = private$.coordiates |> MODULES_SCALE$scale_matrix(FALSE)




                              },

                              #' @details
                              #' Fit the model using MCMC.
                              fit = function() {

                                super$fit()

                                private$.unscaled_centroids  <- MODULES_SCALE$unscale_array_3d(input_array = private$.extracted_data$kernels,
                                                                                               result_list = private$.scaled_coordinates)

                              },

                              get_unscaled_centroids = function() {

                                result <- private$.unscaled_centroids |> apply(c(2, 3), mean)

                                return(result)

                              },

                              plot_map = function(base_shape=NA, coord_system) {


                                if(! all( c("tmap", "sf") %in% rownames(installed.packages())) ) {

                                  stop("you need tmap an sf installed in order to use this functon")

                                }

                                library(sf)
                                library(tmap)

                                m_plot_df <- private$.build_points_df(coord_system)

                                if(!is.na(base_shape)) {

                                  base_shape <-  base_shape + st_transform(base_shape)

                                  mapa_result <-

                                    tm_shape(base_shape) +
                                    tm_grid() +
                                    tm_polygons(col = RColorBrewer::brewer.pal(8, 'Blues')[3],
                                                palette = "-Blues",
                                                title = "", contrast = 0.7, border.col = "grey30", id = "name",
                                                legend.show = FALSE)

                                } else {
                                  mapa_result <- tm_grid()
                                }

                                mapa_result <- mapa_result  +
                                  tm_shape(m_plot_df) +
                                  tm_dots("Type", border.col = "black",
                                          palette = "-RdYlGn", size = 0.03)  +




                                if(!is.na(base_shape)) {

                                  mapa_result <-  mapa_result + tm_borders(lty = "solid", lwd =1.5)
                                }

                                mapa_result <- mapa_result +
                                  tm_compass(type = "8star", position = c("center", "top")) +
                                  tm_scale_bar()


                                return(mapa_result)

                              }
                          ),

                          private = list(
                            .coordiates = NA_real_,
                            .scaled_coordinates = NA_real_,
                            .unscaled_centroids = NA_real_,

                            .get_stan_data = function(event_initial) {

                               event_initial = private$.get_event_initial(event_initial)

                               temp_coordinates_lower = min( private$.scaled_coordinates$scaled_matrix )
                               temp_coordinates_upper = max( private$.scaled_coordinates$scaled_matrix )

                               temp_predefined_var  <-  private$.get_predefined_cov_matrix(event_initial)

                                stan_data = list(
                                  N = private$.N_time,
                                  N_before = event_initial,
                                  K = private$.N_elem,
                                  P = private$.N_pred_var,
                                  Y = self$Y_data[,1,],
                                  X = self$X_data[,1,],
                                  use_predefined_stations_var = private$.use_predefined_stations_var,
                                  predefined_stations_var = temp_predefined_var,

                                  J = self$n_spatial_kernels,
                                  COORDINATES = private$.scaled_coordinates$scaled_matrix,

                                  coordinates_lower = temp_coordinates_lower,
                                  coordinates_upper = temp_coordinates_upper

                                )


                               return(stan_data)

                             },

                            .build_points_df = function(coord_system) {


                              # browser()

                              m_coordinates <-    private$.coordiates |> data.frame()
                              m_coordinates$Type <- "Data"

                              colnames(m_coordinates) <- c("v1", "v2", "Type")

                              m_centroids <- self$get_unscaled_centroids() |> data.frame()
                              m_centroids$Type <- "Kernel"

                              colnames(m_centroids) <- c("v1", "v2", "Type")

                              m_df <- rbind(m_coordinates, m_centroids) |>
                                       sf::st_as_sf(coords = 1:2, crs = coord_system)




                             return(m_df)
                            }


                        )

)
