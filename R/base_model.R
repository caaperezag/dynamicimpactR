BaseImpactModel <- R6::R6Class('BaseImpactModel', public = list(
  name = NA_character_,
  event_initial = NA_integer_,
  X_data = NA_real_,
  Y_data = NA_real_,
  vector_name = NA_character_,
  variables_names = NA_character_,
  confidence_level=0.9,
  initialize = function(name='model impact', event_initial=NULL, X_data, Y_data, vector_name, variables_names, confidence_level) {

    self$name <-  name

    if(is.null(event_initial)) {
      event_initial  <- dim(X_data)[1] # the total time
    }

    self$event_initial <-  event_initial

    # browser()
    #

    if( (confidence_level >= 1) |  (confidence_level <= 0)) {
      stop("Confidence level must be between 0 and 1")
    }

    private$.original_variance <-  NA

    temp_n_dim <- length(dim(Y_data))

    if(temp_n_dim == 3) {
      private$.original_variance <- MODULE_IMPACT$estamate_ml_from_array( Y_data [1:event_initial,,] )$U
    }

    if(length(dim(Y_data)) == 2) {
      Y_data <-  Y_data |> array(dim = c(dim(Y_data)[1], 1, dim(Y_data)[2]))
    }

    if(length(dim(X_data)) == 2) {
      X_data <-  X_data |> array(dim = c(dim(X_data)[1], 1, dim(X_data)[2]))
    }

    if(is.na(private$.original_variance)) {
      private$.original_variance <- MODULE_IMPACT$estamate_ml_from_array( Y_data [1:event_initial,,] )$U
    }


    private$.original_x  <- X_data
    private$.original_y  <- Y_data

    private$.scaled_data_x  <- X_data |> MODULES_SCALE$scale_matrix()
    private$.scaled_data_y  <- Y_data |> MODULES_SCALE$scale_matrix()

    self$X_data <-  private$.scaled_data_x$scaled_matrix
    self$Y_data <-  private$.scaled_data_y$scaled_matrix

    self$vector_name <-  vector_name
    self$variables_names <-  variables_names



  },

  get_end_time = function() {
    return( dim(private$.original_x)[1] )
  },

  fit = function() {

    stop('No implemented')

  },
  predic = function() {

    stop('No implemented')

  },
  plot = function() {
    stop('No implemented')
  }



),
private = list(
  .residuals = NA_real_,
  .original_x = NA_real_,
  .original_y = NA_real_,
  .scaled_data_x = NA_real_,
  .scaled_data_y = NA_real_,
  .original_variance = NA_real_,
  .get_unscaled_Y = function() {
    MODULES_SCALE$unscale_matrix(input_matrix = self$Y_data,
                                 result_list = private$.scaled_data_y)
  },
  .get_unscaled_X = function() {
    MODULES_SCALE$unscale_matrix(input_matrix = self$X_data,
                                 result_list = private$.scaled_data_x)
  }
))


