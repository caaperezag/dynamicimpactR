BaseImpactModel <- R6::R6Class('BaseImpactModel', public = list(
  name = NA_character_,
  event_initial = NA_integer_,
  X_data = NA_real_,
  Y_data = NA_real_,
  log_x = FALSE, 
  log_y = FALSE,
  vector_name = NA_character_,
  variables_names = NA_character_,
  confidence_level=0.9,
  initialize = function(name='model impact', event_initial=NULL, X_data, Y_data, vector_name, variables_names, confidence_level,log_x, log_y, dates=NULL) {

    # browser()

    self$name <-  name

    self$log_x <- log_x
    self$log_y <- log_y

    if(self$log_x) {

      if(any(X_data <= 0) ) {
        stop("when using log elements of X must be positive")
      }

    }

    if(self$log_y) {

      if(any(Y_data <= 0) ) {
        stop("when using log elements of Y must be positive")
      }

    }
    
    if(is.null(dim(Y_data))) {
      stop('Incorrect Y_data dimensions')
    }

    if(is.null(dim(X_data))) {
      stop('Incorrect X_data dimensions')
    }

    if( !(length(dim(Y_data)) %in% c(2, 3)) ) {
      stop('Incorrect Y_data dimensions')

    }

    if( !(length(dim(X_data)) %in% c(2, 3)) ) {
      stop('Incorrect X_data dimensions')

    }


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

    if(any(is.na(Y_data))) {
      stop("There should be no missing values in Y.")
    }

    if(any(is.na(X_data))) {
      stop("There should be no missing values in X.")
    }

     if(dim(Y_data)[1] != dim(X_data)[1]) {
      stop("The length of the time series(first dim of X and Y) is diferent.")
     }

    if(length(dim(Y_data)) == 2) {
      Y_data <-  Y_data |> array(dim = c(dim(Y_data)[1], 1, dim(Y_data)[2]))
    }

    if(length(dim(X_data)) == 2) {
      X_data <-  X_data |> array(dim = c(dim(X_data)[1], 1, dim(X_data)[2]))
    }

    if(dim(Y_data)[2] != dim(X_data)[2]) {
      stop("The number of rows(second dim of X and Y) is diferent.")
    }

    # browser()


    private$.original_x  <- X_data
    private$.original_y  <- Y_data

    if( (dim(Y_data)[2]) == 1) {

      private$.scaled_data_x  <- X_data |> MODULES_SCALE$scale_matrix(self$log_x)
      private$.scaled_data_y  <- Y_data |> MODULES_SCALE$scale_matrix(self$log_y)

    } else {

      private$.scaled_data_x  <- X_data |> MODULES_SCALE$scale_3d_array(self$log_x)
      private$.scaled_data_y  <- Y_data |> MODULES_SCALE$scale_3d_array(self$log_y)

    }

    if(!is.matrix(private$.original_variance)) {
      private$.original_variance <- MODULE_IMPACT$estamate_ml_from_array( Y_data [1:event_initial,,] )$U
      
      #private$.original_variance  <- MODULE_IMPACT$estamate_ml_from_array(private$.scaled_data_y$scaled_matrix[1:event_initial,,])$U
    }

    self$X_data <-  private$.scaled_data_x$scaled_matrix
    self$Y_data <-  private$.scaled_data_y$scaled_matrix

    self$vector_name <-  vector_name
    self$variables_names <-  variables_names

    private$.set_dates_df(X_data, dates)


  },

  get_dates_df = function() {

    if(is.null(private$.dates_df)) {
      return(NULL)
    }


    if( any(is.na(private$.dates_df)) ) {
      return(NULL)
    }

    return(private$.dates_df)

  },

  get_end_time = function() {
    return( dim(private$.original_x)[1] )
  },

  fit = function() {

    stop('No implemented')

  },
  predict = function(event_initial=NULL) {

    stop('No implemented')

  },
  plot = function() {
    stop('No implemented')
  },

  summary = function(dates_list, ci=0.95) {
    stop('No implemented')
  }



),
private = list(
  .dates_df = NA_real_,
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
  },

  .get_event_initial = function(event_initial=NULL) {

    if(is.null(event_initial)) {
      event_initial = self$event_initial
    }

    if(is.null(event_initial)) {
      stop("if event_initial is no defiend when then object is created, it must be given as parameter for this function.")
    }

    return(event_initial)

  },

  .set_dates_df = function(X_data, dates) {

    if(is.null(dates)) {
      return(dates)
    }

    time_index  <- 1:(dim(X_data)[1])

    if(!lubridate::is.Date(dates) ) {
      stop("the dates are not invalid")
    }


    if(length(time_index) != length(time_index)) {

      stop("the legend of dates is different from the first dimension of the X_data and Y_data matrices")

    }

    dates_df  <- data.frame(
      "time_index" = time_index,
      "Date" = dates

    )

    private$.dates_df = dates_df


  }
))


