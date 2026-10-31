#  SL ------------
#' Built-in Model Wrappers
#'
#' Wrappers for linear models, generalized linear models, BART, and
#' Super Learner. Usually called internally by flexPaths using the model
#' specifications supplied to [pathsInfo()].
#'
#' @param data A wide-format data frame with one row per individual.
#' @param X Predictor column names in `data`.
#' @param Y Response column name in `data`.
#' @param fl_call An unevaluated model call, for example
#'   `quote(glm(family = binomial()))`. All arguments in the call must be named.
#' @param object A fitted object returned by the corresponding wrapper.
#' @param newdata A data frame of predictors for prediction.
#' @param ... Additional arguments accepted for S3 compatibility; currently unused.
#'
#' @details
#' The wrappers supply training data, predictors, and the response automatically.
#' For `lm` and `glm`, an additive formula is used unless a named, two-sided
#' `formula` is supplied; its response is replaced by `Y`.
#' For `bart`, `keeptrees = TRUE` is set for prediction.
#' For `SuperLearner`, supply `SL.library` and an appropriate `family`;
#' optional `id` is a column name included in `X`, and that column is excluded
#' from the training predictors.
#'
#' @return A list with the corresponding `fl_*` class, containing `model`
#'   (the fitted model) and `dataX` (the predictors selected by `X`).
#'   `fl_SuperLearner()` also returns `coef` (the ensemble weights).
#'   Prediction methods return predictions for `newdata`.
#'
#' @seealso [pathsInfo()], [stats::lm()], [stats::glm()], [dbarts::bart()],
#'   [SuperLearner::SuperLearner()]
#' @rdname fl_models
#' @import stats
#' @import SuperLearner
#' @export
fl_SuperLearner <- function(data, X = c(""), Y = " ", fl_call) {
  # get the data
  dataX <- data[, X, drop=FALSE]
  Y_train <- data[,Y]

  # cluster data
  if(!is.null(fl_call$id)){
    X_train = dataX[, names(dataX) != fl_call$id]
    X_id = dataX[, fl_call$id]
    fl_call$id = quote(X_id)
  }else{X_train = dataX}

  # write the full formula
  fl_call$X =quote(X_train); fl_call$Y = quote(Y_train)

  # Fit the model
  model <-eval(fl_call)

  # output
  out = list(model = model, coef=model$coef,dataX=dataX) #pred = model$SL.predict
  class(out) = c("fl_SuperLearner")
  return(out)
}


#' @rdname fl_models
#' @export
predict.fl_SuperLearner <- function(object, newdata, ...) {
  # order of X!!!!!
  # Extract the fitted model object from the input
  model <- object$model
  # Predict using the fitted model
  pred <- predict(model, newdata = newdata, onlySL = T)$pred

  return(pred)
}

# revise: do not need to specify the train.X,
# just input whole data, it will directly get the predictors.name

# fl_bart ------------
#' @rdname fl_models
#' @import stats  dbarts
#' @export
fl_bart <- function(data, X = c(""), Y = " ", fl_call) {
  # Extract training data
  X_train <- data[, X, drop = FALSE]
  Y_train <- data[, Y]

  # Constructing the function call string for BART
  fl_call$x.train = quote(X_train); fl_call$y.train = quote(Y_train)
  fl_call$keeptrees = TRUE
  # Fit the model
  model <- eval(fl_call)

  # Output
  out <- list(model = model,dataX=X_train)
  class(out) <- "fl_bart"
  return(out)
}



#' @rdname fl_models
#' @export
predict.fl_bart <- function(object, newdata, ...) {
  # Extract the fitted model object from the input
  model <- object$model

  # Predict using the fitted model
  pred <- colMeans(predict(model, newdata = newdata))

  return(pred)
}


# glm ------------
#' @rdname fl_models
#' @import stats
#' @export
#'
fl_glm <- function(data, X = c(""), Y = " ", fl_call) {

  # Fit the model
  dataX <- data[, X,drop = FALSE]
  data <- data[,c(X,Y)]

  if ("formula" %in% names(fl_call)) {
    # Extract the current formula

    # Modify the formula by replacing the left-hand side with Y
     fl_call$formula <- substitute(Y ~ predictors, list(Y = as.name(Y), predictors = fl_call$formula[[3]]))

  }else{
    # Replace the formula in the original call
    fl_call$formula <- as.formula(paste0(Y, "~ ."))
  }

  fl_call$data = quote(data)
  model <- eval(fl_call)

  # Output
  out <- list(model = model,dataX=dataX) # pred = model$fitted.values
  class(out) <- "fl_glm"
  return(out)
}



# Define the prediction function for fl_glm
#' @rdname fl_models
#' @export
predict.fl_glm <- function(object, newdata, ...) {
  # Extract the fitted model object from the input
  model <- object$model

  # Predict using the fitted model
  pred <- predict(model, newdata = newdata, type = "response")

  return(pred)
}



# lm ------------
#' @rdname fl_models
#' @import stats
#' @export
#'
fl_lm <- function(data, X = c(""), Y = " ", fl_call) {

  # Fit the model
  dataX <- data[, X,drop = FALSE]
  data <- data[,c(X,Y)]

  if ("formula" %in% names(fl_call)) {
    # Extract the current formula

    # Modify the formula by replacing the left-hand side with Y
    fl_call$formula <- substitute(Y ~ predictors, list(Y = as.name(Y), predictors = fl_call$formula[[3]]))

  }else{
    # Replace the formula in the original call
    fl_call$formula <- as.formula(paste0(Y, "~ ."))
  }

  fl_call$data = quote(data)
  model <- eval(fl_call)

  # Output
  out <- list(model = model,dataX=dataX) # pred = model$fitted.values
  class(out) <- "fl_lm"
  return(out)
}



# Define the prediction function for fl_lm
#' @rdname fl_models
#' @export
predict.fl_lm <- function(object, newdata, ...) {
  # Extract the fitted model object from the input
  model <- object$model

  # Predict using the fitted model
  pred <- predict(model, newdata = newdata)

  return(pred)
}
