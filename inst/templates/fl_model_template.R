
# A working linear-model template for user-defined models.
# Copy this file, then rename fl_template, predict.fl_template, and the class
# "fl_template" together (for example, fl_my_model and "fl_my_model").
# Adapt the fitting and prediction code for your model, then source the file
# into the global environment and use model.outcome = ~my_model().
# X: predictor column names; Y: response column name; fl_call: model call.
fl_template <- function(
    data, X = c(""), Y = " ", fl_call  # do not change these inputs
    ) {

  # prepare for dataX
  dataX <- data[, X,drop = FALSE]

  # Fit the model
  # Replace these lines with the arguments needed by your model
  fl_call[[1]] <- quote(stats::lm)
  fl_call$formula <- stats::reformulate(X, response = Y)
  fl_call$data <- quote(data)


  # run the model
  model <- eval(fl_call)

  # Output
  out <- list(model = model,dataX=dataX) # pred = model$fitted.values

  # the same class
  class(out) <- "fl_template"
  return(out)
}


predict.fl_template<- function(object, newdata, ...) {
  # Extract the fitted model object from the input
  model <- object$model

  # Return one conditional mean per row (P(Y = 1) for a binary response).
  # Adapt the prediction call for your model; ... is accepted for S3 compatibility.
  pred <- stats::predict(model, newdata = newdata)

  return(pred)
}
