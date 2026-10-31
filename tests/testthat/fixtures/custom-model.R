# A complete user model, following fl_model_template.R.
fl_testlm <- function(data, X, Y, fl_call) {
  dataX <- data[, X, drop = FALSE]
  fl_call[[1]] <- quote(stats::lm)
  fl_call$formula <- reformulate(X, response = Y)
  fl_call$data <- quote(data)
  out <- list(model = eval(fl_call), dataX = dataX)
  class(out) <- "fl_testlm"
  out
}

predict.fl_testlm <- function(object, newdata, ...) {
  predict(object$model, newdata = newdata)
}
