
# template for the modeling function
# The standard function input for data,
# X: vectors of variable names
# Y: outcome variable
# args: for
fl_template <- function(
    data, X = c(""), Y = " ", fl_call  # do not change these inputs
    ) {

  # prepare for dataX
  dataX <- data[, X,drop = FALSE]

  # Fit the model
  # Write the needed call for your model
  # fl_call$data =  quote(data)
  # fl_call$needed_elments = .....


  # run the model
  model <- eval(fl_call)

  # Output
  out <- list(model = model,dataX=dataX) # pred = model$fitted.values

  # the same class
  class(out) <- "fl_template"
  return(out)
}


predict.fl_template<- function(object, newdata) {
  # Extract the fitted model object from the input
  model <- object$model

  # Predict using the fitted model
  pred <- predict(...)

  return(pred)
}
