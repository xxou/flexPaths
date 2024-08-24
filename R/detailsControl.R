
# Function to reorder the list
reorder_list <- function(x_list) {
  # Extract the names
  names_list <- names(x_list)

  if(!is.null(names_list)){
    # Separate 'cov_x' and 'M' elements
    cov_x_name <- if("cov_x" %in% names_list) "cov_x" else character(0)
    M_names <- names_list[grepl("^M\\d+$", names_list)]

    # Sort 'M' elements
    sorted_M_names <- sort(M_names)

    # Combine in the desired order
    reordered_names <- c(cov_x_name, sorted_M_names)

    # Reorder the list based on the new order
    reordered_list <- x_list[reordered_names]
  }else{
    reordered_list<- x_list
  }

  return(reordered_list)
}


# Function to extract model details
#' @import stats purrr
#'
extract_model_details <- function(model_formula) {
  # Extract the expression from the formula
  expr <- model_formula[[2]]

  # Extract the function name
  fit_name <- as.character(expr[[1]])

  # Check if all elements in the call have names
  all_have_names <- all(names(expr)[-1] != "")

  if (all_have_names==FALSE) {
    stop("Error: All arguments must be named.")
  }
  # Return the details
  return(list(fit_name = fit_name, model_call = expr))
}





# replace family
replace_family <- function(call) {
  if ("family" %in% names(call)) {
    call$family = quote(gaussian())
  }
  return(call)
}


# function to get index in flex path
# index for mediators to be used
create_index <- function(active) {
  # Define the corresponding position vector p
  p <- seq_along(active)
  # Initialize the index vector
  index <- rep(NA, length(active))
  # Set the last value first
  index[length(active)] <- 0

  # Loop through the active vector backwards
  for (i in (length(active) - 1):1) {
    if (active[i] == active[i + 1]) { index[i] <- index[i + 1]
    } else { index[i] <- p[i]}
  }
  return(index)
}


# multiple
mcreate_index <- function(active_matrix) {
  # Get the number of columns in the matrix
  n_cols <- dim(active_matrix)[2]

  # Initialize the index vector
  index <- rep(NA, n_cols)

  # Set the index for the rightmost column to 0
  index[n_cols] <- 0

  # Loop through the columns from right to left
  active_matrix[is.na(active_matrix)] <- -1
  for (i in (n_cols - 1):1) {
    # Compare the entire column vector with the next column vector
    if (all(active_matrix[, i] == active_matrix[, i + 1])) {
      index[i] <- index[i + 1]  # Keep the index same as the right one
    } else {
      index[i] <- i  # Set index to the column number
    }
  }

  return(index)
}


# conbine active value
simplify_active <- function(active_matrix) {
  active_matrix[is.na(active_matrix)] <- -1
  # Initialize the simplified matrix with the first column
  simplified_m <- active_matrix[, 1, drop = FALSE]

  # Loop through columns 2 to n
  for (i in 2:ncol(active_matrix)) {
    # Check if the current column is different from the last column in the simplified matrix
    if (!all(active_matrix[, i] == simplified_m[, ncol(simplified_m)], na.rm = TRUE)) {
      # If different, add the current column to the simplified matrix
      simplified_m <- cbind(simplified_m, active_matrix[, i])
    }
  }
  simplified_m[simplified_m==-1] <- NA
  return(simplified_m)
}


# combine mediators
combine_item <- function(element_list, index) {
  # Filter out elements with index 0
  valid_indices <- which(index != 0)
  filtered_list <- element_list[valid_indices]
  filtered_index <- index[valid_indices]

  # Initialize the result list
  combined_list <- list()

  # Loop through unique indices and combine elements
  unique_indices <- unique(filtered_index)
  for (i in unique_indices) {
    combined_elements <- unlist(filtered_list[filtered_index == i], use.names = FALSE)
    combined_list[[paste0("M", i)]] <- combined_elements
  }

  return(combined_list)
}

# cumulate variables
cumulate_vars<- function(input, first_null = FALSE){
  if(first_null){
    var.list <-c(list(NULL), Reduce(c, c(input), accumulate = TRUE)[-length(input)])
  }else{
    var.list <- Reduce(c, c(input), accumulate = TRUE)
  }
  return(var.list)
}


# product of predicted value from listed model

# model_list<- pro.modelfits[[4]] the model to be predicted
# a_vars <- a.list[[3]] its value
# a_values <- active_s.list[[3]]
#' @import stats purrr
#'
prod_assign_modellist <- function(model_list,a_vars, a_values){
  model_list = model_list[1:length(a_vars)]
  a_values = a_values[1:length(a_vars)]

  Xa_vars = cumulate_vars(a_vars, first_null = T)
  Xa_values = cumulate_vars(a_values, first_null = T)


  # predict the model list
  pred <- pmap(list(model_list, Xa_vars, Xa_values, a_values),
               function(model, Xa_var, Xa_value, Y_value){

                 Ai_data = model$dataX
                 # print( head(Ai_data))
                 if(!is.null(Xa_var)){
                   Ai_data[,Xa_var] = suppressMessages({map_dfc(Xa_value, ~rep(.x, dim(Ai_data)[1]))})
                 }
                 # print( head(Ai_data)); cat(Xa_var,',',Xa_value,",",Y_value,'\n')

                 predY1 <-predict(model, newdata = Ai_data )
                 pred <- predY1*Y_value + (1-predY1)*(1-Y_value)
               })
  prod.modelpred <- reduce(pred, `*`)
  return(prod.modelpred)
}


# boot strap
# object<-out
#' @import stats purrr
#' @import SuperLearner dbarts
#' @importFrom dplyr mutate
one_boot.pathsEffect <- function(object){
  # get the full data
  list2env(object, envir = environment())

  # replaced sample for boot
  index <- sample(1:dim(data)[1], replace = TRUE)
  boot_data = data[index,]

  # call for boot_data
  object$call$data <- substitute(boot_data)
  boot_sim = eval(object$call)
  return(boot_sim$results[,c('Path','Effect')])
}

# object = call_all
# data = emdata
# fit = pathsFit
#' @import stats purrr
#' @import SuperLearner dbarts
one_boot.flexEffect <- function(data, pathsFit, call_all, index_p1, index_p0){
  # replaced sample for boot
  index <- sample(1:dim(data)[1], replace = TRUE)
  boot_data = data[index,]
  pathsFit$data <- boot_data

  # call for boot_data
  object<-map(call_all, ~ {.x$pathsFit <- substitute(pathsFit); return(.x)})

  boot_sim = map_dfr(object, function(x){
    sim = eval(x)
    return(sim$results[,c('active','value')])
  })

  sim_p1 <- boot_sim[index_p1,]
  sim_p0 <- boot_sim[index_p0,]

  out<- data.frame(active = paste(sim_p1$active, "vs" , sim_p0$active),
                   value1 = sim_p1$value,
                   value0 = sim_p0$value)
  return(out)
}





# check SuperLearner
# input a call
#' @import  purrr
#'
check_SuperLearner <- function(call){
  # Extract the model lists from the call
  model_names <- grep("model", names(call), value = TRUE)
  model_lists <- map(model_names, ~ call[[.x]])
  contains_superlearner <- any(map_lgl(model_lists, ~ "SuperLearner" %in% all.names(.x)))
  return(contains_superlearner)
}


