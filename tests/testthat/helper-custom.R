# Source a user model and restore the previous bindings after each test.
local_custom_model <- function(env = globalenv(), .local_envir = parent.frame()) {
  models <- c("fl_testlm", "predict.fl_testlm")
  present <- vapply(models, exists, logical(1), envir = env, inherits = FALSE)
  old <- mget(models[present], envir = env, inherits = FALSE)
  withr::defer({
    for (name in models) {
      if (name %in% names(old)) assign(name, old[[name]], envir = env) else
        if (exists(name, envir = env, inherits = FALSE)) rm(list = name, envir = env)
    }
  }, envir = .local_envir)
  source(test_path("fixtures", "custom-model.R"), local = env)
}
