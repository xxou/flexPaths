###########################################################
#       inner model for single treatment effect
###########################################################



##################################
#   refer 0 decomposition
##################################
#' @import stats purrr
#' @import SuperLearner dbarts
#' @importFrom dplyr mutate
#'
refer0_effect <- function(
    data, A = A, Y = Y, cov_x = cov_x, M.list = M.list,
    CI_level = 0.95, scale = "diff", estimation =estimation,
    Pmodel.lists = Pmodel.lists, Omodel.lists = Omodel.lists, Imodel.lists = Imodel.lists,
    Imodel.source = Imodel.source
){
  cl <- match.call()

  #### 1. data preparation ---------
  # A is a binary;
  I = data[[A]]
  # counterfactual data
  A0_data <- A1_data <- data;
  A0_data[[A]] <- 0; A1_data[[A]] <- 1

  ### cov_x and mediators
  M.list <- reorder_list(M.list)
  cum_mediators <- Reduce(c, c(list(cov_x),M.list), accumulate = TRUE)
  K <- length(cum_mediators)    # have K -1 mediators


  #### 2. propensity model ----------
  ## corresponding model
  if(estimation != "G"){
   # propensity model for IPW or EIF
  pro.details<- purrr::map(Pmodel.lists,extract_model_details)
  pro.fit_names <- purrr::map(pro.details, ~ paste0("fl_", .x$fit_name))
  pro.calls <- purrr::map(pro.details, "model_call")

  pro.modelfits<- purrr::pmap(list(cum_mediators, pro.fit_names, pro.calls),
                       function(cum_mediators, pro.fit_names, pro.calls){
                         model_fun <- match.fun(pro.fit_names)
                         model<-model_fun(data, X = cum_mediators,Y = A, fl_call=pro.calls)
                         return(model)
                       })

  # P(A=1|cov_X); P(A=1|M1,cov_X); P(A=1|M1,...,MK,X)
  gA1_Mk_X.list <- purrr::map2(cum_mediators, pro.modelfits, ~ {predict(.y, newdata = A1_data[,.x, drop=FALSE])})
  # P(A=0|cov_X); P(A=0|M1,cov_X); P(A=0|M1,...,MK,X)
  gA0_Mk_X.list <- purrr::map(gA1_Mk_X.list, ~ {1-.x})
  }

  #### 3. outcome model modelO.list ----------
  if(estimation != "IPW"){
  # outcome regression for G computation or EIF

  out.details<- purrr::map(Omodel.lists,extract_model_details)
  out.fit_names <- purrr::map(out.details, ~ paste0("fl_", .x$fit_name))
  out.calls <- purrr::map(out.details, "model_call")

  ## inner layer outcome regression
  out.modelfits<- purrr::pmap(list(cum_mediators, out.fit_names, out.calls),
                       function(cum_mediators, out.fit_names, out.calls){
                         model_fun <- match.fun(out.fit_names)
                         model<-model_fun(data, X = c(A,cum_mediators),Y = Y, fl_call=out.calls)
                         return(model)
                       })

  EY_A0_Mk_X.list <-purrr::map2(out.modelfits,cum_mediators, ~ {
    pred=predict(.x, A0_data[,c(A,.y)])
    return(pred)
  })
  # E[Y|A=a,X]:
  EY_A1_X <- predict(out.modelfits[[1]], A1_data[,c(A,cum_mediators[[1]])])
  # E[Y|A=a',X]:
  EY_A0_X <- EY_A0_Mk_X.list[[1]]

  ## integrate_mk
  # K-1 iter outcome models
  iter.details<- purrr::map(Imodel.lists,extract_model_details)
  iter.fit_names <- purrr::map(iter.details, ~ paste0("fl_", .x$fit_name))
  iter.calls <- purrr::map(iter.details, "model_call")
  if(Imodel.source == "default"){iter.calls <-purrr::map(iter.calls, replace_family)}


  integrate_mk_1 <- purrr::pmap(list(cum_mediators[-K], iter.fit_names, iter.calls, EY_A0_Mk_X.list[-1]),
                         function(M, fname, icalls, mu){
                           data$mu<- mu

                           model_fun <- match.fun(fname)
                           model<-model_fun(data, X = c(A,M),Y = "mu", fl_call=icalls)

                           pred <- predict(model,  A1_data[,c(A,M)])
                           return(pred)
                         })

  integrate_mk_2 <- purrr::pmap(list(iter.fit_names[1], iter.calls[1],integrate_mk_1[-1]),
                         function(fname, icalls, mu){
                           data$mu<- mu

                           model_fun <- match.fun(fname)
                           model<-model_fun(data, X = c(A,cum_mediators[[1]]),Y = "mu", fl_call=icalls)
                           # print(model)

                           pred <- predict(model, A0_data[,c(A,cum_mediators[[1]])])
                           return(pred)
                         })

  # fit model for direct A->Y; Y(a,Mk(a'))
  EY_A1_MK_X <- predict(out.modelfits[[K]], A1_data[,c(A,cum_mediators[[K]])])

  E_EY_A1_MK_X_given_A0X <- unlist( # E(E(Y|M1,...,Mk, a,X)|a',X)
    purrr::map(list(EY_A1_MK_X), ~ {
      data$mu<- .x; icalls = iter.calls[[1]]; fname = iter.fit_names[[1]]

      model_fun <- match.fun(fname)
      model<-model_fun(data, X = c(A,cum_mediators[[1]]),Y = "mu", fl_call=icalls)
      # print(model)
      pred <- predict(model,  A0_data[,c(A,cum_mediators[[1]])])
      return(pred)
    }))


  }

  #### 4. potential outcome ---------
  if(estimation=="EIF"){
    #---------- (1). A-> M1 ->...-> Y
    Y_via_M1 <-
      ((1-I)/gA0_Mk_X.list[[1]])*(gA1_Mk_X.list[[2]]/gA0_Mk_X.list[[2]])*(gA0_Mk_X.list[[1]]/gA1_Mk_X.list[[1]])*(data[,Y]- EY_A0_Mk_X.list[[2]])+
      (I/gA1_Mk_X.list[[1]])*(EY_A0_Mk_X.list[[2]]-integrate_mk_1[[1]])+
      integrate_mk_1[[1]]

    # ----------- (2). A-> Mi -> ... Y; i=2,3,...,K-1
    Y_via_Mi =data.frame(Y_via_M1)
    if(K>=3){
      for (num in 2:(K-1)) {
        Y_via_Mi[,num] <-
          ((1-I)/gA0_Mk_X.list[[1]])*(gA1_Mk_X.list[[num+1]]/gA0_Mk_X.list[[num+1]])*(gA0_Mk_X.list[[num]]/gA1_Mk_X.list[[num]])*(data[,Y]- EY_A0_Mk_X.list[[num+1]])+
          (I/gA1_Mk_X.list[[1]])*(gA0_Mk_X.list[[num]]/gA1_Mk_X.list[[num]])*(gA1_Mk_X.list[[1]]/gA0_Mk_X.list[[1]])*(EY_A0_Mk_X.list[[num+1]]-integrate_mk_1[[num]])+
          ((1-I)/gA0_Mk_X.list[[1]])*(integrate_mk_1[[num]]-integrate_mk_2[[num-1]])+
          integrate_mk_2[[num-1]] }
    }

    # ----------- (3). A -> Y; Mk, k=K
    Y_direct <-
      (I/gA1_Mk_X.list[[1]])*(gA0_Mk_X.list[[K]]/gA1_Mk_X.list[[K]])*(gA1_Mk_X.list[[1]]/gA0_Mk_X.list[[1]])*(data[,Y]- EY_A1_MK_X)+
      ((1-I)/gA0_Mk_X.list[[1]])*(EY_A1_MK_X-E_EY_A1_MK_X_given_A0X)+
      E_EY_A1_MK_X_given_A0X


    #------------ (4). EY(a) and EY(a')
    # E[Y(A=1)] & E[Y(A=0)]
    Y_A1 =  (I/gA1_Mk_X.list[[1]])*(data[,Y]-EY_A1_X)+EY_A1_X
    Y_A0 =  ((1-I)/gA0_Mk_X.list[[1]])*(data[,Y]-EY_A0_X)+EY_A0_X
  }

  if(estimation == "IPW"){
    #---------- (1). A-> M1 ->...-> Y
    Y_via_M1 <- ((1-I)/gA0_Mk_X.list[[1]])*(gA1_Mk_X.list[[2]]/gA0_Mk_X.list[[2]])*(gA0_Mk_X.list[[1]]/gA1_Mk_X.list[[1]])*data[,Y]

    #------------ (2). A-> Mi -> ... Y; i=2,3,...,K-1
    Y_via_Mi =data.frame(Y_via_M1)
    if(K>=3){
      for (num in 2:(K-1)) {
        Y_via_Mi[,num] <-
          ((1-I)/gA0_Mk_X.list[[1]])*(gA1_Mk_X.list[[num+1]]/gA0_Mk_X.list[[num+1]])*(gA0_Mk_X.list[[num]]/gA1_Mk_X.list[[num]])*data[,Y]
      }
    }
    # ----------- (3). A -> Y; Mk, k=K
    Y_direct <- (I/gA1_Mk_X.list[[1]])*(gA0_Mk_X.list[[K]]/gA1_Mk_X.list[[K]])*(gA1_Mk_X.list[[1]]/gA0_Mk_X.list[[1]])*data[,Y]

    #------------ (4). EY(a) and EY(a')
    # E[Y(A=1)] & E[Y(A=0)]
    Y_A1 =  (I/gA1_Mk_X.list[[1]])*data[,Y]
    Y_A0 =  ((1-I)/gA0_Mk_X.list[[1]])*data[,Y]
  }

  if(estimation=="G"){
    #---------- (1). A-> M1 ->...-> Y
    Y_via_M1 <- integrate_mk_1[[1]]

    # ----------- (2). A-> Mi -> ... Y; i=2,3,...,K-1
    Y_via_Mi =data.frame(Y_via_M1)
    if(K>=3){
      for (num in 2:(K-1)) {
        Y_via_Mi[,num] <- integrate_mk_2[[num-1]] }
    }
    # ----------- (3). A -> Y; Mk, k=K
    Y_direct <- E_EY_A1_MK_X_given_A0X
    #------------ (4). EY(a) and EY(a')
    # E[Y(A=1)] & E[Y(A=0)]
    Y_A1 =  EY_A1_X
    Y_A0 =  EY_A0_X
  }

  ####  5. calculate effects --------------------------
  Y_via_i<-cbind(Y_via_Mi,Y_direct,Y_A1,Y_A0) # head(Y_via_i)
  names(Y_via_i)[1:(K-1)] <- paste0("Y_M", 1:(K - 1))

  if(scale == "diff"){
    out <- data.frame(Path = c(paste0("A->M", 1:(K-1), "->...->Y"), "A->Y", "total effect: A->...->Y")) %>%
      mutate(
        Effect = map_dbl(1:(K+1), ~ mean(Y_via_i[, .x] - Y_via_i[,(K+2)], na.rm = TRUE)),
        SE = map_dbl(1:(K+1), ~ sqrt(var(Y_via_i[, .x] - Y_via_i[,(K+2)], na.rm = TRUE) / nrow(data))),
        CI.lower = Effect + qnorm((1-CI_level)/2)*SE,
        CI.upper = Effect - qnorm((1-CI_level)/2)*SE,
        P.value = round(2 * (1 - pnorm(abs(Effect / SE))), 4)
      )
  }

  if(scale == "risk"){
    out <- data.frame(Path = c(paste0("A->M", 1:(K-1), "->...->Y"), "A->Y", "total effect: A->...->Y")) %>%
      mutate(
        Effect = map_dbl(1:(K+1), ~ mean(Y_via_i[, .x], na.rm = TRUE)/mean(Y_via_i[,(K+2)],na.rm = TRUE )),
      )
  }

  if(scale == "oddsratio"){
    out <- data.frame(Path = c(paste0("A->M", 1:(K-1), "->...->Y"), "A->Y", "total effect: A->...->Y")) %>%
      mutate(
        Effect = map_dbl(1:(K+1), ~ (mean(Y_via_i[, .x],na.rm = TRUE)/(1-mean(Y_via_i[, .x],na.rm = TRUE)))/
                           (mean(Y_via_i[,(K+2)],na.rm = TRUE )/(1-mean(Y_via_i[,(K+2)],na.rm = TRUE)))
        ),
      )
  }

  output <- list(call = cl,
                 decomposition = "0 reference decomposition",
                 results = out,data = data,potential_data = Y_via_i, A = A, Y = Y, cov_x = cov_x, M.list = M.list,
                 estimation =estimation,scale = scale, CI_level=CI_level,
                 Pmodel.lists = Pmodel.lists, Omodel.lists = Omodel.lists, Imodel.lists = Imodel.lists,
                 Imodel.source = Imodel.source)
  class(output) <- "pathsEffect"
  return(output)
}



##################################
#   sequential decomposition
##################################

#' @import stats purrr dplyr
#' @import SuperLearner dbarts
sequential_effect <-  function(
    data, A = A, Y = Y, cov_x = cov_x, M.list = M.list,
    CI_level = 0.95, scale = "diff", estimation =estimation,
    Pmodel.lists = Pmodel.lists, Omodel.lists = Omodel.lists, Imodel.lists = Imodel.lists,
    Imodel.source = Imodel.source
){

  cl <- match.call()
  #### 1. data preparation ---------
  # A is a binary;
  I = data[[A]]
  # counterfactual data
  A0_data <- A1_data <- data;
  A0_data[[A]] <- 0; A1_data[[A]] <- 1

  ### cov_x and mediators
  cum_mediators <- Reduce(c, c(list(cov_x),M.list), accumulate = TRUE)
  K <- length(cum_mediators)    # have K -1 mediators


  #### 2. propensity model ----------
  ## corresponding model
  if(estimation !="G"){
    # propensity model for IPW or EIF
    pro.details<- purrr::map(Pmodel.lists,extract_model_details)
    pro.fit_names <- purrr::map(pro.details, ~ paste0("fl_", .x$fit_name))
    pro.calls <- purrr::map(pro.details, "model_call")

    pro.modelfits<- purrr::pmap(list(cum_mediators, pro.fit_names, pro.calls),
                         function(cum_mediators, pro.fit_names, pro.calls){
                           model_fun <- match.fun(pro.fit_names)
                           model<-model_fun(data, X = cum_mediators,Y = A, fl_call = pro.calls)
                          # print(model)
                           return(model)
                         })

    # P(A=1|cov_X); P(A=1|M1,cov_X); P(A=1|M1,...,MK,X)
    gA1_Mk_X.list <- purrr::map2(cum_mediators, pro.modelfits, ~ {predict(.y, newdata = A1_data[,.x, drop=FALSE])})
    # P(A=0|cov_X); P(A=0|M1,cov_X); P(A=0|M1,...,MK,X)
    gA0_Mk_X.list <- purrr::map(gA1_Mk_X.list, ~ {1-.x})
  }

  #### 3. outcome model modelO.list ----------
  # outcome model

  if(estimation !="IPW"){
  out.details<- purrr::map(Omodel.lists,extract_model_details)
  out.fit_names <- purrr::map(out.details, ~ paste0("fl_", .x$fit_name))
  out.calls <- purrr::map(out.details, "model_call")

  # fit model for Y vs predictors
  out.modelfits<- purrr::pmap(list(cum_mediators, out.fit_names, out.calls),
                       function(cum_mediators, out.fit_names, out.calls){
                         model_fun <- match.fun(out.fit_names)
                         model<-model_fun(data, X = c(A,cum_mediators),Y = Y, fl_call=out.calls)
                         # print(model)
                         return(model)
                       })


  EY_A1_Mk_X.list <-purrr::map2(out.modelfits,cum_mediators, ~ {
    pred=predict(.x, A1_data[,c(A,.y)])
    return(pred)
  })

  # E[Y|A=a,X]:
  EY_A1_X <- EY_A1_Mk_X.list[[1]]
  # E[Y|A=a',X]:
  EY_A0_X <- predict(out.modelfits[[1]], A0_data[,c(A,cum_mediators[[1]])])

  ## integrate_m
  # K-1 iter outcome models
  iter.details<- purrr::map(Imodel.lists,extract_model_details)
  iter.fit_names <- purrr::map(iter.details, ~ paste0("fl_", .x$fit_name))
  iter.calls <- purrr::map(iter.details, "model_call")
  if(Imodel.source == "default"){iter.calls <-purrr::map(iter.calls, ~ replace_family(.x))}

  # E(E(Y|M1, a,X)|a',X); ...; E(E(Y|M1,..., Mk, a,X)|a',X)
  integrate_m <-purrr::pmap(list(iter.fit_names[1], iter.calls[1], EY_A1_Mk_X.list[-1]),
                     function(fname, icalls, mu){
                       data$mu<- mu

                       model_fun <- match.fun(fname)
                       model<-model_fun(data, X = c(A,cum_mediators[[1]]),Y = "mu", fl_call=icalls)
                      # print(model)

                       pred <- predict(model,  A0_data[,c(A,cum_mediators[[1]])])
                       return(pred)}
  )


}
  #### 4. potential outcome ---------------------------
  if(estimation=="EIF"){
  # Y(a,Mk(a'))
  Y_A1_Mk_A0.list <- purrr::pmap(list(gA1_Mk_X.list[-1], gA0_Mk_X.list[-1], EY_A1_Mk_X.list[-1], integrate_m),
                          function(gA1_Mk_X,gA0_Mk_X,bn_Mk_A1_X,integrate_m){
                            (I/gA1_Mk_X.list[[1]])*(gA0_Mk_X/gA1_Mk_X)*(gA1_Mk_X.list[[1]]/gA0_Mk_X.list[[1]])*(data[,Y]- bn_Mk_A1_X)+
                              ((1-I)/gA0_Mk_X.list[[1]])*(bn_Mk_A1_X-integrate_m)+
                              integrate_m
                          })


  # Y(A=1) & Y(A=0)
  Y_A1 =  (I/gA1_Mk_X.list[[1]])*(data[,Y]-EY_A1_X)+EY_A1_X
  Y_A0 =  ((1-I)/gA0_Mk_X.list[[1]])*(data[,Y]-EY_A0_X)+EY_A0_X
  }

  if(estimation == "IPW"){
    Y_A1_Mk_A0.list <- map2(gA1_Mk_X.list[-1], gA0_Mk_X.list[-1],~{
      (I/gA1_Mk_X.list[[1]])*(.y/.x)*(gA1_Mk_X.list[[1]]/gA0_Mk_X.list[[1]])*data[,Y]
    })

    # Y(A=1) & Y(A=0)
    Y_A1 =  (I/gA1_Mk_X.list[[1]])*data[,Y]
    Y_A0 =  ((1-I)/gA0_Mk_X.list[[1]])*data[,Y]
  }

  if(estimation=="G"){
    # Y(a,Mk(a'))
    Y_A1_Mk_A0.list <- integrate_m

    # Y(A=1) & Y(A=0)
    Y_A1 = EY_A1_X
    Y_A0 = EY_A0_X
  }
  ####  5. calculate effects --------------------------
  Y_A1_Mk_A0 <-do.call(cbind, Y_A1_Mk_A0.list)
  potential_outcome <-as.data.frame(cbind(Y_A1, Y_A1_Mk_A0, Y_A0))
  names(potential_outcome)[2:K] <- paste0("Y_M", 1:(K - 1)) # head(potential_outcome)

  if(scale == "diff"){
  out <- data.frame(Path = rep(NA,K+1),Effect = rep(NA,K+1), SE= rep(NA,K+1))
  for (i in 1:K ) {
    effect = potential_outcome[,i]-potential_outcome[,i+1]
    out$Effect[i] = mean(effect,na.rm=T)
    out$SE[i] = sqrt(var(effect,na.rm=T)/nrow(data))
  }
  out$Effect[K+1] = mean(Y_A1 - Y_A0,na.rm=T)
  out$SE[K+1] = sqrt(var(Y_A1 - Y_A0,na.rm=T)/nrow(data))

  out<-out %>% mutate(
    Path = c(paste0("A->M", seq_len(K-1), "->...->Y"), "A->Y", "total effect: A->...->Y"),
    CI.lower = Effect + qnorm((1-CI_level)/2)*SE,
    CI.upper = Effect - qnorm((1-CI_level)/2)*SE,
    P.value = round(2 * (1 - pnorm(abs(Effect / SE))), 4)
  )
  }

  # E(Y(a))/E(Y(a'))
  if(scale == "risk"){
    out <- data.frame(Path = rep(NA,K+1), Effect = rep(NA,K+1))
    for (i in 1:K ) {
      out$Effect[i] = mean(potential_outcome[,i],na.rm=T)/mean(potential_outcome[,i+1],na.rm=T)
    }
    out$Effect[K+1] = mean(Y_A1,na.rm=T)/ mean(Y_A0,na.rm=T)
    out$Path = c(paste0("A->M", seq_len(K-1), "->...->Y"), "A->Y", "total effect: A->...->Y")
  }

  if(scale == "oddsratio"){
    out <- data.frame(Path = rep(NA,K+1), Effect = rep(NA,K+1))
    for (i in 1:K ) {
      out$Effect[i] = (mean(potential_outcome[,i],na.rm=T)/(1-mean(potential_outcome[,i],na.rm=T)))/
                      (mean(potential_outcome[,i+1],na.rm=T)/1-mean(potential_outcome[,i+1],na.rm=T))
    }
    out$Effect[K+1] = (mean(Y_A1,na.rm=T)/(1-mean(Y_A1,na.rm=T)))/
                      (mean(Y_A0,na.rm=T)/(1-mean(Y_A0,na.rm=T)))
    out$Path = c(paste0("A->M", seq_len(K-1), "->...->Y"), "A->Y", "total effect: A->...->Y")
  }

  # output
  output <- list(call = cl,
                 decomposition = "sequential decomposition",
                 results = out,data = data,potential_data = potential_outcome, A = A, Y = Y, cov_x = cov_x, M.list = M.list,
                 estimation =estimation,scale = scale, CI_level=CI_level,
                 Pmodel.lists = Pmodel.lists, Omodel.lists = Omodel.lists, Imodel.lists = Imodel.lists,
                 Imodel.source = Imodel.source)
  class(output) <- "pathsEffect"
  return(output)
}








