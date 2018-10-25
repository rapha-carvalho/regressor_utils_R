library(dplyr)
library(MASS)
library(gtools)
library(imputeTS)


find.best.model <- function(y, x, path, r2=0, n.var){
  
  #Read the data
  read.data <- function (filepath, cols=NULL) {
    df <- read.csv(filepath, sep = ";", stringsAsFactors = FALSE)
    df <- na.replace(df, 0)
    var.names <- colnames(df)
    sel.names <- setdiff(var.names, cols)
    df <- df %>% select(sel.names)
    df
  }

  #Calculate the adjusted R2 of every X variable vs Y variable.
  adj.r2.x <- function(y, x, data, threshold){
    r2adj <- list()
    var.names <- colnames(data)[x]
    y.var <- paste(y, "~")
    for (i in seq_along(var.names)){
      x.var <- var.names[i]
      form <- as.formula(paste(y.var, x.var))
      lm.model <- lm(data = data, formula = form)
      lm.sum_model <- summary(lm.model)
      r2adj[[x.var]] <- lm.sum_model$adj.r.squared 
    }
    unlisted.r2adj <- unlist(r2adj)
    eligible.variables <- unlisted.r2adj[unlisted.r2adj >= threshold]
    eligible.variables <- eligible.variables[sort.list(eligible.variables, decreasing = TRUE)]
    return(eligible.variables)
  }

  #Calculate the R2 for a regression between each X variable vs the remaining X variables.
  vif <- function(data, columns){
    r2adj <- list()
    for (i in columns){
      y.var <- paste(i, "~")
      x.var <- paste0(setdiff(columns, i), collapse = "+")
      form <- as.formula(paste(y.var, x.var))
      lm.model <- lm(data = data, formula = form)
      lm.sum_model <- summary(lm.model)
      r2adj[[i]] <- lm.sum_model$adj.r.squared 
    }
    unlisted.r2adj <- unlist(r2adj)
    rank.var <- unlisted.r2adj[sort.list(unlisted.r2adj)]
    return (rank.var)
  }


  #Select the top N variables.
  var.selection <- function(y, x, data, r2, n.var){
    
    #Get the adjusted R^2 of all variables x against y.
    r2adj.var <- adj.r2.x(y, x, data, r2)
    r2adj.var.names <- names(r2adj.var)
    
    #Calculate the number of x variables to input in the model.
    if (length(r2adj.var) < n.var){
      num.vars <- length(r2adj.var)
    } else {
      num.vars <- n.var
    }
    
    #Get the vif of the selected variables.
    var.vif <- vif(data, r2adj.var.names)  
    
    #Rank the variables 
    var.rank <- (r2adj.var * (1-var.vif))[sort.list((r2adj.var * (1-var.vif)), decreasing = TRUE)]
    
    #Filter the variables by the number of variables wanted.
    var.filtered <- var.rank[1:num.vars]
    
    #Return the filtered variable names.
    return(names(var.filtered))
  }
  
  
  #Return a list of all possible combination of models.
  model.combinations <- function(y, x, data){
    
    #List of all lm models by using 'combn'.
    all.models <- list()
    
    #List of all combinations of models.
    variable.combinations <- list()
    
    print("getting all the combinations and running the lm functions")
    num.models <- 0
    for (i in 1:length(x)) {
      #Get all the formula combinations possible.
      variable.combinations[[i]]<- combn(x, i, simplify = TRUE)
      
      for (j in 1:dim(variable.combinations[[i]])[2]){
        #Run a lm for all the formulas obtained on the previous for loop.
        all.models[[paste0(i, ",", j)]] <- reg.model(y, variable.combinations[[i]][,j], data)
        #Get the number of the current combination.
        num.models <- num.models + 1
        print(paste("running combination number:", num.models))
      }
    }
    print("all lm functions are computed")
    return(all.models)
  }
  
  
  
  #Get the name of the best models per coefficient used.
  get.best.model <- function(models) {
    aic.models <- c()
    
    for (i in 1:length(models)) {
        #Calculate the AIC for each model.
        aic.models[[names(models)[i]]] <- AIC(models[[names(models)[i]]])
    }
    #Find the model with the lowest AIC.
    best <- names(aic.models[which.min(aic.models)])
    
    return(best)
  }
    #Read the data in.
    data <- read.data(path)
    #Get the top N variables.
    variables <- var.selection(y, x, data, r2, n.var)
    #Get all the model combinations.
    all.models <- model.combinations(y, variables, data)
    #Get the model with the lowest AIC.
    best.model <- get.best.model(all.models)
    best.model <- all.models[[best.model]]
    
    return(best.model)
    
    
}
