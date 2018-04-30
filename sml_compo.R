####################################################################################
## Supervised machine learning
####################################################################################

# we permute the farm to be tested, so that we train a model on all the others farms and test performance. 

sml_compo <- function(otu_table, comp, index) {
  
  ## we have to assign everything that ranger will use to the global environment (data.frame, intermediate object like mod...)
  ## even if something is created just before passing it to the function, it won't be found by ranger 
  ## maybe test creating an environment only for "inside" the sml_compo function
  assign("otu_table", otu_table, envir = .GlobalEnv)
  assign("comp", comp, envir = .GlobalEnv)
  assign("index", index, envir = .GlobalEnv)

  ## library
  require(ranger)
  require(irr)
  require(vegan)
  
  OTU <- data.frame(otu_table) #data.frame(otu_aqua_train)
  COMP <- comp # comp_aqua_train

  ## test of right length between data and compo
  if (dim(OTU)[[1]] != length(COMP[,index])) print("No same size between datasets, check it out...")
  
  classification <- F
  
  # get the list of farms to be permuted
  farms <- unique(COMP$Locality)
  
  # prepare the output
  combined1_rf <- c()
  combined2_rf <- c()
  
  
  farm_nam <- c()
  increm <- 0
  
  for (i in farms)
  {
    # subsetting
    assign("nam_t",paste("test_farm_m", i,sep="_"), envir = .GlobalEnv)
    assign("nam_comp_t", paste("test_farm_c_m", i,sep="_"), envir = .GlobalEnv)
    assign(nam_t, subset(OTU, COMP$Locality ==i), envir = .GlobalEnv)
    assign(nam_comp_t, subset(COMP, COMP$Locality ==i), envir = .GlobalEnv)
    
    assign("nam", paste("train_farm_m", i,sep="_"), envir = .GlobalEnv)
    assign("nam_comp", paste("train_farm_c_m", i,sep="_"), envir = .GlobalEnv)
    assign(nam, subset(OTU, COMP$Locality !=i), envir = .GlobalEnv)
    assign(nam_comp, subset(COMP, COMP$Locality !=i), envir = .GlobalEnv)
    
    # ## now the fitting, random forest with Ranger package with default mtry (1/3) for regression according to Breiman
    # similar here the mod objet has to be in the global env. to be fetchable afterwards
    assign("mod", paste("RF_farm_m", i,sep="_"), envir = .GlobalEnv)
    # to control the random process
    set.seed(1)
    assign(mod, ranger(get(nam_comp)[,index] ~ ., data=get(nam), mtry=floor(dim(get(nam))[2]/3), classification = classification, num.trees = 300, importance= "impurity", write.forest = T),envir = .GlobalEnv)
    
    ## prediction for new data
    predict_tr_rf <- predict(get(mod), get(nam_t))
    combined1_rf <- c(combined1_rf,predict_tr_rf$predictions)
    combined2_rf <- c(combined2_rf,get(nam_comp_t)[,index])
    
    # and keep track af farms to samples correspondance
    farm_nam <- c(farm_nam, rep(x = i, dim(get(nam_t))[1]))
    cat("\n", i, ":done - \n\n")
  }
  # output a vector of predictions
  return(combined1_rf)

}


