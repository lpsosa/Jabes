# Auxiliary function cortes(x, places, values) for evaluating step functions
# used in the linear predictor in the Cox regression model

cortes = function(x, places, values){
  # Inputs:
  # places is vector of length k, positive entries
  # values is vector of length k+1 of constants defining the step function
  # x is the vector containing places to evaluate the function
  
  # Outputs:
  # A vector of constants containing step function evaluated at x
  
  temp = cut(x, breaks=c(0, places, Inf), labels=FALSE)
  return(values[temp])
}
