#   ____________________________________________________________________________
#   Define and prepare covariates as functions of T                         ####
require(tidyverse)

plotflag = FALSE
# Turn on if individual plots required

seaHarrison = read.csv("sealevelcurveHarrison.csv")

##  ............................................................................
##  Interpolated sea level                                                  ####

maxT = max(seaHarrison$T)
nivelmar = approxfun(seaHarrison$T, seaHarrison$SL)
if (plotflag == TRUE) {
  plot(0:maxT, nivelmar(0:maxT), type = "l", main = "sea level")
  points(seaHarrison)
}

##  ............................................................................
##  CO2                                                                     ####

CO2table = read.csv("CO2.csv")
maxT = max(CO2table$T)
CO2 = approxfun(CO2table$T, CO2table$CO2)
if (plotflag == TRUE) {
  plot(0:maxT, CO2(0:maxT), type = "l", main = "CO2")
  points(CO2table)
}

##  ............................................................................
##  Sea Surface Temperature SST                                             ####

SSTtable = read.csv("SST.csv")
maxT = max(SSTtable$T)
SST = approxfun(SSTtable$T, SSTtable$SST)
if (plotflag == TRUE) {
  plot(0:maxT, SST(0:maxT), type = "l", main = "Sea surface temperature")
  points(SSTtable)
}

##  ............................................................................
##  Differential of Sea Surface Temperature                                 ####

SSTdiff = approxfun((maxT - 1):0, diff(SST(maxT:0)))
if (plotflag == TRUE) {
  plot(1:maxT, SSTdiff(1:maxT), type = "l",
       main = "Diff Sea surface temperature")
}

##  ............................................................................
##  Sine wave 62 mda                                                        ####

sinusoidal_1 = function(T) {
  20 * sin((2 * pi / 62) * T - 20)
}

##  ............................................................................
##  Sine wave 142 mda                                                       ####

sinusoidal_2 = function(T) {
  20 * sin((2 * pi / 142) * T + 40)
}

##  ............................................................................
##  Fitness function                                                        ####

evolucion = function(T) {
  (T + 250) ^ (-1)
}
if (plotflag == 1) {
  plot(0:maxT, evolucion(0:maxT), type = "l", main = "Fitness")
}

##  ............................................................................
##  Prepare matrix and function for covariates                              ####

COVARIABLES = cbind(
  nivelmar = nivelmar(1:536),
  sin1 = sinusoidal_1(1:536),
  sin2 = sinusoidal_2(1:536),
  evol = evolucion(1:536),
  # co2 = CO2((1:536)),
  sst = SST(1:536),
  sstdiff = SSTdiff(1:536)
)

covariables_temporales = function(T) {
  return(COVARIABLES[T,])
}

remove(plotflag)
remove(seaHarrison, maxT)
remove(SSTtable)
remove(CO2table)
