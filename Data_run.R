
source("Covariates_time.R")
library(survival)

sepkoski00 = readRDS("sepkoski00.rds")



#   ____________________________________________________________________________
#   Create new database with time-dependent covariates                       ####


new_database = function(sepkoski00){
  x = sepkoski00
  x = subset(x, x$grupoid < 8 & x$longevidad >= 1)
  x$longevidad = as.integer(x$longevidad)
  d = sum(x$longevidad)
  
  ##new matrix with time-dependent covariates
  t_dep = matrix(, ncol = 13, nrow = d)
  
  ##Event
  t_dep[1:dim(t_dep)[1], 4] = rep(0, dim(t_dep)[1])
  ii = length(x$longevidad) - 1
  ini = 1
  for(i in 1:ii){
    fin = ini + x$longevidad[i]
    lon = x$longevidad[i] + 1
    
    ##id
    t_dep[ini:fin, 1] = rep(i, lon)
    
    ##genera
    t_dep[ini:fin, 10] = rep(x$grupoid[i],lon)
    
    ##event
    if(x$observed[i] == TRUE) {t_dep[fin-1,4] = 1}
    
    ##First appearance to later increase the longevity and see the covariates.
    t_dep[ini:fin,8]=rep(x$FA[i],lon)  
    
    ##Update
    ini=fin
  }
  t_dep[ini:dim(t_dep)[1], 10] = rep(x$grupoid[length(x$longevidad)], x$longevidad[length(x$longevidad)])
  t_dep[ini:dim(t_dep)[1], 8] = rep(x$FA[length(x$longevidad)], x$longevidad[length(x$longevidad)])
  t_dep[ini:dim(t_dep)[1], 1] = rep(length(x$longevidad), x$longevidad[length(x$longevidad)])
  if(x$observed[length(x$longevidad)] == TRUE) {t_dep[dim(t_dep)[1], 4] = 1}
  
  ##intervals for each genera
  time0 = c()
  time1 = c()
  ind1 = 1
  iii = ii + 1
  for(j in 1:iii){
    ind = ind1 + x$longevidad[j]
    long = x$longevidad[j]
    long1 = x$longevidad[j] + 1
    time0[ind1:ind] = c(0:long)
    time1[ind1:ind] = c(1:long1)
    ind1 = ind
  }
  t_dep[1:dim(t_dep)[1], 2] = time0[1:dim(t_dep)[1]]
  t_dep[1:dim(t_dep)[1], 3] = time1[1:dim(t_dep)[1]]
  
  ##value of the covariates in each life interval by gender
  t_dep[, 8] = t_dep[, 8] - t_dep[, 3]
  t_dep[, 5] = nivelmar(t_dep[, 8])
  t_dep[, 6] = sinusoidal_1(t_dep[, 8])
  t_dep[, 7] = sinusoidal_2(t_dep[, 8])
  t_dep[, 9] = evolucion(t_dep[, 8])
  t_dep[, 11] = SST(t_dep[, 8])
  t_dep[, 12] = CO2(t_dep[, 8])
  t_dep[, 13] = SSTdiff(t_dep[, 8])
  t_dep=data.frame(t_dep)
  colnames(t_dep)=c("id","t_ini","t_fin","event","Sea_level","sinusoid62","sinusoid140","T","fitness","g_id","sst","CO2","SSTdiff")
  return(t_dep)
}


#   ____________________________________________________________________________
#   Estimates, standard errors and hazard multipliers per unit covariate for####

simple_cox = function(t_dep){
  return(coxph(Surv(t_dep$t_ini, t_dep$t_fin, t_dep$event) ~ t_dep$Sea_level + t_dep$sinusoid62 + t_dep$sinusoid140 + t_dep$`fitness`+ t_dep$SSTdiff +relevel(as.factor(t_dep$g_id),ref = "7"), cluster = t_dep$id))
  }


cox_complex = function(t_dep, model_cuts){
  cut_Cox_mar = survSplit(Surv(t_dep$t_ini, t_dep$t_fin, t_dep$event) ~ t_dep$Sea_level + t_dep$sinusoid62 + t_dep$sinusoid140 + t_dep$`fitness`+t_dep$SSTdiff + relevel(as.factor(t_dep$g_id),ref = "7"), data=data.frame(t_dep), cut=model_cuts[[1]],episode= "tgroup", id="id")
  cut_Cox_sin_2 = survSplit(Surv(t_dep$t_ini, t_dep$t_fin, t_dep$event) ~ t_dep$Sea_level + t_dep$sinusoid62 + t_dep$sinusoid140 + t_dep$`fitness`+t_dep$SSTdiff + relevel(as.factor(t_dep$g_id),ref = "7"), data=data.frame(t_dep), cut=model_cuts[[3]],episode= "tgroup9", id="id")
  cut_Cox_evo = survSplit(Surv(t_dep$t_ini, t_dep$t_fin, t_dep$event) ~ t_dep$Sea_level + t_dep$sinusoid62 + t_dep$sinusoid140 + t_dep$`fitness` +t_dep$SSTdiff+ relevel(as.factor(t_dep$g_id),ref = "7"), data=data.frame(t_dep), cut=model_cuts[[4]],episode= "tgroup2", id="id")
  cut_Cox_grupos = survSplit(Surv(t_dep$t_ini, t_dep$t_fin, t_dep$event) ~ t_dep$Sea_level + t_dep$sinusoid62 + t_dep$sinusoid140 + t_dep$`fitness`+t_dep$SSTdiff + relevel(as.factor(t_dep$g_id),ref = "7"), data=data.frame(t_dep), cut=model_cuts[[6]],episode= "tgroup11", id="id")
  cut_Cox_f=cbind(cut_Cox_mar,cut_Cox_sin_2[,11],cut_Cox_evo[,11], cut_Cox_grupos[,11])
  cox=coxph(Surv(t_dep$t_ini, t_dep$t_fin, t_dep$event) ~ t_dep$Sea_level:strata(tgroup) + t_dep$sinusoid62  + t_dep$sinusoid140:strata(cut_Cox_sin_2[,11])+ t_dep$`fitness`:strata(cut_Cox_evo[,11]) +t_dep$SSTdiff + as.factor(t_dep$g_id):strata(cut_Cox_grupos[,11]) , data=cut_Cox_f)
  return(list(cox,model_cuts,cox$coefficients))
}


list_of_coefficients = function(cox,nomb_coef){
  coef = list()
  for (i in 1:length(nomb_coef)) {
    a = which(grepl(nomb_coef[i], names(cox[[1]]$coefficients), ignore.case = T) == TRUE)
    coef[[i]]=cox[[1]]$coefficients[a]
  }
  names(coef)= nomb_coef
  return(coef)
}


base = new_database(sepkoski00)
model_cuts = list(c(9,26,72), c(), c(6,17,42), c(6,17,42,72),c(),c(5,11,20,30,42))
names(model_cuts)= c("Sea_level","sinusoid62","sinusoid140", "fitness","SSTdiff","grupos")


execute_simple_cox = function(base, model_cuts){
  
  cox_sin_estrato = simple_cox(base)
  # Proportional hazards test
  temp = cox.zph(cox_sin_estrato)
  print(temp)
  plot(temp)
}



execute_cox_complex = function(){
  cox = cox_complex(base,model_cuts)
  nomb_coef = c("Sea_level","sinusoid62","sinusoid140", "fitness","SSTdiff","grupos")
  coeficientes = list_of_coefficients(cox,nomb_coef)
  
  for (i in 1:length(coeficientes)) {
    j=length(coeficientes[[i]])
    names(coeficientes[[i]])=rep(names(model_cuts)[i],j)
  }
  
  for (i in 1:6) {
    names(coeficientes[[6]])[(i*7-6):(i*7)]=c("Trilobita","Brachiopoda","Gastropoda","Cephalopoda","Anthozoa","Bivalvia","Otros")
  }
  
  # Proportional hazards test
  temp_1 = cox.zph(cox[[1]])
  
  print(temp_1)
  return(cox)
 # return(list(cox,model_cuts,coeficientes))
}
c=execute_cox_complex()



#   ____________________________________________________________________________
#   Residuals analysis                                                      ####

# Deviance
plot(resid(c[[1]],type="deviance",collapse = base$id),xlab="indice",ylab="residuos(tipo desvio)", main="Residuos (tipo deviance)" )
summary(resid(c[[1]],type="deviance",collapse = base$id))

# Martingale
plot(resid(c[[1]],type="martingale",collapse = base$id),xlab="indice",ylab="residuos(tipo martingala)", main="Residuos (tipo martingala)" )

# Cox-Snell
di=x$observed
mt=resid(c[[1]],type="martingale",collapse = base$id)
cs_res=di-mt
haz_res=survfit(Surv(cs_res,di)~1,type = "fleming-harrington")
ch=haz_res$cumhaz
plot(haz_res$time,ch,xlab = "Cox-Snell residuals", ylab = "Cumulative hazard of residual", main = "Cox-Snell residuals")
lines(c(0,12),c(0,12))

#saveRDS(model_cuts,"model_cuts.rds")
#saveRDS(coeficientes,"coeficientes.rds")
