
source("Covariates_time.R")
#source("Prepara_Sepkoski.R")
library(survival)
#prepara_Sepkoski()
#sepkoski00 = readRDS("sepkoski00.rds")


crea_base_temporalidad = function(sepkoski00){
  
  x = sepkoski00
  x = subset(x, x$grupoid < 8 & x$longevidad >= 1)
  #x = subset(x, x$longevidad >= 1)
  #dim(x)
  
  ## Convenientemente para moverse en la tabla
  x$longevidad = as.integer(x$longevidad)
  
  d = sum(x$longevidad)
  
  ##Crea matriz transformada para tomar en cuenta covariables dependientes del tiempo
  t_dep = matrix(, ncol = 20, nrow = d)
  
  ##Eventos
  t_dep[1:dim(t_dep)[1], 4] = rep(0, dim(t_dep)[1])
  ii = length(x$longevidad) - 1
  ini = 1
  
  
  for(i in 1:ii){
    
    fin = ini + x$longevidad[i]
    lon = x$longevidad[i] + 1
    
    ##id
    t_dep[ini:fin, 1] = rep(i, lon)
    
    ##genero
    t_dep[ini:fin, 10] = rep(x$grupoid[i],lon)
    
    ##eventos
    if(x$observed[i] == TRUE) {t_dep[fin-1,4] = 1}
    
    ##Primera aparición para después ir aumentando la longevidad y ver las covariables.
    t_dep[ini:fin,8]=rep(x$FA[i],lon)  
    
    ##Actualizo
    ini=fin
  }
  
  
  t_dep[ini:dim(t_dep)[1], 10] = rep(x$grupoid[length(x$longevidad)], x$longevidad[length(x$longevidad)])
  t_dep[ini:dim(t_dep)[1], 8] = rep(x$FA[length(x$longevidad)], x$longevidad[length(x$longevidad)])
  t_dep[ini:dim(t_dep)[1], 1] = rep(length(x$longevidad), x$longevidad[length(x$longevidad)])
  if(x$observed[length(x$longevidad)] == TRUE) {t_dep[dim(t_dep)[1], 4] = 1}
  
  ##intervalos para cada género
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
  
  ##valor de covariables en cada t de cada género
  t_dep[, 8] = t_dep[, 8] - t_dep[, 3]
  t_dep[, 5] = nivelmar(t_dep[, 8])
  t_dep[, 6] = sinusoidal_1(t_dep[, 8])
  t_dep[, 7] = sinusoidal_2(t_dep[, 8])
  t_dep[, 9] = evolucion(t_dep[, 8])
  t_dep[, 18] = SST(t_dep[, 8])
  t_dep[, 19] = CO2(t_dep[, 8])
  t_dep[, 20] = SSTdiff(t_dep[, 8])
  
  ##grupo
  dummies = matrix(0, nrow=dim(t_dep)[1], ncol=length(unique(t_dep[,10])))
  
  dim(dummies)
  for(i in 1:dim(t_dep)[1]){
    dummies[i, t_dep[,10][i]] = 1
  }
  
  t_dep[,11]=dummies[,1]
  t_dep[,12]=dummies[,2]
  t_dep[,13]=dummies[,3]
  t_dep[,14]=dummies[,4]
  t_dep[,15]=dummies[,5]
  t_dep[,16]=dummies[,6]
  t_dep[,17]=dummies[,7]
  
  t_dep=data.frame(t_dep)
  colnames(t_dep)=c("id","t_ini","t_fin","event","nivelmar","sinusoidal1","sinusoidal2","T","fitness","g_id","trilobita","brachiopoda","gastropoda","cephalopoda","anthozoa","bivalvia","otros","sst","CO2","SSTdiff")
  return(t_dep)
}

ejecuta_cox_sin_estrato = function(t_dep){
  return(coxph(Surv(t_dep$t_ini, t_dep$t_fin, t_dep$event) ~ t_dep$nivelmar + t_dep$sinusoidal1 + t_dep$sinusoidal2 + t_dep$`fitness`+ t_dep$SSTdiff +relevel(as.factor(t_dep$g_id),ref = "7"), cluster = t_dep$id))
  #return(coxph(Surv(t_dep$t_ini, t_dep$t_fin, t_dep$event) ~ t_dep$nivelmar + t_dep$sinusoidal1 + t_dep$sinusoidal2 + t_dep$`fitness`+ t_dep$SSTdiff+ t_dep$trilobita+ t_dep$brachiopoda+ t_dep$gastropoda+ t_dep$cephalopoda+ t_dep$anthozoa+ t_dep$bivalvia+ t_dep$otros, cluster = t_dep$id))
  
  }


ejecuta_cox_particiones_coeficientes = function(t_dep, modelo){
  cut_Cox_mar = survSplit(Surv(t_dep$t_ini, t_dep$t_fin, t_dep$event) ~ t_dep$nivelmar + t_dep$sinusoidal1 + t_dep$sinusoidal2 + t_dep$`fitness` + t_dep$trilobita+ t_dep$brachiopoda+ t_dep$gastropoda+ t_dep$cephalopoda+ t_dep$anthozoa, data=data.frame(t_dep), cut=modelo[[1]],episode= "tgroup", id="id")
  cut_Cox_sin_1 = survSplit(Surv(t_dep$t_ini, t_dep$t_fin, t_dep$event) ~ t_dep$nivelmar + t_dep$sinusoidal1 + t_dep$sinusoidal2 + t_dep$`fitness` + t_dep$trilobita+ t_dep$brachiopoda+ t_dep$gastropoda+ t_dep$cephalopoda+ t_dep$anthozoa, data=data.frame(t_dep), cut=modelo[[2]],episode= "tgroup1", id="id")
  cut_Cox_sin_2 = survSplit(Surv(t_dep$t_ini, t_dep$t_fin, t_dep$event) ~ t_dep$nivelmar + t_dep$sinusoidal1 + t_dep$sinusoidal2 + t_dep$`fitness` + t_dep$trilobita+ t_dep$brachiopoda+ t_dep$gastropoda+ t_dep$cephalopoda+ t_dep$anthozoa, data=data.frame(t_dep), cut=modelo[[3]],episode= "tgroup9", id="id")
  
  cut_Cox_evo = survSplit(Surv(t_dep$t_ini, t_dep$t_fin, t_dep$event) ~ t_dep$nivelmar + t_dep$sinusoidal1 + t_dep$sinusoidal2 + t_dep$`fitness` + t_dep$trilobita+ t_dep$brachiopoda+ t_dep$gastropoda+ t_dep$cephalopoda+ t_dep$anthozoa, data=data.frame(t_dep), cut=modelo[[4]],episode= "tgroup2", id="id")
  cut_Cox_g1 = survSplit(Surv(t_dep$t_ini, t_dep$t_fin, t_dep$event) ~ t_dep$nivelmar + t_dep$sinusoidal1 + t_dep$sinusoidal2 + t_dep$`fitness` + t_dep$trilobita+ t_dep$brachiopoda+ t_dep$gastropoda+ t_dep$cephalopoda+ t_dep$anthozoa, data=data.frame(t_dep), cut=modelo[[6]],episode= "tgroup3", id="id")
  cut_Cox_g3 = survSplit(Surv(t_dep$t_ini, t_dep$t_fin, t_dep$event) ~ t_dep$nivelmar + t_dep$sinusoidal1 + t_dep$sinusoidal2 + t_dep$`fitness` + t_dep$trilobita+ t_dep$brachiopoda+ t_dep$gastropoda+ t_dep$cephalopoda+ t_dep$anthozoa, data=data.frame(t_dep), cut=modelo[[8]],episode= "tgroup4", id="id")
  cut_Cox_g4 = survSplit(Surv(t_dep$t_ini, t_dep$t_fin, t_dep$event) ~ t_dep$nivelmar + t_dep$sinusoidal1 + t_dep$sinusoidal2 + t_dep$`fitness` + t_dep$trilobita+ t_dep$brachiopoda+ t_dep$gastropoda+ t_dep$cephalopoda+ t_dep$anthozoa, data=data.frame(t_dep), cut=modelo[[9]],episode= "tgroup5", id="id")
  cut_Cox_g5 = survSplit(Surv(t_dep$t_ini, t_dep$t_fin, t_dep$event) ~ t_dep$nivelmar + t_dep$sinusoidal1 + t_dep$sinusoidal2 + t_dep$`fitness` + t_dep$trilobita+ t_dep$brachiopoda+ t_dep$gastropoda+ t_dep$cephalopoda+ t_dep$anthozoa, data=data.frame(t_dep), cut=modelo[[10]],episode= "tgroup6", id="id")
  cut_Cox_CO2 = survSplit(Surv(t_dep$t_ini, t_dep$t_fin, t_dep$event) ~ t_dep$nivelmar + t_dep$sinusoidal1 + t_dep$sinusoidal2 + t_dep$`fitness` + t_dep$trilobita+ t_dep$brachiopoda+ t_dep$gastropoda+ t_dep$cephalopoda+ t_dep$anthozoa, data=data.frame(t_dep), cut=c(9,23,60),episode= "tgroup7", id="id")
  cut_Cox_sst = survSplit(Surv(t_dep$t_ini, t_dep$t_fin, t_dep$event) ~ t_dep$nivelmar + t_dep$sinusoidal1 + t_dep$sinusoidal2 + t_dep$`fitness` + t_dep$trilobita+ t_dep$brachiopoda+ t_dep$gastropoda+ t_dep$cephalopoda+ t_dep$anthozoa, data=data.frame(t_dep), cut=modelo[[5]],episode= "tgroup8", id="id")
  cut_Cox_SSTdiff = survSplit(Surv(t_dep$t_ini, t_dep$t_fin, t_dep$event) ~ t_dep$nivelmar + t_dep$sinusoidal1 + t_dep$sinusoidal2 + t_dep$`fitness` + t_dep$trilobita+ t_dep$brachiopoda+ t_dep$gastropoda+ t_dep$cephalopoda+ t_dep$anthozoa, data=data.frame(t_dep), cut=modelo[[11]],episode= "tgroup10", id="id")
  cut_Cox_grupos = survSplit(Surv(t_dep$t_ini, t_dep$t_fin, t_dep$event) ~ t_dep$nivelmar + t_dep$sinusoidal1 + t_dep$sinusoidal2 + t_dep$`fitness` + t_dep$trilobita+ t_dep$brachiopoda+ t_dep$gastropoda+ t_dep$cephalopoda+ t_dep$anthozoa, data=data.frame(t_dep), cut=modelo[[12]],episode= "tgroup11", id="id")
  
   ##Quitar sst y quitar cortes de co2
  cut_Cox_f=cbind(cut_Cox_mar,cut_Cox_sin_1[,14],cut_Cox_sin_2[,14],cut_Cox_evo[,14],cut_Cox_g1[,14],cut_Cox_g3[,14],cut_Cox_g4[,14],cut_Cox_g5[,14],cut_Cox_CO2[,14],cut_Cox_sst[,14],cut_Cox_SSTdiff[,14], cut_Cox_grupos[,14])
  #cox=coxph(Surv(t_dep$t_ini, t_dep$t_fin, t_dep$event) ~ t_dep$nivelmar:strata(tgroup) + t_dep$sinusoidal1:strata(cut_Cox_sin_1[,14])  + t_dep$sinusoidal2:strata(cut_Cox_sin_2[,14]) + t_dep$`fitness`:strata(cut_Cox_evo[,14]) +t_dep$`sst`:strata(cut_Cox_sst[,14]) + t_dep$trilobita:strata(cut_Cox_g1[,14])+ t_dep$brachiopoda+ t_dep$gastropoda:strata(cut_Cox_g3[,14])+ t_dep$cephalopoda:strata(cut_Cox_g4[,14])+ t_dep$anthozoa:strata(cut_Cox_g5[,14]) , data=cut_Cox_f)
  cox=coxph(Surv(t_dep$t_ini, t_dep$t_fin, t_dep$event) ~ t_dep$nivelmar:strata(tgroup) + t_dep$sinusoidal1  + t_dep$sinusoidal2:strata(cut_Cox_sin_2[,14])+ t_dep$`fitness`:strata(cut_Cox_evo[,14]) +t_dep$SSTdiff + as.factor(t_dep$g_id):strata(cut_Cox_grupos[,14]) , data=cut_Cox_f)
  
  return(list(cox,modelo,cox$coefficients))
}


lista_de_coefficientes = function(cox,nomb_coef){
  coef = list()
  for (i in 1:length(nomb_coef)) {
    a = which(grepl(nomb_coef[i],   names(cox[[1]]$coefficients), ignore.case = T) == TRUE)
    coef[[i]]=cox[[1]]$coefficients[a]
  }
  names(coef)= nomb_coef
  return(coef)
}



base = base_diff
names(base)=c("id","t_ini","t_fin","event","nivelmar","sinusoidal1","sinusoidal2","T","fitness","g_id","trilobita","brachiopoda","gastropoda","cephalopoda","anthozoa","bivalvia","otros","sst","CO2","SSTdiff")
#base = crea_base_temporalidad(sepkoski00)
saveRDS(base,"base_diff.rds")
cox_sin_estrato = ejecuta_cox_sin_estrato(base)

temp = cox.zph(cox_sin_estrato)
print(temp)
plot(temp)


#modelo = list(c(9,30), c(4,10,37), c(6,15,30,40,60), c(9,23,37,50), c(4,8,15,37,60), c(4,6,11), c(), c(37), c(3,5,7,20), c(30))
#names(modelo)= c("nivelmar","sinusoidal1","sinusoidal2", "fitness","sst","trilobita","brachiopoda","gastropoda","cephalopoda","anthozoa")

#para evaluar
modelo = list(c(9,26,72), c(26), c(6,17,42), c(6,17,42,72),c(4,8,15,37,60), c(4,6,11), c(23), c(37), c(3,5,7,20), c(30),c(6,17),c(5,11,20,30,42))
names(modelo)= c("nivelmar","sinusoidal1","sinusoidal2", "fitness","sst","trilobita","brachiopoda","gastropoda","cephalopoda","anthozoa","SSTdiff","grupos")
#para entregar
#modelo = list(c(9,26,72), c(), c(6,17,42), c(6,17,42,72),c(),c(5,11,20,30,42))
#names(modelo)= c("nivelmar","sinusoidal1","sinusoidal2", "fitness","SSTdiff","grupos")


cox = ejecuta_cox_particiones_coeficientes(base,modelo)
nomb_coef = c("nivelmar","sinusoidal1","sinusoidal2", "fitness","SSTdiff","grupos")
coeficientes = lista_de_coefficientes(cox,nomb_coef)

for (i in 1:length(coeficientes)) {
  print(i)
    j=length(coeficientes[[i]])
    names(coeficientes[[i]])=rep(names(modelo)[i],j)
}


names(coeficientes[[6]])[1:7]=c("trilobita","brachiopoda","gastropoda","cephalopoda","anthozoa","bivalvia","otros")
names(coeficientes[[6]])[8:14]=c("trilobita","brachiopoda","gastropoda","cephalopoda","anthozoa","bivalvia","otros")
names(coeficientes[[6]])[15:21]=c("trilobita","brachiopoda","gastropoda","cephalopoda","anthozoa","bivalvia","otros")
names(coeficientes[[6]])[22:28]=c("trilobita","brachiopoda","gastropoda","cephalopoda","anthozoa","bivalvia","otros")
names(coeficientes[[6]])[29:35]=c("trilobita","brachiopoda","gastropoda","cephalopoda","anthozoa","bivalvia","otros")
names(coeficientes[[6]])[36:42]=c("trilobita","brachiopoda","gastropoda","cephalopoda","anthozoa","bivalvia","otros")


temp_1 = cox.zph(cox[[1]])

print(temp_1)
plot(temp_1)

saveRDS(modelo,"modelo.rds")
saveRDS(coeficientes,"coeficientes.rds")
