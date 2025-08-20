#FIT PARAMETRIC MODEL
#Authors: Darlin Soto y Gustavo Soto
#Date: August 20, 2025

#Functions
source("~/error_guggenberger.R")
source("~/tvpar_guggenberger.R")

#Library
library(pso)

#Data
data= as.data.frame(read.table('V445.dat'))
indices_filtro = which(data[,1] < 350 & data[,1] > 185)

t1=data[indices_filtro,1]
y1=data[indices_filtro,2]

#Fit Guggenberger
n = 33
lower_bounds = rep(-2*pi, n)
upper_bounds = rep(2*pi, n)

# Initial values
initial_parameters = runif(n, min = -2*pi, max = 2*pi)

# PSO
initial_time2=proc.time()
fit2=psoptim(fn = error_guggenberger, lower = lower_bounds, upper = upper_bounds,
               control=list(maxit=2200,s=350),par = initial_parameters,t=t1,y=y1)
final_time2 <- proc.time()
final_time2-initial_time2
#    user   system  elapsed 
#1159.249    5.033 1163.707 

fit2_tvpar=tvpar_guggenberger(fit2$par,t1)

saveRDS(fit2_tvpar,"fit2.RDS")


