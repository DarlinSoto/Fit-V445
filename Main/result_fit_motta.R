#FIT SEMI-PARAMETRIC MODEL
#Authors: Darlin Soto y Gustavo Soto
#Date: August 18, 2025

#Functions
source("~ajuste_pspline.R")
source("~bspline.R")

#Data
data= as.data.frame(read.table('V445.dat'))
indices_filtro = which(data[,1] < 350 & data[,1] > 185)

t=data[indices_filtro,1]
y=data[indices_filtro,2]

length(y)
#7325 observation

#Frequencies
f1=(1:9)*1.94903*2*pi

initial_time1 <- proc.time()
fit1=ajuste_pspline(bdeg=3, t=t, knots=20, f=f1, n=length(f1), lam=rep(0,2*length(f1)+1), pord=1, mux=y, boolean_ic=0)
final_time1 <- proc.time()
final_time1-initial_time1
#user   system  elapsed 
#5139.591    6.317 5142.119

saveRDS(fit1,"fit1.RDS")
