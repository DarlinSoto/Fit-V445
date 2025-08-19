dar_tvpar_guggenberger=function(x,t){
  
  aA_B=x[1:2]
  phiA_B=x[3]
  f_B=x[4]
  aF_B=x[5:7]
  phiF_B=x[8:10]
  
  #secondary modulation and parameters
  aA_S=x[11:12]
  phiA_S=x[13]
  f_S=x[14]
  
  #Main frequency and parameters
  a0=x[15]
  a=x[16:24]
  phi=x[25:33]
  f_0=1.94903
  
  p1=(aA_B[1]+aA_B[2]*sin(2*pi*f_B*t+phiA_B))*(aA_S[1]+aA_S[2]*sin(2*pi*f_S*t+phiA_S))
  
  FM_B=0
  for (i in 1:3) {FM_B=FM_B+aF_B[i]*sin(2*pi*i*f_B*t+phiF_B[i])
  
  }
  C=a0
  for (i in 1:9) {C=C+a[i]*sin(2*pi*i*f_0*t+phi[i]+i*FM_B)
  }
  
  k=9
  g1=list()
  mufor=p1*a0
  for (i in 1:k) {
    g1[[i]]=(a[i]*sin(phi[i]+i*FM_B))*p1
    mufor=mufor+g1[[i]]*cos(2*pi*i*f_0*t)
  }
  
  g2=list()
  for (i in 1:k) {
    g2[[i]]=(a[i]*cos(phi[i]+i*FM_B))*p1
    mufor=mufor+g2[[i]]*sin(2*pi*i*f_0*t)
  }
  
  return(list(yhat=mufor,g1hat=g1,g2hat=g2,trend=p1*a0))
}