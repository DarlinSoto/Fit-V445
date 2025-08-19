dar_error_guggenberger=function(x,t,y){
  
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
  
  
  #Model
  AM_B=aA_B[1]+aA_B[2]*sin(2*pi*f_B*t+phiA_B)
  AM_S=aA_S[1]+aA_S[2]*sin(2*pi*f_S*t+phiA_S)
  FM_B=aF_B[1]*sin(2*pi*1*f_B*t+phiF_B[1])+aF_B[2]*sin(2*pi*2*f_B*t +phiF_B[2])+aF_B[3]*sin(2*pi*3*f_B*t+phiF_B[3])
  C=a0+a[1]*sin(2*pi*1*f_0*t+phi[1]+1*FM_B)+a[2]*sin(2*pi*2*f_0*t+phi[2]+2*FM_B)+
    a[3]*sin(2*pi*3*f_0*t+phi[3]+3*FM_B)+a[4]*sin(2*pi*4*f_0*t+phi[4]+4*FM_B)+
    a[5]*sin(2*pi*5*f_0*t+phi[5]+5*FM_B)+a[6]*sin(2*pi*6*f_0*t+phi[6]+6*FM_B)+
    a[7]*sin(2*pi*7*f_0*t+phi[7]+7*FM_B)+a[8]*sin(2*pi*8*f_0*t+phi[8]+8*FM_B)+
    a[9]*sin(2*pi*9*f_0*t+phi[9]+9*FM_B)
  mu=AM_B*AM_S*C

  return(sum(abs(mu-y)^2))
}