#PLOTS
#Authors: Darlin Soto y Gustavo Soto
#Date: August 20, 2025

#Functions
source("~/lomb_scargle_classic.R")

#Fits
fit1 = readRDS("~fit1.RDS")
fit2 = readRDS("~fit2.RDS")

#Data
data= as.data.frame(read.table('~V445.dat'))
indices_filtro = which(data[,1] < 350 & data[,1] > 185)

t=data[indices_filtro,1]
y=data[indices_filtro,2]

sqrt(mean((y-fit1$yhat)^2))
#0.01229696

sqrt(mean(abs(y-fit2$yhat)^2))
#0.0514613


#Graph light curve
pdf("Data.pdf", width=10, height=5)
layout(matrix(1, nrow=1, ncol=2), widths=c(1,1))
par(mar = c(4.5,4.5,1,1))
plot(t,y,t='l',ylim=c(max(y),min(y)),xlab='Time [BJD-2454833]',ylab=bquote(paste(K[p],'[mag]')))
dev.off()


#YHAT
pdf("fit1_fit2_yhat.pdf", width=13, height=9)

layout(matrix(c(1, 2), nrow=2, ncol=1), widths=c(1,1))
par(mar = c(4.5,4.5,1,1))

minn = min(c(fit1$yhat, y, fit2$yhat))
maxx = max(c(fit1$yhat, y, fit2$yhat))

#Fit1
plot(t, y, col=1, lwd=2, pch=".", ylim=c(maxx,minn-0.3),t='l',xlab='Time [BJD-2454833]',
     ylab=bquote(paste(K[p],'[mag]')))
lines(t, fit1$yhat, col=2,lwd=1)
legend("topright",
       legend = c("Observations", "Fit"),
       pch = c(".", "."),
       col = c(1, 2),
       lwd = 2)
text(min(t)+(max(t)-min(t))/2,minn-0.2,cex=1.3,
     expression(paste('Fit obtained with time-varying parameters (Motta et al.',' 2022)')))

#Fit2
plot(t, y, col=1, lwd=2, pch=".", ylim=c(maxx,minn-0.3),t='l',xlab='Time [BJD-2454833]',ylab=bquote(paste(K[p],'[mag]')))
lines(t, fit2$yhat, col=2,lwd=1)
legend("topright",
       legend = c("Observations", "Fit"),
       pch = c(".", "."),
       col = c(1, 2),
       lwd = 2)
text(min(t)+(max(t)-min(t))/2,minn-0.2,cex=1.3,
     expression(paste('Fit obtained with time-invariant parameters (Guggenberger et al.',' 2012)')))

dev.off()


#RESIDUALS

pdf("fit1_fit2_residuals.pdf", width=13, height=9)
layout(matrix(c(1, 2), nrow=2, ncol=1), widths=c(1,1))
par(mar = c(5,5.5,2.5,2))
minn=-0.45
maxx=0.23

plot(t,y-fit1$yhat,"l",ylab=bquote(paste(Residuals,' (',K[p],'[mag])')),ylim=c(maxx,minn),xlab='Time [BJD-2454833]')
text(min(t)+(max(t)-min(t))/2,minn+0.08,cex=1.3,
     expression(paste('Residuals of fit obtained with time-varying parameters (Motta et al.',' 2022)')))

plot(t,y-fit2$yhat,"l",ylab=bquote(paste(Residuals,' (',K[p],'[mag])')),ylim=c(maxx,minn),xlab='Time [BJD-2454833]')
text(min(t)+(max(t)-min(t))/2,minn+0.08,cex=1.3,
     expression(paste('Residuals of fit obtained with time-invariant parameters (Guggenberger et al.',' 2012)')))
dev.off()


#LS PERIODOGRAM

t=data[indices_filtro,1]
Tspan <- max(t) - min(t)
dtmin <- min(diff(sort(unique(t))))
fmin <- 1 / (Tspan * 10)
fmax <- 0.5 / dtmin

freqs <- seq(from = fmin, to = fmax, length.out = length(t))

res_fit1=y-fit1$yhat
ls_res_fit1 <- lomb_scargle_classic(t, res_fit1, freqs)
res_fit2=y-fit2$yhat
ls_res_fit2 <- lomb_scargle_classic(t, res_fit2, freqs)


pdf("fit1_fit2_lsper.pdf", width=13, height=9)

layout(matrix(c(1, 2), nrow=2, ncol=1), widths=c(1,1))
par(mar = c(4.5,4.5,1,1))


plot(freqs,ls_res_fit1$power,t='l',ylim=c(0,1.7),xlab='Frequency',ylab='Lomb-Scargle periodogram')
text(min(freqs)+(max(freqs)-min(freqs))/2,1.6,cex=1.3,
     expression(paste('LS periodogram of the residual obtained with time-varying parameters (Motta et al.',' 2022)')))

plot(freqs,ls_res_fit2$power,t='l',ylim=c(0,1.7),xlab='Frequency',ylab='Lomb-Scargle periodogram')
text(min(freqs)+(max(freqs)-min(freqs))/2,1.6,cex=1.3,
     expression(paste('LS periodogram of the residual obtained with time-invariant parameters (Guggenberger et al.',' 2012)')))

dev.off()

freqs[which.max(ls_res_fit2$power)]
abline(v=3*1.95517,col=2)


#TV PARAMETERS

pdf("fit1_fit2_tvpar.pdf",width=24,height=20)

layout(matrix(c(rep(1,3),2:19),7, 3, byrow = TRUE),widths=c(1,1), heights=c(1,1))
par(mar = c(5.5,5.5,2.5,2))


minn=min(c(fit1$trend,fit2$trend))
maxx=max(c(fit1$trend,fit2$trend))

plot(t,fit1$trend,col=1,ylim=c(minn,maxx),ylab=bquote(paste(hat(m),' and ',hat(u))),
     xlab='Time [BJD-2454833]',t='l',cex.lab = 1.5)
lines(t,fit2$trend,col=2)


for (i in 1:9) {
  minimo=min(c(fit1$g1hat[[i]],fit2$g1hat[[i]]))
  maximo=max(c(fit1$g1hat[[i]],fit2$g1hat[[i]]))
  
  plot(t,fit1$g1hat[[i]],t="l",col=1,ylim=c(maximo,minimo),ylab="",
       xlab='Time [BJD-2454833]',cex.lab = 1.5)
  lines(t,fit2$g1hat[[i]],col=2)
  assay <- as.character(i)
  title(ylab=bquote(paste(hat(g)[1*','*~.(assay)],' and ',hat(h)[1*','*~.(assay)])), 
        line=2.5, font=2,cex.lab = 1.5)
}

for (i in 1:9) {
  minimo=min(c(fit1$g2hat[[i]],fit2$g2hat[[i]]))
  maximo=max(c(fit1$g2hat[[i]],fit2$g2hat[[i]]))
  
  plot(t,fit1$g2hat[[i]],t="l",col=1,ylim=c(maximo,minimo),ylab="",
       xlab='Time [BJD-2454833]',cex.lab = 1.5)
  lines(t,fit2$g2hat[[i]],col=2)
  assay <- as.character(i)
  title(ylab=bquote(paste(hat(g)[2*','*~.(assay)],' and ',hat(h)[2*','*~.(assay)])), 
        line=2.5, font=2,cex.lab = 1.5)
}

dev.off()




