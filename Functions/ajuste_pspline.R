ajuste_pspline <- function(bdeg, t, knots, f, n, lam, pord, mux, boolean_ic) {

  B = bspline(t, min(t) - 0.00001, max(t) + 0.00001, knots, bdeg)
  X = B

  for (i in 1:n) {
      X = cbind(X, B * cos(f[[i]] * t)  )
  }

  for (i in 1:n) {
      X = cbind(X, B * sin( f[[i]] * t ))
  }
  
  nb = ncol(B)
  nX = ncol(X)
  D = diff(diag(nb), diff = pord)
  D=t(D)%*%D
  lambda = diag(lam)
  P = kronecker(sqrt(lambda), D)
  
  X2 = rbind(X, P)
  optimizador=solve(t(X2) %*% X2) %*% t(X2)
  ahat=optimizador%*% c(mux, rep(0, nrow(P)))
  yhat = X %*% ahat
  S=X2%*%optimizador
  
  df = sum(diag(S))
  rss = sum((mux - yhat) ^ 2)
  sigma2hat = rss / (length(t) - df)
  
  deviance = sum((mux - yhat) ^ 2)
  AIC = deviance / length(t) + 2 * df * sigma2hat / length(t)

  g1 = list()
  g2 = list()
  varg1 = list()
  
  for (i in 1:n) {
    g1[[i]] = B %*% ahat[(i * ncol(B) + 1):((i + 1) * ncol(B))]
    xaux = X[, (i * ncol(B) + 1):((i + 1) * ncol(B))]
    paux = B %*% solve(t(xaux) %*% xaux) %*% t(xaux)
    varg1[[i]] = sigma2hat * diag(paux %*% t(paux))
  }
  
  contador = 1
  varg2 = list()
  aux = (n + 1):(2 * n)
  
  for (i1 in aux) {
    g2[[contador]] = B %*% ahat[(i1 * ncol(B) + 1):((i1 + 1) * ncol(B))]
    xaux = X[, (i1 * ncol(B) + 1):((i1 + 1) * ncol(B))]
    paux = B %*% solve(t(xaux) %*% xaux) %*% t(xaux)
    varg2[[contador]] = sigma2hat * diag(paux %*% t(paux))
    contador = contador + 1
  }
  
  trend = B %*% ahat[1:ncol(B)]
  aux = B %*% solve(t(X[, 1:ncol(B)]) %*% X[, 1:ncol(B)]) %*% t(X[, 1:ncol(B)])
  vartrend = sigma2hat * diag(aux %*% t(aux))
  
  varyhat = sigma2hat * diag(S %*% t(S))
  
  return(list(yhat = yhat, trend = trend, g1hat = g1, g2hat = g2, varg1 = varg1, 
              varg2 = varg2, vartrend = vartrend, residuos = mux - yhat, aic = AIC, 
              varyhat = varyhat, coef = ahat))
}


