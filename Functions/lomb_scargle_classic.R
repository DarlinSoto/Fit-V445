lomb_scargle_classic <- function(t, y, freqs) {
  
  t <- as.numeric(t)
  y <- as.numeric(y)
  freqs <- as.numeric(freqs)
  
  M <- length(freqs)
  P <- numeric(M)
  
  for (m in seq_len(M)) {
    f <- freqs[m]
    w <- 2 * pi * f
    S <- sum(sin(2 * w * t))
    C <- sum(cos(2 * w * t))
    tau <- 0.5 * atan2(S, C) / w
    
    ct <- cos(w * (t - tau))
    st <- sin(w * (t - tau))
    
    yc <- sum(y * ct)
    ys <- sum(y * st)
    cc <- sum(ct^2)
    ss <- sum(st^2)
    
    P[m] <- 0.5 *( (yc^2 / cc) + (ys^2 / ss) )
  }
  
  data.frame(freq = freqs, power = P)
}
