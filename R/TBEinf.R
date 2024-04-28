#####################################################################################################
# Link function and its derivative
link_inverse <- function(eta) -log(- eta) # Inverse link function
link_derivative <- function(p) exp(p) # Derivative of the inverse link function, m with respect to eta

f <- function(x,y,phi) {
  exp((1-y)*x-x*log(x))*phi^(-x)*sin(pi*x)
}

f_phi <- function(y,phi) {
  n <- length(y)
  f_phi <- numeric(0)
  for (i in 1:n){
    f_phi[i] <- integrate(function(x) f(x,y[i],phi), lower = 0, upper = Inf)$value / pi
  }
  return(f_phi)
} # Landau distribution with phi

df_phi <- function(y,phi) {
  df <- function(x,y,phi) {
    exp((1-y)*x-x*log(x))*phi^(-x-1)*sin(pi*x)
  }
  n <- length(y)
  df_phi <- numeric(0)
  for (i in 1:n){
    df_phi[i] <- integrate(function(x) df(x,y[i],phi), lower = 0, upper = Inf)$value / pi
  }
  return(df_phi)
} # Derivative of Landau distribution


dTBEinf <- function(y, X, beta, phi, method="saddle", m = NULL) {
  if(is.null(m)){
    eta <- X %*% beta
    m <- -log(- eta)
  }
  if ( method == "real" ) { # Real density
      return(exp(exp(-m)*(m-y+1)/phi)/pi*f_phi(y,phi))
    }
    else if (method == "saddle") { # Saddlepoint density 
      n <- length(y)
      return(exp(-(exp(-y)+exp(-m)*(y-m-1))/phi)/sqrt(2*pi*phi^2*exp(m)))
    }
    else if (method == "finverse") { # Density evaluated by using Fourier inversion
      Fi <- function(x, y, m, phi){
        Re(exp( ( (-exp(-m)+x*1i)*(1-log(exp(-m)-x*1i))+exp(-m)*(1+m)-x*y*1i )/phi ))
      }
      n <- length(y)
      inversion <- numeric(0)
      for (i in 1:n){
        inversion[i] <- quadgk(function(x) Fi(x, y[i], m[i], phi), -1, 1, tol = 1e-12) / (2*pi*phi)
      }
      return(inversion)
    }
    else if (method == "mWtrans") { # Density evaluated by the modified W-transformation which can be applied to integrating oscillating functions
      f_phi_mWtrans <- function(y, phi, max_p = 10, max_s = 100){ # Evaluate f_phi by the modified W-transformation
        n <- length(y)
        f_phi_mW <- numeric(0)
        for (i in 1:n){
        f_phi_j <- function(y,phi,j) {
          integrate(function(x) f(x,y,phi), lower = 0, upper = j)$value / pi
        } # Landau distribution with phi
      
        ###### Split the integral into two parts: 0 to 1 and 1 to infinity
        f_phi_1 <- f_phi_j(y[i],phi,1) # Calculate the integral from 0 to 1 first
        t_zp <- 1:max_s # zero points of the integrand (start from 1)
        f_phi_v <- numeric(0)
        for ( k in 1:(max_s+1) ) {
          f_phi_v[k] <- f_phi_j(y[i],phi,k)
        }
        F_t <- f_phi_v[t_zp] - f_phi_1
        w_t <- f_phi_v[2:(max_s+1)]-f_phi_v[t_zp]
        M <- matrix(data = 0, nrow = max_p, ncol = max_s)
        N <- matrix(data = 0, nrow = max_p, ncol = max_s)
        W <- matrix(data = 0, nrow = max_p, ncol = max_s)
        M[1,] <- F_t/w_t # M_{-1}
        N[1,] <- 1/w_t # N_{-1}
      
        for (p in 2:max_p) {
          for (s in 1:(max_s-1)) {
            M[p,s] <- (M[p-1,s]-M[p-1,s+1])/(1/t_zp[s]-1/t_zp[s+p+1])
            N[p,s] <- (N[p-1,s]-N[p-1,s+1])/(1/t_zp[s]-1/t_zp[s+p+1])
            W[p,s] <- M[p,s]/N[p,s]
            if ( s>2 ){
              if ( (abs(W[p,s]-W[p,s-1]) + abs(W[p,s]-W[p,s-2]))/W[p,s]<1e-10 ){
                break
              }
            }
          }
        }
        f_phi_mW[i] <- W[p,s]
      }
      return(f_phi_mW)
    }
  f_phi_mWtrans(y, phi) * exp( exp(-m)*(m+1-y)/phi )
  }
  else {
    stop("Undefined method!")
  }
}

####### Total deviance
TBEinf.Dev <- function(y, X, beta, m = NULL){
  if(is.null(m)){
    eta <- X %*% beta
    m <- -log(- eta)
  }
  2 * sum( exp(-y) + exp(-m) * (y-m-1) )
}

####### Log-likelihood
logLikTBEinf <- function(y, X, beta, phi, m = NULL){
  if(is.null(m)){
    eta <- X %*% beta
    m <- -log(- eta)
  }
  log( sum( dTBEinf(y, y, phi) ) ) - 1/(2*phi) * TBEinf.Dev(y, X, beta)
}

####### Derivative of the log-likelihood with respect to phi
dTBEinf.dldphi <- function(y, X, beta, phi, m = NULL){
  if(is.null(m)){
    eta <- X %*% beta
    m <- -log(- eta)
  }
  (TBEinf.Dev(y, X, beta)-2*sum(exp(-y)))/(2*phi^2)+sum(df_phi(y, phi)/f_phi(y, phi))
}

####### Derivative of the saddlepoint log-likelihood with respect to phi
dTBEinf.dldphi.saddle <- function(y, X, beta, phi, m = NULL){
  if(is.null(m)){
    eta <- X %*% beta
    m <- -log(- eta)
  }
  n <- length(y)
  -n/phi + 1/(2*phi^2) * TBEinf.Dev(y, X, beta)
}

####### Estimate beta in the glm for the TBE_infinity models based on IRLS (Iteratively Reweighted Least Squares)
TBEinf.glm.fit <- function(y, X, beta0 = NULL, tau = 1e-8, max.iter = 100){
  n <- nrow(X)
  pp <- ncol(X)
  if(is.null(beta0)){
    beta0 <- rep(-0.01, pp)
  }
  beta <- beta0
  for (iter in 1:max.iter) {
    eta <- X %*% beta
    m_temp <- link_inverse(eta)
    
    # Check NaN. If there is NaN, break.
    if ( any(is.nan(m_temp)) ) {
      cat("There is NaN in", iter, "iterations\n")
      break
    }
    m <- m_temp
    
    W <- as.vector( exp(m) ) # Working weights
    z <- eta + (y - m) / link_derivative(m) # Working response
    
    # Update step using weighted least squares
    XT <- t(X)
    WX <- XT
    for (j in 1:n ) {
      WX[,j] <- XT[,j] * W[j]
    }
    beta_new <- solve(WX %*% X) %*% (WX %*% z) # eta (z) regress on X
    
    # Check convergence
    if (max(abs(beta_new - beta)) < tau) {
      # cat("Converged in", iter, "iterations\n")
      break
    }
    
    beta <- beta_new
  }
  
  total_deviance <- TBEinf.Dev(y, X, beta)
  phi <- total_deviance / (n - pp)
  
  val = list(beta     = beta_new, 
             iter     = iter, 
             m        = m,
             phi      = phi)
  
  return(val)
}

####### Predict the value of y~TBE_infinity by its mean=link_inverse(X*beta)
TBEinf.glm.pre <- function(X, beta){
  eta_pre <- X %*% beta
  m_pre <- link_inverse(eta_pre)
  y_pre <- m_pre
  return(y_pre)
}



