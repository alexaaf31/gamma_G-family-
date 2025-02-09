######################################GFW BFGS##########################################################
require(gsl)
# Melhorando a função RGFW com verificações adicionais
RZBMW <- function(n, a, alpha, beta, lambda) {
  u <- runif(n, 0, 1)
  A1 <- qgamma(1 - u, a, 1)
  A2 <- (lambda / beta) * ((A1 / alpha)^(1 / beta))
  Q1 <- (beta / lambda) * lambert_W0(A2)
  return(Q1)
}

# Exemplo de uso
n <- 1000
a <- 2.0
alpha  <- 0.2
beta   <- 0.2
lambda <- 0.6
y <- RZBMW(n, a, alpha, beta, lambda)
hist(y)

# Define the f(x) function
f <- function(x) {
  fdp = (1/gamma(a))*alpha*x^(beta - 1)*(beta + lambda*x)*(alpha*x^(beta)*exp(lambda*x))^(a-1)*exp(lambda*x - alpha*x^(beta)*exp(lambda*x))
  return(fdp)
}

# Plot the PDF
hist(y, probability = TRUE, border = "black", col = "white", main = "")
curve(f(x), add = TRUE, col = 2, lwd = 2)
legend("topright", lwd = c(1), lty = c(1),
       c("Histogram", "ZB-MW PDF"),
       col = c("black", "red"),
       bty = "n")

#fun??o para encontrar o estimador de maxima revossimilha?a via monte carlo
kama.mc = function(nrep=1000, nobs=500, semente=2011, a = 1.9, alpha=0.4, beta=0.8,
                   lambda = 2.5)
{
  
  tempo.inicio = Sys.time()
  logLikk = function(theta){
    a = theta[1]; alpha = theta[2]; beta = theta[3]; lambda = theta[4]; nobs = length(y)
    f = (1/gamma(a))*alpha*y^(beta - 1)*(beta + lambda*y)*(alpha*y^(beta)*exp(lambda*y))^(a-1)*exp(lambda*y - alpha*y^(beta)*exp(lambda*y))
    loglik = log(f)
    minus.logL <- sum(-loglik, na.rm = F)
    return(minus.logL)
  }
  emva      = rep(0, nrep)
  emvalpha  = rep(0, nrep)
  emvbeta   = rep(0,nrep)
  emvlambda = rep(0,nrep)
  set.seed(semente)
  contadorFalhas = 0
  # laco de Monte Carlo
  i = 0
  while(i < nrep){
    y = RZBMW(nobs, a, alpha, beta, lambda)
    ir =(optim(c(a, alpha, beta, lambda),logLikk, method="Nelder-Mead"))
    if(ir$convergence == 0){
      i = i + 1
      emva[i]      = ir$par[1]
      emvalpha[i]  = ir$par[2]
      emvbeta[i]   = ir$par[3]
      emvlambda[i] = ir$par[4]
      
    }
    else{
      contadorFalhas = contadorFalhas + 1
    }
  } 
  # fim do laco de Monte Carlo
  amedio = mean(emva)
  alphamedio = mean(emvalpha)
  betamedio = mean(emvbeta)
  lambdamedio = mean(emvlambda)
  
  avies      = amedio      - a
  alphavies  = alphamedio  - alpha
  betavies   = betamedio   - beta
  lambdavies = lambdamedio - lambda
  
  tempo.fim  = Sys.time()
  tempo.exec = tempo.fim - tempo.inicio
  aEQM       = var(emva)      + (avies)^2
  alphaEQM   = var(emvalpha)  + (alphavies)^2
  betaEQM    = var(emvbeta)   + (betavies)^2
  lambdaEQM  = var(emvlambda) + (lambdavies)^2
  
  resultado = data.frame(nobs=nobs, nrep=nrep, semente=semente,a=a, alpha=alpha,
                         beta=beta,lambda=lambda,amedio = amedio, alphamedio = alphamedio, betamedio=betamedio,
                         lambdamedio=lambdamedio,avies = avies,alphavies=alphavies, betavies=betavies,lambdavies=lambdavies,
                         aEQM= aEQM,alphaEQM=alphaEQM,betaEQM = betaEQM,lambdaEQM=lambdaEQM,
                         falhas=contadorFalhas, tempo=tempo.exec)
  return(resultado)
}
kama.mc(nobs=50)
kama.mc(nobs=100)
kama.mc(nobs=300) 
kama.mc(nobs=500)


################################################################################################################################
####################################################RB-BXII#####################################################################
################################################################################################################################

# Melhorando a função RGFW com verificações adicionais
RRBXII <- function(n, a, c, d, s) {
  u <- runif(n, 0, 1)
  A1 <- qgamma(u, a, 1)
  Q1 <- s*((1 - (1/exp(A1)))^(-(1/d)) - 1)^(1/c)
  return(Q1)
}

# Exemplo de uso
n <- 1000
a <- 7.5
c <- 0.2
d <- 0.8
s <- 0.6
y <- RRBXII(n, a, c, d, s)
hist(y)

# Define the f(x) function
f <- function(x) {
  fdp = ((c*d*x^(c-1))/(s^c * gamma(a)))*(1 + (x/s)^c)^(-(d+1))*(-log(1 - (1 + (x/s)^c)^(-d)))^(a-1)
  return(fdp)
}

# Plot the PDF
hist(y, probability = TRUE, border = "black", col = "white", main = "")
curve(f(x), add = TRUE, col = 2, lwd = 2)
legend("topright", lwd = c(1), lty = c(1),
       c("Histogram", "ISW PDF"),
       col = c("black", "red"),
       bty = "n")

#fun??o para encontrar o estimador de maxima revossimilha?a via monte carlo
kama.mc = function(nrep=1000, nobs=500, semente=2011, a = 1.5, c=3.2, d=0.8,
                   s = 0.1)
{
  
  tempo.inicio = Sys.time()
  logLikk = function(theta){
    a = theta[1]; c = theta[2]; d = theta[3]; s = theta[4]; nobs = length(y)
    f = ((c*d*y^(c-1))/(s^c * gamma(a)))*(1 + (y/s)^c)^(-(d+1))*(-log(1 - (1 + (y/s)^c)^(-d)))^(a-1)
    loglik = log(f)
    minus.logL <- sum(-loglik, na.rm = F)
    return(minus.logL)
  }
  emva = rep(0, nrep)
  emvc = rep(0, nrep)
  emvd = rep(0,nrep)
  emvs = rep(0,nrep)
  set.seed(semente)
  contadorFalhas = 0
  # laco de Monte Carlo
  i = 0
  while(i < nrep){
    y = RRBXII(nobs, a, c, d, s)
    ir =(optim(c(a,c,d,s),logLikk, method="Nelder-Mead"))
    if(ir$convergence == 0){
      i = i + 1
      emva[i]  = ir$par[1]
      emvc[i]  = ir$par[2]
      emvd[i]  = ir$par[3]
      emvs[i]  = ir$par[4]
      
    }
    else{
      contadorFalhas = contadorFalhas + 1
    }
  } 
  # fim do laco de Monte Carlo
  amedio = mean(emva)
  cmedio = mean(emvc)
  dmedio = mean(emvd)
  smedio = mean(emvs)
  
  avies = amedio - a
  cvies = cmedio - c
  dvies = dmedio - d
  svies = smedio - s
  
  tempo.fim = Sys.time()
  tempo.exec = tempo.fim - tempo.inicio
  aEQM = var(emva) + (avies)^2
  cEQM = var(emvc) + (cvies)^2
  dEQM = var(emvd) + (dvies)^2
  sEQM = var(emvs) + (svies)^2
  
  resultado = data.frame(nobs=nobs, nrep=nrep, semente=semente,a=a, c=c,
                         d=d,s=s,amedio = amedio, cmedio = cmedio, dmedio=dmedio,
                         smedio=smedio,avies = avies,cvies=cvies, dvies=dvies,svies=svies,
                         aEQM= aEQM,cEQM=cEQM,dEQM = dEQM,sEQM=sEQM,
                         falhas=contadorFalhas, tempo=tempo.exec)
  return(resultado)
}
kama.mc(nobs=50)
kama.mc(nobs=100)
kama.mc(nobs=300) 
kama.mc(nobs=500)
