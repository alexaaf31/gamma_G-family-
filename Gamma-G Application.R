
require(AdequacyModel)
require(zipfR)
require(survival)

t = c(61, 77, 17, 72, 71, 9, 190, 69, 98, 7, 36, 128, 128, 157, 
      128, 36, 127, 5, 126, 34, 122, 121, 213, 152, 152, 39, 39,
      1, 243, 91, 91, 60, 30, 60, 23, 1, 151, 31, 115, 61, 183, 
      122, 26, 30, 46, 108, 136, 136, 8, 23, 23, 112, 6, 6, 13,
      83, 13, 7, 53, 23, 12, 6, 36, 65, 95, 42, 16, 13, 14, 131,
      41, 48, 137, 75, 9, 14, 14, 8, 14, 14, 7, 7, 7, 21, 14, 9,
      14, 14, 18, 67, 14, 14, 14, 14, 14, 19, 13, 20, 13, 45, 14, 14, 14)

hist(t)
ekm=survfit(Surv(t)~1)
########################################################################################
###########################BurrXII######################################################
########################################################################################

pdf_BXII <- function(par,x){
  c = par[1]
  k = par[2]
  s = par[3]
  c*k*s^(-c)*x^(c-1)*(1 + (x/s)^c)^(-(k+1))
}

# bXII acumulada

cdf_BXII <- function(par,x){
  c = par[1]
  k = par[2]
  s = par[3]
  1 - (1 + (x/s)^c)^(-k)
}

# bXII sobrevivencia

sdf_BXII <- function(par,x){
  c = par[1]
  k = par[2]
  s = par[3]
  (1 + (x/s)^c)^(-k)
}

results_BXII <- goodness.fit(pdf = pdf_BXII, cdf = cdf_BXII, 
                             starts = c(1.113493, 4.535670, 262.586574), data = t, 
                             method = "B", mle = NULL, domain = c(0, Inf),lim_inf = c(0,0,0),
                             lim_sup = c(2,2,2), S = 250, prop=0.1, N=500)

results_BXII

####################################################################################################
#################################=====Gamma_BXII-ZB=====################################################
#####################################################################################################
pdf_GBXII = function(par,x){
  a = par[1]
  c = par[2]
  k = par[3]
  s = par[4]
  g = c*k*s^(-c)*x^(c-1)*(1 + (x/s)^c)^(-(k+1))
  G = 1 - (1 + (x/s)^c)^(-k)
  (g/gamma(a))*(-log(1-G))^(a-1)
}

cdf_GBXII = function(par, x){
  a = par[1]
  c = par[2]
  k = par[3]
  s = par[4]
  G = 1 - (1 + (x/s)^c)^(-k)
  ((Igamma(a,-log(1-G)))/gamma(a))
}

sdf_GBXII = function(par, x){
  a = par[1]
  c = par[2]
  k = par[3]
  s = par[4]
  G = 1 - (1 + (x/s)^c)^(-k)
  1 - (((Igamma(a,-log(1-G)))/gamma(a)))
}

results_GBXII <- goodness.fit(pdf = pdf_GBXII, cdf = cdf_GBXII, 
                              starts = c(16.2316952, 0.2792122, 30.8685800, 131.7939758), data = t, 
                              method = "B", domain = c(0, Inf), lim_inf = c(0,0,0,0),
                              lim_sup = c(2,2,2,2), mle = NULL, S = 250, prop=0.1,N=50)

results_GBXII

####################################################################################################
#################################=====Gamma_BXII-RB=====################################################
#####################################################################################################
pdf_GBXII_RB = function(par,x){
  a = par[1]
  c = par[2]
  k = par[3]
  s = par[4]
  g = c*k*s^(-c)*x^(c-1)*(1 + (x/s)^c)^(-(k+1))
  G = 1 - (1 + (x/s)^c)^(-k)
  (g/gamma(a))*(-log(G))^(a-1)
}

cdf_GBXII_RB = function(par, x){
  a = par[1]
  c = par[2]
  k = par[3]
  s = par[4]
  G = 1 - (1 + (x/s)^c)^(-k)
  1 - ((Igamma(a,-log(G)))/gamma(a))
}

sdf_GBXII_RB = function(par, x){
  a = par[1]
  c = par[2]
  k = par[3]
  s = par[4]
  G = 1 - (1 + (x/s)^c)^(-k)
  1 - (1 - ((Igamma(a,-log(G)))/gamma(a)))
}

results_GBXII_RB <- goodness.fit(pdf = pdf_GBXII_RB, cdf = cdf_GBXII_RB, 
                                 starts = c(0.1063074,   0.9846576,  63.3638482, 315.0819893), data = t, 
                                 method = "B", domain = c(0, Inf), lim_inf = c(0,0,0,0),
                                 lim_sup = c(2,2,2,2), mle = NULL, S = 250, prop=0.1,N=50)

results_GBXII_RB

########################################################################################
######################gamma_Weibull_Modificada_ZB###########################################
########################################################################################

pdf_GMW <- function(par,x){
  a = par[1]
  alpha = par[2]
  lambda = par[3]
  theta  = par[4]
  g <- alpha*x^(theta-1)*(theta+lambda*x)*exp(lambda*x - alpha*x^(theta)*exp(lambda*x))
  G <- 1 - exp(-alpha*x^theta*exp(lambda*x))
  (g/gamma(a))*(-log(1-G))^(a-1)
  
}

# bmw acumulada

cdf_GMW <- function(par,x){
  a = par[1]
  alpha = par[2]
  lambda = par[3]
  theta  = par[4]
  G <- 1 - exp(-alpha*x^theta*exp(lambda*x))
  (Igamma(a,-log(1-G)))/gamma(a)
}

# bmw sobrevivencia

sdf_GMW <- function(par,x){
  a = par[1]
  alpha = par[2]
  lambda = par[3]
  theta  = par[4]
  G <- 1 - exp(-alpha*x^theta*exp(lambda*x))
  fd <- (Igamma(a,-log(1-G)))/gamma(a)
  1-fd
}
results_GMW <- goodness.fit(pdf = pdf_GMW, cdf = cdf_GMW, 
                            starts = c(1.28610948, 0.05161952, 0.00099830, 0.79366544), data = t, 
                            method = "B", mle = NULL, domain = c(0, Inf),lim_inf = c(0,0,0,0),
                            lim_sup = c(2,2,2,2), S = 250, prop=0.1, N=50)

results_GMW

########################################################################################
######################gamma_Weibull_Modificada_RB###########################################
########################################################################################

pdf_GMW_RB <- function(par,x){
  a = par[1]
  alpha = par[2]
  lambda = par[3]
  theta  = par[4]
  g <- alpha*x^(theta-1)*(theta+lambda*x) * exp(lambda*x - alpha*x^(theta)*exp(lambda*x))
  G <- 1 - exp(-alpha*x^theta*exp(lambda*x))
  (g/gamma(a))*(-log(G))^(a-1)
  
}

# bmw acumulada

cdf_GMW_RB <- function(par,x){
  a = par[1]
  alpha = par[2]
  lambda = par[3]
  theta  = par[4]
  G <- 1 - exp(-alpha*x^theta*exp(lambda*x))
  1 - ((Igamma(a,-log(G)))/gamma(a))
}

# bmw sobrevivencia

sdf_GMW_RB <- function(par,x){
  a = par[1]
  alpha = par[2]
  lambda = par[3]
  theta  = par[4]
  G <- 1 - exp(-alpha*x^theta*exp(lambda*x))
  fd <- 1 - ((Igamma(a,-log(1-G)))/gamma(a))
  1-fd
}
results_GMW_RB <- goodness.fit(pdf = pdf_GMW_RB, cdf = cdf_GMW_RB, 
                               starts = c(1.0433814330, 0.0186780121, 0.0006464495, 0.9632677587), data = t, 
                               method = "B", mle = NULL, domain = c(0, Inf),lim_inf = c(0,0,0,0),
                               lim_sup = c(2,2,2,2), S = 250, prop=0.1, N=50)

results_GMW_RB
########################################################################################
######################Weibull_Modificada###########################################
########################################################################################

pdf_MW <- function(par,x){
  alpha = par[1]
  lambda = par[2]
  theta  = par[3]
  g <- alpha*x^(theta-1)*(theta+lambda*x) * exp(lambda*x - alpha*x^(theta)*exp(lambda*x))
  
}

# bmw acumulada

cdf_MW <- function(par,x){
  alpha = par[1]
  lambda = par[2]
  theta  = par[3]
  G <- 1 - exp(-alpha*x^theta*exp(lambda*x))
}

# bmw sobrevivencia

sdf_MW <- function(par,x){
  alpha = par[1]
  lambda = par[2]
  theta  = par[3]
  G <- 1 - exp(-alpha*x^theta*exp(lambda*x))
  1-G
}
results_MW <- goodness.fit(pdf = pdf_MW, cdf = cdf_MW, 
                           starts = c(0.01460810, 0.01089828, 0.26776126), data = t, 
                           method = "B", mle = NULL, domain = c(0, Inf),lim_inf = c(0,0,0,0),
                           lim_sup = c(2,2,2,2), S = 250, prop=0.1, N=50)

results_MW

####################################################################################################
################################ ========Log-logistica=========#####################################
####################################################################################################
pdf_LL = function(par,x){
  c = par[1]
  s = par[2]
  g = c*s^(-c)*x^(c-1)*(1 + (x/s)^c)^(-2)
  g
  
}

cdf_LL = function(par, x){
  c = par[1]
  s = par[2]
  G = 1 - (1 + (x/s)^c)^(-1)
  G
}

sdf_LL = function(par, x){
  c = par[1]
  s = par[2]
  G = 1 - (1 + (x/s)^c)^(-1)
  1-G
}

results_LL <- goodness.fit(pdf = pdf_LL, cdf = cdf_LL, 
                           starts = c(1.456074, 32.076042), data = t, 
                           method = "B", domain = c(0, Inf),mle = NULL, lim_inf = c(0,0),
                           lim_sup = c(2,2), S = 250, prop=0.1, N=500)

results_LL

#####################################################################################################
#################################=====G_LL_RB=====######################################################
#####################################################################################################
pdf_GLL = function(par,x){
  a     = par[1]
  alpha = par[2]
  beta  = par[3]
  g = beta*alpha^(-beta)*x^(beta-1)*(1 + (x/alpha)^beta)^(-2)
  G = 1 - (1 + (x/alpha)^beta)^(-1)
  (g/gamma(a))*(-log(G))^(a-1)
}

cdf_GLL = function(par, x){
  a     = par[1]
  alpha = par[2]
  beta  = par[3]
  G = 1 - (1 + (x/alpha)^beta)^(-1)
  1 - ((Igamma(a,-log(G)))/gamma(a))
}

sdf_GLL = function(par, x){
  a     = par[1]
  alpha = par[2]
  beta  = par[3]
  G = 1 - (1 + (x/alpha)^beta)^(-1)
  1 - (1 - ((Igamma(a,-log(G)))/gamma(a)))
}

results_GLL <- goodness.fit(pdf = pdf_GLL, cdf = cdf_GLL, 
                            starts = c(9.174254, 1116.631669, 2.573863), data = t, 
                            method = "B", domain = c(0, Inf), lim_inf = c(0,0,0,0),
                            lim_sup = c(2,2,2,2), mle = NULL, S = 250, prop=0.1,N=50)

results_GLL

#####################################################################################################
#################################=====G_LL_ZB=====######################################################
#####################################################################################################
pdf_GLL_ZB = function(par,x){
  a     = par[1]
  alpha = par[2]
  beta  = par[3]
  g = beta*alpha^(-beta)*x^(beta-1)*(1 + (x/alpha)^beta)^(-2)
  G = 1 - (1 + (x/alpha)^beta)^(-1)
  (g/gamma(a))*(-log(1-G))^(a-1)
}

cdf_GLL_ZB = function(par, x){
  a     = par[1]
  alpha = par[2]
  beta  = par[3]
  G = 1 - (1 + (x/alpha)^beta)^(-1)
  ((Igamma(a,-log(1-G)))/gamma(a))
}

sdf_GLL_ZB = function(par, x){
  a     = par[1]
  alpha = par[2]
  beta  = par[3]
  G = 1 - (1 + (x/alpha)^beta)^(-1)
  1 - (((Igamma(a,-log(1-G)))/gamma(a)))
}

results_GLL_ZB <- goodness.fit(pdf = pdf_GLL_ZB, cdf = cdf_GLL_ZB, 
                               starts = c(1.081941, 28.403620,  1.431716), data = t, 
                               method = "B", domain = c(0, Inf), lim_inf = c(0,0,0,0),
                               lim_sup = c(2,2,2,2), mle = NULL, S = 250, prop=0.1,N=50)

results_GLL_ZB


##grafico do histograma da distribuição

x = seq(0,300, by = 0.01)
hist(t, probability = TRUE, xlab = "t", ylab = "f(t)", main = "",ylim = c(0,0.020))
lines(x, pdf_GBXII(par = results_GBXII$mle, x),col = 1, lwd = 1, lty = 1)
lines(x, pdf_GBXII_RB(par = results_GBXII_RB$mle, x),col = 1, lwd = 1, lty = 2)
lines(x, pdf_LL(par = results_LL$mle, x), col = 1,lwd = 1, lty = 4)
lines(x, pdf_GLL(par = results_GLL$mle, x), col = 1,lwd = 1,lty = 3)
legend("topright",lwd=c(1,1,1,1),lty=c(1,2,3,4),
       c("ZB-GBXII","RB-GBXII","RB-GLL","LL"),
       col=c(1,1,1,1),
       bty="n")

##gráfico da função empírica
x = seq(0,350, by = 0.01)
plot(ekm$time, 1-ekm$surv, xlab = "t", ylab = "F(t)",type = "l", lwd=1, col=1, lty=1)
lines(x, cdf_GBXII(par = results_GBXII$mle, x),col = 1, lwd = 1, lty = 1)
lines(x, cdf_GBXII_RB(par = results_GBXII_RB$mle, x),col = 1, lwd = 1, lty = 2)
lines(x, cdf_GLL(par = results_GLL$mle, x), col = 1, lwd = 1, lty = 3)
lines(x, cdf_LL(par = results_LL$mle, x), col = 1, lwd = 1, lty = 4)
legend("bottomright",lwd=c(1,1,1,1),lty=c(1,2,3,4), yjust = 0.7,
       c("ZB-GBXII","RB-GBXII","RB-GLL","LL"),
       col=c(1,1,1,1),
       bty="n")

# tabela de estimativas e erro padrao

GBXII_estimativas = c(results_GBXII$mle)
GBXII_erro = c(results_GBXII$Erro)
GBXII_RB_estimativas = c(results_GBXII_RB$mle)
GBXII_RB_erro = c(results_GBXII_RB$Erro)
BXII_estimativas = c(results_BXII$mle)
BXII_erro = c(results_BXII$Erro)
GLL_estimativa = results_GLL$mle
GLL_erro = results_GLL$Erro
GLL_ZB_estimativas = c(results_GLL_ZB$mle)
GLL_ZB_erro = c(results_GLL_ZB$Erro)
GMW_ZB_estimativa = c(results_GMW_RB$mle)
GMW_ZB_erro = c(results_GMW_RB$Erro)
GMW_estimativa = c(results_GMW$mle)
GMW_erro = c(results_GMW$Erro)
MW_estimavia = c(results_MW$mle)
MW_erro = c(results_MW$Erro)
LL_estimativas = c(results_LL$mle)
LL_erro = c(results_LL$Erro)


tabela1 =  rbind(GBXII_estimativas, GBXII_erro, BXII_estimativas,BXII_erro,GBXII_RB_estimativas,
                 GBXII_RB_erro,GLL_estimativa,GLL_erro, GLL_ZB_estimativas,
                 GLL_ZB_erro, LL_estimativas, LL_erro, GMW_ZB_estimativa, GMW_ZB_erro,GMW_estimativa,
                 GMW_erro,MW_estimavia,MW_erro)                  
colnames(tabela1) = c("a","b","alpha","beta","gama")
stargazer::stargazer(round(tabela1,6))
tabela1

A.est <- c(results_GBXII$A, results_GBXII_RB$A, results_GMW$A, results_GMW_RB$A, results_GLL$A,
           results_GLL_ZB$A, results_BXII$A, results_MW$A,results_LL$A)
W.est <- c(results_GBXII$W, results_GBXII_RB$W, results_GMW$W, results_GMW_RB$W, results_GLL$W,
           results_GLL_ZB$W, results_BXII$W, results_MW$W,results_LL$W)
CAIC.est <- c(results_GBXII$`CAIC `, results_GBXII_RB$`CAIC `, results_GMW$`CAIC `, results_GMW_RB$`CAIC `,
              results_GLL$`CAIC `, results_GLL_ZB$`CAIC `,results_BXII$`CAIC `, results_MW$`CAIC `,
              results_LL$`CAIC `) 
AIC.est <- c(results_GBXII$AIC, results_GBXII_RB$AIC, results_GMW$AIC, results_GMW_RB$AIC, results_GLL$AIC,
             results_GLL_ZB$AIC,results_BXII$AIC, results_MW$AIC,results_LL$AIC) 
BIC.est <- c(results_GBXII$BIC, results_GBXII_RB$BIC, results_GMW$BIC, results_GMW_RB$BIC, results_GLL$BIC,
             results_GLL_ZB$BIC,results_BXII$BIC, results_MW$BIC,results_LL$BIC)
HQIC.est <- c(results_GBXII$HQIC, results_GBXII_RB$HQIC, results_GMW$HQIC, results_GMW_RB$HQIC, results_GLL$HQIC,
              results_GLL_ZB$HQIC,results_BXII$HQIC, results_MW$HQIC,results_LL$HQIC)
KS.est <- c(results_GBXII$KS, results_GBXII_RB$KS, results_GMW$KS, results_GMW_RB$KS, results_GLL$KS,
            results_GLL_ZB$KS,results_BXII$KS, results_MW$KS,results_LL$KS)
tabela <- cbind(W.est, A.est, AIC.est, CAIC.est, BIC.est, HQIC.est)
rownames(tabela) <- c("ZB-GBXII","RB-GBXII","ZB-GMW","RB-GMW","RB-GLL", "ZB-GLL","BXII", "MW", "LL")
tabela
stargazer::stargazer(round(tabela,6))
