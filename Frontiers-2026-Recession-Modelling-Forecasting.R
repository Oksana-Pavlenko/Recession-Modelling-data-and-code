################################################################################
# Modelling and Forecasting Recessions in Latvia: 
# Classical and Regime-Dependent Logistic Models
################################################################################

library(readxl)
library(pscl)
library(pROC)
library(caret)
library(forecast)
library(lmtest)
################################################################################
# CORRELATION ANALYSIS
################################################################################

setwd(...)
Dati<-read_excel("DataF2026.xlsx",sheet="1",range='A1:J75',col_names = TRUE) # 
# Data for the period from 2007 in which all variables are available
head(Dati)
#
E<-Dati$Export
H<-Dati$HousingCost
C<-Dati$Confidence
I<-Dati$Inflation
R3<-Dati$R3m
Tn<-Dati$Tension
U<-Dati$Unemployment
Cov<-Dati$Cov
R<-Dati$Recesion
#
D<-data.frame(E,H,C,I,R3,Tn,U,Cov,R)
head(D)  
tail(D)

# CORELATION
r<-cor(data.frame(R,E,H,C,I,R3,Tn,U,Cov)) # H,C,I,R3 signif.
round(r,4)

################################################################################
# MODELLING WITH DATA TILL THE END 2025
################################################################################

Dati <- read_excel("DataF2026.xlsx",sheet = "bezH",range="A1:AC83",col_names = TRUE)
# Data available for wider period excluding H
head(Dati)
E<-Dati$Export; C<-Dati$Confidence; I<-Dati$Inflation; dR3<-Dati$dR3
R3<-Dati$R3m; Tn<-Dati$Tension; U<-Dati$Unemployment; Cov<-Dati$Cov
GDP<-Dati$RealGDP; R<-Dati$Recesion 
R1<-Dati$R1; E1<-Dati$Export1; C1<-Dati$C1; I1<-Dati$I1
R31<-Dati$R3m1; Tn1<-Dati$T1; U1<-Dati$U1; Cov1<-Dati$Cov1; GDP1<-Dati$GDP1
E2<-Dati$E2; C2<-Dati$C2; I2<-Dati$I2; R32<-Dati$R3m2
Tn2<-Dati$T2; U2<-Dati$U2; Cov2<-Dati$Cov2; GDP2<-Dati$GDP2; R2<-Dati$R2
#
factors<-data.frame(E,C,I,R3,Tn,U,Cov,GDP,E1,C1,I1,R31,Tn1,U1,Cov1,
                    GDP1,E2,C2,I2,R32,Tn2,U2,Cov2,GDP2)
factors<-scale(factors)
D<-data.frame(factors,R,R1,R2,dR3)
#
################################################################################
# LOGISTIC REGRESSION  - LARGE MODEL
################################################################################
# 
Mod<-glm(R~I+R1+C1+U1+E1+Tn2,     #    all signif.! 
         family=binomial,data=D)
summary(Mod)$coefficients # sign
ggtsdisplay(residuals(Mod))
pR2(Mod) #                              McFaden 0.6545707 USE IT IN PAPER!
car::vif(Mod) # Multicollinearity!

################################################################################
# LOGISTIC REGRESSION  - CHOSEN MODEL
################################################################################
#  
Mod<-glm(R~R1+I+C,family=binomial,data=D)     #   RIC model
summary(Mod)$coefficients # sign
ggtsdisplay(residuals(Mod))
pR2(Mod) #                              McFaden 0.4247203 USE IT IN PAPER!
car::vif(Mod) # NO Multicollinearity!
AIC(Mod)

# AUC
p1 <- predict(Mod, type = "response")
rocLogit<-roc(R, p1)
auc(rocLogit) # AUC (Area under the curve)

################################################################################
# THRESHOLD REGRESSION  - CHOSEN MODEL
################################################################################

C<-scale(C); C1<-scale(C1); I<-scale(I)
Y <- R
X <- as.matrix(cbind(R1,I,C)) 
X1<-R1
X2<-I 
X3<-C
Z<-C1 # threshold variable
#
thresholds <- quantile(C1,probs = seq(0.1, 0.9, 0.01),na.rm = TRUE)

# Log-likelihood function:
loglik_fun <- function(gamma) {
  X_temp <- X
  X_temp$regime <- ifelse(C1 <= gamma, 0, 1)
  model <- glm(Y ~ (1+X1+X2+X3) * regime, 
               data = X_temp,
               family = binomial)
  as.numeric(logLik(model))
}

# Threshold optimization:
ll_values <- suppressWarnings(sapply(thresholds, loglik_fun))
best_gamma <- thresholds[which.max(ll_values)]
cat(best_gamma,"- best threshold")

# Final model:
X$regime <- ifelse(C1 <= best_gamma, 0, 1)
final_model <- glm(Y ~ (1+X1+X2+X3) * regime,
                   data = X,
                   family = binomial)
summary(final_model)


X$pred_prob <- predict(final_model,type = "response")
X$pred_class <- ifelse(X$pred_prob > 0.5, 1, 0)
table(True = Y, Predicted = X$pred_class)

rocThresh <- roc(Y, X$pred_prob)
auc(rocThresh) # AUC

LL_model <- as.numeric(logLik(final_model))
null_model <- glm(Y ~ 1, family = binomial)
LL_null <- as.numeric(logLik(null_model))
pseudo_R2 <- 1 - (LL_model / LL_null)
pseudo_R2 # 0.5626605

################################################################################
# STL  - CHOSEN MODEL
################################################################################

R  <- Dati$Recesion
R1 <- Dati$R1
I1 <- Dati$Inflation
C1 <- Dati$C1
C <- Dati$Confidence

################################################################################
# STANDARDIZATION
################################################################################

X <- scale(data.frame(I1, C1,C))
X <- data.frame(X, R = R, R1=R1)

################################################################################
# TRANSITION FUNCTION
################################################################################

G_fun <- function(gamma, c0, C1){
  1 / (1 + exp(-gamma * (C1 - c0)))
}

################################################################################
# LOG-LIKELIHOOD FUNCTION
################################################################################

loglik_STR <- function(par){
  
  b0 <- par[1]
  b1 <- par[2]
  b2 <- par[3]
  b3 <- par[4]
  
  d0 <- par[5]
  d1 <- par[6]
  d2 <- par[7]
  d3 <- par[8]
  
  gamma <- par[9]
  c0    <- par[10]
  
  G <- G_fun(gamma, c0, X$C1)
  
  eta <- (b0 + d0*G) +
    (b1 + d1*G)*X$R1 +
    (b2 + d2*G)*X$I1 +
    (b3 + d3*G)*X$C
  
  p <- 1/(1 + exp(-eta))
  p <- pmin(pmax(p, 1e-6), 1 - 1e-6)
  
  -sum(X$R * log(p) + (1 - X$R) * log(1 - p))
}

# STARTING VALUES
start_par <- c(
  0, 0, 0, 0,   # b0 b1 b2 b3
  0, 0, 0, 0,   # d0 d1 d2 d3
  10,            # gamma
  0.3  # c
)

# ESTIMATION:
fit <- optim(
  par = start_par,
  fn = loglik_STR,
  method = "BFGS",
  hessian = TRUE,
  control = list(maxit = 5000)
)


# PARAMETERS:
b0 <- fit$par[1]; b1 <- fit$par[2]; b2 <- fit$par[3]; b3 <- fit$par[4]
d0 <- fit$par[5]; d1 <- fit$par[6]; d2 <- fit$par[7]; d3 <- fit$par[8]
gamma_hat <- fit$par[9]
c_hat     <- fit$par[10]

# STANDARD ERRORS:
vcov_mat <- solve(fit$hessian)
se <- sqrt(diag(vcov_mat))

results <- data.frame(
  Parameter = c("b0","b1","b2","b3",
                "d0","d1","d2","d3",
                "gamma","c"),
  Estimate = fit$par,
  StdError = se,
  z = fit$par / se,
  p = 2*(1 - pnorm(abs(fit$par / se)))
)

print(results)

# PREDICTIONS:
G <- G_fun(gamma_hat, c_hat, X$C1)

eta <- (b0 + d0*G) +
  (b1 + d1*G)*X$R1 +
  (b2 + d2*G)*X$I1 +
  (b3 + d3*G)*X$C

prob <- 1/(1 + exp(-eta))

pred <- ifelse(prob > 0.5, 1, 0)

table(Actual = X$R, Predicted = pred)

# ROC:
# First we draw previous!
plot(rocLogit, col = 2, lwd = 2, main=c("ROC curve analysis"))
df <- coords(rocLogit, x = "all", ret = c("specificity", "sensitivity"))
x <- df$specificity
y <- df$sensitivity
o <- order(x)
polygon(c(0, x[o], 1),
        c(0, y[o], 0),
        col = rgb(1, 0, 0, 0.2),
        border = NA)
roc_objSTL <- roc(X$R, prob)
plot(roc_objSTL,add=TRUE,col=4)
plot(rocThresh,add=TRUE,col=1)
#
legend("bottomright",
       legend = c( "Threshold","Smooth Transition", "Logit"),
       col = c("black", "blue", "red"),
       lty = 1, lwd = 1, bty = "n",inset = c(0.05, 0.05))
x2 <- roc_objSTL$specificities
y2 <- roc_objSTL$sensitivities
o2 <- order(x2)
polygon(c(0, x2[o2], 1),
        c(0, y2[o2], 0),
        col = rgb(0, 0, 1, 0.15),
        border = NA)

auc(roc_objSTL)

# PSEUDO R2 (McFadden):
loglik_null <- function(par0){
  b0 <- par0[1]
  eta <- b0
  p <- 1/(1 + exp(-eta))
  p <- pmin(pmax(p, 1e-6), 1 - 1e-6)
  -sum(X$R * log(p) + (1 - X$R) * log(1 - p))
}

null_fit <- optim(par = c(0),
                  fn = loglik_null,
                  method = "BFGS")

LL_null <- -null_fit$value
LL_full <- -fit$value

pseudo_R2 <- 1 - (LL_full / LL_null)
cat("McFadden pseudo R2:", pseudo_R2)

# AIC:
k <- length(fit$par)     # number of parameters
AIC_value <- 2*k + 2*fit$value
cat("AIC =", AIC_value) # AIC




################################################################################
###
###               FORECAST COMPARISON OF THE THREE MODELS                   ####
###
################################################################################

# Training period for model
train.period <- 2:70
center <- apply(factors[train.period, ],2,mean,na.rm=TRUE)
scalev <- apply(factors[train.period, ],2,sd, na.rm=TRUE)
factors.sc <- scale(factors,center=center,scale=scalev)
D <- data.frame(factors.sc, R=Dati$Recesion,R1=Dati$R1,R2=Dati$R2,dR3=Dati$dR3)

################################################################################
# LOGIT MODEL
################################################################################

# In training sample
Mod <- glm(R ~ R1 + I + C,family = binomial,data = D[2:70,])

# FORECAST FUNCTION:
forecast_logit <- function(start.period, end.period){
  n <- end.period-start.period+1
  I.f <- D$I[start.period:end.period]
  C.f <- D$C[start.period:end.period]
  R1.f <- D$R[start.period-1]  # last known recession value
  P <- numeric(n)
  for(i in 1:n){newdata <- data.frame(R1=R1.f,I=I.f[i],C=C.f[i])
        P[i] <- predict(Mod,newdata=newdata,type="response")
        R1.f <- ifelse(P[i]>=0.5,1,0) # recursive forecast
  }
  return(P)}

# LAST 12 PERIODS  (LOGISTIC REGRESSION FORECAST OF RECESSION)
test12 <- 71:82
P12 <- forecast_logit(71,82)

Result12 <- data.frame(
  Period=test12,
  Actual=D$R[test12],
  Probability=P12,
  Forecast=ifelse(P12>=0.5,1,0))
Result12


################################################################################
# THRESHOLD LOGISTIC REGRESSION FORECAST
################################################################################

train <- 2:70
test  <- 71:82

Y_train <- D$R[train]

X_train <- data.frame(R1 = D$R1[train],I  = D$I[train],
  C  = D$C[train],  C1 = D$C1[train])

# THRESHOLD SEARCH:

thresholds <- quantile( X_train$C1,probs = seq(0.5,0.8,0.05),na.rm=TRUE)

loglik_fun <- function(gamma){X_temp <- X_train
    X_temp$regime <- ifelse(
    X_temp$C1 <= gamma,0,1  )
  
  model <- glm(Y_train ~ (R1 + I + C)*regime,
    data=X_temp,family=binomial  )
    as.numeric(logLik(model))
}

ll_values <- suppressWarnings(sapply(thresholds,loglik_fun))
 
best_gamma <- thresholds[which.max(ll_values)]


# FINAL THRESHOLD MODEL (TRAINING ONLY):

X_train$regime <- ifelse(X_train$C1 <= best_gamma,  0,  1)

threshold_model <- glm(Y_train ~ (R1 + I + C)*regime,data=X_train,
                       family=binomial)

summary(threshold_model)

# RECURSIVE FORECAST FUNCTION:

forecast_threshold <- function(start.period,end.period){
      n <- end.period-start.period+1
      I.f  <- D$I[start.period:end.period]
  C.f  <- D$C[start.period:end.period]
  C1.f <- D$C1[start.period:end.period]
      R1.f <- D$R[start.period-1]
  
  P <- numeric(n)
  
    for(i in 1:n){
    regime <- ifelse(C1.f[i] <= best_gamma, 0,1 )
    newdata <- data.frame(R1=R1.f,I=I.f[i],C=C.f[i],regime=regime   )
    P[i] <- predict(threshold_model,newdata=newdata,type="response"    )
    
    # recursive update
    R1.f <- ifelse(P[i]>=0.5, 1,0)
  }
  
  return(P)
}

# 12 PERIOD FORECAST: 

P_threshold12 <- forecast_threshold(  71,  82)

Threshold_Result12 <- data.frame( Period=test, Actual=D$R[test],
    Probability=P_threshold12, 
    Forecast=ifelse(
    P_threshold12>=0.5,1,0  )
)

#THRESHOLD MODEL CONFUSION MATRIX
# table(Actual=Threshold_Result12$Actual,Forecast=Threshold_Result12$Forecast) 
# mean(Threshold_Result12$Actual == Threshold_Result12$Forecast)

roc_threshold <- roc(Threshold_Result12$Actual,
  Threshold_Result12$Probability)
auc(roc_threshold)


#
# SMOOTH TRANSITION LOGISTIC MODEL:

# Using already prepared D object from LOGIT part

train <- 2:70
test  <- 71:82

X_train <- data.frame(
    R  = D$R[train],
    R1 = D$R1[train],
    I1 = D$I1[train],
    C1 = D$C1[train],
    C  = D$C[train]
)


# STANDARDIZATION IS ALREADY DONE IN D

# TRANSITION FUNCTION
G_fun <- function(gamma,c0,z){
    1/(1+exp(-gamma*(z-c0)))
  
}

# STL LOG-LIKELIHOOD:

loglik_STR <- function(par){
  
  b0 <- par[1]
  b1 <- par[2]
  b2 <- par[3]
  b3 <- par[4]
  
  d0 <- par[5]
  d1 <- par[6]
  d2 <- par[7]
  d3 <- par[8]
  
  gamma <- par[9]
  
  c0 <- par[10]
  
  G <- G_fun(gamma, c0, X_train$C1)
  
  eta <-
    
    (b0+d0*G)+
    (b1+d1*G)*X_train$R1+
    (b2+d2*G)*X_train$I1+
    (b3+d3*G)*X_train$C
  
  p <- 1/(1+exp(-eta))
  
  p <- pmin(pmax(p,1e-6),1-1e-6)
  
  
  -sum(X_train$R*log(p)+  (1-X_train$R)*log(1-p)  )
  
}

# ESTIMATION
################################################################################

start_par <- c(
    0,0,0,0,
    0,0,0,0,
    10,
    0 )

fit_STL <- optim(par=start_par,fn=loglik_STR,method="BFGS",hessian=TRUE,
    control=list(maxit=5000))

# PARAMETERS

par_STL <- fit_STL$par

names(par_STL) <- c(
  "b0","b1","b2","b3",
  "d0","d1","d2","d3",
  "gamma","c")

print(par_STL)

# TRAINING PERFORMANCE

predict_STL <- function(par,newdata,R1_value){
  
  b0 <- par[1]
  b1 <- par[2]
  b2 <- par[3]
  b3 <- par[4]
  
  d0 <- par[5]
  d1 <- par[6]
  d2 <- par[7]
  d3 <- par[8]
  
  gamma <- par[9]
  c0 <- par[10]
  
  G <- G_fun(gamma,c0,newdata$C1)
  
    eta <-
    (b0+d0*G)+
    (b1+d1*G)*R1_value+
    (b2+d2*G)*newdata$I1+
    (b3+d3*G)*newdata$C
  
  
  1/(1+exp(-eta))
  
}

# RECURSIVE 12 PERIOD FORECAST
################################################################################

forecast_STL <- function(start.period,end.period){
      n <- end.period-start.period+1
    I1.f <- D$I1[start.period:end.period]
    C1.f <- D$C1[start.period:end.period]
    C.f  <- D$C[start.period:end.period]
    R1.f <- D$R[start.period-1]
    P <- numeric(n)
  
    for(i in 1:n){
        newdata <- data.frame(I1=I1.f[i],C1=C1.f[i],C=C.f[i])
        P[i] <- predict_STL(par_STL,newdata,R1.f)
    R1.f <- ifelse(P[i]>=0.5,1,0)
      }
    return(P)
  }

# STL FORECAST 12 PERIODS
#
P_STL12 <- forecast_STL(71,82)

STL_Result12 <- data.frame(Period=test,Actual=D$R[test],Probability=P_STL12,
  Forecast=ifelse(P_STL12>=0.5,1,0))
STL_Result12

# MODEL QUALITY

# STL CONFUSION MATRIX:
table(Actual=STL_Result12$Actual,  Forecast=STL_Result12$Forecast)

mean( STL_Result12$Actual==STL_Result12$Forecast)

roc_STL12 <- roc(  STL_Result12$Actual, STL_Result12$Probability)

auc(roc_STL12) # AUC


#  COMMON GRAPH OF PROBABILITIES:

Comparison_all <- data.frame(Period=test,Actual=D$R[test],Logit=P12,
  Threshold=P_threshold12,STL=P_STL12)
Comparison_all

plot(Comparison_all$Period,  Comparison_all$Actual,  type="o",  pch=16,
  lwd=2,  ylim=c(0,1),  xlab="Period",  ylab="Probability / Recession",
    main="Recession probability: LOGIT vs THRESHOLD vs STL",cex.main=1.1)

grid()

lines(  Comparison_all$Period,  Comparison_all$Logit,  type="o",
  pch=16,col=2,  lwd=2,  lty=1)

lines(  Comparison_all$Period,  Comparison_all$Threshold,
  type="o",  pch=16,col=3,  lwd=2,  lty=1)

lines(  Comparison_all$Period,  Comparison_all$STL,  type="o",
  pch=16,col=4,  lwd=2,  lty=1)

legend( "topright",  legend=c("Actual","Logit","Threshold","STL"),
          pch=c(16,17,18,15),lty=1,col=c(1,2,3,4),lwd=2)

#  TABLE OF FORECAST PROBABILITIES:
#
Final_Probabilities <- data.frame(Period = 71:82,Actual_Recession = D$R[71:82],
    Logit_Probability = round(P12,4),Threshold_Probability = round(P_threshold12,4),
    STL_Probability = round(P_STL12,4))

print(Final_Probabilities)

################################################################################
#                     FORECAST ERROR COMPARISON
################################################################################

actual <- Final_Probabilities$Actual_Recession
error_logit <- 
  (actual - Final_Probabilities$Logit_Probability)^2
error_threshold <- 
  (actual - Final_Probabilities$Threshold_Probability)^2
error_STL <- 
  (actual - Final_Probabilities$STL_Probability)^2

# Mean Brier Score
Brier_table <- data.frame(
  Model=c("Logit","Threshold","STL"),
  Brier=c(mean(error_logit),mean(error_threshold),mean(error_STL)))

Brier_table # The smaller, the better

# Diebold–Mariano:
dm.test(error_logit,error_threshold,alternative="two.sided")
dm.test(error_logit,error_STL,alternative="two.sided")
dm.test(error_threshold,error_STL,alternative="two.sided")

DM_results <- data.frame(
    Comparison=c("Logit vs Threshold","Logit vs STL", "Threshold vs STL"),
    DM_statistic=c(
    dm.test(error_logit,error_threshold)$statistic,
    dm.test(error_logit,error_STL)$statistic,
    dm.test(error_threshold,error_STL)$statistic  ),
    p_value=c(
    dm.test(error_logit,error_threshold)$p.value,
    dm.test(error_logit,error_STL)$p.value,
    dm.test(error_threshold,error_STL)$p.value  ))
DM_results

# MAE
MAE_table <- data.frame(
  Model=c("Logit","Threshold","STL"),
  MAE=c(mean(abs(actual-
               Final_Probabilities$Logit_Probability)),
        mean(abs(actual-
               Final_Probabilities$Threshold_Probability)),
        mean(abs(actual-
               Final_Probabilities$STL_Probability))   ))
MAE_table
