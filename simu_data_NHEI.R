library(remotes)
library(sampling)
library(survey)
library(nonprobsvy)
library(ggplot2)
library(pps)
library(compositions)
library(robustbase)
library(MatchIt)
library(dummies)
library(MASS)
library(plyr)
library(KWML)
library(lme4)
library(depmixS4)
library(snowfall)
library(NonProbEst)
library(sandwich)
library(jmuOutlier)
library(tidyverse)

#批注
##############################################################################################
###########################################Generate Common Variables
##############################################################################################
CVdatagenerate <- function(N_LCi,p_LCi,l_LCi,u_LCi,lamda_LCi,df_LCi,m_LCi,sd_LCi){
  X1 <- rbinom(n=N_LCi, size=1, prob=p_LCi) 
  X2 <- runif(n=N_LCi, min=l_LCi, max=u_LCi)
  X3 <- rexp(n=N_LCi, rate=lamda_LCi)
  X4 <- rchisq(n=N_LCi, df=df_LCi)
  X5 <- rnorm(n=N_LCi, mean=m_LCi, sd=sd_LCi)
  
  var_X1 <- p_LCi*(1-p_LCi)
  var_X2 <- ((u_LCi-l_LCi)^2)/12
  var_X3 <- 1/(lamda_LCi^2)
  var_X4 <- 2*df_LCi
  var_X5 <- sd_LCi^2
  
  m_X1 <- p_LCi
  m_X2 <- (l_LCi+u_LCi)/2
  m_X3 <- 1/lamda_LCi
  m_X4 <- df_LCi
  m_X5 <- m_LCi
  
  var_X <- var_X1+var_X2+var_X3+var_X4+var_X5
  
  out <- list(X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,
              var_X1=var_X1,var_X2=var_X2,var_X3=var_X3,var_X4=var_X4,var_X5=var_X5,
              m_X1=m_X1,m_X2=m_X2,m_X3=m_X3,m_X4=m_X4,m_X5=m_X5,
              var_X=var_X)
  return(out)
}
##############################################################################################
###########################################Sub-function: Used to calculate var_XW and coef_WtoY.
##############################################################################################
Y_generate <- function(X,W1,W2,W3,W4,W5,var_X,var_W,R2_WtoY,coef_ZXtoY){
  EX <- mean(X)
  EW1 <- mean(W1);EW2 <- mean(W2);EW3 <- mean(W3);EW4 <- mean(W4);EW5 <- mean(W5)
  COV1 <- cov(X,W1);COV2 <- cov(X,W2);COV3 <- cov(X,W3);COV4 <- cov(X,W4);COV5 <- cov(X,W5)
  var_XW <- (coef_ZXtoY^2)*(var_X*var_W+(EX^2)*var_W+var_X*(EW1^2+EW2^2+EW3^2+EW4^2+EW5^2)+COV1^2+COV2^2+COV3^2+COV4^2+COV5^2)
  coef_WtoY <- sqrt(R2_WtoY/(1-R2_WtoY)*(var_X+var_XW+1)/var_W)
  out <- list(var_XZ=var_XW,coef_ZtoY=coef_WtoY)
  return(out)
}

###########################################################################################################################################################
###########################################Main function:generate the nonprobability sample and reference sample
###########################################################################################################################################################
# Parameter Description                                                                                                                                   #
# N: Size of the finite population                                                                                                                        #
# n_A: size of the nonprobability sample                                                                                                                  #
# n_B: size of the reference sample                                                                                                                       #
# n_LC: number of heterogeneous subpopulations                                                                                                            #
# icc: degree of heterogeneity                                                                                                                            #
# n_Z: the number of covariates whose effects differ in sign across subpopulations, possible values range from 1 to 5                                     #
# n_P: the number of subpopulations whose covariate-effect directions are opposite to those of the remaining subpopulations for the n_Z covariates        #
#      for n_Z taking values of 2 or 3, n_P can take values of 1; for n_Z taking values of 4 or 5, n_P can take values of 1 and 2.                        #
# Ytype: type of target outcome, possible value include 'gaussian' and 'binomial'                                                                         #
# dist: type of parameterizations, possible values include 1 (parameters increasing monotonically) and 2 (parameters vary irregularly).                   #
# path_scenario: DAG for X, Y, Z, and nonprobability sampling, possible values range from 1 to 12. DAGs 9-12 correspond to scenarios being affected.      #
# R2_ZtoR: proportion of variance of nonprobability sampling (R) explained by Z, correlations of Z and R.                                                 #
# R2_ZtoX: proportion of variance of X explained by Z, correlations of Z and X.                                                                           #
# R2_ZtoY: proportion of variance of Y explained by Z, correlations of Z and Y.                                                                           #
# coef_ZXtoY: coefficient of interaction term ZX and Y.                                                                                                   #
# coef_XtoR: coefficient of X and R.                                                                                                                      #
# coef_YtoR: coefficient of Y and R.                                                                                                                      #
###########################################################################################################################################################
simu_data_NHEI <- function(N,n_A,n_B,n_LC,icc,n_Z,n_P,Ytype,dist,path_scenario,R2_ZtoR,R2_ZtoX,R2_ZtoY,coef_ZXtoY,coef_XtoR,coef_YtoR){
  
  
  #设定不同子人群的分布参数
  if(dist==1){
    PV <- 
      rbind(    
        data.frame(p_LCi=0.1,l_LCi=0, u_LCi=3, lamda_LCi=0.5,df_LCi=1.0, m_LCi=1, sd_LCi=1),   #1
        data.frame(p_LCi=0.2,l_LCi=5, u_LCi=8, lamda_LCi=1.0,df_LCi=2.0, m_LCi=5, sd_LCi=1),   #2
        data.frame(p_LCi=0.3,l_LCi=10,u_LCi=13,lamda_LCi=1.5,df_LCi=3.0, m_LCi=9, sd_LCi=1),   #3
        data.frame(p_LCi=0.4,l_LCi=15,u_LCi=18,lamda_LCi=2.0,df_LCi=4.0, m_LCi=13,sd_LCi=1),   #4
        data.frame(p_LCi=0.5,l_LCi=20,u_LCi=23,lamda_LCi=2.5,df_LCi=5.0, m_LCi=17,sd_LCi=1))   #5
  }
  
  if(dist==2){
    PV <- 
      rbind(    
        data.frame(p_LCi=0.4,l_LCi=15,u_LCi=18,lamda_LCi=1.5,df_LCi=3, m_LCi=5, sd_LCi=1),   #1
        data.frame(p_LCi=0.1,l_LCi=0, u_LCi=3, lamda_LCi=0.5,df_LCi=4, m_LCi=17,sd_LCi=1),   #2
        data.frame(p_LCi=0.3,l_LCi=5, u_LCi=8, lamda_LCi=2.0,df_LCi=5, m_LCi=1, sd_LCi=1),   #3
        data.frame(p_LCi=0.2,l_LCi=10,u_LCi=13,lamda_LCi=1.0,df_LCi=1, m_LCi=13,sd_LCi=1),   #4
        data.frame(p_LCi=0.5,l_LCi=20,u_LCi=23,lamda_LCi=2.5,df_LCi=2, m_LCi=9, sd_LCi=1))   #5
  }
  
  
  #First, calculate the subpopulation size for each heterogeneity scenario.----
  P_LCgroup <- 1/n_LC
  if(n_LC==1){
    N1<-N;n_A1<-n_A}
  if(n_LC==2){
    N1<-round(N*P_LCgroup,0);N2<-N-N1
    n_A1<-round(n_A*P_LCgroup,0);n_A2<-n_A-n_A1}
  if(n_LC==3){
    N1<-N2<-round(N*P_LCgroup,0);N3<-N-N1-N2
    n_A1<-n_A2<-round(n_A*P_LCgroup,0);n_A3<-n_A-n_A1-n_A2}
  if(n_LC==4){
    N1<-N2<-N3<-round(N*P_LCgroup,0);N4<-N-N1-N2-N3
    n_A1<-n_A2<-n_A3<-round(n_A*P_LCgroup,0);n_A4<-n_A-n_A1-n_A2-n_A3}
  if(n_LC==5){
    N1<-N2<-N3<-N4<-round(N*P_LCgroup,0);N5<-N-N1-N2-N3-N4
    n_A1<-n_A2<-n_A3<-n_A4<-round(n_A*P_LCgroup,0);n_A5<-n_A-n_A1-n_A2-n_A3-n_A4}
  
  S1 <- S2 <- S3 <- S4 <- S5 <- 1:5
  S12 <- S22 <- S32 <- S42 <- S52 <- 1:5
  
  #The signs of the coefficients for Z1–Z5 across different subpopulations.
  if(n_LC==1){
    PS <- matrix(1,nrow=n_LC,ncol=5) }
  if(n_LC==2 & n_Z==1 & n_P==1){
    PS <-
      matrix(c(-1,1,1,1,1,  #1
               1,1,1,1,1)   #2
             ,nrow=2,byrow = T)    }
  if(n_LC==2 & n_Z==2 & n_P==1){
    PS <-
      matrix(c(-1,1,1,1,1,  #1
               1,-1,1,1,1)   #2
             ,nrow=2,byrow = T)    }
  if(n_LC==2 & n_Z==3 & n_P==1){
    PS <-
      matrix(c(-1,1,-1,1,1,  #1
               1,-1,1,1,1)   #2
             ,nrow=2,byrow = T)    }
  if(n_LC==2 & n_Z==4 & n_P==1){
    PS <-
      matrix(c(-1,1,-1,1,1,  #1
               1,-1,1,-1,1)   #2
             ,nrow=2,byrow = T)    }
  if(n_LC==2 & n_Z==5 & n_P==1){
    PS <-
      matrix(c(-1,1,-1,1,-1,  #1
               1,-1,1,-1,1)   #2
             ,nrow=2,byrow = T)    }
  if(n_LC==3 & n_Z==2 & n_P==1){
    PS <-
      matrix(c(-1,1,1,1,1,   #1
               1,-1,1,1,1,   #2
               1,1,1,1,1)   #3
             ,nrow=3,byrow = T)    }
  if(n_LC==3 & n_Z==3 & n_P==1){
    PS <-
      matrix(c(-1,1,1,1,1,   #1
               1,-1,1,1,1,   #2
               1,1,-1,1,1)   #3
             ,nrow=3,byrow = T)    }
  if(n_LC==3 & n_Z==4 & n_P==1){
    PS <-
      matrix(c(-1,1,1,-1,1,   #1
               1,-1,1,1,1,   #2
               1,1,-1,1,1)   #3
             ,nrow=3,byrow = T)    }
  if(n_LC==3 & n_Z==5 & n_P==1){
    PS <-
      matrix(c(-1,1,1,-1,1,   #1
               1,-1,1,1,-1,   #2
               1,1,-1,1,1)   #3
             ,nrow=3,byrow = T)    }
  if(n_LC==4 & n_Z==3 & n_P==1){
    PS <-
      matrix(c(-1,1,1,1,1,   #1
               1,-1,1,1,1,   #2
               1,1,-1,1,1,   #3
               1,1,1,1,1)    #4
             ,nrow=4,byrow = T)    }
  if(n_LC==4 & n_Z==4 & n_P==1){
    PS <-
      matrix(c(-1,1,1,1,1,   #1
               1,-1,1,1,1,   #2
               1,1,-1,1,1,   #3
               1,1,1,-1,1)   #4
             ,nrow=4,byrow = T)    }
  if(n_LC==4 & n_Z==5 & n_P==1){
    PS <-
      matrix(c(-1,1,1,1,-1,   #1
               1,-1,1,1,1,   #2
               1,1,-1,1,1,   #3
               1,1,1,-1,1)   #4
             ,nrow=4,byrow = T)    }
  if(n_LC==4 & n_Z==3 & n_P==2){
    PS <-
      matrix(c(-1,1,1,1,1,   #1
               -1,-1,1,1,1,   #2
               1,-1,-1,1,1,   #3
               1,1,-1,1,1)    #4
             ,nrow=4,byrow = T)    }
  if(n_LC==4 & n_Z==4 & n_P==2){
    PS <-
      matrix(c(-1,1,1,-1,1,   #1
               -1,-1,1,1,1,   #2
               1,-1,-1,1,1,   #3
               1,1,-1,-1,1)   #4
             ,nrow=4,byrow = T)    }
  if(n_LC==4 & n_Z==5 & n_P==2){
    PS <-
      matrix(c(-1,1,1,-1,-1,   #1
               -1,-1,1,1,-1,   #2
               1,-1,-1,1,1,   #3
               1,1,-1,-1,1)   #4
             ,nrow=4,byrow = T)    }
  if(n_LC==5 & n_Z==4 & n_P==1){
    PS <-
      matrix(c(-1,1,1,1,1,   #1
               1,-1,1,1,1,   #2
               1,1,-1,1,1,   #3
               1,1,1,-1,1,   #4
               1,1,1,1,1)    #5
             ,nrow=5,byrow = T)    }
  if(n_LC==5 & n_Z==5 & n_P==1){
    PS <-
      matrix(c(-1,1,1,1,1,   #1
               1,-1,1,1,1,   #2
               1,1,-1,1,1,   #3
               1,1,1,-1,1,   #4
               1,1,1,1,-1)    #5
             ,nrow=5,byrow = T)    }
  if(n_LC==5 & n_Z==4 & n_P==2){
    PS <-
      matrix(c(-1,1,1,1,1,   #1
               -1,-1,1,1,1,   #2
               1,-1,-1,1,1,   #3
               1,1,-1,-1,1,   #4
               1,1,1,-1,1)    #5
             ,nrow=5,byrow = T)    }
  if(n_LC==5 & n_Z==5 & n_P==2){
    PS <-
      matrix(c(-1,1,1,1,-1,   #1
               -1,-1,1,1,1,   #2
               1,-1,-1,1,1,   #3
               1,1,-1,-1,1,   #4
               1,1,1,-1,-1)   #5
             ,nrow=5,byrow = T)    }
  
  
  #The number of heterogeneous populations is【1】----------
  if(n_LC %in% c(1,2,3,4,5)){
    
    ########################################################Generate correction variable Z
    #The first independent variable of the latent subclass
    CVdata_LC1 <- CVdatagenerate(N_LCi=N1,p_LCi=PV[S1[1],1],l_LCi=PV[S2[1],2],u_LCi=PV[S2[1],3],lamda_LCi=PV[S3[1],4],df_LCi=PV[S4[1],5],m_LCi=PV[S5[1],6],sd_LCi=PV[S5[1],7])
    Z11 <- CVdata_LC1$X1;var_Z11 <- CVdata_LC1$var_X1;m_Z11 <- CVdata_LC1$m_X1
    Z12 <- CVdata_LC1$X2;var_Z12 <- CVdata_LC1$var_X2;m_Z12 <- CVdata_LC1$m_X2
    Z13 <- CVdata_LC1$X3;var_Z13 <- CVdata_LC1$var_X3;m_Z13 <- CVdata_LC1$m_X3
    Z14 <- CVdata_LC1$X4;var_Z14 <- CVdata_LC1$var_X4;m_Z14 <- CVdata_LC1$m_X4
    Z15 <- CVdata_LC1$X5;var_Z15 <- CVdata_LC1$var_X5;m_Z15 <- CVdata_LC1$m_X5
    var_Z1 <- (PS[S12[1],1]^2)*var_Z11+(PS[S22[1],2]^2)*var_Z12+(PS[S32[1],3]^2)*var_Z13+(PS[S42[1],4]^2)*var_Z14+(PS[S52[1],5]^2)*var_Z15
    miu1 <- PS[S12[1],1]*m_Z11+PS[S22[1],2]*m_Z12+PS[S32[1],3]*m_Z13+PS[S42[1],4]*m_Z14+PS[S52[1],5]*m_Z15
    
    ###################################################################################Generating X
    ########################################################Scenarios 1–4, 9–12: Generating X using Z
    if(path_scenario %in% c(1,2,3,4,9,10,11,12)){
      coef_Z1toX1 <- sqrt(R2_ZtoX/(1-R2_ZtoX)*1/var_Z1)
      X1 <- coef_Z1toX1*(PS[S12[1],1]*Z11+PS[S22[1],2]*Z12+PS[S32[1],3]*Z13+PS[S42[1],4]*Z14+PS[S52[1],5]*Z15)+rnorm(N1)
      var_X1 <- (coef_Z1toX1^2)*var_Z1+1
      m_X1 <- coef_Z1toX1*miu1}
    ########################################################Scenarios 5-8: Generating X Independently
    if(path_scenario %in% c(5,6,7,8)){
      X1 <- rnorm(n=N1, mean=5, sd=1)
      var_X1 <- 1
      m_X1 <- 5}
    
    ###################################################################################Generating Y
    ########################################################Scenarios 1234: Generating Y using X
    if(path_scenario %in% c(1,2,3,4)){
      if(Ytype=='gaussian'){
        Y1 <- 2*X1+rnorm(N1)
        var_Y1 <- 4*var_X1+1
        m_Y1 <- 2*m_X1
      }#Ytype=='gaussian'
      
      if(Ytype=='binomial'){
        Y1s <- 2*X1+rnorm(N1)
        Y1c <- Y1s-min(Y1s)+2;pi_Y1 <- inclusionprobabilities(Y1c, N1/2);Y1 <- UPpivotal(pi_Y1)
        m_Y1 <- 1/2;var_Y1 <- 1/4
      }#Ytype=='binomial'
      
    }#path_scenario %in% c(1,2,3,4)
    ########################################################Scenarios 5-12: Generating Y Using X and Z
    if(path_scenario %in% c(5,6,7,8,9,10,11,12)){
      
      OUT1 <- Y_generate(X1,Z11,Z12,Z13,Z14,Z15,var_X1,var_Z1,R2_ZtoY,coef_ZXtoY)
      var_X1Z1 <-OUT1$var_XZ
      coef_Z1toY1 <-OUT1$coef_ZtoY 
      
      if(Ytype=='gaussian'){
        Y1 <- 2*X1+coef_Z1toY1*(PS[S12[1],1]*Z11+PS[S22[1],2]*Z12+PS[S32[1],3]*Z13+PS[S42[1],4]*Z14+PS[S52[1],5]*Z15)+coef_ZXtoY*X1*(PS[S12[1],1]*Z11+PS[S22[1],2]*Z12+PS[S32[1],3]*Z13+PS[S42[1],4]*Z14+PS[S52[1],5]*Z15)+rnorm(N1)
        var_Y1 <- 4*var_X1+(coef_Z1toY1^2)*var_Z1+var_X1Z1+1
        m_Y1 <- 2*m_X1+coef_Z1toY1*miu1+coef_ZXtoY*m_X1*miu1
      }#Ytype=='gaussian'
      
      if(Ytype=='binomial'){
        Y1s <- 2*X1+coef_Z1toY1*(PS[S12[1],1]*Z11+PS[S22[1],2]*Z12+PS[S32[1],3]*Z13+PS[S42[1],4]*Z14+PS[S52[1],5]*Z15)+coef_ZXtoY*X1*(PS[S12[1],1]*Z11+PS[S22[1],2]*Z12+PS[S32[1],3]*Z13+PS[S42[1],4]*Z14+PS[S52[1],5]*Z15)+rnorm(N1)
        Y1c <- Y1s-min(Y1s)+2;pi_Y1 <- inclusionprobabilities(Y1c, N1/2);Y1 <- UPpivotal(pi_Y1)
        m_Y1 <- 1/2;var_Y1 <- 1/4
      }#Ytype=='binomial'
      
    }#path_scenario %in% c(5,6,7,8,9,10,11,12)
    
    ###################################################################################Mean and Variance of the Non-Probability Inclusion Probability Model
    ########################################################Scenarios 1, 5, 9: Generating R using Z
    if(path_scenario %in% c(1,5,9)){
      coef_Z1toR1 <- sqrt(R2_ZtoR/(1-R2_ZtoR)*(1)/var_Z1)
      miu_A1 <- coef_Z1toR1*miu1
      var_A1 <- (coef_Z1toR1^2)*var_Z1+1
    }
    ########################################################Scenarios 2, 6, 10: Generate R using X and Z.
    if(path_scenario %in% c(2,6,10)){
      coef_Z1toR1 <- sqrt(R2_ZtoR/(1-R2_ZtoR)*(1+(coef_XtoR^2)*var_X1)/var_Z1)
      miu_A1 <- coef_Z1toR1*miu1+coef_XtoR*m_X1
      var_A1 <- (coef_Z1toR1^2)*var_Z1+(coef_XtoR^2)*var_X1+1
    }
    ########################################################Scenarios 3, 7, 11: Generate R using Y and Z.
    if(path_scenario %in% c(3,7,11)){
      coef_Z1toR1 <- sqrt(R2_ZtoR/(1-R2_ZtoR)*((coef_YtoR^2)*var_Y1+1)/var_Z1)
      miu_A1 <- coef_Z1toR1*miu1+coef_YtoR*m_Y1
      var_A1 <- (coef_Z1toR1^2)*var_Z1+(coef_YtoR^2)*var_Y1+1
    }
    ########################################################Scenarios 4, 8, 12: Generating R using X, Y, and Z
    if(path_scenario %in% c(4,8,12)){
      coef_Z1toR1 <- sqrt(R2_ZtoR/(1-R2_ZtoR)*((coef_XtoR^2)*var_X1+1+(coef_YtoR^2)*var_Y1)/var_Z1)
      miu_A1 <- coef_Z1toR1*miu1+coef_XtoR*m_X1+coef_YtoR*m_Y1
      var_A1 <- (coef_Z1toR1^2)*var_Z1+(coef_XtoR^2)*var_X1+(coef_YtoR^2)*var_Y1+1
    }
  }#n_LC %in% c(1,2,3,4,5)
  
  if(n_LC == 1){
    ###################################################################################Generate indicator variables for non-probability sampling.
    ########################################################Scenarios 1, 5, 9: Generate R with Z.
    if(path_scenario %in% c(1,5,9)){
      linear_A1 <- coef_Z1toR1*(PS[S12[1],1]*Z11+PS[S22[1],2]*Z12+PS[S32[1],3]*Z13+PS[S42[1],4]*Z14+PS[S52[1],5]*Z15)+rnorm(N1)
    }
    ########################################################Scenarios 2, 6, 10: Generate R using X and Z.
    if(path_scenario %in% c(2,6,10)){
      linear_A1 <- coef_Z1toR1*(PS[S12[1],1]*Z11+PS[S22[1],2]*Z12+PS[S32[1],3]*Z13+PS[S42[1],4]*Z14+PS[S52[1],5]*Z15)+coef_XtoR*X1+rnorm(N1)
    }
    ########################################################Scenarios 3, 7, 11: Generate R using Y and Z.
    if(path_scenario %in% c(3,7,11)){
      linear_A1 <- coef_Z1toR1*(PS[S12[1],1]*Z11+PS[S22[1],2]*Z12+PS[S32[1],3]*Z13+PS[S42[1],4]*Z14+PS[S52[1],5]*Z15)+coef_YtoR*Y1+rnorm(N1)
    }
    ########################################################Scenarios 4, 8, 12: Generating R using X, Y, and Z
    if(path_scenario %in% c(4,8,12)){
      linear_A1 <- coef_Z1toR1*(PS[S12[1],1]*Z11+PS[S22[1],2]*Z12+PS[S32[1],3]*Z13+PS[S42[1],4]*Z14+PS[S52[1],5]*Z15)+coef_XtoR*X1+coef_YtoR*Y1+rnorm(N1)
    }
    ########################################################Generate the sampling probabilities and indicator variables based on the values ​​generated above.
    linear_A1_offset <- linear_A1-min(linear_A1)+2
    pi_A1 <- inclusionprobabilities(linear_A1_offset, n_A1)
    ind_A1 <- UPpoisson(pi_A1)
  }#LC=1
  
  #The number of heterogeneous populations is【2】----------
  if(n_LC %in% c(2,3,4,5)){
    ########################################################
    #The second independent variable of the latent subclass
    CVdata_LC2 <- CVdatagenerate(N_LCi=N2,p_LCi=PV[S1[2],1],l_LCi=PV[S2[2],2],u_LCi=PV[S2[2],3],lamda_LCi=PV[S3[2],4],df_LCi=PV[S4[2],5],m_LCi=PV[S5[2],6],sd_LCi=PV[S5[2],7])
    Z21 <- CVdata_LC2$X1;var_Z21 <- CVdata_LC2$var_X1;m_Z21 <- CVdata_LC2$m_X1
    Z22 <- CVdata_LC2$X2;var_Z22 <- CVdata_LC2$var_X2;m_Z22 <- CVdata_LC2$m_X2
    Z23 <- CVdata_LC2$X3;var_Z23 <- CVdata_LC2$var_X3;m_Z23 <- CVdata_LC2$m_X3
    Z24 <- CVdata_LC2$X4;var_Z24 <- CVdata_LC2$var_X4;m_Z24 <- CVdata_LC2$m_X4
    Z25 <- CVdata_LC2$X5;var_Z25 <- CVdata_LC2$var_X5;m_Z25 <- CVdata_LC2$m_X5
    var_Z2 <- (PS[S12[2],1]^2)*var_Z21+(PS[S22[2],2]^2)*var_Z22+(PS[S32[2],3]^2)*var_Z23+(PS[S42[2],4]^2)*var_Z24+(PS[S52[2],5]^2)*var_Z25
    miu2 <- PS[S12[2],1]*m_Z21+PS[S22[2],2]*m_Z22+PS[S32[2],3]*m_Z23+PS[S42[2],4]*m_Z24+PS[S52[2],5]*m_Z25
    
    ###################################################################################
    ########################################################
    if(path_scenario %in% c(1,2,3,4,9,10,11,12)){
      coef_Z2toX2 <- sqrt(R2_ZtoX/(1-R2_ZtoX)*1/var_Z2)
      X2 <- coef_Z2toX2*(PS[S12[2],1]*Z21+PS[S22[2],2]*Z22+PS[S32[2],3]*Z23+PS[S42[2],4]*Z24+PS[S52[2],5]*Z25)+rnorm(N2)
      var_X2 <- (coef_Z2toX2^2)*var_Z2+1
      m_X2 <- coef_Z2toX2*miu2}
    ########################################################
    if(path_scenario %in% c(5,6,7,8)){
      X2 <- rnorm(n=N2, mean=5, sd=1)
      var_X2 <- 1
      m_X2 <- 5}
    ###################################################################################
    ########################################################
    if(path_scenario %in% c(1,2,3,4)){
      
      if(Ytype=='gaussian'){
        Y2 <- 2*X2+rnorm(N2)
        var_Y2 <- 4*var_X2+1
        m_Y2 <- 2*m_X2
      }#Ytype=='gaussian'
      
      if(Ytype=='binomial'){
        Y2s <- 2*X2+rnorm(N2)
        Y2c <- Y1s-min(Y2s)+2;pi_Y2 <- inclusionprobabilities(Y2c, N2/2);Y2 <- UPpivotal(pi_Y2)
        m_Y2 <- 1/2;var_Y2 <- 1/4
      }#Ytype=='binomial'
      
    }#path_scenario %in% c(1,2,3,4)
    ########################################################
    if(path_scenario %in% c(5,6,7,8,9,10,11,12)){
      
      OUT2 <- Y_generate(X2,Z21,Z22,Z23,Z24,Z25,var_X2,var_Z2,R2_ZtoY,coef_ZXtoY)
      var_X2Z2 <-OUT2$var_XZ
      coef_Z2toY2 <-OUT2$coef_ZtoY 
      
      if(Ytype=='gaussian'){
        Y2 <- 2*X2+coef_Z2toY2*(PS[S12[2],1]*Z21+PS[S22[2],2]*Z22+PS[S32[2],3]*Z23+PS[S42[2],4]*Z24+PS[S52[2],5]*Z25)+coef_ZXtoY*X2*(PS[S12[2],1]*Z21+PS[S22[2],2]*Z22+PS[S32[2],3]*Z23+PS[S42[2],4]*Z24+PS[S52[2],5]*Z25)+rnorm(N2)
        var_Y2 <- 4*var_X2+(coef_Z2toY2^2)*var_Z2+var_X2Z2+1
        m_Y2 <- 2*m_X2+coef_Z2toY2*miu2+coef_ZXtoY*m_X2*miu2
      }#Ytype=='gaussian'
      
      if(Ytype=='binomial'){
        Y2s <- 2*X2+coef_Z2toY2*(PS[S12[2],1]*Z21+PS[S22[2],2]*Z22+PS[S32[2],3]*Z23+PS[S42[2],4]*Z24+PS[S52[2],5]*Z25)+coef_ZXtoY*X2*(PS[S12[2],1]*Z21+PS[S22[2],2]*Z22+PS[S32[2],3]*Z23+PS[S42[2],4]*Z24+PS[S52[2],5]*Z25)+rnorm(N2)
        Y2c <- Y1s-min(Y2s)+2;pi_Y2 <- inclusionprobabilities(Y2c, N2/2);Y2 <- UPpivotal(pi_Y2)
        m_Y2 <- 1/2;var_Y2 <- 1/4
      }#Ytype=='binomial'
      
    }#path_scenario %in% c(5,6,7,8,9,10,11,12)
    ###################################################################################
    ########################################################
    if(path_scenario %in% c(1,5,9)){
      coef_Z2toR2 <- sqrt(R2_ZtoR/(1-R2_ZtoR)*(1)/var_Z2)
      miu_A2 <- coef_Z2toR2*miu2
      var_A2 <- (coef_Z2toR2^2)*var_Z2+1
    }
    ########################################################
    if(path_scenario %in% c(2,6,10)){
      coef_Z2toR2 <- sqrt(R2_ZtoR/(1-R2_ZtoR)*(1+(coef_XtoR^2)*var_X2)/var_Z2)
      miu_A2 <- coef_Z2toR2*miu2+coef_XtoR*m_X2
      var_A2 <- (coef_Z2toR2^2)*var_Z2+(coef_XtoR^2)*var_X2+1
    }
    ########################################################
    if(path_scenario %in% c(3,7,11)){
      coef_Z2toR2 <- sqrt(R2_ZtoR/(1-R2_ZtoR)*((coef_YtoR^2)*var_Y2+1)/var_Z2)
      miu_A2 <- coef_Z2toR2*miu2+coef_YtoR*m_Y2
      var_A2 <- (coef_Z2toR2^2)*var_Z2+(coef_YtoR^2)*var_Y2+1
    }
    ########################################################
    if(path_scenario %in% c(4,8,12)){
      coef_Z2toR2 <- sqrt(R2_ZtoR/(1-R2_ZtoR)*((coef_XtoR^2)*var_X2+1+(coef_YtoR^2)*var_Y2)/var_Z2)
      miu_A2 <- coef_Z2toR2*miu2+coef_XtoR*m_X2+coef_YtoR*m_Y2
      var_A2 <- (coef_Z2toR2^2)*var_Z2+(coef_XtoR^2)*var_X2+(coef_YtoR^2)*var_Y2+1
    }
    
    
  }#n_LC %in% c(2,3,4,5)
  
  if(n_LC == 2){
    ###################################################################################
    #Within-group variance
    var_I <- 0.5*var_A1+0.5*var_A2
    var_B <- var_I/(1-icc)*icc
    a <- sqrt(var_B)
    
    b10 <- -1*a-miu_A1
    b20 <- a-miu_A2
    ########################################################
    if(path_scenario %in% c(1,5,9)){
      linear_A1 <- coef_Z1toR1*(PS[S12[1],1]*Z11+PS[S22[1],2]*Z12+PS[S32[1],3]*Z13+PS[S42[1],4]*Z14+PS[S52[1],5]*Z15)+rnorm(N1)+b10
      linear_A2 <- coef_Z2toR2*(PS[S12[2],1]*Z21+PS[S22[2],2]*Z22+PS[S32[2],3]*Z23+PS[S42[2],4]*Z24+PS[S52[2],5]*Z25)+rnorm(N2)+b20
    }
    ########################################################
    if(path_scenario %in% c(2,6,10)){
      linear_A1 <- coef_Z1toR1*(PS[S12[1],1]*Z11+PS[S22[1],2]*Z12+PS[S32[1],3]*Z13+PS[S42[1],4]*Z14+PS[S52[1],5]*Z15)+coef_XtoR*X1+rnorm(N1)+b10
      linear_A2 <- coef_Z2toR2*(PS[S12[2],1]*Z21+PS[S22[2],2]*Z22+PS[S32[2],3]*Z23+PS[S42[2],4]*Z24+PS[S52[2],5]*Z25)+coef_XtoR*X2+rnorm(N2)+b20
    }
    ########################################################
    if(path_scenario %in% c(3,7,11)){
      linear_A1 <- coef_Z1toR1*(PS[S12[1],1]*Z11+PS[S22[1],2]*Z12+PS[S32[1],3]*Z13+PS[S42[1],4]*Z14+PS[S52[1],5]*Z15)+coef_YtoR*Y1+rnorm(N1)+b10
      linear_A2 <- coef_Z2toR2*(PS[S12[2],1]*Z21+PS[S22[2],2]*Z22+PS[S32[2],3]*Z23+PS[S42[2],4]*Z24+PS[S52[2],5]*Z25)+coef_YtoR*Y2+rnorm(N2)+b20
    }
    ########################################################
    if(path_scenario %in% c(4,8,12)){
      linear_A1 <- coef_Z1toR1*(PS[S12[1],1]*Z11+PS[S22[1],2]*Z12+PS[S32[1],3]*Z13+PS[S42[1],4]*Z14+PS[S52[1],5]*Z15)+coef_XtoR*X1+coef_YtoR*Y1+rnorm(N1)+b10
      linear_A2 <- coef_Z2toR2*(PS[S12[2],1]*Z21+PS[S22[2],2]*Z22+PS[S32[2],3]*Z23+PS[S42[2],4]*Z24+PS[S52[2],5]*Z25)+coef_XtoR*X2+coef_YtoR*Y2+rnorm(N2)+b20
    }
    ########################################################
    linear_A1_offset <- linear_A1-min(linear_A1)+2
    pi_A1 <- inclusionprobabilities(linear_A1_offset, n_A1)
    ind_A1 <- UPpoisson(pi_A1)
    
    linear_A2_offset <- linear_A2-min(linear_A2)+2
    pi_A2 <- inclusionprobabilities(linear_A2_offset, n_A2)
    ind_A2 <- UPpoisson(pi_A2)
    
  }#LC=2
  
  
  #The number of heterogeneous populations is【3】----------
  if(n_LC %in% c(3,4,5)){
    #The third independent variable of the latent subclass
    CVdata_LC3 <- CVdatagenerate(N_LCi=N3,p_LCi=PV[S1[3],1],l_LCi=PV[S2[3],2],u_LCi=PV[S2[3],3],lamda_LCi=PV[S3[3],4],df_LCi=PV[S4[3],5],m_LCi=PV[S5[3],6],sd_LCi=PV[S5[3],7])
    Z31 <- CVdata_LC3$X1;var_Z31 <- CVdata_LC3$var_X1;m_Z31 <- CVdata_LC3$m_X1
    Z32 <- CVdata_LC3$X2;var_Z32 <- CVdata_LC3$var_X2;m_Z32 <- CVdata_LC3$m_X2
    Z33 <- CVdata_LC3$X3;var_Z33 <- CVdata_LC3$var_X3;m_Z33 <- CVdata_LC3$m_X3
    Z34 <- CVdata_LC3$X4;var_Z34 <- CVdata_LC3$var_X4;m_Z34 <- CVdata_LC3$m_X4
    Z35 <- CVdata_LC3$X5;var_Z35 <- CVdata_LC3$var_X5;m_Z35 <- CVdata_LC3$m_X5
    var_Z3 <- (PS[S12[3],1]^2)*var_Z31+(PS[S22[3],2]^2)*var_Z32+(PS[S32[3],3]^2)*var_Z33+(PS[S42[3],4]^2)*var_Z34+(PS[S52[3],5]^2)*var_Z35
    miu3 <- PS[S12[3],1]*m_Z31+PS[S22[3],2]*m_Z32+PS[S32[3],3]*m_Z33+PS[S42[3],4]*m_Z34+PS[S52[3],5]*m_Z35
    
    ###################################################################################
    ########################################################
    if(path_scenario %in% c(1,2,3,4,9,10,11,12)){
      coef_Z3toX3 <- sqrt(R2_ZtoX/(1-R2_ZtoX)*1/var_Z3)
      X3 <- coef_Z3toX3*(PS[S12[3],1]*Z31+PS[S22[3],2]*Z32+PS[S32[3],3]*Z33+PS[S42[3],4]*Z34+PS[S52[3],5]*Z35)+rnorm(N3)
      var_X3 <- (coef_Z3toX3^2)*var_Z3+1
      m_X3 <- coef_Z3toX3*miu3}
    ########################################################
    if(path_scenario %in% c(5,6,7,8)){
      X3 <- rnorm(n=N3, mean=5, sd=1)
      var_X3 <- 1
      m_X3 <- 5}
    
    ###################################################################################
    ########################################################
    if(path_scenario %in% c(1,2,3,4)){
      
      if(Ytype=='gaussian'){
        Y3 <- 2*X3+rnorm(N3)
        var_Y3 <- 4*var_X3+1
        m_Y3 <- 2*m_X3
      }#Ytype=='gaussian'
      
      if(Ytype=='binomial'){
        Y3s <- 2*X3+rnorm(N3)
        Y3c <- Y3s-min(Y3s)+2;pi_Y3 <- inclusionprobabilities(Y3c, N3/2);Y3 <- UPpivotal(pi_Y3)
        m_Y3 <- 1/2;var_Y3 <- 1/4
      }#Ytype=='binomial'
      
    }#path_scenario %in% c(1,2,3,4)
    ########################################################
    if(path_scenario %in% c(5,6,7,8,9,10,11,12)){
      
      OUT3 <- Y_generate(X3,Z31,Z32,Z33,Z34,Z35,var_X3,var_Z3,R2_ZtoY,coef_ZXtoY)
      var_X3Z3 <-OUT3$var_XZ
      coef_Z3toY3 <-OUT3$coef_ZtoY 
      
      if(Ytype=='gaussian'){
        Y3 <- 2*X3+coef_Z3toY3*(PS[S12[3],1]*Z31+PS[S22[3],2]*Z32+PS[S32[3],3]*Z33+PS[S42[3],4]*Z34+PS[S52[3],5]*Z35)+coef_ZXtoY*X3*(PS[S12[3],1]*Z31+PS[S22[3],2]*Z32+PS[S32[3],3]*Z33+PS[S42[3],4]*Z34+PS[S52[3],5]*Z35)+rnorm(N3)
        var_Y3 <- 4*var_X3+(coef_Z3toY3^2)*var_Z3+var_X3Z3+1
        m_Y3 <- 2*m_X3+coef_Z3toY3*miu3+coef_ZXtoY*m_X3*miu3
      }#Ytype=='gaussian'
      
      if(Ytype=='binomial'){
        Y3s <- 2*X3+coef_Z3toY3*(PS[S12[3],1]*Z31+PS[S22[3],2]*Z32+PS[S32[3],3]*Z33+PS[S42[3],4]*Z34+PS[S52[3],5]*Z35)+coef_ZXtoY*X3*(PS[S12[3],1]*Z31+PS[S22[3],2]*Z32+PS[S32[3],3]*Z33+PS[S42[3],4]*Z34+PS[S52[3],5]*Z35)+rnorm(N3)
        Y3c <- Y3s-min(Y3s)+2;pi_Y3 <- inclusionprobabilities(Y3c, N3/2);Y3 <- UPpivotal(pi_Y3)
        m_Y3 <- 1/2;var_Y3 <- 1/4
      }#Ytype=='binomial'
      
    }#path_scenario %in% c(5,6,7,8,9,10,11,12)
    
    ###################################################################################
    ########################################################
    if(path_scenario %in% c(1,5,9)){
      coef_Z3toR3 <- sqrt(R2_ZtoR/(1-R2_ZtoR)*(1)/var_Z3)
      miu_A3 <- coef_Z3toR3*miu2
      var_A3 <- (coef_Z3toR3^2)*var_Z3+1
    }
    ########################################################
    if(path_scenario %in% c(2,6,10)){
      coef_Z3toR3 <- sqrt(R2_ZtoR/(1-R2_ZtoR)*(1+(coef_XtoR^2)*var_X3)/var_Z3)
      miu_A3 <- coef_Z3toR3*miu2+coef_XtoR*m_X3
      var_A3 <- (coef_Z3toR3^2)*var_Z3+(coef_XtoR^2)*var_X3+1
    }
    ########################################################
    if(path_scenario %in% c(3,7,11)){
      coef_Z3toR3 <- sqrt(R2_ZtoR/(1-R2_ZtoR)*((coef_YtoR^2)*var_Y3+1)/var_Z3)
      miu_A3 <- coef_Z3toR3*miu2+coef_YtoR*m_Y3
      var_A3 <- (coef_Z3toR3^2)*var_Z3+(coef_YtoR^2)*var_Y3+1
    }
    ########################################################
    if(path_scenario %in% c(4,8,12)){
      coef_Z3toR3 <- sqrt(R2_ZtoR/(1-R2_ZtoR)*((coef_XtoR^2)*var_X3+1+(coef_YtoR^2)*var_Y3)/var_Z3)
      miu_A3 <- coef_Z3toR3*miu2+coef_XtoR*m_X3+coef_YtoR*m_Y3
      var_A3 <- (coef_Z3toR3^2)*var_Z3+(coef_XtoR^2)*var_X3+(coef_YtoR^2)*var_Y3+1
    }
  }#n_LC %in% c(3,4,5)
  
  if(n_LC == 3){
    ###################################################################################
    var_I <- 1/3*var_A1+1/3*var_A2+1/3*var_A3
    var_B <- var_I/(1-icc)*icc
    a <- sqrt(3/2*var_B)
    
    b10 <- -1*a-miu_A1
    b20 <- -miu2
    b30 <- a-miu_A3
    ########################################################
    if(path_scenario %in% c(1,5,9)){
      linear_A1 <- coef_Z1toR1*(PS[S12[1],1]*Z11+PS[S22[1],2]*Z12+PS[S32[1],3]*Z13+PS[S42[1],4]*Z14+PS[S52[1],5]*Z15)+rnorm(N1)+b10
      linear_A2 <- coef_Z2toR2*(PS[S12[2],1]*Z21+PS[S22[2],2]*Z22+PS[S32[2],3]*Z23+PS[S42[2],4]*Z24+PS[S52[2],5]*Z25)+rnorm(N2)+b20
      linear_A3 <- coef_Z3toR3*(PS[S12[3],1]*Z31+PS[S22[3],2]*Z32+PS[S32[3],3]*Z33+PS[S42[3],4]*Z34+PS[S52[3],5]*Z35)+rnorm(N3)+b30
    }
    ########################################################
    if(path_scenario %in% c(2,6,10)){
      linear_A1 <- coef_Z1toR1*(PS[S12[1],1]*Z11+PS[S22[1],2]*Z12+PS[S32[1],3]*Z13+PS[S42[1],4]*Z14+PS[S52[1],5]*Z15)+coef_XtoR*X1+rnorm(N1)+b10
      linear_A2 <- coef_Z2toR2*(PS[S12[2],1]*Z21+PS[S22[2],2]*Z22+PS[S32[2],3]*Z23+PS[S42[2],4]*Z24+PS[S52[2],5]*Z25)+coef_XtoR*X2+rnorm(N2)+b20
      linear_A3 <- coef_Z3toR3*(PS[S12[3],1]*Z31+PS[S22[3],2]*Z32+PS[S32[3],3]*Z33+PS[S42[3],4]*Z34+PS[S52[3],5]*Z35)+coef_XtoR*X3+rnorm(N3)+b30
    }
    ########################################################
    if(path_scenario %in% c(3,7,11)){
      linear_A1 <- coef_Z1toR1*(PS[S12[1],1]*Z11+PS[S22[1],2]*Z12+PS[S32[1],3]*Z13+PS[S42[1],4]*Z14+PS[S52[1],5]*Z15)+coef_YtoR*Y1+rnorm(N1)+b10
      linear_A2 <- coef_Z2toR2*(PS[S12[2],1]*Z21+PS[S22[2],2]*Z22+PS[S32[2],3]*Z23+PS[S42[2],4]*Z24+PS[S52[2],5]*Z25)+coef_YtoR*Y2+rnorm(N2)+b20
      linear_A3 <- coef_Z3toR3*(PS[S12[3],1]*Z31+PS[S22[3],2]*Z32+PS[S32[3],3]*Z33+PS[S42[3],4]*Z34+PS[S52[3],5]*Z35)+coef_YtoR*Y3+rnorm(N3)+b30
    }
    ########################################################
    if(path_scenario %in% c(4,8,12)){
      linear_A1 <- coef_Z1toR1*(PS[S12[1],1]*Z11+PS[S22[1],2]*Z12+PS[S32[1],3]*Z13+PS[S42[1],4]*Z14+PS[S52[1],5]*Z15)+coef_XtoR*X1+coef_YtoR*Y1+rnorm(N1)+b10
      linear_A2 <- coef_Z2toR2*(PS[S12[2],1]*Z21+PS[S22[2],2]*Z22+PS[S32[2],3]*Z23+PS[S42[2],4]*Z24+PS[S52[2],5]*Z25)+coef_XtoR*X2+coef_YtoR*Y2+rnorm(N2)+b20
      linear_A3 <- coef_Z3toR3*(PS[S12[3],1]*Z31+PS[S22[3],2]*Z32+PS[S32[3],3]*Z33+PS[S42[3],4]*Z34+PS[S52[3],5]*Z35)+coef_XtoR*X3+coef_YtoR*Y3+rnorm(N3)+b30
    }
    ########################################################
    linear_A1_offset <- linear_A1-min(linear_A1)+2
    pi_A1 <- inclusionprobabilities(linear_A1_offset, n_A1)
    ind_A1 <- UPpoisson(pi_A1)
    
    linear_A2_offset <- linear_A2-min(linear_A2)+2
    pi_A2 <- inclusionprobabilities(linear_A2_offset, n_A2)
    ind_A2 <- UPpoisson(pi_A2)
    
    linear_A3_offset <- linear_A3-min(linear_A3)+2
    pi_A3 <- inclusionprobabilities(linear_A3_offset, n_A3)
    ind_A3 <- UPpoisson(pi_A3)
  }#LC=3
  
  #The number of heterogeneous populations is【4】----------
  if(n_LC %in% c(4,5)){
    ########################################################
    #The fourth independent variable of the latent subclass
    CVdata_LC4 <- CVdatagenerate(N_LCi=N4,p_LCi=PV[S1[4],1],l_LCi=PV[S2[4],2],u_LCi=PV[S2[4],3],lamda_LCi=PV[S3[4],4],df_LCi=PV[S4[4],5],m_LCi=PV[S5[4],6],sd_LCi=PV[S5[4],7])
    Z41 <- CVdata_LC4$X1;var_Z41 <- CVdata_LC4$var_X1;m_Z41 <- CVdata_LC4$m_X1
    Z42 <- CVdata_LC4$X2;var_Z42 <- CVdata_LC4$var_X2;m_Z42 <- CVdata_LC4$m_X2
    Z43 <- CVdata_LC4$X3;var_Z43 <- CVdata_LC4$var_X3;m_Z43 <- CVdata_LC4$m_X3
    Z44 <- CVdata_LC4$X4;var_Z44 <- CVdata_LC4$var_X4;m_Z44 <- CVdata_LC4$m_X4
    Z45 <- CVdata_LC4$X5;var_Z45 <- CVdata_LC4$var_X5;m_Z45 <- CVdata_LC4$m_X5
    var_Z4 <- (PS[S12[4],1]^2)*var_Z41+(PS[S22[4],2]^2)*var_Z42+(PS[S32[4],3]^2)*var_Z43+(PS[S42[4],4]^2)*var_Z44+(PS[S52[4],5]^2)*var_Z45
    miu4 <- PS[S12[4],1]*m_Z41+PS[S22[4],2]*m_Z42+PS[S32[4],3]*m_Z43+PS[S42[4],4]*m_Z44+PS[S52[4],5]*m_Z45
    
    ###################################################################################
    ########################################################
    if(path_scenario %in% c(1,2,3,4,9,10,11,12)){
      coef_Z4toX4 <- sqrt(R2_ZtoX/(1-R2_ZtoX)*1/var_Z4)
      X4 <- coef_Z4toX4*(PS[S12[4],1]*Z41+PS[S22[4],2]*Z42+PS[S32[4],3]*Z43+PS[S42[4],4]*Z44+PS[S52[4],5]*Z45)+rnorm(N4)
      var_X4 <- (coef_Z4toX4^2)*var_Z4+1
      m_X4 <- coef_Z4toX4*miu4}
    
    ########################################################
    if(path_scenario %in% c(5,6,7,8)){
      X4 <- rnorm(n=N4, mean=5, sd=1)
      var_X4 <- 1
      m_X4 <- 5}
    
    ###################################################################################
    ########################################################
    if(path_scenario %in% c(1,2,3,4)){
      
      if(Ytype=='gaussian'){
        Y4 <- 2*X4+rnorm(N4)
        var_Y4 <- 4*var_X4+1
        m_Y4 <- 2*m_X4
      }#Ytype=='gaussian'
      
      if(Ytype=='binomial'){
        Y4s <- 2*X4+rnorm(N4)
        Y4c <- Y4s-min(Y4s)+2;pi_Y4 <- inclusionprobabilities(Y4c, N4/2);Y4 <- UPpivotal(pi_Y4)
        m_Y4 <- 1/2;var_Y4 <- 1/4
      }#Ytype=='binomial'
      
    }#path_scenario %in% c(1,2,3,4)
    ########################################################
    if(path_scenario %in% c(5,6,7,8,9,10,11,12)){
      
      OUT4 <- Y_generate(X4,Z41,Z42,Z43,Z44,Z45,var_X4,var_Z4,R2_ZtoY,coef_ZXtoY)
      var_X4Z4 <-OUT4$var_XZ
      coef_Z4toY4 <-OUT4$coef_ZtoY 
      
      if(Ytype=='gaussian'){
        Y4 <- 2*X4+coef_Z4toY4*(PS[S12[4],1]*Z41+PS[S22[4],2]*Z42+PS[S32[4],3]*Z43+PS[S42[4],4]*Z44+PS[S52[4],5]*Z45)+coef_ZXtoY*X4*(PS[S12[4],1]*Z41+PS[S22[4],2]*Z42+PS[S32[4],3]*Z43+PS[S42[4],4]*Z44+PS[S52[4],5]*Z45)+rnorm(N4)
        var_Y4 <- 4*var_X4+(coef_Z4toY4^2)*var_Z4+var_X4Z4+1
        m_Y4 <- 2*m_X4+coef_Z4toY4*miu4+coef_ZXtoY*m_X4*miu4
      }#Ytype=='gaussian'
      
      if(Ytype=='binomial'){
        Y4s <- 2*X4+coef_Z4toY4*(PS[S12[4],1]*Z41+PS[S22[4],2]*Z42+PS[S32[4],3]*Z43+PS[S42[4],4]*Z44+PS[S52[4],5]*Z45)+coef_ZXtoY*X4*(PS[S12[4],1]*Z41+PS[S22[4],2]*Z42+PS[S32[4],3]*Z43+PS[S42[4],4]*Z44+PS[S52[4],5]*Z45)+rnorm(N4)
        Y4c <- Y4s-min(Y4s)+2;pi_Y4 <- inclusionprobabilities(Y4c, N4/2);Y4 <- UPpivotal(pi_Y4)
        m_Y4 <- 1/2;var_Y4 <- 1/4
      }#Ytype=='binomial'
      
    }#path_scenario %in% c(5,6,7,8,9,10,11,12)
    
    ###################################################################################
    ########################################################
    if(path_scenario %in% c(1,5,9)){
      coef_Z4toR4 <- sqrt(R2_ZtoR/(1-R2_ZtoR)*(1)/var_Z4)
      miu_A4 <- coef_Z4toR4*miu2
      var_A4 <- (coef_Z4toR4^2)*var_Z4+1
    }
    ########################################################
    if(path_scenario %in% c(2,6,10)){
      coef_Z4toR4 <- sqrt(R2_ZtoR/(1-R2_ZtoR)*(1+(coef_XtoR^2)*var_X4)/var_Z4)
      miu_A4 <- coef_Z4toR4*miu2+coef_XtoR*m_X4
      var_A4 <- (coef_Z4toR4^2)*var_Z4+(coef_XtoR^2)*var_X4+1
    }
    ########################################################
    if(path_scenario %in% c(3,7,11)){
      coef_Z4toR4 <- sqrt(R2_ZtoR/(1-R2_ZtoR)*((coef_YtoR^2)*var_Y4+1)/var_Z4)
      miu_A4 <- coef_Z4toR4*miu2+coef_YtoR*m_Y4
      var_A4 <- (coef_Z4toR4^2)*var_Z4+(coef_YtoR^2)*var_Y4+1
    }
    ########################################################
    if(path_scenario %in% c(4,8,12)){
      coef_Z4toR4 <- sqrt(R2_ZtoR/(1-R2_ZtoR)*((coef_XtoR^2)*var_X4+1+(coef_YtoR^2)*var_Y4)/var_Z4)
      miu_A4 <- coef_Z4toR4*miu2+coef_XtoR*m_X4+coef_YtoR*m_Y4
      var_A4 <- (coef_Z4toR4^2)*var_Z4+(coef_XtoR^2)*var_X4+(coef_YtoR^2)*var_Y4+1
    }
    
  }#n_LC %in% c(4,5)
  
  if(n_LC == 4){
    ###################################################################################
    var_I <- 1/4*var_A1+1/4*var_A2+1/4*var_A3+1/4*var_A4
    var_B <- var_I/(1-icc)*icc
    a <- sqrt(2/5*var_B)
    
    b10 <- -2*a-miu_A1
    b20 <- -1*a-miu_A2
    b30 <- a-miu_A3
    b40 <- 2*a-miu_A4
    ########################################################
    if(path_scenario %in% c(1,5,9)){
      linear_A1 <- coef_Z1toR1*(PS[S12[1],1]*Z11+PS[S22[1],2]*Z12+PS[S32[1],3]*Z13+PS[S42[1],4]*Z14+PS[S52[1],5]*Z15)+rnorm(N1)+b10
      linear_A2 <- coef_Z2toR2*(PS[S12[2],1]*Z21+PS[S22[2],2]*Z22+PS[S32[2],3]*Z23+PS[S42[2],4]*Z24+PS[S52[2],5]*Z25)+rnorm(N2)+b20
      linear_A3 <- coef_Z3toR3*(PS[S12[3],1]*Z31+PS[S22[3],2]*Z32+PS[S32[3],3]*Z33+PS[S42[3],4]*Z34+PS[S52[3],5]*Z35)+rnorm(N3)+b30
      linear_A4 <- coef_Z4toR4*(PS[S12[4],1]*Z41+PS[S22[4],2]*Z42+PS[S32[4],3]*Z43+PS[S42[4],4]*Z44+PS[S52[4],5]*Z45)+rnorm(N4)+b40
    }
    ########################################################
    if(path_scenario %in% c(2,6,10)){
      linear_A1 <- coef_Z1toR1*(PS[S12[1],1]*Z11+PS[S22[1],2]*Z12+PS[S32[1],3]*Z13+PS[S42[1],4]*Z14+PS[S52[1],5]*Z15)+coef_XtoR*X1+rnorm(N1)+b10
      linear_A2 <- coef_Z2toR2*(PS[S12[2],1]*Z21+PS[S22[2],2]*Z22+PS[S32[2],3]*Z23+PS[S42[2],4]*Z24+PS[S52[2],5]*Z25)+coef_XtoR*X2+rnorm(N2)+b20
      linear_A3 <- coef_Z3toR3*(PS[S12[3],1]*Z31+PS[S22[3],2]*Z32+PS[S32[3],3]*Z33+PS[S42[3],4]*Z34+PS[S52[3],5]*Z35)+coef_XtoR*X3+rnorm(N3)+b30
      linear_A4 <- coef_Z4toR4*(PS[S12[4],1]*Z41+PS[S22[4],2]*Z42+PS[S32[4],3]*Z43+PS[S42[4],4]*Z44+PS[S52[4],5]*Z45)+coef_XtoR*X4+rnorm(N4)+b40
    }
    ########################################################
    if(path_scenario %in% c(3,7,11)){
      linear_A1 <- coef_Z1toR1*(PS[S12[1],1]*Z11+PS[S22[1],2]*Z12+PS[S32[1],3]*Z13+PS[S42[1],4]*Z14+PS[S52[1],5]*Z15)+coef_YtoR*Y1+rnorm(N1)+b10
      linear_A2 <- coef_Z2toR2*(PS[S12[2],1]*Z21+PS[S22[2],2]*Z22+PS[S32[2],3]*Z23+PS[S42[2],4]*Z24+PS[S52[2],5]*Z25)+coef_YtoR*Y2+rnorm(N2)+b20
      linear_A3 <- coef_Z3toR3*(PS[S12[3],1]*Z31+PS[S22[3],2]*Z32+PS[S32[3],3]*Z33+PS[S42[3],4]*Z34+PS[S52[3],5]*Z35)+coef_YtoR*Y3+rnorm(N3)+b30
      linear_A4 <- coef_Z4toR4*(PS[S12[4],1]*Z41+PS[S22[4],2]*Z42+PS[S32[4],3]*Z43+PS[S42[4],4]*Z44+PS[S52[4],5]*Z45)+coef_YtoR*Y4+rnorm(N4)+b40
    }
    ########################################################
    if(path_scenario %in% c(4,8,12)){
      linear_A1 <- coef_Z1toR1*(PS[S12[1],1]*Z11+PS[S22[1],2]*Z12+PS[S32[1],3]*Z13+PS[S42[1],4]*Z14+PS[S52[1],5]*Z15)+coef_XtoR*X1+coef_YtoR*Y1+rnorm(N1)+b10
      linear_A2 <- coef_Z2toR2*(PS[S12[2],1]*Z21+PS[S22[2],2]*Z22+PS[S32[2],3]*Z23+PS[S42[2],4]*Z24+PS[S52[2],5]*Z25)+coef_XtoR*X2+coef_YtoR*Y2+rnorm(N2)+b20
      linear_A3 <- coef_Z3toR3*(PS[S12[3],1]*Z31+PS[S22[3],2]*Z32+PS[S32[3],3]*Z33+PS[S42[3],4]*Z34+PS[S52[3],5]*Z35)+coef_XtoR*X3+coef_YtoR*Y3+rnorm(N3)+b30
      linear_A4 <- coef_Z4toR4*(PS[S12[4],1]*Z41+PS[S22[4],2]*Z42+PS[S32[4],3]*Z43+PS[S42[4],4]*Z44+PS[S52[4],5]*Z45)+coef_XtoR*X4+coef_YtoR*Y4+rnorm(N4)+b40
    }
    ########################################################
    linear_A1_offset <- linear_A1-min(linear_A1)+2
    pi_A1 <- inclusionprobabilities(linear_A1_offset, n_A1)
    ind_A1 <- UPpoisson(pi_A1)
    
    linear_A2_offset <- linear_A2-min(linear_A2)+2
    pi_A2 <- inclusionprobabilities(linear_A2_offset, n_A2)
    ind_A2 <- UPpoisson(pi_A2)
    
    linear_A3_offset <- linear_A3-min(linear_A3)+2
    pi_A3 <- inclusionprobabilities(linear_A3_offset, n_A3)
    ind_A3 <- UPpoisson(pi_A3)
    
    linear_A4_offset <- linear_A4-min(linear_A4)+2
    pi_A4 <- inclusionprobabilities(linear_A4_offset, n_A4)
    ind_A4 <- UPpoisson(pi_A4)
  }#LC=3
  
  #The number of heterogeneous populations is【5】----------
  if(n_LC == 5){
    #The fifth independent variable of the latent subclass
    ########################################################
    CVdata_LC5 <- CVdatagenerate(N_LCi=N5,p_LCi=PV[S1[5],1],l_LCi=PV[S2[5],2],u_LCi=PV[S2[5],3],lamda_LCi=PV[S3[5],4],df_LCi=PV[S4[5],5],m_LCi=PV[S5[5],6],sd_LCi=PV[S5[5],7])
    Z51 <- CVdata_LC5$X1;var_Z51 <- CVdata_LC5$var_X1;m_Z51 <- CVdata_LC5$m_X1
    Z52 <- CVdata_LC5$X2;var_Z52 <- CVdata_LC5$var_X2;m_Z52 <- CVdata_LC5$m_X2
    Z53 <- CVdata_LC5$X3;var_Z53 <- CVdata_LC5$var_X3;m_Z53 <- CVdata_LC5$m_X3
    Z54 <- CVdata_LC5$X4;var_Z54 <- CVdata_LC5$var_X4;m_Z54 <- CVdata_LC5$m_X4
    Z55 <- CVdata_LC5$X5;var_Z55 <- CVdata_LC5$var_X5;m_Z55 <- CVdata_LC5$m_X5
    var_Z5 <- (PS[S12[5],1]^2)*var_Z51+(PS[S22[5],2]^2)*var_Z52+(PS[S32[5],3]^2)*var_Z53+(PS[S42[5],4]^2)*var_Z54+(PS[S52[5],5]^2)*var_Z55
    miu5 <- PS[S12[5],1]*m_Z51+PS[S22[5],2]*m_Z52+PS[S32[5],3]*m_Z53+PS[S42[5],4]*m_Z54+PS[S52[5],5]*m_Z55
    
    ###################################################################################
    ########################################################
    if(path_scenario %in% c(1,2,3,4,9,10,11,12)){
      coef_Z5toX5 <- sqrt(R2_ZtoX/(1-R2_ZtoX)*1/var_Z5)
      X5 <- coef_Z5toX5*(PS[S12[5],1]*Z51+PS[S22[5],2]*Z52+PS[S32[5],3]*Z53+PS[S42[5],4]*Z54+PS[S52[5],5]*Z55)+rnorm(N5)
      var_X5 <- (coef_Z5toX5^2)*var_Z5+1
      m_X5 <- coef_Z5toX5*miu5}
    ########################################################
    if(path_scenario %in% c(5,6,7,8)){
      X5 <- rnorm(n=N5, mean=5, sd=1)
      var_X5 <- 1
      m_X5 <- 5}
    
    ###################################################################################
    ########################################################
    if(path_scenario %in% c(1,2,3,4)){
      
      if(Ytype=='gaussian'){
        Y5 <- 2*X5+rnorm(N5)
        var_Y5 <- 4*var_X5+1
        m_Y5 <- 2*m_X5
      }#Ytype=='gaussian'
      
      if(Ytype=='binomial'){
        Y5s <- 2*X5+rnorm(N5)
        Y5c <- Y5s-min(Y5s)+2;pi_Y5 <- inclusionprobabilities(Y5c, N5/2);Y5 <- UPpivotal(pi_Y5)
        m_Y5 <- 1/2;var_Y5 <- 1/4
      }#Ytype=='binomial'
      
    }#path_scenario %in% c(5,6,7,8,9,10,11,12)
    ########################################################
    if(path_scenario %in% c(5,6,7,8,9,10,11,12)){
      
      OUT5 <- Y_generate(X5,Z51,Z52,Z53,Z54,Z55,var_X5,var_Z5,R2_ZtoY,coef_ZXtoY)
      var_X5Z5 <-OUT5$var_XZ
      coef_Z5toY5 <-OUT5$coef_ZtoY 
      
      if(Ytype=='gaussian'){
        Y5 <- 2*X5+coef_Z5toY5*(PS[S12[5],1]*Z51+PS[S22[5],2]*Z52+PS[S32[5],3]*Z53+PS[S42[5],4]*Z54+PS[S52[5],5]*Z55)+coef_ZXtoY*X5*(PS[S12[5],1]*Z51+PS[S22[5],2]*Z52+PS[S32[5],3]*Z53+PS[S42[5],4]*Z54+PS[S52[5],5]*Z55)+rnorm(N5)
        var_Y5 <- 4*var_X5+(coef_Z5toY5^2)*var_Z5+var_X5Z5+1
        m_Y5 <- 2*m_X5+coef_Z5toY5*miu5+coef_ZXtoY*m_X5*miu5
      }#Ytype=='gaussian'
      
      if(Ytype=='binomial'){
        Y5s <- 2*X5+coef_Z5toY5*(PS[S12[5],1]*Z51+PS[S22[5],2]*Z52+PS[S32[5],3]*Z53+PS[S42[5],4]*Z54+PS[S52[5],5]*Z55)+coef_ZXtoY*X5*(PS[S12[5],1]*Z51+PS[S22[5],2]*Z52+PS[S32[5],3]*Z53+PS[S42[5],4]*Z54+PS[S52[5],5]*Z55)+rnorm(N5)
        Y5c <- Y5s-min(Y5s)+2;pi_Y5 <- inclusionprobabilities(Y5c, N5/2);Y5 <- UPpivotal(pi_Y5)
        m_Y5 <- 1/2;var_Y5 <- 1/4
      }#Ytype=='binomial'
      
    }
    
    ###################################################################################
    ########################################################
    if(path_scenario %in% c(1,5,9)){
      coef_Z5toR5 <- sqrt(R2_ZtoR/(1-R2_ZtoR)*(1)/var_Z5)
      miu_A5 <- coef_Z5toR5*miu2
      var_A5 <- (coef_Z5toR5^2)*var_Z5+1
    }
    ########################################################
    if(path_scenario %in% c(2,6,10)){
      coef_Z5toR5 <- sqrt(R2_ZtoR/(1-R2_ZtoR)*(1+(coef_XtoR^2)*var_X5)/var_Z5)
      miu_A5 <- coef_Z5toR5*miu2+coef_XtoR*m_X5
      var_A5 <- (coef_Z5toR5^2)*var_Z5+(coef_XtoR^2)*var_X5+1
    }
    ########################################################
    if(path_scenario %in% c(3,7,11)){
      coef_Z5toR5 <- sqrt(R2_ZtoR/(1-R2_ZtoR)*((coef_YtoR^2)*var_Y5+1)/var_Z5)
      miu_A5 <- coef_Z5toR5*miu2+coef_YtoR*m_Y5
      var_A5 <- (coef_Z5toR5^2)*var_Z5+(coef_YtoR^2)*var_Y5+1
    }
    ########################################################
    if(path_scenario %in% c(4,8,12)){
      coef_Z5toR5 <- sqrt(R2_ZtoR/(1-R2_ZtoR)*((coef_XtoR^2)*var_X5+1+(coef_YtoR^2)*var_Y5)/var_Z5)
      miu_A5 <- coef_Z5toR5*miu2+coef_XtoR*m_X5+coef_YtoR*m_Y5
      var_A5 <- (coef_Z5toR5^2)*var_Z5+(coef_XtoR^2)*var_X5+(coef_YtoR^2)*var_Y5+1
    }
    
    ###################################################################################
    var_I <- 1/5*var_A1+1/5*var_A2+1/5*var_A3+1/5*var_A4+1/5*var_A5
    var_B <- var_I/(1-icc)*icc
    a <- sqrt(1/2*var_B)
    
    b10 <- -2*a-miu_A1
    b20 <- -1*a-miu_A2
    b30 <- -miu_A3
    b40 <- a-miu_A4
    b50 <- 2*a-miu_A5
    ########################################################
    if(path_scenario %in% c(1,5,9)){
      linear_A1 <- coef_Z1toR1*(PS[S12[1],1]*Z11+PS[S22[1],2]*Z12+PS[S32[1],3]*Z13+PS[S42[1],4]*Z14+PS[S52[1],5]*Z15)+rnorm(N1)+b10
      linear_A2 <- coef_Z2toR2*(PS[S12[2],1]*Z21+PS[S22[2],2]*Z22+PS[S32[2],3]*Z23+PS[S42[2],4]*Z24+PS[S52[2],5]*Z25)+rnorm(N2)+b20
      linear_A3 <- coef_Z3toR3*(PS[S12[3],1]*Z31+PS[S22[3],2]*Z32+PS[S32[3],3]*Z33+PS[S42[3],4]*Z34+PS[S52[3],5]*Z35)+rnorm(N3)+b30
      linear_A4 <- coef_Z4toR4*(PS[S12[4],1]*Z41+PS[S22[4],2]*Z42+PS[S32[4],3]*Z43+PS[S42[4],4]*Z44+PS[S52[4],5]*Z45)+rnorm(N4)+b40
      linear_A5 <- coef_Z5toR5*(PS[S12[5],1]*Z51+PS[S22[5],2]*Z52+PS[S32[5],3]*Z53+PS[S42[5],4]*Z54+PS[S52[5],5]*Z55)+rnorm(N5)+b50
    }
    ########################################################
    if(path_scenario %in% c(2,6,10)){
      linear_A1 <- coef_Z1toR1*(PS[S12[1],1]*Z11+PS[S22[1],2]*Z12+PS[S32[1],3]*Z13+PS[S42[1],4]*Z14+PS[S52[1],5]*Z15)+coef_XtoR*X1+rnorm(N1)+b10
      linear_A2 <- coef_Z2toR2*(PS[S12[2],1]*Z21+PS[S22[2],2]*Z22+PS[S32[2],3]*Z23+PS[S42[2],4]*Z24+PS[S52[2],5]*Z25)+coef_XtoR*X2+rnorm(N2)+b20
      linear_A3 <- coef_Z3toR3*(PS[S12[3],1]*Z31+PS[S22[3],2]*Z32+PS[S32[3],3]*Z33+PS[S42[3],4]*Z34+PS[S52[3],5]*Z35)+coef_XtoR*X3+rnorm(N3)+b30
      linear_A4 <- coef_Z4toR4*(PS[S12[4],1]*Z41+PS[S22[4],2]*Z42+PS[S32[4],3]*Z43+PS[S42[4],4]*Z44+PS[S52[4],5]*Z45)+coef_XtoR*X4+rnorm(N4)+b40
      linear_A5 <- coef_Z5toR5*(PS[S12[5],1]*Z51+PS[S22[5],2]*Z52+PS[S32[5],3]*Z53+PS[S42[5],4]*Z54+PS[S52[5],5]*Z55)+coef_XtoR*X5+rnorm(N5)+b50
    }
    ########################################################
    if(path_scenario %in% c(3,7,11)){
      linear_A1 <- coef_Z1toR1*(PS[S12[1],1]*Z11+PS[S22[1],2]*Z12+PS[S32[1],3]*Z13+PS[S42[1],4]*Z14+PS[S52[1],5]*Z15)+coef_YtoR*Y1+rnorm(N1)+b10
      linear_A2 <- coef_Z2toR2*(PS[S12[2],1]*Z21+PS[S22[2],2]*Z22+PS[S32[2],3]*Z23+PS[S42[2],4]*Z24+PS[S52[2],5]*Z25)+coef_YtoR*Y2+rnorm(N2)+b20
      linear_A3 <- coef_Z3toR3*(PS[S12[3],1]*Z31+PS[S22[3],2]*Z32+PS[S32[3],3]*Z33+PS[S42[3],4]*Z34+PS[S52[3],5]*Z35)+coef_YtoR*Y3+rnorm(N3)+b30
      linear_A4 <- coef_Z4toR4*(PS[S12[4],1]*Z41+PS[S22[4],2]*Z42+PS[S32[4],3]*Z43+PS[S42[4],4]*Z44+PS[S52[4],5]*Z45)+coef_YtoR*Y4+rnorm(N4)+b40
      linear_A5 <- coef_Z5toR5*(PS[S12[5],1]*Z51+PS[S22[5],2]*Z52+PS[S32[5],3]*Z53+PS[S42[5],4]*Z54+PS[S52[5],5]*Z55)+coef_YtoR*Y5+rnorm(N5)+b50
    }
    ########################################################
    if(path_scenario %in% c(4,8,12)){
      linear_A1 <- coef_Z1toR1*(PS[S12[1],1]*Z11+PS[S22[1],2]*Z12+PS[S32[1],3]*Z13+PS[S42[1],4]*Z14+PS[S52[1],5]*Z15)+coef_XtoR*X1+coef_YtoR*Y1+rnorm(N1)+b10
      linear_A2 <- coef_Z2toR2*(PS[S12[2],1]*Z21+PS[S22[2],2]*Z22+PS[S32[2],3]*Z23+PS[S42[2],4]*Z24+PS[S52[2],5]*Z25)+coef_XtoR*X2+coef_YtoR*Y2+rnorm(N2)+b20
      linear_A3 <- coef_Z3toR3*(PS[S12[3],1]*Z31+PS[S22[3],2]*Z32+PS[S32[3],3]*Z33+PS[S42[3],4]*Z34+PS[S52[3],5]*Z35)+coef_XtoR*X3+coef_YtoR*Y3+rnorm(N3)+b30
      linear_A4 <- coef_Z4toR4*(PS[S12[4],1]*Z41+PS[S22[4],2]*Z42+PS[S32[4],3]*Z43+PS[S42[4],4]*Z44+PS[S52[4],5]*Z45)+coef_XtoR*X4+coef_YtoR*Y4+rnorm(N4)+b40
      linear_A5 <- coef_Z5toR5*(PS[S12[5],1]*Z51+PS[S22[5],2]*Z52+PS[S32[5],3]*Z53+PS[S42[5],4]*Z54+PS[S52[5],5]*Z55)+coef_XtoR*X5+coef_YtoR*Y5+rnorm(N5)+b50
    }
    ########################################################
    linear_A1_offset <- linear_A1-min(linear_A1)+2
    pi_A1 <- inclusionprobabilities(linear_A1_offset, n_A1)
    ind_A1 <- UPpoisson(pi_A1)
    
    linear_A2_offset <- linear_A2-min(linear_A2)+2
    pi_A2 <- inclusionprobabilities(linear_A2_offset, n_A2)
    ind_A2 <- UPpoisson(pi_A2)
    
    linear_A3_offset <- linear_A3-min(linear_A3)+2
    pi_A3 <- inclusionprobabilities(linear_A3_offset, n_A3)
    ind_A3 <- UPpoisson(pi_A3)
    
    linear_A4_offset <- linear_A4-min(linear_A4)+2
    pi_A4 <- inclusionprobabilities(linear_A4_offset, n_A4)
    ind_A4 <- UPpoisson(pi_A4)
    
    linear_A5_offset <- linear_A5-min(linear_A5)+2
    pi_A5 <- inclusionprobabilities(linear_A5_offset, n_A5)
    ind_A5 <- UPpoisson(pi_A5)
    
  }#LC=5
  
  
  if(n_LC==1){
    Z1 <- c(Z11)
    Z2 <- c(Z12)
    Z3 <- c(Z13)
    Z4 <- c(Z14)
    Z5 <- c(Z15)
    X <- c(X1)
    Y <- c(Y1)
    pi_A <- c(pi_A1)
    ind_A <- c(ind_A1)}
  
  if(n_LC==2){
    Z1 <- c(Z11,Z21)
    Z2 <- c(Z12,Z22)
    Z3 <- c(Z13,Z23)
    Z4 <- c(Z14,Z24)
    Z5 <- c(Z15,Z25)
    X <- c(X1,X2)
    Y <- c(Y1,Y2)
    pi_A <- c(pi_A1,pi_A2)
    ind_A <- c(ind_A1,ind_A2)}
  
  if(n_LC==3){
    Z1 <- c(Z11,Z21,Z31)
    Z2 <- c(Z12,Z22,Z32)
    Z3 <- c(Z13,Z23,Z33)
    Z4 <- c(Z14,Z24,Z34)
    Z5 <- c(Z15,Z25,Z35)
    X <- c(X1,X2,X3)
    Y <- c(Y1,Y2,Y3)
    pi_A <- c(pi_A1,pi_A2,pi_A3)
    ind_A <- c(ind_A1,ind_A2,ind_A3)}
  
  if(n_LC==4){
    Z1 <- c(Z11,Z21,Z31,Z41)
    Z2 <- c(Z12,Z22,Z32,Z42)
    Z3 <- c(Z13,Z23,Z33,Z43)
    Z4 <- c(Z14,Z24,Z34,Z44)
    Z5 <- c(Z15,Z25,Z35,Z45)
    X <- c(X1,X2,X3,X4)
    Y <- c(Y1,Y2,Y3,Y4)
    pi_A <- c(pi_A1,pi_A2,pi_A3,pi_A4)
    ind_A <- c(ind_A1,ind_A2,ind_A3,ind_A4)}
  
  if(n_LC==5){
    Z1 <- c(Z11,Z21,Z31,Z41,Z51)
    Z2 <- c(Z12,Z22,Z32,Z42,Z52)
    Z3 <- c(Z13,Z23,Z33,Z43,Z53)
    Z4 <- c(Z14,Z24,Z34,Z44,Z54)
    Z5 <- c(Z15,Z25,Z35,Z45,Z55)
    X <- c(X1,X2,X3,X4,X5)
    Y <- c(Y1,Y2,Y3,Y4,Y5)
    pi_A <- c(pi_A1,pi_A2,pi_A3,pi_A4,pi_A5)
    ind_A <- c(ind_A1,ind_A2,ind_A3,ind_A4,ind_A5)}
  
  ###################################################################################Generates probability samples via random sampling
  pi_B <- rep(n_B/N,times=N)
  sample_B <- sample(1:N,n_B)
  ind_B <- rep(0,times=N)
  ind_B[sample_B] <- 1
  
  ###################################################################################Organize all variables into a single dataset
  pop_data <- data.frame(Z1,Z2,Z3,Z4,Z5,X,Y,pi_A,pi_B,ind_A,ind_B)
  
  
  #把非概率样本、参考样本分别抠出来，再合并成一个大样本
  data_A0 <- pop_data[ind_A==1,]
  data_B0 <- pop_data[ind_B==1,]
  
  data_B1 <- data_B0 %>% mutate(sample_weight=1/pi_B)

  output <- list('nonprobability_sample'=data_A0,
                 'reference_sample'=data_B1)
  return(output)
  
}#simu_data_NHEI

