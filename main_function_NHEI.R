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

###########################################################################################################################################################
###########################################Main function
###########################################################################################################################################################
# Parameter Description                                                                                                                                   #
# data_NP: Nonprobability sample                                                                                                                          #
# data_P: Reference sample (Probability sample)                                                                                                           #
# CIV: Names of common variables between nonprobability sample and reference sample                                                                       #
# type_CIV:Types of common variables between nonprobability sample and reference sample, possible values include 'gaussian','binomial', and 'multinomial' #
# sample_weight: Sample weights for reference sample                                                                                                      #
# maxclassnum:Maximum number of latent classes                                                                                                            #
###########################################################################################################################################################
NPasso_ALC <-
    function(data_NP=data_A,data_P=data_B,
             CIV=c('W1','W2'),type_CIV=c('gaussian','gaussian'),
             sample_weight='sample_weight',maxclassnum=5)
  {
    #################################################################################################
    ##################################################################################Preprocessing the data-----
    #################################################################################################
    
    #计算校正因子
    a <- 1-nrow(data_P)/sum(data_P[,sample_weight])
    
    #计算参考样本校正后的抽样权重
    SW_P <- data_P[,sample_weight]*a
    OSW_P <- data_P[,sample_weight]
    data_P <- cbind(data_P,SW=SW_P,OSW=OSW_P)
    
    #定义非概率样本的两个权重变量
    SW_NP <- rep(1,times=nrow(data_NP))
    OSW_NP <- rep(1,times=nrow(data_NP))
    data_NP <- cbind(data_NP,SW=SW_NP,OSW=OSW_NP)
    
    #生成表示id的变量
    id <- 1:length(nrow(data_NP)+nrow(data_P))
    #合并两个样本以及id列,生成表示样本来源的字符型变量datatype和数值型变量group(1-NP;0-P)
    data_com <- cbind(rbind.fill(cbind(data_NP,group=1),cbind(data_P,group=0)),id)
    
    #################################################################################################
    ################################################################Indentifying the latent class-------
    #################################################################################################
    data_com1 <- data_com
    for (CIV_i in 1:length(CIV)) {
      if(type_CIV[CIV_i]=='binomial'){ data_com1[,CIV[CIV_i]] <- as.numeric(factor(data_com1[,CIV[CIV_i]]))-1 }
      
      if(type_CIV[CIV_i]=='multinomial'){ data_com1[,CIV[CIV_i]] <- as.numeric(factor(data_com1[,CIV[CIV_i]])) }
      
      if(type_CIV[CIV_i]=='gaussian'){ data_com1[,CIV[CIV_i]] <- data_com1[,CIV[CIV_i]] }
    }
    
    #用于潜分类分析的表达式
    expr_variable <- paste('list(',paste(paste(CIV,'~1',sep=''),collapse = ','),')',sep='')
    expr_type <- paste('list(',paste(paste(type_CIV,'()',sep=''),collapse = ','),')',sep='')
    
    #Identification of latent class: the BIC vector is used to store the BIC values for each model with specific latent class number
    BIC <- NULL;results_FM <- NULL
    for (i_class in 1:maxclassnum) {
      lca <- depmixS4::mix(eval(str2lang(expr_variable)),
                           data=data_com1,
                           nstates = i_class,
                           family = eval(str2lang(expr_type)))
      Ind_lca1 <- try(fit <- fit(lca),silent=TRUE)
      
      T <- 0
      while(methods::is(Ind_lca1,"try-error")==TRUE & T<=10) {
        lca <- depmixS4::mix(eval(str2lang(expr_variable)),
                             data=data_com1,
                             nstates = i_class,
                             family = eval(str2lang(expr_type)))
        Ind_lca1 <- try(fit <- fit(lca),silent=TRUE)
        T <- T+1
      }
      if(methods::is(Ind_lca1,"try-error")==TRUE){stop('The finite mixture model cannot be fitted')}
      
      BIC[i_class] <- BIC(fit)
      results_FM <- cbind(results_FM,fit@posterior$state)
    }
    bestnum_vec <- which(BIC==min(BIC))
    if(length(bestnum_vec)>=1){bestnum <- bestnum_vec[1]}
    if(length(bestnum_vec)==1){bestnum <- bestnum_vec}
    LC <- results_FM[,bestnum]
    
    #################################################################################################
    ################################################################Construct a model based on the identified latent subclasses.-------
    #################################################################################################
    
    #Convert multi-class categories into dummy variables.
    data_com2 <- cbind(data_com,LC)
    CIV_2 <- NULL
    
    for (CIV_i in 1:length(CIV)){
      
      if(type_CIV[CIV_i]=='binomial'){ 
        data_com2[,CIV[CIV_i]] <- as.numeric(factor(data_com2[,CIV[CIV_i]]))-1 
        CIV_2 <- c(CIV_2,CIV[CIV_i])
        }
      
      if(type_CIV[CIV_i]=='multinomial'){ 
        dummy_var <- model.matrix(~ data_com2[,CIV[CIV_i]] - 1, data = data_com2)
        dummy_num <- length(unique(data_com2[,CIV[CIV_i]]))
        dummy_name <- paste(CIV[CIV_i],1:dummy_num,sep='')
        colnames(dummy_var) <- dummy_name
        data_com2 <- cbind(data_com2,dummy_var)
        CIV_2 <- c(CIV_2,dummy_name)
        }
      
      if(type_CIV[CIV_i]=='gaussian'){ 
        data_com2[,CIV[CIV_i]] <- data_com2[,CIV[CIV_i]] 
        CIV_2 <- c(CIV_2,CIV[CIV_i])
        }
      
    }
    
    
    #When the optimal number of latent classes is 1, there is no need to account for potential heterogeneity; one can simply use standard regression.
    if(bestnum==1){
      data_NP2 <- data_com2[data_com2$group==1,]
      data_P2 <- data_com2[data_com2$group==0,]
      
      #Model Expression and Construction of Inverse Probability Weighting Models
      formula_myipw <- paste('group~',paste(CIV_2,collapse = '+'),sep='')
      fit_toNP <- glm(eval(str2lang(formula_myipw)), data=data_com2, family='binomial', weights = SW)
    }#bestnum==1
    
    #If the optimal number of latent subgroups exceeds one, a multilevel model is employed to account for heterogeneity; finally, `fit_toNP` and `fit_Y` are generated for subsequent calculations.
    if(bestnum!=1){
      data_NP2 <- data_com2[data_com2$group==1,]
      data_P2 <- data_com2[data_com2$group==0,]
      
      #Model Expression for Inverse Probability Weighting
      formula_myipw <- paste('group~(',paste(CIV,collapse = '+'),'|LC)',sep='')
      suppressWarnings(fit_toNP <- glmer(eval(str2lang(formula_myipw)), data=data_com2, family=binomial, weights = SW))
    }#bestnum!=1
    
    ########################################################################################################################
    ################################################################Improving Weights in Inverse Probability Weighting------
    ########################################################################################################################
    
    #First, calculate the propensity score for inclusion in the sample for each individual in the combined sample
    ps_BOTH <- predict(fit_toNP, newdata=data_com2, type="link")
    #Extracting Propensity Scores for Non-Probability Samples
    ps_NP <- ps_BOTH[data_com2$group==1]
    #Extracting Propensity Scores for Reference Samples
    ps_P <- ps_BOTH[data_com2$group==0]
    #Sampling Weights of Reference Samples (after adjustment)
    sweight_P <- data_com2[data_com2$group==0,'SW']
    #Calculating Weights for Non-Probability Samples After Kernel Smoothing: Standard Normal Transformation
    fit_kernal1 <- kw.wt(p_score.c = ps_NP,p_score.s = ps_P,svy.wt = sweight_P,krn = "dnorm") 
    kwweight_dorm <- fit_kernal1$pswt

    
    ################################################################将所计算的点估计值和方差输出
    output <- kwweight_dorm
    return(output)
    
  }#function:NPasso_ALC

