---
# title: "GLM_A&C"
#author: "Rafael S. de Souza, ..."
#date: "27 de junho de 2014"
#output: pdf_document
---
# 18/08/2014
# Update by Madhura Killedar

# Required libraries

require(MASS)
require(grDevices)
require(plyr)
require(msme)
require(ggplot2)
require(boot)
require(ggthemes)
library(scales)
require(MCMCpack)
require(arm)
library(pROC)
require(plyr)
require(reshape)
library(COUNT)
library(SDMTools)
library(CosmoGLM)
require(stargazer)

#Reading the data 

#Set working directory to the data folder (replace by your own directory)

#data_path<-"/Users/killedar/Documents/Research/projects/CosmoStatsIAA/GLM/PopIIIstars/"
data_path<-"/Users/rafael/Dropbox/artigos/Meusartigos/IAA-WGC/GLMs/Simulation/data/"


Biffi_data<-read.table(file=paste(data_path,"Biffi2014.csv",sep=""),
                       header=TRUE)



#Biffi_data_original<-Biffi_data
# Problem 1: xmol, Z, SFR. SFR is the response variable
Biffi_data<-Biffi_data[,c("SFR","Xmol","Z")]





#Transforming variable into numeric (required by GLM packages) 
Biffi_data$SFR[which(Biffi_data$SFR==0)]<-0
Biffi_data$SFR[which(Biffi_data$SFR>0)]<-1






Biffi_data$SFR<-as.numeric(Biffi_data$SFR)
#Biffi_data$Z<-scale(Biffi_data$Z)#Scaling and centering
#Biffi_data$Xmol<-scale(Biffi_data$Xmol)#Scaling and centering

#------------------------------------------------------#

Fglm <-glm(SFR~Z+Xmol,family=binomial(link="logit"),
           data = Biffi_data)


new.data <- expand.grid(Xmol = seq(min(Biffi_data$Xmol),max(Biffi_data$Xmol), length = 10000), 
                        Z = mean(Biffi_data$Z))


preds <- predict(Fglm, newdata = new.data, type = 'response',se = TRUE)
new.data$SFR <- preds$fit

new.data$ymin <- new.data$SFR - 2*preds$se.fit
new.data$ymax <- new.data$SFR + 2*preds$se.fit 

ggplot(new.data ,aes(x = Xmol, y = SFR)) + 
  geom_point() + 
  geom_line(data = new.data,aes(y = SFR),colour = "blue")


dataset <- expand.grid(Temp = rnorm(30), Age = runif(10))
pretty(dataset$Age, 10)


#Maximum Likelihood Logit 

Fglm <-glm(SFR~Z+Xmol,family=binomial(link="logit"),
                 data = Biffi_data)
#pdf("FLogitSF.pdf",height=8,width=9)
#plot_Prob(Fglm)+ylab("Predicted Probabilities for SF")+xlab("")+
#ggtitle("Frequentist Logistic Regression")+coord_cartesian(ylim = c(0,1.05))
#dev.off()




with(Fglm, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = FALSE))

#------------------------------------------------------#

#Bayesian logit 
# Cauchy prior is default
#pdf("BLogitSF1.pdf",height=8,width=9)
Bglm <- bayesglm(SFR~Z+Xmol,family=binomial(link="logit"),
                 data = Biffi_data)



plot_Prob(Bglm)+ylab("Predicted Probabilities for SF")+xlab("")+
    ggtitle("Bayesian GLM with Logistic Link")+coord_cartesian(ylim = c(0,1.05))
summary(Bglm) 

with(Bglm, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = FALSE))
#dev.off()






X <- model.matrix(~Z+Xmol,
                  data = Biffi_data)
eta<-X%*%coef(Bglm)
Biffi_data$pi <- exp(eta) / (1 + exp(eta))







#Bayesian probit 
# Cauchy prior is default
#pdf("BLogitSF1.pdf",height=8,width=9)
Bprobitglm <- bayesglm(SFR~Xmol,family=binomial(link="probit"),
                 data = Biffi_data)
plot_Prob(Bprobitglm)+ylab("Predicted Probabilities for SF")+xlab("")+
  ggtitle("Bayesian GLM with Probit Link")+coord_cartesian(ylim = c(0,1.05))
with(Bglm, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = FALSE))
#dev.off()












#------------------------------------------------------#
# Bayesian GLM
# Bayesian with Cauchy prior plot posterior
logpriorfun <- function(beta){
  sum(dcauchy(beta, log=TRUE))
}
MCMC.glm<-MCMClogit(SFR~Xmol+Z,
                    data = Biffi_data,user.prior.density=logpriorfun,
                    logfun=TRUE,burnin=100,mcmc = 10000)

#pdf("BLogitSF.pdf",height=8,width=9)
plotProbmcmc(MCMC.glm)+ylab("Predicted Probabilities for SF")+xlab("")+
  ggtitle("Bayesian Logistic Regression 2")+coord_cartesian(ylim = c(0,1.05))
#annotate("text", x =1.5, y=0.15, label=paste("AIC: ",round(modelfit(Fglm)$AIC,2),sep=""),size=6)+  
#annotate("text", x =1.5, y=0.05, label=paste("AICn: ",round(modelfit(Fglm)$AICn,2),sep=""),size=6)+  
# annotate("text", x =0.5, y=0.15, label=paste("BIC: ",round(modelfit(Fglm)$BIC,2),sep=""),size=6)+
#annotate("text", x =0.5, y=0.05, label=paste("BICqh: ",round(modelfit(Fglm)$BICqh,2),sep=""),size=6)  
#dev.off()
MCMC.glm<-MCMCprobit(SFR~Xmol,
                    data = Biffi_data,user.prior.density=logpriorfun,
                    logfun=TRUE,burnin=100,mcmc = 10000)


#------------------------------------------------------#
# ROC curve
#Frequentist




tc <- trainControl("cv", 10, savePredictions=T, classProbs = TRUE, summaryFunction = twoClassSummary)  #"cv" = cross-validation, 10-fold

fit <- train(SFR~Xmol+Z,
             data      = Biffi_data   ,
             method    = "glm"    ,
             family    = binomial(link="logit") ,
             metric = "ROC",
             trControl = tc)

Fold1<-subset(fit$pred,Resample=="Fold05")
F1 <-roc(Fold1$obs,Fold1[,4])
F1


coords(F1,x="best")[1]
confusion.matrix(Fold1$obs,Fold1[,4], threshold = coords(F1,x="best")[1])




inTraining <- createDataPartition(Biffi_data$SFR, p = 0.75,list=FALSE)
training <- Biffi_data[inTraining, ]
testing <- Biffi_data[-inTraining, ]
Fglm <-glm(SFR~Z+Xmol,family=binomial(link="logit"),
           data = training)

ROCF<- data.frame(True=training$SFR,predicted=predict(Fglm, newdata=training,type = "response"))

F1 <-roc(ROCF$True,ROCF$predicted)
coords(F1,x="best")[1]
confusion.matrix(ROCF$True, ROCF$predicted, threshold = coords(F1,x="best")[1])

ROCF2<- data.frame(True=testing$SFR,predicted=predict(Fglm, newdata=testing,type = "response"))

F2 <-roc(ROCF2$True,ROCF2$predicted)
coords(F1,x="best")[1]
cm2<-confusion.matrix(ROCF2$True, ROCF2$predicted, threshold = coords(F1,x="best")[1])

ROCF2$class<-ROCF2$predicted
ROCF2$class[which(ROCF2$class>=coords(F1,x="best")[1])]<-1
ROCF2$class[which(ROCF2$class<coords(F1,x="best")[1])]<-0
#ROCF2$class[which(ROCF2$class==1)]<-"SF"
#ROCF2$class[which(ROCF2$class==0)]<-"No SF"

omission(cm2)
sensitivity(cm2)
specificity(cm2)
prop.correct(cm2)

require(mlearning)
cm2<-confusion(factor(ROCF2$True),factor(ROCF2$class),labels=c("Actual","Predicted"),names=c("No SF","SF"))
plot(cm2)





F1gg<-data.frame(Sensitivity=F1$sensitivities,Specificity=1-F1$specificities)
F2gg<-data.frame(Sensitivity=F2$sensitivities,Specificity=1-F2$specificities)
#pdf("ROC_GLM.pdf")
ggplot(data=F1gg,aes(x=Specificity,y=Sensitivity))+geom_line(size=1.5,color="blue4")+
  geom_line(data=F2gg,aes(x=Specificity,y=Sensitivity),size=1.5,color="green4")+
  theme_economist_white(gray_bg = F, base_size = 11, base_family = "sans")+
  theme(plot.title = element_text(hjust=0.5),axis.title.y=element_text(vjust=0.75),axis.title.x=element_text(vjust=-0.25),
        text = element_text(size=20))+ggtitle("ROC Curve GLM")+
  annotate("text",size=8, x = 0.6, y = 0.4,color="green4", label = paste("AUC: ",round(F1$auc[1],2),sep="") ) +
  annotate("text",size=8, x = 0.6, y = 0.3,color="blue4", label = paste("AUC: ",round(F2$auc[1],2),sep="") ) +
  geom_segment(x=0,y=0,xend=1,yend=1,colour="red",size=1)+
  xlab("1-Specificity")
#dev.off()


     


#------------------------------------------------------#












#myplot<-ggplot(posterior, aes(x=value),group=X2,fill=X2,color=X2)+geom_density(data=posterior,aes(fill=Var2))+
 # scale_color_stata(guide="none")+scale_fill_stata(guide="none")+
#  theme_economist_white(gray_bg = F, base_size = 11, base_family = "sans")+
 # xlab("Posterior for the Odds Ratio")+ylab("Density")+
  #theme(plot.title = element_text(hjust=0.5),axis.title.y=element_text(vjust=0.75),axis.title.x=element_text(vjust=-0.25),
   #     text = element_text(size=20))+
  #facet_wrap(~Var2,ncol=1,scales="free")
#plot2<-facet_wrap_labeller(myplot, labels = c(expression(x[mol1]), expression(Z[1]), expression(Z[2])))

#pdf("Odss_MCMCglm.pdf",width=8,height=9) 
#plot2
#dev.off() 

#------------------------------------------------------#


#------------------------------------------------------#



ggdata<-Biffi_data
ggdata$SFR<-as.factor(ggdata$SFR)
ggdata$SFR<-revalue(ggdata$SFR, c("0"="No SF", "1"="SF"))


# Scatterplot
#pdf("Scatter-xmolZ.pdf", width=4, height=3.5)
#scat <- ggplot(data=ggdata,aes(x=Xmol,y=1e-9+Z))

asinh_trans <- function(){
  trans_new(name = 'asinh', transform = function(x) asinh(x), 
            inverse = function(x) sinh(x))
}

cairo_pdf("Scatter-xmolZ.pdf", height=8,width=9)
scat <- ggplot(data=ggdata,aes(x=Xmol,y=asinh(1e8*Z)))
scat <- scat + xlab(expression(X[mol])) + 
  ylab( expression(ArcSinh(Z%*%10^8/Z['\u0298'])))

scat <- scat + geom_point(size=3,aes(color=SFR,shape=SFR))
#scat <- scat + geom_point(aes(color=SFR,shape=factor(SFR))
#scat <- scat + scale_color_discrete(name="SFR",labels=c("SFR=0","SFR>0"))
scat <- scat + theme(legend.position=c(0.85,0.15)) 
scat <- scat + theme(legend.background=element_rect(fill="white"))
scat <- scat + theme(panel.background =element_rect(fill="grey75"))

scat+theme_economist_white(gray_bg = F, base_size = 11, base_family = "sans")+
  theme(plot.title = element_text(hjust=0.5),axis.title.y=element_text(vjust=0.75),axis.title.x=element_text(vjust=-0.25),
        text = element_text(size=20))+scale_colour_solarized(name="")+scale_shape_tremmel(name="")+
  scale_x_continuous(trans = 'log10',
                     breaks=trans_breaks("log10",function(x) 10^x),
                     labels=trans_format("log10",math_format(10^.x)))
#  scale_y_continuous(expand=c(0.1,0),trans = 'log10',
 #                    breaks=trans_breaks("log10",function(x) 10^x),
  #                   labels=trans_format("log10",math_format(10^.x)))
dev.off()
ggsave("Scatter-xmolZ.pdf", height=8,width=9)















