library(ggplot2)
library(tidyr)
library(devtools)
library(dplyr)
library(party)
library(survival)

# target genes selected
genes=c('ENO2',
        'NCAM1',
        'SYP',
        'CHGA',
        'FOLH1',
        'KLK3 ', "SLC2A1","SLC2A2"
        ,"SLC2A3",
        "SLC2A4","SLC2A5","SLC2A6","SLC2A7","SLC2A8","SLC2A9","SLC2A10","SLC2A11",
        "SLC2A12","SLC2A13","SLC2A14","GCK", "OAX1",  "HK1" ,  "HK2",  "HK3")


###Extracting, creating and cleaning three datasets###

## Extracting and cleaning Cambridge dataset to get a cleaned data frame -- "datacambridge"（Some of the code below has been adapted from https://bioinformatics.cruk.cam.ac.uk/apps/camcAPP/）

# Extracting Gene expressions, survival and demographic variables

library(prostateCancerCamcap)
data(camcap,package = 'prostateCancerCamcap')

pd_camcap <- tbl_df(pData(camcap)) 
fd_camcap <- tbl_df(fData(camcap)) 
exp_camcap <- tbl_df(data.frame(ID = as.character(featureNames(camcap)),exprs(camcap)))
head(exp_camcap)
probes <- fd_camcap %>% 
          filter(Symbol %in% genes)%>%
          dplyr::select(ID)%>% 
          unique %>% 
          as.matrix %>%   
          as.character

data<- exp_camcap %>% 
      filter(ID %in% probes) %>%
      gather(geo_accession,Expression,-ID)  #7761 obs
fd <- fd_camcap
pd=pd_camcap


# Extracting only the gene version (probe) with the largest variability (Most variable probe)

require(dplyr) 
summary_stats <- data %>% 
                 group_by(ID) %>% 
                summarise(mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE),iqr=IQR(Expression,na.rm=TRUE)) #39 obs meaning 39 distinct IDs

data <- left_join(data,summary_stats) %>% 
         mutate(Z = (Expression - mean) / sd) #join by ID,obs become 7761(vars in summary_stats duplicate the same value by ID)

mostVarProbes <- left_join(summary_stats,fd) %>% 
                 arrange(Symbol,desc(iqr)) %>% 
                 distinct(Symbol,.keep_all=TRUE) %>% 
                 dplyr::select(ID) %>%  
                 as.matrix %>%  
                 as.character

data <- filter(data, ID %in% mostVarProbes)


# Join the data sets and define time to death and death indicators and drop where time of death indicator is missing

data <- left_join(data, dplyr::select(fd, ID, Symbol))

data <- left_join(data, pd)%>%
        mutate( Time = as.numeric(as.character(FollowUpTime)), Event = ifelse(BCR=='Y',1,0))

data <- data %>% 
        filter(!is.na(Time) & !is.na(Event))


# Drop the derived variables when we were computing the the most variable probes and make the data horizontal.In other words, make rows by patient record and columns as the variables


datacambridge=data%>%dplyr::select(-Z,-mean,-sd,-iqr,-ID)%>%tidyr::spread(Symbol,Expression) # 111 obs, 40 vars


##Extracting and cleaning Taylor dataset to get a cleaned data frame -- "datataylor"(Same procedure)

if(!require(prostateCancerTaylor)) {
  source('http://www.bioconductor.org/biocLite.R') 
  biocLite('prostateCancerTaylor')
}

#Convert into data convenient for dplyr
library(prostateCancerTaylor)
data(taylor,package = 'prostateCancerTaylor')
pd_taylor <- tbl_df(pData(taylor))
fd_taylor <- tbl_df(fData(taylor))
exp_taylor <- tbl_df(data.frame(ID = as.character(featureNames(taylor)),log2(exprs(taylor))))
probes <- fd_taylor %>% filter(Gene %in% genes) %>% dplyr::select(ID) %>% unique %>% as.matrix %>%  as.character
data <- exp_taylor  %>% filter(ID %in% probes) %>% 
  gather(geo_accession,Expression,-ID)
fd <- fd_taylor %>% mutate(Symbol = Gene)
pd <- pd_taylor
summary_stats <- data %>% group_by(ID) %>% 
  summarise(mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE),iqr=IQR(Expression,na.rm=TRUE))
data <- left_join(data,summary_stats) %>% mutate(Z = (Expression - mean) / sd)
mostVarProbes <- left_join(summary_stats,fd) %>% 
  arrange(Symbol,desc(iqr)) %>% 
  distinct(Symbol,.keep_all=TRUE) %>% 
  dplyr::select(ID) %>%  as.matrix %>%  as.character

data <- filter(data, ID %in% mostVarProbes)
data <- left_join(data, dplyr::select(fd, ID, Symbol))
data <- left_join(data, pd)
data <- data %>% filter(!is.na(Time) & !is.na(Event))

datataylor =data%>%dplyr::select(-Z,-mean,-sd,-iqr,-ID)%>%tidyr::spread(Symbol,Expression) #140 obs, 33 vars



## Extracting and cleaning stockholm dataset to get a cleaned data frame -- "datastockholm"(Same procedure)

if(!require(prostateCancerstockholm)) {
  source('http://www.bioconductor.org/biocLite.R') 
  biocLite('prostateCancerStockholm')
}

# Convert into data convenient for dplyr
library(prostateCancerStockholm)
data(stockholm,package = 'prostateCancerStockholm')
pd_stockholm <- tbl_df(pData(stockholm))
fd_stockholm <- tbl_df(fData(stockholm))
exp_stockholm <- tbl_df(data.frame(ID = as.character(featureNames(stockholm)),log2(exprs(stockholm))))

probes <- fd_stockholm %>% 
  filter(Symbol %in% genes) %>% dplyr::select(ID) %>% unique %>% as.matrix %>%  as.character

data <- exp_stockholm  %>% filter(ID %in% probes) %>% 
  gather(geo_accession,Expression,-ID)
fd <- fd_stockholm %>% mutate(Gene=Symbol)


pd <- mutate(pd_stockholm, Time = as.numeric(FollowUpTime), Event = ifelse(BCR=='Y',1,0))


summary_stats <- data %>% group_by(ID) %>% 
  summarise(mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE),iqr=IQR(Expression,na.rm=TRUE))
data <- left_join(data,summary_stats) %>% mutate(Z = (Expression - mean) / sd)
mostVarProbes <- left_join(summary_stats,fd) %>% 
  arrange(Symbol,desc(iqr)) %>% 
  distinct(Symbol,.keep_all=TRUE) %>% 
  dplyr::select(ID) %>%  as.matrix %>%  as.character

data <- filter(data, ID %in% mostVarProbes)
data <- left_join(data, dplyr::select(fd, ID, Symbol))
data <- left_join(data, pd)
data <- data %>% filter(!is.na(Time) & !is.na(Event))

datastockholm =data%>%dplyr::select(-Z,-mean,-sd,-iqr,-ID)%>%tidyr::spread(Symbol,Expression) #92 obs, 36 vars


## joining three datasets 

library(tcltk)
library(tcltk2)
library(forcats)
datacambridge=datacambridge%>%mutate(dataset="Cam") #add a variable that has value--"Cam" for every observation
datataylor=datataylor%>%mutate(dataset="taylor")
datastockholm=datastockholm%>%mutate(dataset="stoch")

# keep the demographic vars and genes we want only
datacambridge1=datacambridge%>%mutate(Stage=ClinicalStage, iCluster=as.factor(iCluster))%>%select(geo_accession, Age, Time, PSA, Event, Stage, iCluster, Gleason, dataset,ECE, names(datacambridge)[18:40] )

datataylor1=datataylor%>%mutate(iCluster=as.factor(Copy.number.Cluster))%>%
  select(geo_accession,  Time,   Event, Stage, iCluster, Gleason, dataset,  names(datataylor)[11:33] )

datastockholm1=datastockholm%>%mutate(Stage=ClinicalStage, iCluster=as.factor(iCluster))%>%select(geo_accession, Time, PSA, Event, Stage, iCluster, Gleason, dataset,ECE, names(datastockholm)[14:36] )


datacamtaylor=full_join(datacambridge1,datataylor1) #join horizontally
datacamtaylorstock=full_join(datacamtaylor,datastockholm1)

datacamtaylorstock=datacamtaylorstock%>%mutate(PSA=as.numeric(PSA))
datacamtaylorstock=datacamtaylorstock%>%mutate(Age=as.integer(Age))

# clean iCluster and Gleason, combining values for the same levels 
datacamtaylorstock=datacamtaylorstock%>%mutate(iCluster=fct_collapse(as.factor(iCluster),"1"="clust1","2"="clust2","3"="clust3","4"="clust4","5"="clust5"), Gleason=fct_collapse(as.factor(Gleason),"5"="5=3+2","6"="6=3+3","6"="3+3","6"="6=2+4","7"="7=3+4","7"="7=4+3","7"="3+4","7"="4+3","8"="8=3+5","8"="4+4","8"="8=4+4","8"="3+5","8"="5+3","9"="9=4+5","9"="4+5","9"="9=5+4","10"="10=5+5"))
datacamtaylorstock$ECE=as.factor(datacamtaylorstock$ECE)
datacamtaylorstock$dataset=as.factor(datacamtaylorstock$dataset)
datacamtaylorstock$geo_accession=as.factor(datacamtaylorstock$geo_accession)

#Clean Stage, combining values for the same levels
datacamtaylorstock=datacamtaylorstock%>%mutate(Stage=fct_collapse((Stage),"1c"="T1c","1c"="T1C","1c"="T1c N0M0","1c"="T1c N0Mx","1c"="T1c NxMx","1c"="T1c NXMX","2"="T2 N0M0","2"="T2 N0Mx","2"="T2 NxMx","2"="T2 NxMX","2"="T2N0M0",    "2"="T2N0Mx","2a"="T2a","2a"="T2A","2a"="T2a N0M0","2a"="T2a N1M0","2a"= "T2a NxMX","2b"="T2b","2b"="T2B","2b"="T2b N0M0","2b"="T2b N0Mx","2b"="T2b N1Mx","2b"="T2b Nxmx","2b"="T2b NxMx","2c"="T2c","2c"="T2C","2c"="T2c NxMx","2c"="T2cN0Mx","2c"="T2N0M0","2c"="T2N0Mx", "3"="T3","3"="T3 N0M0","3"="T3 N0Mx","3"="T3 NxMx","3a"="T3a","3a"="T3A","3a"="T3a N0M0","3a"="T3a N0Mx","3a"="T3a NxMx","3b"="T3b N0M0","3b"="T3b NxMx","1c"="pT1c NxMx","2"="T1/T2"))
datacamtaylorstock$Stage[datacamtaylorstock$Stage=="NA"]=NA
datacamtaylorstock$Stage[datacamtaylorstock$Stage==""]=NA

### section 1: we conduct a univariate analysis based on three datasets###
library(compareGroups)
library(tcltk2)
data1=as.data.frame(datacamtaylorstock)
cGroupsGUI(data1) # output is Table 1


### section 2: KM Curves###
#KM curves among different Gleason Scores(Graph 2)
library(rms)
library(survival)
library(npsurv)
prostate.out1=rms::npsurv(formula=Surv(Time,Event)~Gleason, data=datacamtaylorstock)
survplot(fit=prostate.out1,conf = c("none","bands","bars")[2], xlab = "Time To BCR", ylab = "Probability of Freedom from BCR",xlim=c(0,150),ylim = c(0,1),label.curves = TRUE, levels.only = FALSE, abbrev.label = FALSE, loglog = FALSE, logt = FALSE,dots = TRUE, 
         n.risk = TRUE)

# KM curves among three datasets (Graph 1)
library(rms)
library(survival)
library(npsurv)
prostate.out2=rms::npsurv(formula=Surv(Time,Event)~dataset, data=datacamtaylorstock)
survplot(prostate.out2,conf = c("none","bands","bars")[2], xlab = "Time to BCR", ylab = "Probability of Freedom from BCR",xlim=c(0,150),ylim = c(0,1),label.curves = TRUE, levels.only = FALSE, abbrev.label = FALSE, loglog = FALSE, logt = FALSE,dots = TRUE, 
         n.risk = TRUE, main="Graph 1  KM curves among three datasets")

###Section 3:Stratified Logrank Test among three datasets ###

#Logrank Test among three datasets 
library(survival)
prostate.out3 =survdiff( Surv(Time,Event) ~ dataset,  rho=0, data=datacamtaylorstock) 
prostate.out3

prostate.out4 =survdiff( Surv(Time,Event) ~ dataset,  rho=1, data=datacamtaylorstock) 
prostate.out4


#found "iCluster" is a cofounder. so use stratified logrank test.

prostate.out5 =survdiff( Surv(Time,Event) ~ iCluster,  rho=0, data=datacamtaylorstock) 
prostate.out5
chisq.test( table( datacamtaylorstock$dataset, datacamtaylorstock$iCluster ) )

prostate.out6 =survdiff( Surv(Time,Event) ~ dataset+strata(iCluster) ,  rho=0, data=datacamtaylorstock) 
prostate.out6

###Section 4:Cox’s proportional hazards regression model###

#Delete Age, PSA and ECE since not all the dataset have these three variables.

library(survival)
library(MASS)
datacamtaylorstock1=datacamtaylorstock%>%select(-Age, -PSA,-ECE)
datacamtaylorstock2=na.omit(datacamtaylorstock1)


#imputation for missing values. 20 NAs. impute 18 in "iCluster" and "Gleason", and delete 2 in"Stage". datacamtaylorstock6 has 341 obs and 30 variables, which is the final dataset for cox model

library(mice)
library(randomForest)
md.pattern(datacamtaylorstock1)
names(datacamtaylorstock1)
datacamtaylorstock3=datacamtaylorstock1%>%select(5:20)
imp <- mice(datacamtaylorstock3,maxit = 5,m=5,method = "rf",seed=500)
datacamtaylorstock4<- complete(imp,1) 
sum(is.na(datacamtaylorstock4))
datacamtaylorstock5=cbind(datacamtaylorstock1[1:4],datacamtaylorstock4,datacamtaylorstock1[21:30])
sum(is.na(datacamtaylorstock5))
datacamtaylorstock6=na.omit(datacamtaylorstock5)


# cox model using full imputed dataset without "Age, PSA and ECE"

fitall=coxph(formula = Surv(Time,Event,type = "right") ~ Stage+iCluster+Gleason+dataset+ENO2+NCAM1+SYP+CHGA+FOLH1+SLC2A1+SLC2A2+SLC2A3+SLC2A4+SLC2A5+SLC2A6+SLC2A7+SLC2A8+SLC2A9+SLC2A10+SLC2A11+SLC2A12+SLC2A13+SLC2A14+GCK+HK1+HK2+HK3,data = datacamtaylorstock6)
aic1=stepAIC(fitall,scope =list(upper = ~Stage+iCluster+Gleason+dataset+ENO2+NCAM1+SYP+CHGA+FOLH1+SLC2A1+SLC2A2+SLC2A3+SLC2A4+SLC2A5+SLC2A6+SLC2A7+SLC2A8+SLC2A9+SLC2A10+SLC2A11+SLC2A12+SLC2A13+SLC2A14+GCK+HK1+HK2+HK3 , lower = ~1,
                                direction="both"),steps=100000,k=2)
aic1$anova
aic1


#cox model with only cambridge dataset plus deduct NAs

library(survival)
library(MASS)
datacamtaylorstock22=na.omit(datacamtaylorstock)
fitCam=coxph(formula = Surv(Time,Event,type = "right") ~ 
               Age+PSA+Stage+iCluster+Gleason+dataset+ECE+ENO2+NCAM1+SYP+CHGA+FOLH1+SLC2A1+SLC2A2+SLC2A3+SLC2A4+SLC2A5+SLC2A6+SLC2A7+SLC2A8+SLC2A9+SLC2A10+SLC2A11+SLC2A12+SLC2A13+SLC2A14+GCK+HK1+HK2+HK3,data = datacamtaylorstock22)
aicCam=stepAIC(fitCam,scope =list(upper = ~Age+PSA+Stage+iCluster+Gleason+dataset+ECE+ENO2+NCAM1+SYP+CHGA+FOLH1+SLC2A1+SLC2A2+SLC2A3+SLC2A4+SLC2A5+SLC2A6+SLC2A7+SLC2A8+SLC2A9+SLC2A10+SLC2A11+SLC2A12+SLC2A13+SLC2A14+GCK+HK1+HK2+HK3 , lower = ~1,
                                  direction="both"),steps=100000,k=2)
aicCam$anova
aicCam




#cox model using datacamtaylorstock6 with interactions 

library(survival)
library(MASS)
fitall_intersection=coxph(formula = Surv(Time,Event,type = "right") ~ 
                            Stage+iCluster+Gleason+dataset+ENO2+NCAM1+SYP+CHGA+FOLH1+SLC2A1+SLC2A2+SLC2A3+SLC2A4+SLC2A5+SLC2A6+SLC2A7+SLC2A8+SLC2A9+SLC2A10+SLC2A11+SLC2A12+SLC2A13+SLC2A14+GCK+HK1+HK2+HK3+ENO2*dataset+NCAM1*dataset+SYP*dataset+CHGA*dataset+FOLH1*dataset+SLC2A1*dataset+SLC2A2*dataset+SLC2A3*dataset+SLC2A4*dataset+SLC2A5*dataset+SLC2A6*dataset+SLC2A7*dataset+SLC2A8*dataset+SLC2A9*dataset+SLC2A10*dataset+SLC2A11*dataset+SLC2A12*dataset+SLC2A13*dataset+SLC2A14*dataset+GCK*dataset+HK1*dataset+HK2*dataset+HK3*dataset,data = datacamtaylorstock6)

aic1_interaction=stepAIC(fitall_intersection,scope =list(upper = ~Stage+iCluster+Gleason+dataset+ENO2+NCAM1+SYP+CHGA+FOLH1+SLC2A1+SLC2A2+SLC2A3+SLC2A4+SLC2A5+SLC2A6+SLC2A7+SLC2A8+SLC2A9+SLC2A10+SLC2A11+SLC2A12+SLC2A13+SLC2A14+GCK+HK1+HK2+HK3+ENO2*dataset+NCAM1*dataset+SYP*dataset+CHGA*dataset+FOLH1*dataset+SLC2A1*dataset+SLC2A2*dataset+SLC2A3*dataset+SLC2A4*dataset+SLC2A5*dataset+SLC2A6*dataset+SLC2A7*dataset+SLC2A8*dataset+SLC2A9*dataset+SLC2A10*dataset+SLC2A11*dataset+SLC2A12*dataset+SLC2A13*dataset+SLC2A14*dataset+GCK*dataset+HK1*dataset+HK2*dataset+HK3*dataset , lower = ~1,
                                                         direction="both"),steps=100000,k=2)

aic1_interaction$anova
aic1_interaction



#cox model with full data no imputation (simply delete NAs)!

library(survival)
library(MASS)
fitall_intersection2=coxph(formula = Surv(Time,Event,type = "right") ~ 
                             Stage+iCluster+Gleason+dataset+ENO2+NCAM1+SYP+CHGA+FOLH1+SLC2A1+SLC2A2+SLC2A3+SLC2A4+SLC2A5+SLC2A6+SLC2A7+SLC2A8+SLC2A9+SLC2A10+SLC2A11+SLC2A12+SLC2A13+SLC2A14+GCK+HK1+HK2+HK3+ENO2*dataset+NCAM1*dataset+SYP*dataset+CHGA*dataset+FOLH1*dataset+SLC2A1*dataset+SLC2A2*dataset+SLC2A3*dataset+SLC2A4*dataset+SLC2A5*dataset+SLC2A6*dataset+SLC2A7*dataset+SLC2A8*dataset+SLC2A9*dataset+SLC2A10*dataset+SLC2A11*dataset+SLC2A12*dataset+SLC2A13*dataset+SLC2A14*dataset+GCK*dataset+HK1*dataset+HK2*dataset+HK3*dataset,data = datacamtaylorstock2)

aic1_interaction2=stepAIC(fitall_intersection2,scope =list(upper = ~Stage+iCluster+Gleason+dataset+ENO2+NCAM1+SYP+CHGA+FOLH1+SLC2A1+SLC2A2+SLC2A3+SLC2A4+SLC2A5+SLC2A6+SLC2A7+SLC2A8+SLC2A9+SLC2A10+SLC2A11+SLC2A12+SLC2A13+SLC2A14+GCK+HK1+HK2+HK3+ENO2*dataset+NCAM1*dataset+SYP*dataset+CHGA*dataset+FOLH1*dataset+SLC2A1*dataset+SLC2A2*dataset+SLC2A3*dataset+SLC2A4*dataset+SLC2A5*dataset+SLC2A6*dataset+SLC2A7*dataset+SLC2A8*dataset+SLC2A9*dataset+SLC2A10*dataset+SLC2A11*dataset+SLC2A12*dataset+SLC2A13*dataset+SLC2A14*dataset+GCK*dataset+HK1*dataset+HK2*dataset+HK3*dataset , lower = ~1,
                                                           direction="both"),k=2)

aic1_interaction2$anova
aic1_interaction2



###Section 5  Check proportionality assumption ###

library(survival)
library(survminer)
#check PH assumption of aic1_interaction
check1=cox.zph(aic1_interaction,"log")	
check1
# par(mfrow=c(4,4))
# plot(check1)
#ggcoxzph(check1)




library(randomForestSRC)

fit.rf <- rfsrc(Surv(Time, Event) ~Stage+iCluster+Gleason+dataset+ENO2+NCAM1+SYP+CHGA+FOLH1+SLC2A1+SLC2A2+SLC2A3+SLC2A4+SLC2A5+SLC2A6+SLC2A7+SLC2A8+SLC2A9+SLC2A10+SLC2A11+SLC2A12+SLC2A13+SLC2A14+GCK+HK1+HK2+HK3+ENO2*dataset+NCAM1*dataset+SYP*dataset+CHGA*dataset+FOLH1*dataset+SLC2A1*dataset+SLC2A2*dataset+SLC2A3*dataset+SLC2A4*dataset+SLC2A5*dataset+SLC2A6*dataset+SLC2A7*dataset+SLC2A8*dataset+SLC2A9*dataset+SLC2A10*dataset+SLC2A11*dataset+SLC2A12*dataset+SLC2A13*dataset+SLC2A14*dataset+GCK*dataset+HK1*dataset+HK2*dataset+HK3*dataset,data = datacamtaylorstock6, ntree = 400, block.size = 1,importance =TRUE, Splitrules="logrank")
print(fit.rf)
plot(fit.rf)


#just for cam rf

fit.rf.Cam <- rfsrc(Surv(Time, Event) ~ENO2+NCAM1+SYP+CHGA+FOLH1+SLC2A1+SLC2A2+SLC2A3+SLC2A4+SLC2A5+SLC2A6+SLC2A7+SLC2A8+SLC2A9+SLC2A10+SLC2A11+SLC2A12+SLC2A13+SLC2A14+GCK+HK1+HK2+HK3,data =datacamtaylorstock22, ntree = 400, block.size = 1,importance =TRUE, Splitrules="logrank")
print(fit.rf.Cam)
plot(fit.rf.Cam)
