
library(riskRegression, verbose = FALSE, quietly = TRUE)
library(survival)
library(Hmisc)
library(rms)
library(splines)
library(ggplot2)
#install.packages("splitstackshape")
library(splitstackshape) ##expandrows function
#install.packages("gee")
library(gee) # version 4.13-25


############### Step 1. define the intervention population
###Extract baseline data, including baseline disease states and covariates

ARIC_base=select(mesa_all, ID_C, age_base, S0, S1, S2, S3, edu, fam_chd, obesity_v1, gender, race)

summary(ARIC_base)


##restrict to those with S1=1, starting from at-risk state

ARIC_base=ARIC_base[which(ARIC_base$S1==1),]


##exclude missing values

ARIC_base=ARIC_base[which(is.na(ARIC_base$edu)==F & 
is.na(ARIC_base$fam_chd)==F & is.na(ARIC_base$obesity_v1)==F
& is.na(ARIC_base$gender)==F & is.na(ARIC_base$race)==F),]

nrow(ARIC_base)


#check, no participants start with S2 or s3

table(ARIC_base$S1, ARIC_base$S2, useNA = "ifany")

table(ARIC_base$S1, ARIC_base$S3, useNA = "ifany")


###change the projection population to long format

# summary(ARIC_base$age_base)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  44.00   49.00   54.00   54.23   59.00   66.00 


####Specify parameters

var_loop=1000   ##bootstrap loops 100

#tm=110-44  ###expand each id to tm observations, we allow the youngest to live to 110 years

tm=30

max_i=nrow(ARIC_base)

#max_i=10

num_rg=2 #num_rg is the number of exposures

define_p_death=1  ##no need to have restrictions, as this is not estimating life expectancy


###create long-format data and time-varying age

ARIC_base=ARIC_base

ARIC_long <- expandRows(ARIC_base, count=tm, count.is.col=F)
nrow(ARIC_long)

ARIC_long$agecon <- as.numeric(ARIC_long$age_base) + rep(0:(tm - 1), length.out = tm)

summary(ARIC_long$agecon)


#############Step 2. Define exposure regimes


ARIC_S1_rg1 = ARIC_long
ARIC_S2_rg1 = ARIC_long
ARIC_S3_rg1 = ARIC_long

ARIC_S1_rg2 = ARIC_long
ARIC_S2_rg2 = ARIC_long
ARIC_S3_rg2 = ARIC_long


ARIC_S1_rg1$smkcon=0
ARIC_S2_rg1$smkcon=0
ARIC_S3_rg1$smkcon=0

ARIC_S1_rg2$smkcon=1
ARIC_S2_rg2$smkcon=1
ARIC_S3_rg2$smkcon=1



##############Step 3. Construct age-specific TR and estimate MALY and disease path for this population


####################################################
##########for each bootstrap population#############
####################################################


#####################################
#######mediation analysis############
#####################################


##survival risks

NDE_marginal_rg1 <- matrix(0, nrow=var_loop,ncol=1) 
NIE_m1_marginal_rg1 <- matrix(0, nrow=var_loop,ncol=1) 
NIE_m2_marginal_rg1 <- matrix(0, nrow=var_loop,ncol=1) 
NIE_m1_m2_marginal_rg1 <- matrix(0, nrow=var_loop,ncol=1) 
TE_marginal_rg1 <- matrix(0, nrow=var_loop,ncol=1) 


NDE_marginal_rg2 <- matrix(0, nrow=var_loop,ncol=1) 
NIE_m1_marginal_rg2 <- matrix(0, nrow=var_loop,ncol=1) 
NIE_m2_marginal_rg2 <- matrix(0, nrow=var_loop,ncol=1) 
NIE_m1_m2_marginal_rg2 <- matrix(0, nrow=var_loop,ncol=1) 
TE_marginal_rg2 <- matrix(0, nrow=var_loop,ncol=1) 


###

PIE_m1_marginal_rg2 <- matrix(0, nrow=var_loop,ncol=1) 
PIE_m2_marginal_rg2 <- matrix(0, nrow=var_loop,ncol=1) 
PIE_m1_m2_marginal_rg2 <- matrix(0, nrow=var_loop,ncol=1) 


CIE_m1_marginal_rg2 <- matrix(0, nrow=var_loop,ncol=1) 
CIE_m2_marginal_rg2 <- matrix(0, nrow=var_loop,ncol=1) 
CIE_m1_m2_marginal_rg2 <- matrix(0, nrow=var_loop,ncol=1) 


##mediation effect

TE_marginal_rg21 <- matrix(0, nrow=var_loop,ncol=1) 
NDE_marginal_rg21 <- matrix(0, nrow=var_loop,ncol=1) 

NIE_marginal_rg21 <- matrix(0, nrow=var_loop,ncol=1) 
NIE_m1_marginal_rg21 <- matrix(0, nrow=var_loop,ncol=1) 
NIE_m2_marginal_rg21 <- matrix(0, nrow=var_loop,ncol=1) 
NIE_m1_m2_marginal_rg21 <- matrix(0, nrow=var_loop,ncol=1) 

PIE_marginal_rg21 <- matrix(0, nrow=var_loop,ncol=1) 
PIE_m1_marginal_rg21 <- matrix(0, nrow=var_loop,ncol=1) 
PIE_m2_marginal_rg21 <- matrix(0, nrow=var_loop,ncol=1) 
PIE_m1_m2_marginal_rg21 <- matrix(0, nrow=var_loop,ncol=1) 

CIE_marginal_rg21 <- matrix(0, nrow=var_loop,ncol=1) 
CIE_m1_marginal_rg21 <- matrix(0, nrow=var_loop,ncol=1) 
CIE_m2_marginal_rg21 <- matrix(0, nrow=var_loop,ncol=1) 
CIE_m1_m2_marginal_rg21 <- matrix(0, nrow=var_loop,ncol=1) 


##percentage

NDE_100_marginal_rg21 <- matrix(0, nrow=var_loop,ncol=1) 

NIE_100_marginal_rg21 <- matrix(0, nrow=var_loop,ncol=1) 
NIE_100_m1_marginal_rg21<- matrix(0, nrow=var_loop,ncol=1) 
NIE_100_m2_marginal_rg21<- matrix(0, nrow=var_loop,ncol=1) 
NIE_100_m1_m2_marginal_rg21<- matrix(0, nrow=var_loop,ncol=1) 

PIE_100_marginal_rg21 <- matrix(0, nrow=var_loop,ncol=1) 
PIE_100_m1_marginal_rg21<- matrix(0, nrow=var_loop,ncol=1) 
PIE_100_m2_marginal_rg21<- matrix(0, nrow=var_loop,ncol=1) 
PIE_100_m1_m2_marginal_rg21<- matrix(0, nrow=var_loop,ncol=1) 

CIE_100_marginal_rg21 <- matrix(0, nrow=var_loop,ncol=1) 
CIE_100_m1_marginal_rg21<- matrix(0, nrow=var_loop,ncol=1) 
CIE_100_m2_marginal_rg21<- matrix(0, nrow=var_loop,ncol=1) 
CIE_100_m1_m2_marginal_rg21<- matrix(0, nrow=var_loop,ncol=1) 


########start bootstrap


for (simu in 1:var_loop){


set.seed(simu)  


mesa_all_boot<-mesa_all[sample(1:nrow(mesa_all), size=nrow(mesa_all), replace=TRUE),]


data1=mesa_all_boot[which(mesa_all_boot$S1==1),]   ##start from at risk
data2=mesa_all_boot[which(mesa_all_boot$S2==1),]   ##start from CHD
data3=mesa_all_boot[which(mesa_all_boot$S3==1),]   ##start from HF


data1$tte=data1$S1_end-data1$S1_start
data2$tte=data2$S2_end-data2$S2_start
data3$tte=data3$S3_end-data3$S3_start


################################
##change data from wide to long#
################################


##data1 long format

data1<-data1[which(data1$tte>0),]    ###exclude missing tte or tte=0

summary(data1$tte)

data1$S1_ylong<-0
data1[which(data1$S1_y>0),]$S1_ylong<-1

addmargins(table(data1$S1_y))
table(data1$S1_ylong)

library(survival)

data1_long <- survSplit(Surv(tte, S1_ylong)~.,
                        data = data1,
                        cut = c(0:(trunc(max(data1$tte))+1)), ##cut follow up year into intervals
                        episode="interval",
                        start="tstart",
                        end = "tstop")

##check if the number of events matches

table(data1_long$S1_ylong)

data1_long$S1_mylong<-0
data1_long[which(data1_long$S1_y==2 & data1_long$S1_ylong==1),]$S1_mylong<-2
data1_long[which(data1_long$S1_y==3 & data1_long$S1_ylong==1),]$S1_mylong<-3
data1_long[which(data1_long$S1_y==4 & data1_long$S1_ylong==1),]$S1_mylong<-4


table(data1_long$S1_mylong)
table(data1$S1_y)


##define age 

data1_long$agecon=data1_long$age_base+data1_long$interval-1


##define time-varying variables


data1_long <- data1_long %>% 
  
  mutate (
    obesitycon = case_when (
   tstart<V2DAYS_assign ~ obesity_v1,
   V2DAYS_assign<=tstart & tstart<V3DAYS_assign ~ obesity_v2,
   V3DAYS_assign<=tstart & tstart<V4DAYS_assign ~ obesity_v3,
   V4DAYS_assign<=tstart & tstart<V5DATE51_DAYS_assign~ obesity_v4,
   V5DATE51_DAYS_assign<=tstart & tstart<V6DATE61_DAYS_assign~ obesity_v5,
   V6DATE61_DAYS_assign<=tstart ~ obesity_v6),

    smkcon = case_when (
   tstart<V2DAYS_assign ~ smk1,
   V2DAYS_assign<=tstart & tstart<V3DAYS_assign ~ smk2,
   V3DAYS_assign<=tstart & tstart<V4DAYS_assign ~ smk3,
   V4DAYS_assign<=tstart & tstart<V5DATE51_DAYS_assign~ smk4,
   V5DATE51_DAYS_assign<=tstart & tstart<V6DATE61_DAYS_assign~ smk5,
   V6DATE61_DAYS_assign<=tstart ~ smk6)
       )


##data2 long format

data2<-data2[which(data2$tte>0),]    ###exclude missing tte or tte=0

summary(data2$tte)

data2$S2_ylong<-0
data2[which(data2$S2_y>0),]$S2_ylong<-1

addmargins(table(data2$S2_y))
table(data2$S2_ylong)

library(survival)

data2_long <- survSplit(Surv(tte, S2_ylong)~.,
                        data = data2,
                        cut = c(0:(trunc(max(data2$tte))+1)), ##cut follow up year into intervals
                        episode="interval",
                        start="tstart",
                        end = "tstop")

##check if the number of events matches

table(data2_long$S2_ylong)

data2_long$S2_mylong<-0

data2_long[which(data2_long$S2_y==3 & data2_long$S2_ylong==1),]$S2_mylong<-3
data2_long[which(data2_long$S2_y==4 & data2_long$S2_ylong==1),]$S2_mylong<-4

table(data2_long$S2_mylong)
table(data2$S2_y)

prop.table(table(data2_long$S2_mylong,exclude=NULL))  ##Incident  Source: MI or FATCHD

##define age 

data2_long$agecon=data2_long$age_base+data2_long$interval-1

##define time-varying variables


data2_long <- data2_long %>% 
  
  mutate (
   obesitycon = case_when (
   tstart<V2DAYS_assign ~ obesity_v1,
   V2DAYS_assign<=tstart & tstart<V3DAYS_assign ~ obesity_v2,
   V3DAYS_assign<=tstart & tstart<V4DAYS_assign ~ obesity_v3,
   V4DAYS_assign<=tstart & tstart<V5DATE51_DAYS_assign ~ obesity_v4,
   V5DATE51_DAYS_assign<=tstart & tstart<V6DATE61_DAYS_assign ~ obesity_v5,
   V6DATE61_DAYS_assign<=tstart ~ obesity_v6),

    smkcon = case_when (
   tstart<V2DAYS_assign ~ smk1,
   V2DAYS_assign<=tstart & tstart<V3DAYS_assign ~ smk2,
   V3DAYS_assign<=tstart & tstart<V4DAYS_assign ~ smk3,
   V4DAYS_assign<=tstart & tstart<V5DATE51_DAYS_assign~ smk4,
   V5DATE51_DAYS_assign<=tstart & tstart<V6DATE61_DAYS_assign~ smk5,
   V6DATE61_DAYS_assign<=tstart ~ smk6)

       )


####data3 long format

data3<-data3[which(data3$tte>0),]    ###exclude missing tte or tte=0

summary(data3$tte)

data3$S3_ylong<-0
data3[which(data3$S3_y>0),]$S3_ylong<-1

addmargins(table(data3$S3_y))
table(data3$S3_ylong)

library(survival)

data3_long <- survSplit(Surv(tte, S3_ylong)~.,
                        data = data3,
                        cut = c(0:(trunc(max(data3$tte))+1)),
                        episode="interval",
                        start="tstart",
                        end = "tstop")

##check if the number of events matches

table(data3_long$S3_ylong)

data3_long$S3_mylong<-0

data3_long[which(data3_long$S3_y==4 & data3_long$S3_ylong==1),]$S3_mylong<-4

table(data3_long$S3_mylong)
table(data3$S3_y)

prop.table(table(data3_long$S3_mylong,exclude=NULL))  ##Incident  Source: MI or FATCHD

##define age 

data3_long$agecon=data3_long$age_base+data3_long$interval-1

##define time-varying variables


data3_long <- data3_long %>% 
  
  mutate (
    obesitycon = case_when (
   tstart<V2DAYS_assign ~ obesity_v1,
   V2DAYS_assign<=tstart & tstart<V3DAYS_assign ~ obesity_v2,
   V3DAYS_assign<=tstart & tstart<V4DAYS_assign ~ obesity_v3,
   V4DAYS_assign<=tstart & tstart<V5DATE51_DAYS_assign~ obesity_v4,
   V5DATE51_DAYS_assign<=tstart & tstart<V6DATE61_DAYS_assign~ obesity_v5,
   V6DATE61_DAYS_assign<=tstart ~ obesity_v6),

    smkcon = case_when (
   tstart<V2DAYS_assign ~ smk1,
   V2DAYS_assign<=tstart & tstart<V3DAYS_assign ~ smk2,
   V3DAYS_assign<=tstart & tstart<V4DAYS_assign ~ smk3,
   V4DAYS_assign<=tstart & tstart<V5DATE51_DAYS_assign~ smk4,
   V5DATE51_DAYS_assign<=tstart & tstart<V6DATE61_DAYS_assign~ smk5,
   V6DATE61_DAYS_assign<=tstart ~ smk6)

       )


##########Fit CSC models

##Dataset1: Estimate transition rates from state 1 to states 2-4 

data1_long$tte_cmr=data1_long$tstop-data1_long$tstart
summary(data1_long$tte_cmr)

##before running bootstrap, test the markov assumption in the original population

cfit1 <- CSC(formula = list(
Hist(tte_cmr,S1_mylong) ~ agecon+as.factor(smkcon)+as.factor(edu)+as.factor(fam_chd)+as.factor(obesity_v1)+as.factor(gender) +as.factor(race),
Hist(tte_cmr,S1_mylong) ~ agecon+as.factor(smkcon)+as.factor(edu)+as.factor(fam_chd)+as.factor(obesity_v1)+as.factor(gender) +as.factor(race),
Hist(tte_cmr,S1_mylong) ~ agecon+as.factor(smkcon)+as.factor(edu)+as.factor(fam_chd)+as.factor(obesity_v1)+as.factor(gender) +as.factor(race)), data = data1_long)

summary(data1_long$tte_cmr)


##Dataset2: Estimate transition rates from state 2 to states 3-4 

data2_long$tte_cmr=data2_long$tstop-data2_long$tstart
summary(data2_long$tte_cmr)


cfit2 <- CSC(formula = list(
Hist(tte_cmr, S2_mylong) ~ agecon+as.factor(smkcon)+as.factor(edu)+as.factor(fam_chd)+as.factor(obesity_v1)+as.factor(gender) +as.factor(race),
Hist(tte_cmr, S2_mylong) ~ agecon+as.factor(smkcon)+as.factor(edu)+as.factor(fam_chd)+as.factor(obesity_v1)+as.factor(gender) +as.factor(race)),
             data = data2_long)


##Dataset3: Estimate transition rates from state 3 to state 4 

data3_long<-data3_long %>% 
  mutate (
    
    sub3=case_when (
      S3_history=="bhtn-chd-chf"   ~ "htn-chd-chf",
      S3_history=="bhtn-chf"  ~ "htn-chf"),
    
    S3_mylong_new=case_when(
      S3_mylong==4 ~ 1,
      TRUE~0)
  )

table(data3_long$S3_mylong_new,data3_long$S3_mylong)
table(data3_long$S3_mylong)

data3_long$tte_cmr=data3_long$tstop-data3_long$tstart
summary(data3_long$tte_cmr)

cfit3 <- coxph(
Surv(tte_cmr,S3_mylong_new) ~as.factor(sub3)+agecon+as.factor(smkcon)+as.factor(edu)+as.factor(fam_chd)+as.factor(obesity_v1)+as.factor(gender) +as.factor(race),
    data = data3_long, x = TRUE)


############Step 3b. Regime 1----Estimate transition rates (TR) 

##prediction in dataset1

pfit12_rg1 <- predict(cfit1, newdata = ARIC_S1_rg1, cause = 2, times =0.5, se = F, keep.newdata = FALSE)
pfit13_rg1 <- predict(cfit1, newdata = ARIC_S1_rg1, cause = 3, times =0.5, se = F, keep.newdata = FALSE)
pfit14_rg1 <- predict(cfit1, newdata = ARIC_S1_rg1, cause = 4, times =0.5, se = F, keep.newdata = FALSE)


predict1_rg1<-data.frame(
                pfit12_rg1$absRisk,  
                pfit13_rg1$absRisk,  
                pfit14_rg1$absRisk)

colnames(predict1_rg1)=c("r12","r13","r14")


##predict in dataset 2


pfit23_rg1 <- predict(cfit2, newdata = ARIC_S2_rg1, cause = 3, times =0.5, se = F, keep.newdata = FALSE)
pfit24_rg1 <- predict(cfit2, newdata = ARIC_S2_rg1, cause = 4, times =0.5, se = F, keep.newdata = FALSE)


predict2_rg1<-data.frame(pfit23_rg1$absRisk,pfit24_rg1$absRisk)
colnames(predict2_rg1)=c("r23","r24")


##predict in dataset 3


new_data3_3_rg1<-ARIC_S3_rg1
new_data3_3_rg1$sub3="htn-chd-chf"

new_data3_4_rg1<-ARIC_S3_rg1
new_data3_4_rg1$sub3="htn-chf"


predict3_012_rg1 <- data.frame(1-predictCox(cfit3, newdata= new_data3_3_rg1, times=0.5, se = F)$survival)
colnames(predict3_012_rg1)=c("r34")

predict3_01_rg1 <- data.frame(1-predictCox(cfit3, newdata= new_data3_4_rg1, times=0.5, se = F)$survival)
colnames(predict3_01_rg1)=c("r34")


############Step 3b. Regime 2----Estimate transition rates (TR) 


##prediction in dataset1

pfit12_rg2 <- predict(cfit1, newdata = ARIC_S1_rg2, cause = 2, times =0.5, se = F, keep.newdata = FALSE)
pfit13_rg2 <- predict(cfit1, newdata = ARIC_S1_rg2, cause = 3, times =0.5, se = F, keep.newdata = FALSE)
pfit14_rg2 <- predict(cfit1, newdata = ARIC_S1_rg2, cause = 4, times =0.5, se = F, keep.newdata = FALSE)


predict1_rg2<-data.frame(
                pfit12_rg2$absRisk,  
                pfit13_rg2$absRisk,  
                pfit14_rg2$absRisk)

colnames(predict1_rg2)=c("r12","r13","r14")


##predict in dataset 2


pfit23_rg2 <- predict(cfit2, newdata = ARIC_S2_rg2, cause = 3, times =0.5, se = F, keep.newdata = FALSE)
pfit24_rg2 <- predict(cfit2, newdata = ARIC_S2_rg2, cause = 4, times =0.5, se = F, keep.newdata = FALSE)


predict2_rg2<-data.frame(pfit23_rg2$absRisk,pfit24_rg2$absRisk)
colnames(predict2_rg2)=c("r23","r24")


##predict in dataset 3


new_data3_3_rg2<-ARIC_S3_rg2
new_data3_3_rg2$sub3="htn-chd-chf"

new_data3_4_rg2<-ARIC_S3_rg2
new_data3_4_rg2$sub3="htn-chf"


predict3_012_rg2 <- data.frame(1-predictCox(cfit3, newdata= new_data3_3_rg2, times=0.5, se = F)$survival)
colnames(predict3_012_rg2)=c("r34")

predict3_01_rg2 <- data.frame(1-predictCox(cfit3, newdata= new_data3_4_rg2, times=0.5, se = F)$survival)
colnames(predict3_01_rg2)=c("r34")


######################################
##########for each person#############
######################################


SSLY_person_rg1<- matrix(0, nrow=nrow(ARIC_base),ncol=8) 
SOP_person_rg1<- matrix(0, nrow=nrow(ARIC_base),ncol=8) 

SSLY_person_rg2<- matrix(0, nrow=nrow(ARIC_base),ncol=8) 
SOP_person_rg2<- matrix(0, nrow=nrow(ARIC_base),ncol=8) 
 

CIE_SSLY_person_rg2_M1_M2<- matrix(0, nrow=nrow(ARIC_base),ncol=8) 
CIE_SOP_person_rg2_M1_M2<- matrix(0, nrow=nrow(ARIC_base),ncol=8) 
 
CIE_SSLY_person_rg2_M1<- matrix(0, nrow=nrow(ARIC_base),ncol=8) 
CIE_SOP_person_rg2_M1<- matrix(0, nrow=nrow(ARIC_base),ncol=8) 
 
CIE_SSLY_person_rg2_M2<- matrix(0, nrow=nrow(ARIC_base),ncol=8) 
CIE_SOP_person_rg2_M2<- matrix(0, nrow=nrow(ARIC_base),ncol=8) 
 

PIE_SSLY_person_rg2_M1_M2<- matrix(0, nrow=nrow(ARIC_base),ncol=8) 
PIE_SOP_person_rg2_M1_M2<- matrix(0, nrow=nrow(ARIC_base),ncol=8) 
 
PIE_SSLY_person_rg2_M1<- matrix(0, nrow=nrow(ARIC_base),ncol=8) 
PIE_SOP_person_rg2_M1<- matrix(0, nrow=nrow(ARIC_base),ncol=8) 
 
PIE_SSLY_person_rg2_M2<- matrix(0, nrow=nrow(ARIC_base),ncol=8) 
PIE_SOP_person_rg2_M2<- matrix(0, nrow=nrow(ARIC_base),ncol=8) 
 


########Step 3b2. Regime 1-3. simulate TR, TP, SOP, MALY, disease path


for (i in 1:max_i) {   ##each individual in the data


################################################
##########each person across time###############
################################################

TR_var <- array(0, dim = c(8, 8, tm, num_rg))  

# Vectorized tracking index for the current person 'i' across all time points

person_time_indices <- (1:tm) + (i - 1) * tm


for (rg in 1: num_rg) {
  
  p1 <- get(paste0("predict1_rg", rg))
  p2 <- get(paste0("predict2_rg", rg))
  p3_01 <- get(paste0("predict3_01_rg", rg))
  p3_012 <- get(paste0("predict3_012_rg", rg))
  
  TR_var[1, 2, , rg] <- p1$r12[person_time_indices]
  TR_var[1, 3, , rg] <- p1$r13[person_time_indices]
  TR_var[1, 5, , rg] <- p1$r14[person_time_indices]
  
  TR_var[2, 4, , rg] <- p2$r23[person_time_indices]
  TR_var[2, 6, , rg] <- p2$r24[person_time_indices]
  
  TR_var[3, 7, , rg] <- p3_01$r34[person_time_indices]
  
  TR_var[4, 8, , rg] <- p3_012$r34[person_time_indices]
    
  for (t in 1:tm) {
    for (s_tr in 1:7) {
      TR_var[s_tr, s_tr, t, rg] <- -sum(TR_var[s_tr, (s_tr + 1):8, t, rg])
    }
  }

}



##############################
###########CIE################
##############################


CIE_TR_var_M1_M2 <- array(0, dim = c(8, 8, tm, 1))  

CIE_TR_var_M1 <- array(0, dim = c(8, 8, tm, 1))  

CIE_TR_var_M2 <- array(0, dim = c(8, 8, tm, 1))  


##rg2 vs. rg1, through M1 and M2

##choose to use E=(0,0,1,1)

  CIE_TR_var_M1_M2[1, 2, ,1] <- predict1_rg1$r12[person_time_indices] 
  CIE_TR_var_M1_M2[1, 3, ,1] <- predict1_rg1$r13[person_time_indices] 
  CIE_TR_var_M1_M2[1, 5, ,1] <- predict1_rg1$r14[person_time_indices] 
  
  CIE_TR_var_M1_M2[2, 4, ,1] <- predict2_rg1$r23[person_time_indices] 
  CIE_TR_var_M1_M2[2, 6, ,1] <- predict2_rg1$r24[person_time_indices] 
  
  CIE_TR_var_M1_M2[3, 7, ,1] <- predict3_01_rg2$r34[person_time_indices]   
  CIE_TR_var_M1_M2[4, 8, ,1] <- predict3_012_rg2$r34[person_time_indices]  
   
 
  for (t in 1:tm) {
    for (s_tr in 1:7) {
      CIE_TR_var_M1_M2[s_tr, s_tr, t,1] <- -sum(CIE_TR_var_M1_M2[s_tr, (s_tr + 1):8, t,1])
    }
  }




##rg2 vs. rg1, through M1

##choose to use E=(0,1,1,1)

  CIE_TR_var_M1[1, 2, ,1] <- predict1_rg1$r12[person_time_indices] 
  CIE_TR_var_M1[1, 3, ,1] <- predict1_rg1$r13[person_time_indices]
  CIE_TR_var_M1[1, 5, ,1] <- predict1_rg1$r14[person_time_indices]
  
  CIE_TR_var_M1[2, 4, ,1] <- predict2_rg2$r23[person_time_indices]
  CIE_TR_var_M1[2, 6, ,1] <- predict2_rg2$r24[person_time_indices]
  
  CIE_TR_var_M1[3, 7, ,1] <- predict3_01_rg2$r34[person_time_indices]
  
  CIE_TR_var_M1[4, 8, ,1] <- predict3_012_rg2$r34[person_time_indices]
    
  for (t in 1:tm) {
    for (s_tr in 1:7) {
      CIE_TR_var_M1[s_tr, s_tr, t,1] <- -sum(CIE_TR_var_M1[s_tr, (s_tr + 1):8, t,1])
    }
  }



##rg2 vs. rg1, through M2

##choose to use E=(0,1,1,1)

  CIE_TR_var_M2[1, 2, ,1] <- predict1_rg1$r12[person_time_indices]
  CIE_TR_var_M2[1, 3, ,1] <- predict1_rg1$r13[person_time_indices]  
  CIE_TR_var_M2[1, 5, ,1] <- predict1_rg1$r14[person_time_indices]
  
  CIE_TR_var_M2[2, 4, ,1] <- predict2_rg2$r23[person_time_indices]
  CIE_TR_var_M2[2, 6, ,1] <- predict2_rg2$r24[person_time_indices]
  
  CIE_TR_var_M2[3, 7, ,1] <- predict3_01_rg2$r34[person_time_indices]
  
  CIE_TR_var_M2[4, 8, ,1] <- predict3_012_rg2$r34[person_time_indices]
    
  for (t in 1:tm) {
    for (s_tr in 1:7) {
      CIE_TR_var_M2[s_tr, s_tr, t,1] <- -sum(CIE_TR_var_M2[s_tr, (s_tr + 1):8, t,1])
    }
  }



##############################
###########PIE################
##############################


PIE_TR_var_M1_M2 <- array(0, dim = c(8, 8, tm, 1))  

PIE_TR_var_M1 <- array(0, dim = c(8, 8, tm, 1))  

PIE_TR_var_M2 <- array(0, dim = c(8, 8, tm, 1))  


##rg2 vs. rg1 controlling for M1 and M2

##E=(1,0,0,0)

  PIE_TR_var_M1_M2[1, 2, ,1] <- predict1_rg2$r12[person_time_indices] 
  PIE_TR_var_M1_M2[1, 3, ,1] <- predict1_rg2$r13[person_time_indices] 
  PIE_TR_var_M1_M2[1, 5, ,1] <- predict1_rg2$r14[person_time_indices] 
  
  PIE_TR_var_M1_M2[2, 4, ,1] <- predict2_rg1$r23[person_time_indices] 
  PIE_TR_var_M1_M2[2, 6, ,1] <- predict2_rg1$r24[person_time_indices] 
  
  PIE_TR_var_M1_M2[3, 7, ,1] <- predict3_01_rg1$r34[person_time_indices]   
  PIE_TR_var_M1_M2[4, 8, ,1] <- predict3_012_rg1$r34[person_time_indices]  
    

  for (t in 1:tm) {
    for (s_tr in 1:7) {
      PIE_TR_var_M1_M2[s_tr, s_tr, t,1] <- -sum(PIE_TR_var_M1_M2[s_tr, (s_tr + 1):8, t,1])
    }
  }



##rg2 vs. rg1 controlling for M1

##E=(1,0,0,0)


  PIE_TR_var_M1[1, 2, ,1] <- predict1_rg2$r12[person_time_indices] 
  PIE_TR_var_M1[1, 3, ,1] <- predict1_rg2$r13[person_time_indices]
  PIE_TR_var_M1[1, 5, ,1] <- predict1_rg2$r14[person_time_indices]
  
  PIE_TR_var_M1[2, 4, ,1] <- predict2_rg1$r23[person_time_indices]
  PIE_TR_var_M1[2, 6, ,1] <- predict2_rg1$r24[person_time_indices]
  
  PIE_TR_var_M1[3, 7, ,1] <- predict3_01_rg1$r34[person_time_indices]
  
  PIE_TR_var_M1[4, 8, ,1] <- predict3_012_rg1$r34[person_time_indices]
    
  for (t in 1:tm) {
    for (s_tr in 1:7) {
      PIE_TR_var_M1[s_tr, s_tr, t,1] <- -sum(PIE_TR_var_M1[s_tr, (s_tr + 1):8, t,1])
    }
  }



##rg2 vs. rg1 controlling for M2

##E=(1,0,0,0)


  PIE_TR_var_M2[1, 2, ,1] <- predict1_rg2$r12[person_time_indices]
  PIE_TR_var_M2[1, 3, ,1] <- predict1_rg2$r13[person_time_indices]  ##change to rg1
  PIE_TR_var_M2[1, 5, ,1] <- predict1_rg2$r14[person_time_indices]
  
  PIE_TR_var_M2[2, 4, ,1] <- predict2_rg1$r23[person_time_indices]
  PIE_TR_var_M2[2, 6, ,1] <- predict2_rg1$r24[person_time_indices]
  
  PIE_TR_var_M2[3, 7, ,1] <- predict3_01_rg1$r34[person_time_indices]
  
  PIE_TR_var_M2[4, 8, ,1] <- predict3_012_rg1$r34[person_time_indices]
    
  for (t in 1:tm) {
    for (s_tr in 1:7) {
      PIE_TR_var_M2[s_tr, s_tr, t,1] <- -sum(PIE_TR_var_M2[s_tr, (s_tr + 1):8, t,1])
    }
  }


############################################################
###########Survival risk and survival time################
############################################################


SOP_var    <- array(0, dim = c(tm, 8, num_rg))

SSLY       <- array(0, dim = c(tm, 8, num_rg))


#for each individual, SOP at t=0 all start from S1 state

SOP_var_0 <- c(1,0,0,0,0,0,0,0)


##estimate SOP, SSLY in year 1

for (rg in 1: num_rg) {

SOP_var[1, ,rg]<-SOP_var_0 + SOP_var_0%*%TR_var[, ,1, rg]

  SSLY[1, , rg]<- SOP_var[1, , rg]
}


for (t in 2:tm) {
  for (rg in 1: num_rg) {
    
##once death probability> define_p_death, assume dead

##If flagged dead at t-1, skip recalculation, all following SOP and SSLY remains 0

   if (t > 2 && sum(SOP_var[t-1, 5:8, rg]) >= define_p_death) next
    
     SOP_var[t, , rg] <- SOP_var[t-1, , rg]+ SOP_var[t-1, ,rg]%*%TR_var[, ,t,rg]

    SSLY[t, , rg] <- SSLY[t-1, , rg] + SOP_var[t, , rg]

  }
}


  
valid_years_rg1 <- which(rowSums(SSLY[, , 1]) > 0)   ##identify years before death
last_active_t_rg1 <- max(valid_years_rg1)            ##this is the last year before death
SSLY_person_rg1[i, ] <- SSLY[last_active_t_rg1, , 1] ##extract this last year to estimate SSLY
SOP_person_rg1[i, ] <- SOP_var[last_active_t_rg1, , 1] ##extract this last year to estimate SSLY


valid_years_rg2 <- which(rowSums(SSLY[, , 2]) > 0)   ##identify years before death
last_active_t_rg2 <- max(valid_years_rg2)            ##this is the last year before death
SSLY_person_rg2[i, ] <- SSLY[last_active_t_rg2, , 2] ##extract this last year to estimate SSLY
SOP_person_rg2[i, ] <- SOP_var[last_active_t_rg2, , 2] ##extract this last year to estimate SSLY



#######CIE through M1 and M2


CIE_SOP_var_M1_M2    <- array(0, dim = c(tm, 8, 1))

CIE_SSLY_M1_M2       <- array(0, dim = c(tm, 8, 1))


#for each individual, SOP at t=0 all start from S1 state

CIE_SOP_var_0_M1_M2 <- c(1,0,0,0,0,0,0,0)


##estimate SOP, SSLY in year 1

for (rg in 1:1) {

CIE_SOP_var_M1_M2[1, ,rg] <- CIE_SOP_var_0_M1_M2 + CIE_SOP_var_0_M1_M2%*%CIE_TR_var_M1_M2[, ,1,rg]

  CIE_SSLY_M1_M2[1, , rg]<- CIE_SOP_var_M1_M2[1, , rg]
}


for (t in 2:tm) {
  for (rg in 1:1) {
    
##once death probability> define_p_death, assume dead

##If flagged dead at t-1, skip recalculation, all following SOP and SSLY remains 0

   if (t > 2 && sum(CIE_SOP_var_M1_M2[t-1, 5:8, rg]) >= define_p_death) next
    
     CIE_SOP_var_M1_M2[t, , rg] <- CIE_SOP_var_M1_M2[t-1, , rg]+ CIE_SOP_var_M1_M2[t-1, ,rg]%*%CIE_TR_var_M1_M2[, ,t,rg]

    CIE_SSLY_M1_M2[t, , rg] <- CIE_SSLY_M1_M2[t-1, , rg] + CIE_SOP_var_M1_M2[t, , rg]

  }
}

  
CIE_valid_years_rg2_M1_M2 <- which(rowSums(CIE_SSLY_M1_M2[, , 1]) > 0)   ##identify years before death
CIE_last_active_t_rg2_M1_M2 <- max(CIE_valid_years_rg2_M1_M2)            ##this is the last year before death
CIE_SSLY_person_rg2_M1_M2[i, ] <- CIE_SSLY_M1_M2[CIE_last_active_t_rg2_M1_M2, , 1] ##extract this last year to estimate SSLY
CIE_SOP_person_rg2_M1_M2[i, ] <- CIE_SOP_var_M1_M2[CIE_last_active_t_rg2_M1_M2, , 1] ##extract this last year to estimate SSLY



####CIE through M1


CIE_SOP_var_M1    <- array(0, dim = c(tm, 8, 1))

CIE_SSLY_M1       <- array(0, dim = c(tm, 8, 1))


#for each individual, SOP at t=0 all start from S1 state

CIE_SOP_var_0_M1 <- c(1,0,0,0,0,0,0,0)


##estimate SOP, SSLY in year 1

for (rg in 1:1) {

CIE_SOP_var_M1[1, ,rg] <- CIE_SOP_var_0_M1 + CIE_SOP_var_0_M1%*%CIE_TR_var_M1[, ,1,rg]

  CIE_SSLY_M1[1, , rg]<- CIE_SOP_var_M1[1, , rg]
}


for (t in 2:tm) {
  for (rg in 1:1) {
    
##once death probability> define_p_death, assume dead

##If flagged dead at t-1, skip recalculation, all following SOP and SSLY remains 0

   if (t > 2 && sum(CIE_SOP_var_M1[t-1, 5:8, rg]) >= define_p_death) next
    
     CIE_SOP_var_M1[t, , rg] <- CIE_SOP_var_M1[t-1, , rg]+ CIE_SOP_var_M1[t-1, ,rg]%*%CIE_TR_var_M1[, ,t,rg]

    CIE_SSLY_M1[t, , rg] <- CIE_SSLY_M1[t-1, , rg] + CIE_SOP_var_M1[t, , rg]

  }
}

  
CIE_valid_years_rg2_M1 <- which(rowSums(CIE_SSLY_M1[, , 1]) > 0)   ##identify years before death
CIE_last_active_t_rg2_M1 <- max(CIE_valid_years_rg2_M1)            ##this is the last year before death
CIE_SSLY_person_rg2_M1[i, ] <- CIE_SSLY_M1[CIE_last_active_t_rg2_M1, , 1] ##extract this last year to estimate SSLY
CIE_SOP_person_rg2_M1[i, ] <- CIE_SOP_var_M1[CIE_last_active_t_rg2_M1, , 1] ##extract this last year to estimate SSLY



###CIE through M2



CIE_SOP_var_M2    <- array(0, dim = c(tm, 8, 1))

CIE_SSLY_M2       <- array(0, dim = c(tm, 8, 1))


#for each individual, SOP at t=0 all start from S1 state

CIE_SOP_var_0_M2 <- c(1,0,0,0,0,0,0,0)


##estimate SOP, SSLY in year 1

for (rg in 1:1) {

CIE_SOP_var_M2[1, ,rg] <- CIE_SOP_var_0_M2 + CIE_SOP_var_0_M2%*%CIE_TR_var_M2[, ,1,rg]

  CIE_SSLY_M2[1, , rg]<- CIE_SOP_var_M2[1, , rg]
}


for (t in 2:tm) {
  for (rg in 1:1) {
    
##once death probability> define_p_death, assume dead

##If flagged dead at t-1, skip recalculation, all following SOP and SSLY remains 0

   if (t > 2 && sum(CIE_SOP_var_M2[t-1, 5:8, rg]) >= define_p_death) next
    
     CIE_SOP_var_M2[t, , rg] <- CIE_SOP_var_M2[t-1, , rg]+ CIE_SOP_var_M2[t-1, ,rg]%*%CIE_TR_var_M2[, ,t,rg]

    CIE_SSLY_M2[t, , rg] <- CIE_SSLY_M2[t-1, , rg] + CIE_SOP_var_M2[t, , rg]

  }
}

  
CIE_valid_years_rg2_M2 <- which(rowSums(CIE_SSLY_M2[, , 1]) > 0)   ##identify years before death
CIE_last_active_t_rg2_M2 <- max(CIE_valid_years_rg2_M2)            ##this is the last year before death
CIE_SSLY_person_rg2_M2[i, ] <- CIE_SSLY_M2[CIE_last_active_t_rg2_M2, , 1] ##extract this last year to estimate SSLY
CIE_SOP_person_rg2_M2[i, ] <- CIE_SOP_var_M2[CIE_last_active_t_rg2_M2, , 1] ##extract this last year to estimate SSLY



#######PIE through M1 and M2


PIE_SOP_var_M1_M2    <- array(0, dim = c(tm, 8, 1))

PIE_SSLY_M1_M2       <- array(0, dim = c(tm, 8, 1))


#for each individual, SOP at t=0 all start from S1 state

PIE_SOP_var_0_M1_M2 <- c(1,0,0,0,0,0,0,0)


##estimate SOP, SSLY in year 1

for (rg in 1:1) {

PIE_SOP_var_M1_M2[1, ,rg] <- PIE_SOP_var_0_M1_M2 + PIE_SOP_var_0_M1_M2%*%PIE_TR_var_M1_M2[, ,1,rg]

  PIE_SSLY_M1_M2[1, , rg]<- PIE_SOP_var_M1_M2[1, , rg]
}


for (t in 2:tm) {
  for (rg in 1:1) {
    
##once death probability> define_p_death, assume dead

##If flagged dead at t-1, skip recalculation, all following SOP and SSLY remains 0

   if (t > 2 && sum(PIE_SOP_var_M1_M2[t-1, 5:8, rg]) >= define_p_death) next
    
     PIE_SOP_var_M1_M2[t, , rg] <- PIE_SOP_var_M1_M2[t-1, , rg]+ PIE_SOP_var_M1_M2[t-1, ,rg]%*%PIE_TR_var_M1_M2[, ,t,rg]

    PIE_SSLY_M1_M2[t, , rg] <- PIE_SSLY_M1_M2[t-1, , rg] + PIE_SOP_var_M1_M2[t, , rg]

  }
}

  
PIE_valid_years_rg2_M1_M2 <- which(rowSums(PIE_SSLY_M1_M2[, , 1]) > 0)   ##identify years before death
PIE_last_active_t_rg2_M1_M2 <- max(PIE_valid_years_rg2_M1_M2)            ##this is the last year before death
PIE_SSLY_person_rg2_M1_M2[i, ] <- PIE_SSLY_M1_M2[PIE_last_active_t_rg2_M1_M2, , 1] ##extract this last year to estimate SSLY
PIE_SOP_person_rg2_M1_M2[i, ] <- PIE_SOP_var_M1_M2[PIE_last_active_t_rg2_M1_M2, , 1] ##extract this last year to estimate SSLY



####PIE through M1


PIE_SOP_var_M1    <- array(0, dim = c(tm, 8, 1))

PIE_SSLY_M1       <- array(0, dim = c(tm, 8, 1))


#for each individual, SOP at t=0 all start from S1 state

PIE_SOP_var_0_M1 <- c(1,0,0,0,0,0,0,0)


##estimate SOP, SSLY in year 1

for (rg in 1:1) {

PIE_SOP_var_M1[1, ,rg] <- PIE_SOP_var_0_M1 + PIE_SOP_var_0_M1%*%PIE_TR_var_M1[, ,1,rg]

  PIE_SSLY_M1[1, , rg]<- PIE_SOP_var_M1[1, , rg]
}


for (t in 2:tm) {
  for (rg in 1:1) {
    
##once death probability> define_p_death, assume dead

##If flagged dead at t-1, skip recalculation, all following SOP and SSLY remains 0

   if (t > 2 && sum(PIE_SOP_var_M1[t-1, 5:8, rg]) >= define_p_death) next
    
     PIE_SOP_var_M1[t, , rg] <- PIE_SOP_var_M1[t-1, , rg]+ PIE_SOP_var_M1[t-1, ,rg]%*%PIE_TR_var_M1[, ,t,rg]

    PIE_SSLY_M1[t, , rg] <- PIE_SSLY_M1[t-1, , rg] + PIE_SOP_var_M1[t, , rg]

  }
}

  
PIE_valid_years_rg2_M1 <- which(rowSums(PIE_SSLY_M1[, , 1]) > 0)   ##identify years before death
PIE_last_active_t_rg2_M1 <- max(PIE_valid_years_rg2_M1)            ##this is the last year before death
PIE_SSLY_person_rg2_M1[i, ] <- PIE_SSLY_M1[PIE_last_active_t_rg2_M1, , 1] ##extract this last year to estimate SSLY
PIE_SOP_person_rg2_M1[i, ] <- PIE_SOP_var_M1[PIE_last_active_t_rg2_M1, , 1] ##extract this last year to estimate SSLY



###PIE through M2



PIE_SOP_var_M2    <- array(0, dim = c(tm, 8, 1))

PIE_SSLY_M2       <- array(0, dim = c(tm, 8, 1))


#for each individual, SOP at t=0 all start from S1 state

PIE_SOP_var_0_M2 <- c(1,0,0,0,0,0,0,0)


##estimate SOP, SSLY in year 1

for (rg in 1:1) {

PIE_SOP_var_M2[1, ,rg] <- PIE_SOP_var_0_M2 + PIE_SOP_var_0_M2%*%PIE_TR_var_M2[, ,1,rg]

  PIE_SSLY_M2[1, , rg]<- PIE_SOP_var_M2[1, , rg]
}


for (t in 2:tm) {
  for (rg in 1:1) {
    
##once death probability> define_p_death, assume dead

##If flagged dead at t-1, skip recalculation, all following SOP and SSLY remains 0

   if (t > 2 && sum(PIE_SOP_var_M2[t-1, 5:8, rg]) >= define_p_death) next
    
     PIE_SOP_var_M2[t, , rg] <- PIE_SOP_var_M2[t-1, , rg]+ PIE_SOP_var_M2[t-1, ,rg]%*%PIE_TR_var_M2[, ,t,rg]

    PIE_SSLY_M2[t, , rg] <- PIE_SSLY_M2[t-1, , rg] + PIE_SOP_var_M2[t, , rg]

  }
}

  
PIE_valid_years_rg2_M2 <- which(rowSums(PIE_SSLY_M2[, , 1]) > 0)   ##identify years before death
PIE_last_active_t_rg2_M2 <- max(PIE_valid_years_rg2_M2)            ##this is the last year before death
PIE_SSLY_person_rg2_M2[i, ] <- PIE_SSLY_M2[PIE_last_active_t_rg2_M2, , 1] ##extract this last year to estimate SSLY
PIE_SOP_person_rg2_M2[i, ] <- PIE_SOP_var_M2[PIE_last_active_t_rg2_M2, , 1] ##extract this last year to estimate SSLY



}     ####### end of loop for each individual




##########################################
###########mediation analysis#############
##########################################


##mediation effect for each person

NDE_person_rg1<-as.matrix(SOP_person_rg1[,5])
TE_person_rg1 <- as.matrix(rowSums(SOP_person_rg1[, 5:8]))

NIE_m1_person_rg1<-as.matrix(SOP_person_rg1[,6])
NIE_m2_person_rg1<-as.matrix(SOP_person_rg1[,7])
NIE_m1_m2_person_rg1 <-as.matrix(SOP_person_rg1[,8])


NDE_person_rg2<-as.matrix(SOP_person_rg2[,5])
TE_person_rg2 <- as.matrix(rowSums(SOP_person_rg2[, 5:8]))

NIE_m1_person_rg2<-as.matrix(SOP_person_rg2[,6])
NIE_m2_person_rg2<-as.matrix(SOP_person_rg2[,7])
NIE_m1_m2_person_rg2 <-as.matrix(SOP_person_rg2[,8])


NDE_marginal_rg1[simu]<- colMeans(NDE_person_rg1) 
TE_marginal_rg1[simu]<- colMeans(TE_person_rg1) 

NIE_m1_marginal_rg1[simu]<- colMeans(NIE_m1_person_rg1) 
NIE_m2_marginal_rg1[simu]<- colMeans(NIE_m2_person_rg1) 
NIE_m1_m2_marginal_rg1[simu]<- colMeans(NIE_m1_m2_person_rg1) 


NDE_marginal_rg2[simu]<- colMeans(NDE_person_rg2) 
TE_marginal_rg2[simu]<- colMeans(TE_person_rg2) 

NIE_m1_marginal_rg2[simu]<- colMeans(NIE_m1_person_rg2) 
NIE_m2_marginal_rg2[simu]<- colMeans(NIE_m2_person_rg2) 
NIE_m1_m2_marginal_rg2[simu]<- colMeans(NIE_m1_m2_person_rg2) 



############################
###CIE combined for M1--M2
############################


CIE_m1_person_rg2 <- as.matrix(CIE_SOP_person_rg2_M1[, 6]+ CIE_SOP_person_rg2_M1[, 8])
CIE_m2_person_rg2 <- as.matrix(CIE_SOP_person_rg2_M2[, 7])
CIE_m1_m2_person_rg2 <- as.matrix(CIE_SOP_person_rg2_M1_M2[, 8])


CIE_m1_marginal_rg2[simu]<- colMeans(CIE_m1_person_rg2) 
CIE_m2_marginal_rg2[simu]<- colMeans(CIE_m2_person_rg2) 
CIE_m1_m2_marginal_rg2[simu]<- colMeans(CIE_m1_m2_person_rg2) 


PIE_m1_person_rg2 <- as.matrix(PIE_SOP_person_rg2_M1[, 6])
PIE_m2_person_rg2 <- as.matrix(PIE_SOP_person_rg2_M2[, 7])
PIE_m1_m2_person_rg2 <- as.matrix(PIE_SOP_person_rg2_M1_M2[, 8])


PIE_m1_marginal_rg2[simu]<- colMeans(PIE_m1_person_rg2) 
PIE_m2_marginal_rg2[simu]<- colMeans(PIE_m2_person_rg2) 
PIE_m1_m2_marginal_rg2[simu]<- colMeans(PIE_m1_m2_person_rg2) 



##calculate mediation effect


TE_person_rg21 <- as.matrix(TE_person_rg2-TE_person_rg1)


NDE_person_rg21 <- as.matrix(NDE_person_rg2-NDE_person_rg1)


NIE_m1_person_rg21 <- as.matrix(NIE_m1_person_rg2-NIE_m1_person_rg1)
NIE_m2_person_rg21 <- as.matrix(NIE_m2_person_rg2-NIE_m2_person_rg1)
NIE_m1_m2_person_rg21 <- as.matrix(NIE_m1_m2_person_rg2-NIE_m1_m2_person_rg1)


############################
###CIE combined for M1--M2
############################

CIE_m1_person_rg21 <- as.matrix(NIE_m1_person_rg2+NIE_m1_m2_person_rg2-CIE_m1_person_rg2)
CIE_m2_person_rg21 <- as.matrix(NIE_m2_person_rg2-CIE_m2_person_rg2)
CIE_m1_m2_person_rg21 <- as.matrix(NIE_m1_m2_person_rg2-CIE_m1_m2_person_rg2)


PIE_m1_person_rg21 <- as.matrix(PIE_m1_person_rg2-NIE_m1_person_rg1)
PIE_m2_person_rg21 <- as.matrix(PIE_m2_person_rg2-NIE_m2_person_rg1)
PIE_m1_m2_person_rg21 <- as.matrix(PIE_m1_m2_person_rg2-NIE_m1_m2_person_rg1)


TE_marginal_rg21[simu]<- colMeans(TE_person_rg21) 
NDE_marginal_rg21[simu]<- colMeans(NDE_person_rg21) 

NIE_m1_marginal_rg21[simu]<- colMeans(NIE_m1_person_rg21) 
NIE_m2_marginal_rg21[simu]<- colMeans(NIE_m2_person_rg21) 
NIE_m1_m2_marginal_rg21[simu]<- colMeans(NIE_m1_m2_person_rg21) 

NIE_marginal_rg21[simu]<- NIE_m1_marginal_rg21[simu]+NIE_m2_marginal_rg21[simu]+NIE_m1_m2_marginal_rg21[simu]


CIE_m1_marginal_rg21[simu]<- colMeans(CIE_m1_person_rg21) 
CIE_m2_marginal_rg21[simu]<- colMeans(CIE_m2_person_rg21) 
CIE_m1_m2_marginal_rg21[simu]<- colMeans(CIE_m1_m2_person_rg21) 

CIE_marginal_rg21[simu]<- CIE_m1_marginal_rg21[simu]+CIE_m2_marginal_rg21[simu]+CIE_m1_m2_marginal_rg21[simu]


PIE_m1_marginal_rg21[simu]<- colMeans(PIE_m1_person_rg21) 
PIE_m2_marginal_rg21[simu]<- colMeans(PIE_m2_person_rg21) 
PIE_m1_m2_marginal_rg21[simu]<- colMeans(PIE_m1_m2_person_rg21) 

PIE_marginal_rg21[simu]<- PIE_m1_marginal_rg21[simu]+PIE_m2_marginal_rg21[simu]+PIE_m1_m2_marginal_rg21[simu]


##percentage, %

NDE_100_person_rg21 <- as.matrix(100*NDE_person_rg21/TE_person_rg21)

NIE_100_m1_person_rg21 <- as.matrix(100*NIE_m1_person_rg21/TE_person_rg21)
NIE_100_m2_person_rg21 <- as.matrix(100*NIE_m2_person_rg21/TE_person_rg21)
NIE_100_m1_m2_person_rg21 <- as.matrix(100*NIE_m1_m2_person_rg21/TE_person_rg21)


PIE_100_m1_person_rg21 <- as.matrix(100*PIE_m1_person_rg21/TE_person_rg21)
PIE_100_m2_person_rg21 <- as.matrix(100*PIE_m2_person_rg21/TE_person_rg21)
PIE_100_m1_m2_person_rg21 <- as.matrix(100*PIE_m1_m2_person_rg21/TE_person_rg21)


CIE_100_m1_person_rg21 <- as.matrix(100*CIE_m1_person_rg21/TE_person_rg21)
CIE_100_m2_person_rg21 <- as.matrix(100*CIE_m2_person_rg21/TE_person_rg21)
CIE_100_m1_m2_person_rg21 <- as.matrix(100*CIE_m1_m2_person_rg21/TE_person_rg21)


#allow the TE to be 0

NDE_100_marginal_rg21[simu]<- colMeans(NDE_100_person_rg21, na.rm=T) 

NIE_100_m1_marginal_rg21[simu]<- colMeans(NIE_100_m1_person_rg21, na.rm=T) 
NIE_100_m2_marginal_rg21[simu]<- colMeans(NIE_100_m2_person_rg21, na.rm=T) 
NIE_100_m1_m2_marginal_rg21[simu]<- colMeans(NIE_100_m1_m2_person_rg21, na.rm=T) 

NIE_100_marginal_rg21[simu]= NIE_100_m1_marginal_rg21[simu]+ NIE_100_m2_marginal_rg21[simu]+ NIE_100_m1_m2_marginal_rg21[simu]


PIE_100_m1_marginal_rg21[simu]<- colMeans(PIE_100_m1_person_rg21, na.rm=T) 
PIE_100_m2_marginal_rg21[simu]<- colMeans(PIE_100_m2_person_rg21, na.rm=T) 
PIE_100_m1_m2_marginal_rg21[simu]<- colMeans(PIE_100_m1_m2_person_rg21, na.rm=T) 

PIE_100_marginal_rg21[simu]= PIE_100_m1_marginal_rg21[simu]+ PIE_100_m2_marginal_rg21[simu]+ PIE_100_m1_m2_marginal_rg21[simu]

CIE_100_m1_marginal_rg21[simu]<- colMeans(CIE_100_m1_person_rg21, na.rm=T) 
CIE_100_m2_marginal_rg21[simu]<- colMeans(CIE_100_m2_person_rg21, na.rm=T) 
CIE_100_m1_m2_marginal_rg21[simu]<- colMeans(CIE_100_m1_m2_person_rg21, na.rm=T) 

CIE_100_marginal_rg21[simu]= CIE_100_m1_marginal_rg21[simu]+ CIE_100_m2_marginal_rg21[simu]+ CIE_100_m1_m2_marginal_rg21[simu]


#####################################################
###########add subgroup analysis here if needed#####
#####################################################

## to evaluate the effect modification by gender or race
## restrict LE_marginal_rg1 among subgroups


  print(paste("Bootstrap iteration:", simu))


}     ###End of bootstrap 



##########################################
###########mediation analysis#############
##########################################


####mediation analysis


##survival risk

NDE_final_rg1 <- quantile(NDE_marginal_rg1, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
NIE_m1_final_rg1 <- quantile(NIE_m1_marginal_rg1, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
NIE_m2_final_rg1 <- quantile(NIE_m2_marginal_rg1, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
NIE_m1_m2_final_rg1 <- quantile(NIE_m1_m2_marginal_rg1, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
TE_final_rg1 <- quantile(TE_marginal_rg1, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)

NDE_final_rg2 <- quantile(NDE_marginal_rg2, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
NIE_m1_final_rg2 <- quantile(NIE_m1_marginal_rg2, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
NIE_m2_final_rg2 <- quantile(NIE_m2_marginal_rg2, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
NIE_m1_m2_final_rg2 <- quantile(NIE_m1_m2_marginal_rg2, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
TE_final_rg2 <- quantile(TE_marginal_rg2, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)


##mediation effect

NDE_final_rg21 <- quantile(NDE_marginal_rg21, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
TE_final_rg21 <- quantile(TE_marginal_rg21, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)

NIE_final_rg21 <- quantile(NIE_marginal_rg21, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)

NIE_m1_final_rg21 <- quantile(NIE_m1_marginal_rg21, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
NIE_m2_final_rg21 <- quantile(NIE_m2_marginal_rg21, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
NIE_m1_m2_final_rg21 <- quantile(NIE_m1_m2_marginal_rg21, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)


CIE_final_rg21 <- quantile(CIE_marginal_rg21, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)

CIE_m1_final_rg21 <- quantile(CIE_m1_marginal_rg21, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
CIE_m2_final_rg21 <- quantile(CIE_m2_marginal_rg21, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
CIE_m1_m2_final_rg21 <- quantile(CIE_m1_m2_marginal_rg21, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)


PIE_final_rg21 <- quantile(PIE_marginal_rg21, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)

PIE_m1_final_rg21 <- quantile(PIE_m1_marginal_rg21, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
PIE_m2_final_rg21 <- quantile(PIE_m2_marginal_rg21, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
PIE_m1_m2_final_rg21 <- quantile(PIE_m1_m2_marginal_rg21, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)


##proportion

NDE_100_final_rg21 <- quantile(NDE_100_marginal_rg21, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)

NIE_100_final_rg21 <- quantile(NIE_100_marginal_rg21, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)

NIE_100_m1_final_rg21 <- quantile(NIE_100_m1_marginal_rg21, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
NIE_100_m2_final_rg21 <- quantile(NIE_100_m2_marginal_rg21, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
NIE_100_m1_m2_final_rg21 <- quantile(NIE_100_m1_m2_marginal_rg21, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)


PIE_100_final_rg21 <- quantile(PIE_100_marginal_rg21, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)

PIE_100_m1_final_rg21 <- quantile(PIE_100_m1_marginal_rg21, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
PIE_100_m2_final_rg21 <- quantile(PIE_100_m2_marginal_rg21, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
PIE_100_m1_m2_final_rg21 <- quantile(PIE_100_m1_m2_marginal_rg21, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)


CIE_100_final_rg21 <- quantile(CIE_100_marginal_rg21, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)

CIE_100_m1_final_rg21 <- quantile(CIE_100_m1_marginal_rg21, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
CIE_100_m2_final_rg21 <- quantile(CIE_100_m2_marginal_rg21, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
CIE_100_m1_m2_final_rg21 <- quantile(CIE_100_m1_m2_marginal_rg21, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)



####organize into tables

##TE, NDE, NIE

survival_risk=100*rbind(
TE_final_rg1, NDE_final_rg1, NIE_m1_final_rg1, NIE_m2_final_rg1, NIE_m1_m2_final_rg1,
TE_final_rg2, NDE_final_rg2, NIE_m1_final_rg2, NIE_m2_final_rg2, NIE_m1_m2_final_rg2)


NIE_table =100*rbind(
TE_final_rg21, NDE_final_rg21, NIE_final_rg21, 
NIE_m1_final_rg21, NIE_m2_final_rg21, NIE_m1_m2_final_rg21)


NIE_100 =rbind(
NDE_100_final_rg21, NIE_100_final_rg21,
NIE_100_m1_final_rg21, NIE_100_m2_final_rg21, NIE_100_m1_m2_final_rg21)



##PIE

PIE_table=100*rbind(PIE_final_rg21, PIE_m1_final_rg21, PIE_m2_final_rg21, PIE_m1_m2_final_rg21)

PIE_100=rbind(PIE_100_final_rg21, PIE_100_m1_final_rg21, PIE_100_m2_final_rg21, PIE_100_m1_m2_final_rg21)


##CIE

CIE_table=100*rbind(CIE_m1_final_rg21, CIE_m2_final_rg21, CIE_m1_m2_final_rg21)

CIE_100=rbind(CIE_100_m1_final_rg21, CIE_100_m2_final_rg21, CIE_100_m1_m2_final_rg21)


##output table


survival_risk

NIE_table

PIE_table

CIE_table


NIE_100

PIE_100

CIE_100



