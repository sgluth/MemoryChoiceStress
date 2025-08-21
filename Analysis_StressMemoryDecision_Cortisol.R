################################################
### Project: Stress and memory-based decisions #
### (with Lars Schwabe) #
####################################
### Statistical analyses #
####################################
### sebastian.gluth@uni-hamburg.de #
####################################


### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### ---
### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- 
### --- Preparations
### --- ### --- ### --- ###

#clear working environment
rm(list=ls())

#clear all plots
if(!is.null(dev.list())) dev.off()

#load required libraries
library(psych)
library(readxl)
library(DescTools)
library(R2jags)
library(rtdists)
#parallel computing stuff
library(parallel)
library(doParallel)
library(foreach)
numCores <- detectCores()
registerDoParallel(cores=numCores)

#set and go to working direction
#setwd("/Users/Gluth/Arbeit/Projects/StressMemoryChoice/Analysis/")
setwd("/home/bal1607/Stress/Analysis/")

### --- ### --- ### --- ###
### --- Load files (preprocessed in Matlab)
### --- ### --- ### --- ###

#Rating
rating <- read.csv(file='RATING.csv')

#Encoding
encoding <- read.csv(file='ENCODING.csv')

#Decisions
decision <- read.csv(file='DECISION.csv')

#Cued Recall
cued_recall <- read.csv(file='CUED_RECALL.csv')


### --- ### --- ### --- ###
### --- Get some basic information (e.g., # of subjects)
### --- ### --- ### --- ###
subject_ID <- unique(encoding$SubjectID)
nsubj <- length(subject_ID)
ntrials_R <- nrow(rating)/nsubj
ntrials_E <- nrow(encoding)/nsubj
ntrials_D <- nrow(decision)/nsubj
ntrials_C <- nrow(cued_recall)/nsubj


### --- ### --- ### --- ###
### --- For exclusion criteria: choice consistency
### --- ### --- ### --- ###
consistency <- data.frame(matrix(nrow = nsubj,ncol = 3))
colnames(consistency) <- c('overall','memory','control')
above_chance <- data.frame(matrix(nrow = nsubj,ncol = 3))
colnames(above_chance) <- c('overall','memory','control')
below_chance <- data.frame(matrix(nrow = nsubj,ncol = 3))
colnames(below_chance) <- c('overall','memory','control')
for (s in 1:nsubj){
  #get decision data of current subject
  cdata <- decision[decision$SubjectID==subject_ID[s],]
  #calculate whether decisions are correct (=higher rated item chosen)
  # vL_vs_vR <- (cdata$LeftAvgRatingR1R2R3>cdata$RightAvgRatingR1R2R3)-(cdata$LeftAvgRatingR1R2R3<cdata$RightAvgRatingR1R2R3)
  vL_vs_vR <- (cdata$LeftAvgRatingR1R2>cdata$RightAvgRatingR1R2)-(cdata$LeftAvgRatingR1R2<cdata$RightAvgRatingR1R2)
  choice_correct <- ((cdata$Choice==1)&(vL_vs_vR==1))|((cdata$Choice==2)&(vL_vs_vR==-1))
  #calculate consistency over all trials
  elig_trials <- vL_vs_vR!=0 #only trials with actual rating difference are eligible
  trial_acc_overall <- choice_correct[elig_trials]
  consistency$overall[s] <- mean(trial_acc_overall)
  #calculate consistency per condition (memory, control)
  trial_acc_memory <- choice_correct[elig_trials&(cdata$IsControlBlock==0)]
  trial_acc_control <- choice_correct[elig_trials&(cdata$IsControlBlock==1)]
  consistency$memory[s] <- mean(trial_acc_memory)
  consistency$control[s] <- mean(trial_acc_control)
  #test whether individual data is above chance
  bt_acc_overall <- binom.test(sum(trial_acc_overall),sum(elig_trials),alternative = "greater")
  bt_acc_memory <- binom.test(sum(trial_acc_memory),sum(elig_trials&(cdata$IsControlBlock==0)),alternative = "greater")
  bt_acc_control <- binom.test(sum(trial_acc_control),sum(elig_trials&(cdata$IsControlBlock==1)),alternative = "greater")
  above_chance$overall[s] <- bt_acc_overall$p.value<.05
  above_chance$memory[s] <- bt_acc_memory$p.value<.05
  above_chance$control[s] <- bt_acc_control$p.value<.05
  #test whether individual data is below chance
  bt_acc_overall <- binom.test(sum(trial_acc_overall),sum(elig_trials),alternative = "less")
  bt_acc_memory <- binom.test(sum(trial_acc_memory),sum(elig_trials&(cdata$IsControlBlock==0)),alternative = "less")
  bt_acc_control <- binom.test(sum(trial_acc_control),sum(elig_trials&(cdata$IsControlBlock==1)),alternative = "less")
  below_chance$overall[s] <- bt_acc_overall$p.value<(.05/nsubj) #Bonferroni corrected for the total number of to-be-tested subjects
  below_chance$memory[s] <- bt_acc_memory$p.value<(.05/nsubj) #Bonferroni corrected for the total number of to-be-tested subjects
  below_chance$control[s] <- bt_acc_control$p.value<(.05/nsubj) #Bonferroni corrected for the total number of to-be-tested subjects
}


### --- ### --- ### --- ###
### --- For exclusion criteria: checking the proportion of min and max ratings
### --- ### --- ### --- ###
min_max_rating = c(0,1000) #lowest and highest possible rating (according to internal scale)
extremeRatings_R1R2 <- data.frame(matrix(nrow = nsubj,ncol = 2))
colnames(extremeRatings_R1R2) <- c('min','max')
for (s in 1:nsubj){
  #get rating data of current subject
  cdata <- rating[rating$SubjectID==subject_ID[s],]
  extremeRatings_R1R2$min[s] <- mean(cdata$AvgRatingR1R2==min_max_rating[1])
  extremeRatings_R1R2$max[s] <- mean(cdata$AvgRatingR1R2==min_max_rating[2])
}


### --- ### --- ### --- ###
### --- For exclusion criteria: checking response times
### --- ### --- ### --- ###
RTcrit <- .25
tooFast <- data.frame(matrix(nrow = nsubj, ncol = 3))
colnames(tooFast) <- c('Rating','Decision','Recall')
for (s in 1:nsubj){
  #get rating data of current subject
  cdata <- rating[rating$SubjectID==subject_ID[s],]
  tooFast$Rating[s] <- mean(c(cdata$RespTimeR1,cdata$RespTimeR2,cdata$RespTimeR3)<RTcrit)
  #get decision data of current subject
  cdata <- decision[decision$SubjectID==subject_ID[s],]
  tooFast$Decision[s] <- mean(cdata$RespTime<RTcrit)
  #get recall data of current subject
  cdata <- cued_recall[cued_recall$SubjectID==subject_ID[s],]
  tooFast$Recall[s] <- mean(cdata$RespTimeRating<RTcrit)
}

### --- ### --- ### --- ###
### --- For exclusion criteria: get cortisol data
### --- ### --- ### --- ###
cortDat <- as.data.frame(read_excel("cortMBD.xlsx"))
cortDat <- data.frame(cbind(as.numeric(cortDat$Subject),as.numeric(cortDat$time),as.numeric(cortDat$CORT)))
colnames(cortDat) <- c('Subject','Time','Cortisol')
cortDat$CortLog <- log(cortDat$Cortisol) #data seems to be heavily right-skewed, so log-transform (suggested by Lars) seems plausible

### --- ### --- ### --- ###
### --- For exclusion criteria: get physiological data and subjective ratings (MDBF, STAI)
### --- ### --- ### --- ###

#read in data about heart rate and blood pressure 
GaB <- read_excel("SubjInfoCortisol_AUG25.xlsx")

#read in data from questionnaires
Qd <- read_excel("Vor-Ort-Erhebung_AUG25.xlsx")

#match with other dataset (and stress group assignment)
Qd_m <- c()
Qd_ID <- c()
for (s in 1:length(GaB$Token)){
  m <- which(GaB$Token[s]==Qd$Zugangsschlüssel)
  if (length(m)>0){
    Qd_m[s] <- m
    Qd_ID[s] <- s
  }
}
Qd <- Qd[Qd_m,]
Qd$VpNr <- GaB$VpNr[Qd_ID]
Qd$Group <- GaB$Group[Qd_ID]
Qd_strGr <- as.factor(Qd$Group) #0 = no stress, 1 = stress
Qd_strGr <- Qd_strGr[is.na(Qd_strGr)==FALSE]

### --- ### --- ### --- ###
### --- Exclusion of participants
### --- ### --- ### --- ###

#REASON 1: too many extreme ratings (minimum or maximum possible rating > 40%)
exclusionExtremeRatings <- subject_ID[(extremeRatings_R1R2$min>.4)|(extremeRatings_R1R2$max>.4)]
#REASON 2: extremely low decision consistency in control or memory decisions
#(note: although this exclusion criterion wasn't pre-registered, these participants simply refused to do the task,
# so it should be okay to exclude them)
exclusionPoorConsistency <- subject_ID[(below_chance$control)|(below_chance$memory)]
#REASON 3: too fast decisions in any of the (not directly supervised) tasks
exclusionTooFast <- subject_ID[(tooFast$Rating>.2)|(tooFast$Decision>.2)|(tooFast$Recall>.2)]
#REASON 4: no cortisol data
xCort <- unique(cortDat$Subject[((is.na(cortDat$Cortisol)==TRUE)&(cortDat$Time==1))|((is.na(cortDat$Cortisol)==TRUE)&(cortDat$Time==4))])
exclusionCortisol <- c(xCort[xCort %in% subject_ID],subject_ID[(subject_ID %in% cortDat$Subject)==FALSE])
#REASON 5: no questionnaire data
exclusionQuest <- subject_ID[(subject_ID %in% Qd$VpNr)==FALSE]

#Proceed with data from included participants only (and sort them)
excludedSubj <- unique(c(exclusionExtremeRatings,exclusionPoorConsistency,exclusionTooFast,exclusionCortisol,exclusionQuest))
subject_ID <- sort(subject_ID[(subject_ID %in% excludedSubj)==0])
nsubj <- length(subject_ID)
rating <- rating[(rating$SubjectID %in% excludedSubj)==0,]
encoding <- encoding[(encoding$SubjectID %in% excludedSubj)==0,]
decision <- decision[(decision$SubjectID %in% excludedSubj)==0,]
cued_recall <- cued_recall[(cued_recall$SubjectID %in% excludedSubj)==0,]


### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### ---
### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- 
### --- Stress conditions, demographics
### --- ### --- ### --- ###

#Information about group assignment
stressGroup <- as.factor(GaB$Group) #0 = no stress, 1 = stress
nsubj_Stress <- sum(stressGroup==1)
nsubj_noStress <- sum(stressGroup==0)

#demographics incl. comparison of groups
#age
ageInfo.all <- describe(GaB$Age)
ageInfo.Stress <- describe(GaB$Age[stressGroup==1])
ageInfo.noStress <- describe(GaB$Age[stressGroup==0])
ageGroupComp <- t.test(GaB$Age[stressGroup==1],GaB$Age[stressGroup==0]) #n.s.

#gender
genderInfo <- matrix(data = c(sum(GaB$Gender[stressGroup==1]=="0"),sum(GaB$Gender[stressGroup==1]=="1"),
                              sum(GaB$Gender[stressGroup==0]=="0"),sum(GaB$Gender[stressGroup==0]=="1")),nrow=2,ncol=2)
dimnames(genderInfo) <- list(gender = c("f","m"),condition = c("Stress","noStress"))
genderGroupComp <- chisq.test(genderInfo) #n.s. (obviously)


### --- ### --- ### --- ###
### --- Blood pressure, heart rate, cortisol
### --- ### --- ### --- ###
#blood pressure and heart rate
Sys <- cbind(GaB$T1_Sys,GaB$T2_Sys,GaB$T3_Sys,GaB$T4_Sys,GaB$T5_Sys);colnames(Sys) <- c('T1','T2','T3','T4','T5')
statsSys <- matrix(data=NA,nrow = 5,ncol = 3);colnames(statsSys) <- c('t-value','df','p-value')
ciSys <- matrix(data=NA,nrow = 5,ncol = 2);colnames(ciSys) <- c("Stress","noStress")
Dia <- cbind(GaB$T1_Dia,GaB$T2_Dia,GaB$T3_Dia,GaB$T4_Dia,GaB$T5_Dia);colnames(Dia) <- c('T1','T2','T3','T4','T5')
statsDia <- matrix(data=NA,nrow = 5,ncol = 3);colnames(statsDia) <- c('t-value','df','p-value')
ciDia <- matrix(data=NA,nrow = 5,ncol = 2);colnames(ciDia) <- c("Stress","noStress")
heartRate <- cbind(GaB$T1_Puls,GaB$T2_Puls,GaB$T3_Puls,GaB$T4_Puls,GaB$T5_Puls);colnames(heartRate) <- c('T1','T2','T3','T4','T5')
statsHR <- matrix(data=NA,nrow = 5,ncol = 3);colnames(statsHR) <- c('t-value','df','p-value')
ciHR <- matrix(data=NA,nrow = 5,ncol = 2);colnames(ciHR) <- c("Stress","noStress")
for (T in 1:5){
  stats <- t.test(Sys[stressGroup==1,T],Sys[stressGroup==0,T])
  statsSys[T,] <- c(stats$statistic,stats$parameter,stats$p.value) #significant difference only at T3
  ciSys[T,] <- c(qnorm(.975)*sd(Sys[stressGroup==1,T])/sqrt(nsubj_Stress),
                qnorm(.975)*sd(Sys[stressGroup==0,T])/sqrt(nsubj_noStress))
  stats <- t.test(Dia[stressGroup==1,T],Dia[stressGroup==0,T])
  statsDia[T,] <- c(stats$statistic,stats$parameter,stats$p.value) #significant difference only at T3
  ciDia[T,] <- c(qnorm(.975)*sd(Dia[stressGroup==1,T])/sqrt(nsubj_Stress),
                qnorm(.975)*sd(Dia[stressGroup==0,T])/sqrt(nsubj_noStress))
  stats <- t.test(heartRate[stressGroup==1,T],heartRate[stressGroup==0,T])
  statsHR[T,] <- c(stats$statistic,stats$parameter,stats$p.value) #significant difference only at T3
  ciHR[T,] <- c(qnorm(.975)*sd(heartRate[stressGroup==1,T])/sqrt(nsubj_Stress),
                qnorm(.975)*sd(heartRate[stressGroup==0,T])/sqrt(nsubj_noStress))
}
#Cortisol
CSG <- cortDat$CortLog[which(cortDat$Subject %in% subject_ID[stressGroup==1])]
CCG <- cortDat$CortLog[which(cortDat$Subject %in% subject_ID[stressGroup==0])]
CortStress <- matrix(CSG,ncol = 5,byrow = T);colnames(CortStress) <- c('T1','T2','T3','T4','T5')
CortNostress <- matrix(CCG,ncol = 5,byrow = T);colnames(CortNostress) <- c('T1','T2','T3','T4','T5')
statsCort <- matrix(data=NA,nrow = 5,ncol = 3);colnames(statsCort) <- c('t-value','df','p-value')
ciCort <- matrix(data=NA,nrow = 5,ncol = 2);colnames(ciCort) <- c("Stress","noStress")
for (T in 1:5){
  stats <- t.test(CortStress[,T],CortNostress[,T])
  statsCort[T,] <- c(stats$statistic,stats$parameter,stats$p.value) #significant difference only at T3
  ciCort[T,] <- c(qnorm(.975)*sd(CortStress[,T],na.rm = TRUE)/sqrt(sum(is.na(CortStress[,T])==FALSE)),
                  qnorm(.975)*sd(CortNostress[,T],na.rm = TRUE)/sqrt(sum(is.na(CortNostress[,T])==FALSE)))
}

#plot physiology with 95% CIs
#Systolic
png('Figures/Systolic',width = 600, height = 600)
plot(.95:4.95,colMeans(Sys[stressGroup==1,]),type = 'b',col='red',
     xlim = c(.5,5.5),ylim = c(110,145),xaxt="n",xlab = 'Time point',ylab = 'Systolic blood pressure',lwd=2)
lines(1.05:5.05,colMeans(Sys[stressGroup==0,]),type = 'b',lwd=2)
for (T in 1:5){
  lines(c(T,T)-.05,mean(Sys[stressGroup==1,T])+ciSys[T,1]*c(-1,1),col='red')
  lines(c(T,T)+.05,mean(Sys[stressGroup==0,T])+ciSys[T,2]*c(-1,1))
}
text(3,142,'***',cex=2);text(3,140,'(p < .001)')
axis(1,at = 1:5,labels = c('T1','T2','T3','T4','T5'))
title('Systolic blood pressure')
legend(4,140,c('Stress','No stress'),fill=c('red','black'),bty='n')
dev.off()

#Diastolic
png('Figures/Diastolic',width = 600, height = 600)
plot(.95:4.95,colMeans(Dia[stressGroup==1,]),type = 'b',col='red',
     xlim = c(.5,5.5),ylim = c(60,95),xaxt="n",xlab = 'Time point',ylab = 'Diastolic blood pressure',lwd=2)
lines(1.05:5.05,colMeans(Dia[stressGroup==0,]),type = 'b',lwd=2)
for (T in 1:5){
  lines(c(T,T)-.05,mean(Dia[stressGroup==1,T])+ciDia[T,1]*c(-1,1),col='red')
  lines(c(T,T)+.05,mean(Dia[stressGroup==0,T])+ciDia[T,2]*c(-1,1))
}
text(3,90,'***',cex=2);text(3,88,'(p < .001)')
axis(1,at = 1:5,labels = c('T1','T2','T3','T4','T5'))
title('Diastolic blood pressure')
legend(4,90,c('Stress','No stress'),fill=c('red','black'),bty='n')
dev.off()

#Heart rate
png('Figures/HeartRate',width = 600, height = 600)
plot(.95:4.95,colMeans(heartRate[stressGroup==1,]),type = 'b',col='red',
     xlim = c(.5,5.5),ylim = c(60,80),xaxt="n",xlab = 'Time point',ylab = 'Heart rate',lwd=2)
lines(1.05:5.05,colMeans(heartRate[stressGroup==0,]),type = 'b',lwd=2)
for (T in 1:5){
  lines(c(T,T)-.05,mean(heartRate[stressGroup==1,T])+ciHR[T,1]*c(-1,1),col='red')
  lines(c(T,T)+.05,mean(heartRate[stressGroup==0,T])+ciHR[T,2]*c(-1,1))
}
text(3,78,'*',cex=2);text(3,77,'(p = .026)')
axis(1,at = 1:5,labels = c('T1','T2','T3','T4','T5'))
title('Heart rate')
legend(4,78,c('Stress','No stress'),fill=c('red','black'),bty='n')
dev.off()

#Cortisol
png('Figures/Cortisol',width = 600, height = 600)
plot(.95:4.95,colMeans(CortStress,na.rm = TRUE),type = 'b',col='red',
     xlim = c(.5,5.5),ylim = c(0,2),xaxt="n",xlab = 'Time point',ylab = 'Salivary cortisol (log-transformed)',lwd=2)
lines(1.05:5.05,colMeans(CortNostress,na.rm = TRUE),type = 'b',lwd=2)
for (T in 1:5){
  lines(c(T,T)-.05,mean(CortStress[,T],na.rm = TRUE)+ciCort[T,1]*c(-1,1),col='red')
  lines(c(T,T)+.05,mean(CortNostress[,T],na.rm = TRUE)+ciCort[T,2]*c(-1,1))
}
text(4,1.9,'***',cex=2);text(4,1.8,'(p < .001)')
text(5,1.35,'*',cex=2);text(5,1.25,'(p = .015)')
axis(1,at = 1:5,labels = c('T1','T2','T3','T4','T5'))
title('Cortisol')
legend(4.55,2,c('Stress','No stress'),fill=c('red','black'),bty='n')
dev.off()

### --- ### --- ### --- ###
### --- Subjective ratings (MDBF, STAI)
### --- ### --- ### --- ###

#calculate sum-scores for MDBF and STAI
STAI <- cbind(rowMeans(Qd[,3:22]),rowMeans(Qd[,47:66]),rowMeans(Qd[,91:110]),rowMeans(Qd[,135:154]))
MDBF <- cbind(rowMeans(Qd[,23:46]),rowMeans(Qd[,67:90]),rowMeans(Qd[,111:134]),rowMeans(Qd[,155:178]))

#stats
statsSTAI <- matrix(data=NA,nrow = 4,ncol = 3);colnames(statsSys) <- c('t-value','df','p-value')
ciSTAI <- matrix(data=NA,nrow = 4,ncol = 2);colnames(ciSys) <- c("Stress","noStress")
statsMDBF <- matrix(data=NA,nrow = 4,ncol = 3);colnames(statsSys) <- c('t-value','df','p-value')
ciMDBF <- matrix(data=NA,nrow = 4,ncol = 2);colnames(ciSys) <- c("Stress","noStress")
for (T in 1:4){
  stats <- t.test(STAI[Qd_strGr==1,T],STAI[Qd_strGr==0,T])
  statsSTAI[T,] <- c(stats$statistic,stats$parameter,stats$p.value) #significant difference only at T3
  ciSTAI[T,] <- c(qnorm(.975)*sd(STAI[Qd_strGr==1,T],na.rm = T)/sqrt(sum(Qd_strGr==1,na.rm = T)),
                 qnorm(.975)*sd(STAI[Qd_strGr==0,T],na.rm = T)/sqrt(sum(Qd_strGr==0,na.rm = T)))
  stats <- t.test(MDBF[Qd_strGr==1,T],MDBF[Qd_strGr==0,T])
  statsMDBF[T,] <- c(stats$statistic,stats$parameter,stats$p.value) #significant difference only at T3
  ciMDBF[T,] <- c(qnorm(.975)*sd(MDBF[Qd_strGr==1,T],na.rm = T)/sqrt(sum(Qd_strGr==1,na.rm = T)),
                  qnorm(.975)*sd(MDBF[Qd_strGr==0,T],na.rm = T)/sqrt(sum(Qd_strGr==0,na.rm = T)))
}

#figures
png('Figures/STAI',width = 600, height = 600)
plot(.95:3.95,colMeans(STAI[Qd_strGr==1,],na.rm = TRUE),type = 'b',col='red',
     xlim = c(.5,4.5),ylim = c(1.6,2.4),xaxt="n",xlab = 'Time point',ylab = 'STAI-S (mean)',lwd=2)
lines(1.05:4.05,colMeans(STAI[Qd_strGr==0,],na.rm = TRUE),type = 'b',lwd=2)
for (T in 1:4){
  lines(c(T,T)-.05,mean(STAI[Qd_strGr==1,T],na.rm = TRUE)+ciSTAI[T,1]*c(-1,1),col='red')
  lines(c(T,T)+.05,mean(STAI[Qd_strGr==0,T],na.rm = TRUE)+ciSTAI[T,2]*c(-1,1))
}
text(3,2.35,'***',cex=2);text(3,2.32,'(p < .001)')
text(4,2.35,'*',cex=2);text(4,2.32,'(p = .040)')
axis(1,at = 1:4,labels = c('T1','T2','T3','T4'))
title('STAI-S')
legend(0.5,2.4,c('Stress','No stress'),fill=c('red','black'),bty='n')
dev.off()

png('Figures/MDBF',width = 600, height = 600)
plot(.95:3.95,colMeans(MDBF[Qd_strGr==1,],na.rm = TRUE),type = 'b',col='red',
     xlim = c(.5,4.5),ylim = c(3.3,4.2),xaxt="n",xlab = 'Time point',ylab = 'MDBF (mean)',lwd=2)
lines(1.05:4.05,colMeans(MDBF[Qd_strGr==0,],na.rm = TRUE),type = 'b',lwd=2)
for (T in 1:4){
  lines(c(T,T)-.05,mean(MDBF[Qd_strGr==1,T],na.rm = TRUE)+ciMDBF[T,1]*c(-1,1),col='red')
  lines(c(T,T)+.05,mean(MDBF[Qd_strGr==0,T],na.rm = TRUE)+ciMDBF[T,2]*c(-1,1))
}
text(3,4.16,'***',cex=2);text(3,4.12,'(p < .001)')
axis(1,at = 1:4,labels = c('T1','T2','T3','T4'))
title('MDBF')
legend(3.5,4.2,c('Stress','No stress'),fill=c('red','black'),bty='n')
dev.off()
    
### --- ### --- ### --- ###
### --- Rating consistency before and after stress induction (mentioned in preregistration protocol)
### --- ### --- ### --- ###
corrR12R3 <- c()
for (s in 1:nsubj) {
  corrR12R3[s] <- min(cor(rating[rating$SubjectID==subject_ID[s],8:9]))
}
ratingCorStatsGroupComp <- t.test(FisherZ(corrR12R3[stressGroup==1]),FisherZ(corrR12R3[stressGroup==0])) #n.s.


### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### ---
### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- 
### --- Choice consistency (pre-registered)
### --- ### --- ### --- ###

#add group assignment to datasets (incl. recall) and order by group
decision$StressGroup <- 0
decision$StressGroup[which(decision$SubjectID %in% GaB$VpNr[stressGroup==1])] = 1
decision <- decision[order(decision$StressGroup,decreasing=TRUE),] #put all subjects from stress group first
cued_recall$StressGroup <- 0
cued_recall$StressGroup[which(cued_recall$SubjectID %in% GaB$VpNr[stressGroup==1])] = 1
cued_recall <- cued_recall[order(cued_recall$StressGroup,decreasing=TRUE),] #put all subjects from stress group first
stressGroup <- c(rep(1,sum(stressGroup==1)),rep(0,sum(stressGroup==0))) #needs to be reordered as well
subject_ID <- unique(decision$SubjectID)

#add some important variables to the decision data.frame
decision$vL_vs_vR <- (decision$LeftAvgRatingR1R2>decision$RightAvgRatingR1R2)-(decision$LeftAvgRatingR1R2<decision$RightAvgRatingR1R2)
decision$choiceConsistent <- ((decision$Choice==1)&(decision$vL_vs_vR==1))|((decision$Choice==2)&(decision$vL_vs_vR==-1))
decision$eligTrials <- decision$LeftAvgRatingR1R2!=decision$RightAvgRatingR1R2

#get consistency per subject
consistency <- data.frame(matrix(nrow = nsubj,ncol = 3))
colnames(consistency) <- c('overall','memory','control')
for (s in 1:nsubj){
  consistency$overall[s] <- mean(decision$choiceConsistent[decision$eligTrials
                                                            &(decision$SubjectID==subject_ID[s])])
  consistency$memory[s] <- mean(decision$choiceConsistent[decision$eligTrials
                                                      &(decision$SubjectID==subject_ID[s])&(decision$IsControlBlock==0)])
  consistency$control[s] <- mean(decision$choiceConsistent[decision$eligTrials
                                                       &(decision$SubjectID==subject_ID[s])&(decision$IsControlBlock==1)])
}

#simple frequentist statistics
choiceAllGroupComp <- t.test(consistency$overall[stressGroup==1],consistency$overall[stressGroup==0])
choiceControlGroupComp <- t.test(consistency$control[stressGroup==1],consistency$control[stressGroup==0])
choiceMemGroupComp <- t.test(consistency$memory[stressGroup==1],consistency$memory[stressGroup==0])
#looks like stress was generally more consistent, with effect more driven by control trials -> unexpected

### --- ### Bayesian analysis ### --- ###
# preparations
subjTrialIdx <- c()
stressTrialIdx <- c()
for (t in 1:nrow(decision)){
  subjTrialIdx[t] <- which(subject_ID==decision$SubjectID[t]) #counts trial-wise subjects from 1 to n
  stressTrialIdx[t] <- stressGroup[decision$SubjectID[t]==subject_ID] #counts trial-wise group assignment
}

#data definition
priorSD <- 1 #SD of group-level (half) normal distributions
N1 <- sum(decision$eligTrials[stressTrialIdx==1]) #n of trials, stress group
N2 <- sum(decision$eligTrials[stressTrialIdx==0]) #n of trials, no stress group
S1 <- subjTrialIdx[decision$eligTrials&(stressTrialIdx==1)] #subjects of stress group
S2 <- subjTrialIdx[decision$eligTrials&(stressTrialIdx==0)]-max(S1) #subjects of no stress group (starting from 1)
M1 <- decision$IsControlBlock[decision$eligTrials&(stressTrialIdx==1)]==0 #memory trials, stress group
M2 <- decision$IsControlBlock[decision$eligTrials&(stressTrialIdx==0)]==0 #memory trials, no stress group
C1 <- decision$choiceConsistent[decision$eligTrials&(stressTrialIdx==1)] #consistency, stress group
C2 <- decision$choiceConsistent[decision$eligTrials&(stressTrialIdx==0)] #consistency, no stress group

# get data and initial values together and specify the model
choiceData <- list('N1'=N1,'N2'=N2,'S1'=S1,'S2'=S2,'M1'=M1,'M2'=M2,'C1'=C1,'C2'=C2,
                   'priorSD'=priorSD,'nsubj_Stress'=nsubj_Stress,'nsubj_noStress'=nsubj_noStress)

choiceJAGS <- jags.parallel(choiceData,
                            parameters.to.save = c('mu_Baseline_Stress','sigma_Baseline_Stress','mu_Memory_Stress','sigma_Memory_Stress',
                                                   'mu_Baseline_noStress','sigma_Baseline_noStress','mu_Memory_noStress','sigma_Memory_noStress',
                                                   's_Baseline_Stress','s_Memory_Stress','s_Baseline_noStress','s_Memory_noStress'),
                            model.file = "BayesModel_ChoiceConsistency.txt",working.directory = 'BayesModels',
                            n.chains = 8, n.iter = 100000, n.burnin = 60000,n.thin = 20, n.cluster= 8, DIC = TRUE)

# check convergence
choiceRhats <- choiceJAGS$BUGSoutput$summary[,8] #check with max(Rhats), which should ideally be < 1.01 (1.05 would also be okay)
max(choiceRhats)

# plot performance in all 4 conditions (ChoiceType by Stress)
consistency11 <- pnorm(choiceJAGS$BUGSoutput$sims.list$mu_Baseline_Stress+choiceJAGS$BUGSoutput$sims.list$mu_Memory_Stress)
consistency12 <- pnorm(choiceJAGS$BUGSoutput$sims.list$mu_Baseline_noStress+choiceJAGS$BUGSoutput$sims.list$mu_Memory_noStress)
consistency21 <- pnorm(choiceJAGS$BUGSoutput$sims.list$mu_Baseline_Stress)
consistency22 <- pnorm(choiceJAGS$BUGSoutput$sims.list$mu_Baseline_noStress)

png('Figures/Choice',width = 600, height = 700)
histGran <- (0:200)/200
par(fig=c(0,1,0,1))
hist(consistency21,histGran,xlim=c(.48,.92),ylim=c(0,3500),xaxt='n',yaxt='n',
     col=rgb(1,0,.5,1/2),xlab = 'p(consistent choice)',main = 'Choice consistency')
hist(consistency11,histGran,add=T,col=rgb(1,0,0,1/2))
hist(consistency12,histGran,add=T,col=rgb(0,0,1,1/2))
hist(consistency22,histGran,add=T,col=rgb(0,.5,1,1/2))
axis(1,at = (0:10)/10)
axis(2,at = seq(0,3000,1000))
legend(0.5,3200,legend = c('no Stress/Memory','Stress/Memory','no Stress/Control','Stress/Control'),
       fill = c(rgb(0,0,1,1/2),rgb(1,0,0,1/2),rgb(0,.5,1,1/2),rgb(1,0,.5,1/2)),bty='n')
dev.off()

#plot effects of memory, stress and memory x stress
choiceMemory <- choiceJAGS$BUGSoutput$sims.list$mu_Memory_Stress + choiceJAGS$BUGSoutput$sims.list$mu_Memory_noStress
choiceStressBaseline <- choiceJAGS$BUGSoutput$sims.list$mu_Baseline_Stress - choiceJAGS$BUGSoutput$sims.list$mu_Baseline_noStress
choiceStressMemory <- choiceJAGS$BUGSoutput$sims.list$mu_Memory_Stress - choiceJAGS$BUGSoutput$sims.list$mu_Memory_noStress

png('Figures/ChoiceEffects',width = 600, height = 700)
histGran <- ((0:200)-100)/50
par(fig=c(0,1,.6,1))
h1 <- hist(choiceMemory,histGran,xlim = c(-2,2),main='Memory effect',xlab = '')
lines(c(0,0),c(0,max(h1$counts)),col='red',lwd=2)
HDI <- round(quantile(choiceMemory,c(.025,.975)),2)
text(1.5,max(h1$counts)/2,paste('95% HDI = [',HDI[1],',',HDI[2],']'))
par(fig=c(0,1,.3,.7),new=T)
h1 <- hist(choiceStressBaseline,histGran,xlim = c(-2,2),main='Stress effect on choices (in general)',xlab = '')
lines(c(0,0),c(0,max(h1$counts)),col='red',lwd=2)
HDI <- round(quantile(choiceStressBaseline,c(.025,.975)),2)
text(1.5,max(h1$counts)/2,paste('95% HDI = [',HDI[1],',',HDI[2],']'))
par(fig=c(0,1,0,.4),new=T)
h1 <- hist(choiceStressMemory,histGran,xlim = c(-2,2),main='Stress effect specific to memory-based decisions',xlab = '')
lines(c(0,0),c(0,max(h1$counts)),col='red',lwd=2)
HDI <- round(quantile(choiceStressMemory,c(.025,.975)),2)
text(1.5,max(h1$counts)/2,paste('95% HDI = [',HDI[1],',',HDI[2],']'))
dev.off()

### --- ### Exploratory / Post-hoc analyses ### --- ###

#within stress group: higher T3-T1 change in Blood Pressure associated with higher (baseline) choice consistency: r = .264, p = .024
cs<-match(subject_ID[1:nsubj_Stress],Qd$VpNr)
consistency_and_blood <- cor.test(Sys[cs,3]+Dia[cs,3]-Sys[cs,1]-Dia[cs,1],choiceJAGS$BUGSoutput$mean$s_Baseline_Stress)

#within stress group: higher T3-T1 change in STAI-S trend-associated with higher (baseline) choice consistency: r = .216, p = .067
consistency_and_STAI <- cor.test(STAI[cs,3]-STAI[cs,1],choiceJAGS$BUGSoutput$mean$s_Baseline_Stress)

#preparation for cortisol (matrix)
dm <- matrix(cortDat$CortLog[which(cortDat$Subject %in% subject_ID)],ncol=5,byrow=T)
sm <- matrix(cortDat$Subject[which(cortDat$Subject %in% subject_ID)],ncol=5,byrow=T)
cortMat <- data.frame(cbind(dm,sm[,1]));colnames(cortMat) <- c('T1','T2','T3','T4','T5','sID')
cs <- match(subject_ID[1:nsubj_Stress],cortMat$sID) #cs <- match(subject_ID[(nsubj_Stress+1):(nsubj)],cortMat$sID)
#within stress group: higher T4-T1 change in Cortisol trend-associated with higher (baseline) choice consistency: r = .213, p = .068
consistency_and_Cortisol <- cor.test(cortMat[cs,4]-cortMat[cs,1],choiceJAGS$BUGSoutput$mean$s_Baseline_Stress)


### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### ---
### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- 
### --- Memory performance (pre-registered)
### --- ### --- ### --- ###

#extend cued recall matrix and get average recall performances
cued_recall_answers <- read_excel("CUED_RECALL_ANSWERS.xlsx")
cued_recall$eligTrials <- NA
cued_recall$Correct <- NA
recallPerf <- data.frame(matrix(data=NA,nrow = nsubj,ncol = 3))
colnames(recallPerf) <- c('overall','memory','control')
for (s in 1:nsubj){
  recallCorr <- as.data.frame(cued_recall[cued_recall$SubjectID==subject_ID[s],])
  recallAnsw <- as.data.frame(cued_recall_answers[cued_recall_answers$VpNo==subject_ID[s],])
  cued_recall$eligTrials[cued_recall$SubjectID==subject_ID[s]] = recallAnsw$SnackNo!="NA"
  cued_recall$Correct[cued_recall$SubjectID==subject_ID[s]] = recallCorr$ImageNumber==recallAnsw$SnackNo
  recallPerf$overall[s] <- mean(cued_recall$Correct[(cued_recall$SubjectID==subject_ID[s])&cued_recall$eligTrials])
  recallPerf$memory[s] <- mean(cued_recall$Correct[(cued_recall$SubjectID==subject_ID[s])&(cued_recall$IsControlBlock==0)&cued_recall$eligTrials])
  recallPerf$control[s] <- mean(cued_recall$Correct[(cued_recall$SubjectID==subject_ID[s])&(cued_recall$IsControlBlock==1)&cued_recall$eligTrials])
}

#simple frequentist statistics
recallAllGroupComp <- t.test(recallPerf$overall[stressGroup==1],recallPerf$overall[stressGroup==0])
recallControlGroupComp <- t.test(recallPerf$control[stressGroup==1],recallPerf$control[stressGroup==0])
recallMemGroupComp <- t.test(recallPerf$memory[stressGroup==1],recallPerf$memory[stressGroup==0]) 
#looks like stress recalled better in memory trials, which kind of drives the trend effect overall -> unexpected

### --- ### Bayesian analysis ### --- ###
# preparations
subjTrialIdx <- c()
stressTrialIdx <- c()
for (t in 1:nrow(cued_recall)){
  subjTrialIdx[t] <- which(subject_ID==cued_recall$SubjectID[t]) #counts trial-wise subjects from 1 to n
  stressTrialIdx[t] <- stressGroup[cued_recall$SubjectID[t]==subject_ID] #counts trial-wise group assignment
}

# data definition
#data definition
N1 <- sum(cued_recall$eligTrials[stressTrialIdx==1]) #n of trials, stress group
N2 <- sum(cued_recall$eligTrials[stressTrialIdx==0]) #n of trials, no stress group
S1 <- subjTrialIdx[cued_recall$eligTrials&(stressTrialIdx==1)] #subjects of stress group
S2 <- subjTrialIdx[cued_recall$eligTrials&(stressTrialIdx==0)]-max(S1) #subjects of no stress group (starting from 1)
M1 <- cued_recall$IsControlBlock[cued_recall$eligTrials&(stressTrialIdx==1)]==0 #memory trials, stress group
M2 <- cued_recall$IsControlBlock[cued_recall$eligTrials&(stressTrialIdx==0)]==0 #memory trials, no stress group
C1 <- cued_recall$Correct[cued_recall$eligTrials&(stressTrialIdx==1)] #accuracy, stress group
C2 <- cued_recall$Correct[cued_recall$eligTrials&(stressTrialIdx==0)] #accuracy, no stress group

# get data and initial values together and specify the model
recallData <- list('N1'=N1,'N2'=N2,'S1'=S1,'S2'=S2,'M1'=M1,'M2'=M2,'C1'=C1,'C2'=C2,
                   'priorSD'=priorSD,'nsubj_Stress'=nsubj_Stress,'nsubj_noStress'=nsubj_noStress)

recallJAGS <- jags.parallel(recallData,
                            parameters.to.save = c('mu_Baseline_Stress','sigma_Baseline_Stress','mu_Memory_Stress','sigma_Memory_Stress',
                                                   'mu_Baseline_noStress','sigma_Baseline_noStress','mu_Memory_noStress','sigma_Memory_noStress',
                                                   's_Baseline_Stress','s_Memory_Stress','s_Baseline_noStress','s_Memory_noStress'),
                            model.file = "BayesModel_RecallPerformance.txt",working.directory = 'BayesModels',
                            n.chains = 8, n.iter = 100000, n.burnin = 60000,n.thin = 20, n.cluster= 8, DIC = TRUE)

# check convergence
recallRhats <- recallJAGS$BUGSoutput$summary[,8]
max(recallRhats)

# plot performance in all 4 conditions (ChoiceType by Stress)
accuracy11 <- pnorm(recallJAGS$BUGSoutput$sims.list$mu_Baseline_Stress+recallJAGS$BUGSoutput$sims.list$mu_Memory_Stress)
accuracy12 <- pnorm(recallJAGS$BUGSoutput$sims.list$mu_Baseline_noStress+recallJAGS$BUGSoutput$sims.list$mu_Memory_noStress)
accuracy21 <- pnorm(recallJAGS$BUGSoutput$sims.list$mu_Baseline_Stress)
accuracy22 <- pnorm(recallJAGS$BUGSoutput$sims.list$mu_Baseline_noStress)

png('Figures/Recall',width = 600, height = 700)
histGran <- (0:75)/75
par(fig=c(0,1,0,1))
hist(accuracy21,histGran,xlim=c(.2,.8),ylim=c(0,2500),xaxt='n',yaxt='n',
     col=rgb(1,0,.5,1/2),xlab = 'p(recall correct)',main = 'Recall performance')
hist(accuracy11,histGran,add=T,col=rgb(1,0,0,1/2))
hist(accuracy12,histGran,add=T,col=rgb(0,0,1,1/2))
hist(accuracy22,histGran,add=T,col=rgb(0,.5,1,1/2))
axis(1,at = (0:10)/10)
axis(2,at = seq(0,3000,1000))
legend(0.2,2300,legend = c('no Stress/Memory','Stress/Memory','no Stress/Control','Stress/Control'),
       fill = c(rgb(0,0,1,1/2),rgb(1,0,0,1/2),rgb(0,.5,1,1/2),rgb(1,0,.5,1/2)),bty='n')
dev.off()

#plot effects of memory, stress and memory x stress
recallMemory <- recallJAGS$BUGSoutput$sims.list$mu_Memory_Stress + recallJAGS$BUGSoutput$sims.list$mu_Memory_noStress
recallStress <- recallJAGS$BUGSoutput$sims.list$mu_Baseline_Stress - recallJAGS$BUGSoutput$sims.list$mu_Baseline_noStress
recallStressMemory <- recallJAGS$BUGSoutput$sims.list$mu_Memory_Stress - recallJAGS$BUGSoutput$sims.list$mu_Memory_noStress

png('Figures/RecallEffects',width = 600, height = 700)
histGran <- ((0:100)-50)/50
par(fig=c(0,1,.6,1))
h1 <- hist(recallMemory,histGran,xlim = c(-1,1),main='Memory effect',xlab = '')
lines(c(0,0),c(0,max(h1$counts)),col='red',lwd=2)
HDI <- round(quantile(recallMemory,c(.025,.975)),2)
text(0.75,max(h1$counts)/2,paste('95% HDI = [',HDI[1],',',HDI[2],']'))
par(fig=c(0,1,.3,.7),new=T)
h1 <- hist(recallStress,histGran,xlim = c(-1,1),main='Stress effect on recall (in general)',xlab = '')
lines(c(0,0),c(0,max(h1$counts)),col='red',lwd=2)
HDI <- round(quantile(recallStress,c(.025,.975)),2)
text(0.75,max(h1$counts)/2,paste('95% HDI = [',HDI[1],',',HDI[2],']'))
par(fig=c(0,1,0,.4),new=T)
h1 <- hist(recallStressMemory,histGran,xlim = c(-1,1),main='Stress effect specific to recall in memory-based trials',xlab = '')
lines(c(0,0),c(0,max(h1$counts)),col='red',lwd=2)
HDI <- round(quantile(recallStressMemory,c(.025,.975)),2)
text(0.75,max(h1$counts)/2,paste('95% HDI = [',HDI[1],',',HDI[2],']'))
dev.off()

### --- ### Exploratory / Post-hoc analyses ### --- ###

#within no-stress group: lower T3-T1 change in STAI-S associated with higher (memory-specific) memory accuracy r = -.254, p = .028
cs <- match(subject_ID[(nsubj_Stress+1):nsubj],Qd$VpNr)
accuracyMem_and_STAI <- cor.test(STAI[cs,3]-STAI[cs,1],recallJAGS$BUGSoutput$mean$s_Memory_noStress)
#within no-stress group: lower T3-T1 change in MDBF associated with higher (baseline) memory accuracy r = .222, p = .055
accuracy_and_MDBF <- cor.test(MDBF[cs,3]-MDBF[cs,1],recallJAGS$BUGSoutput$mean$s_Baseline_noStress)


### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### ---
### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- 
### --- Memory bias (mentioned in the exploratory part of preregistration)
### --- ### --- ### --- ###

decision$LeftMem <- NA #whether left item was (later) recalled
decision$RightMem <- NA #whether right item was (later) recalled
decision$MemLevel <- NA #0 = none remembered, 1 = one remembered, 2 = both remembered
decision$zValLeft <- NA #standardized value of left option
decision$zValRight <- NA #standardized value of right option
decision$zValMem <- NA #standardized value of remembered option (if MemLevel == 1)
coeffMemBias <- data.frame(matrix(nrow = nsubj,ncol = 2))
colnames(coeffMemBias) <- c('MemoryBias','ValueRemembered')
for (s in 1:nsubj){
  #get decision and recall data of current subject
  cdata <- decision[decision$SubjectID==subject_ID[s],]
  rdata <- cued_recall[cued_recall$SubjectID==subject_ID[s],]
  
  #go into each (memory) block
  LeftMem <- c()
  RightMem <- c()
  for (b in 1:3){
    cb <- cdata[cdata$BlockNumber==b,]
    rb <- rdata[rdata$BlockNumber==b,]
    LeftMem <- c(LeftMem,rb$Correct[(match(cb$LeftImageNumber,rb$ImageNumber))])
    RightMem <- c(RightMem,rb$Correct[(match(cb$RightImageNumber,rb$ImageNumber))])
  }
  
  #save memory info
  decision$LeftMem[decision$SubjectID==subject_ID[s]] <- LeftMem
  decision$RightMem[decision$SubjectID==subject_ID[s]] <- RightMem
  decision$MemLevel[decision$SubjectID==subject_ID[s]] <- LeftMem+RightMem
  
  #save value info
  M_SD <- c(mean(c(cdata$LeftAvgRatingR1R2,cdata$RightAvgRatingR1R2)),sd(c(cdata$LeftAvgRatingR1R2,cdata$RightAvgRatingR1R2)))
  decision$zValLeft[decision$SubjectID==subject_ID[s]] <- (cdata$LeftAvgRatingR1R2-M_SD[1])/M_SD[2]
  decision$zValRight[decision$SubjectID==subject_ID[s]] <- (cdata$RightAvgRatingR1R2-M_SD[1])/M_SD[2]
  zValMem <- decision$zValLeft[decision$SubjectID==subject_ID[s]]*LeftMem*(LeftMem!=RightMem)+
             decision$zValRight[decision$SubjectID==subject_ID[s]]*RightMem*(LeftMem!=RightMem)
  decision$zValMem[decision$SubjectID==subject_ID[s]] <- zValMem
  
  #perform subject-wise logistic regression for memory bias
  critTrials <- (LeftMem!=RightMem)&(cdata$IsControlBlock==0)&((cdata$Choice==1)|(cdata$Choice==2))
  X <- zValMem[critTrials]
  Y <- ((cdata$Choice[critTrials]==1)&(LeftMem[critTrials]==1))|((cdata$Choice[critTrials]==2)&(RightMem[critTrials]==1))
  if ((sum(Y==1)>=5) && (sum(Y==0)>=5)){ #according to preregistration, only run analysis if at least 5 decisions each were made
    logReg <- glm(Y ~ X, family = binomial())
    coeffMemBias[s,] <- logReg$coefficients
  }
}

#mean choice consistency per subject, memory level and condition
consistentMemLevel <- tapply(decision$choiceConsistent,
                             list(decision$SubjectID,decision$MemLevel,decision$IsControlBlock),mean)

#simple frequentist statistics
MemBiasIntercept <- t.test(coeffMemBias$MemoryBias) #this is the actual memory bias (p = .007; but effect size rather low: d = 0.25)
MemBiasSlope <- t.test(coeffMemBias$ValueRemembered) #this is (just checking) the influence of the remembered option's value (p < .001)
MemBiasInterceptGroupComp <- t.test(coeffMemBias$MemoryBias[stressGroup==1],coeffMemBias$MemoryBias[stressGroup==0]) # n.s.
MemBiasSlopeGroupComp <- t.test(coeffMemBias$ValueRemembered[stressGroup==1],coeffMemBias$ValueRemembered[stressGroup==0]) # n.s.
#memory bias successfully replicated but no influence of stress

### --- ### Bayesian analysis ### --- ###
# preparations
decisionMB <- decision[which(decision$SubjectID %in% subject_ID[is.na(coeffMemBias[,1])==F]),] #dataset only with usable subjects
subject_ID_MB <- unique(decisionMB$SubjectID)
nsubjMB_Stress <- length(unique(decisionMB$SubjectID[decisionMB$StressGroup==1]))
nsubjMB_noStress <- length(unique(decisionMB$SubjectID[decisionMB$StressGroup==0]))
subjTrialIdx <- c()
stressTrialIdx <- c()
for (t in 1:nrow(decisionMB)){
  subjTrialIdx[t] <- which(subject_ID_MB==decisionMB$SubjectID[t]) #counts trial-wise subjects from 1 to n
  stressTrialIdx[t] <- stressGroup[decisionMB$SubjectID[t]==subject_ID] #counts trial-wise group assignment
}
critTrials <- (decisionMB$LeftMem!=decisionMB$RightMem)&(decisionMB$IsControlBlock==0)&((decisionMB$Choice==1)|(decisionMB$Choice==2))

# data definition
decisionMB$MemChosen <- ((decisionMB$Choice==1)&(decisionMB$LeftMem==1))|((decisionMB$Choice==2)&(decisionMB$RightMem==1))
N1 <- sum(critTrials[(stressTrialIdx==1)]) #n of critical trials, stress group
N2 <- sum(critTrials[(stressTrialIdx==0)]) #n of critical trials, no stress group
S1 <- subjTrialIdx[critTrials&(stressTrialIdx==1)] #subjects in critical trials, stress group
S2 <- subjTrialIdx[critTrials&(stressTrialIdx==0)]-max(S1) #subjects in critical trials, no stress group (starting from 1)
V1 <- decisionMB$zValMem[critTrials&(stressTrialIdx==1)] #value of remembered option, stress group
V2 <- decisionMB$zValMem[critTrials&(stressTrialIdx==0)] #value of remembered option, no stress group
C1 <- decisionMB$MemChosen[critTrials&(stressTrialIdx==1)] #remembered option chosen, stress group
C2 <- decisionMB$MemChosen[critTrials&(stressTrialIdx==0)] #remembered option chosen, no stress group

membiasData <- list('N1'=N1,'N2'=N2,'S1'=S1,'S2'=S2,'V1'=V1,'V2'=V2,'C1'=C1,'C2'=C2,
                    'priorSD'=priorSD,'nsubj_Stress'=nsubjMB_Stress,'nsubj_noStress'=nsubjMB_noStress)

membiasJAGS <- jags.parallel(membiasData,
                             parameters.to.save = c('mu_Bias_Stress','sigma_Bias_Stress','mu_Value_Stress','sigma_Value_Stress',
                                                    'mu_Bias_noStress','sigma_Bias_noStress','mu_Value_noStress','sigma_Value_noStress',
                                                    's_Bias_Stress','s_Value_Stress','s_Bias_noStress','s_Value_noStress'),
                             model.file = "BayesModel_MemoryBias.txt",working.directory = 'BayesModels',
                             n.chains = 8, n.iter = 100000, n.burnin = 60000,n.thin = 20, n.cluster= 8, DIC = TRUE)

# THIS WOULD BE AN ALTERNATIVE WAY OF MODELING IT, BUT THE RESULTS ARE PRETTY MUCH THE SAME:
# membiasJAGS_v2 <- jags.parallel(membiasData,
#                              parameters.to.save = c('mu_Bias','sigma_Bias','mu_Value','sigma_Value',
#                                'mu_Bias_Stress','sigma_Bias_Stress','mu_Value_Stress','sigma_Value_Stress',
#                                's_Bias_SG','s_Value_SG','s_Bias_Stress_SG','s_Value_Stress_SG',
#                                's_Bias_NG','s_Value_NG','s_Bias_Stress_NG','s_Value_Stress_NG'),
#                              model.file = "BayesModel_MemoryBias_v2.txt",working.directory = 'BayesModels',
#                              n.chains = 8, n.iter = 100000, n.burnin = 50000,n.thin = 25,n.cluster= 8, DIC = TRUE)

# check convergence
membiasRhats <- membiasJAGS$BUGSoutput$summary[,8] #check with max(Rhats), which should ideally be < 1.01 (1.05 would also be okay)
max(membiasRhats)

#plot effects of bias, value and stress
membiasBias <- membiasJAGS$BUGSoutput$sims.list$mu_Bias_Stress+membiasJAGS$BUGSoutput$sims.list$mu_Bias_noStress
membiasValue <- membiasJAGS$BUGSoutput$sims.list$mu_Value_Stress+membiasJAGS$BUGSoutput$sims.list$mu_Value_noStress
membiasStressOnBias <- membiasJAGS$BUGSoutput$sims.list$mu_Bias_Stress-membiasJAGS$BUGSoutput$sims.list$mu_Bias_noStress
membiasStressOnValue <- membiasJAGS$BUGSoutput$sims.list$mu_Value_Stress-membiasJAGS$BUGSoutput$sims.list$mu_Value_noStress

png('Figures/MemBiasEffectsBias',width = 600, height = 700)
histGran <- ((0:200)-100)/100
par(fig=c(0,1,.45,1))
h1 <- hist(membiasBias,histGran,xlim = c(-.7,.7),main='Memory bias',xlab = '')
lines(c(0,0),c(0,max(h1$counts)),col='red',lwd=2)
HDI <- round(quantile(membiasBias,c(.025,.975)),2)
text(.5,max(h1$counts)/2,paste('95% HDI = [',HDI[1],',',HDI[2],']'))
par(fig=c(0,1,0,.55),new=T)
h1 <- hist(membiasStressOnBias,histGran,xlim = c(-.7,.7),main='Stress effect on memory bias',xlab = '')
lines(c(0,0),c(0,max(h1$counts)),col='red',lwd=2)
HDI <- round(quantile(membiasStressOnBias,c(.025,.975)),2)
text(.5,max(h1$counts)/2,paste('95% HDI = [',HDI[1],',',HDI[2],']'))
dev.off()

png('Figures/MemBiasEffectsValue',width = 600, height = 700)
par(fig=c(0,1,.45,1))
h1 <- hist(membiasValue,histGran,xlim = c(-.7,.7),main='Value of remembered option',xlab = '')
lines(c(0,0),c(0,max(h1$counts)),col='red',lwd=2)
HDI <- round(quantile(membiasValue,c(.025,.975)),2)
text(-.5,max(h1$counts)/2,paste('95% HDI = [',HDI[1],',',HDI[2],']'))
par(fig=c(0,1,0,.55),new=T)
h1 <- hist(membiasStressOnValue,histGran,xlim = c(-.7,.7),main='Stress effect on value',xlab = '')
lines(c(0,0),c(0,max(h1$counts)),col='red',lwd=2)
HDI <- round(quantile(membiasStressOnValue,c(.025,.975)),2)
text(-.5,max(h1$counts)/2,paste('95% HDI = [',HDI[1],',',HDI[2],']'))
dev.off()

#get posterior predictives for the memory bias
nPP <- 1000
zValRange <- (-20:20)/10
pChoiceRem_Stress <- data.frame(matrix(NA,nrow = nPP,ncol = length(zValRange)))
pChoiceRem_noStress <- pChoiceRem_Stress
nSamples <- length(membiasBias)
for (p in 1:nPP){
  cs <- sample(nSamples,1) #current sample
  params <- matrix(c(membiasJAGS$BUGSoutput$sims.list$mu_Bias_Stress[cs],
              membiasJAGS$BUGSoutput$sims.list$mu_Value_Stress[cs],
              membiasJAGS$BUGSoutput$sims.list$mu_Bias_noStress[cs],
              membiasJAGS$BUGSoutput$sims.list$mu_Value_noStress[cs]),ncol = 4)
  pChoiceRem_Stress[p,] <- pnorm(params[1]+params[2]*zValRange)
  pChoiceRem_noStress[p,] <- pnorm(params[3]+params[4]*zValRange)
}
qChoiceRem_Stress <- data.frame(matrix(NA,nrow = 2,ncol = length(zValRange)))
qChoiceRem_noStress <- qChoiceRem_Stress
for (z in 1:length(zValRange)){
  qChoiceRem_Stress[,z] <- quantile(pChoiceRem_Stress[,z],c(.025,.975))
  qChoiceRem_noStress[,z] <- quantile(pChoiceRem_noStress[,z],c(.025,.975))
}

#plot posterior predictives
png('Figures/MemBiasPP',width = 600, height = 600)
par(fig=c(0,1,0,1))
plot(zValRange,rep(.5,length(zValRange)),type='l',lty=2,ylim=c(0,1),bty='l',
     xlab='Value of remembered option',ylab='P(choose remembered)')
title('Memory bias: posterior predictives')
lines(c(0,0),c(0,1),lty=2)
lines(zValRange,qChoiceRem_noStress[1,],col='blue')
lines(zValRange,qChoiceRem_noStress[2,],col='blue')
lines(zValRange,colMeans(qChoiceRem_noStress),lwd=3,col='blue')
lines(zValRange,qChoiceRem_Stress[1,],col='red')
lines(zValRange,qChoiceRem_Stress[2,],col='red')
lines(zValRange,colMeans(qChoiceRem_Stress),lwd=3,col='red')
legend(-2,1,legend = c('Stress','no Stress'),fill = c('red','blue'),bty='n')
dev.off()


### --- ### Exploratory / Post-hoc analyses ### --- ###

#within no-stress group: higher T4-T1 change in cortisol trend-associated with lower memory bias  r = .228, p = .082
cs <- match(subject_ID_MB[(nsubjMB_Stress+1):(nsubjMB_Stress+nsubjMB_noStress)],Qd$VpNr)
memBias_and_Cort2 <- cor.test(cortMat[cs,4]-cortMat[cs,1],membiasJAGS$BUGSoutput$mean$s_Bias_noStress)

### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### ---
### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- 
### --- Response times of decisions (mentioned in the pre-reg protocol)
### --- ### --- ### --- ###

#data cleaning group level (too fast / too slow RTs)
decision$eligTrials[decision$RespTime<.25] = FALSE #basically consistent with exclusion criterion
decision$eligTrials[decision$RespTime>10] = FALSE #these are trials without button presses

#data cleaning individual level
decision$rtIndivOut <- NA
for (s in 1:nsubj){
  eligRT <- decision$RespTime[decision$eligTrials&(decision$SubjectID==subject_ID[s])]
  rtCutoffs <- mean(eligRT)+c(1,-1)*4*sd(eligRT)
  decision$rtIndivOut[decision$SubjectID==subject_ID[s]] = (decision$RespTime[decision$SubjectID==subject_ID[s]]>rtCutoffs[1]) | 
                                                           (decision$RespTime[decision$SubjectID==subject_ID[s]]<rtCutoffs[2])
}
decision$eligTrials[decision$rtIndivOut==TRUE] = FALSE #this seems to exclude trials at block start when participant take very long

#get (log) RT per subject
respTime <- data.frame(matrix(nrow = nsubj,ncol = 3))
colnames(respTime) <- c('overall','memory','control')
for (s in 1:nsubj){
  respTime$overall[s] <- mean(log(decision$RespTime[decision$eligTrials&(decision$SubjectID==subject_ID[s])]))
  respTime$memory[s] <- mean(log(decision$RespTime[decision$eligTrials&(decision$SubjectID==subject_ID[s])&(decision$IsControlBlock==0)]))
  respTime$control[s] <- mean(log(decision$RespTime[decision$eligTrials&(decision$SubjectID==subject_ID[s])&(decision$IsControlBlock==1)]))
}

#simple frequentist statistics
rtAllGroupComp <- t.test(respTime$overall[stressGroup==1],respTime$overall[stressGroup==0]) # n.s.
rtControlGroupComp <- t.test(respTime$control[stressGroup==1],respTime$control[stressGroup==0]) # n.s.
rtMemGroupComp <- t.test(respTime$memory[stressGroup==1],respTime$memory[stressGroup==0]) # n.s.
rtMemControl <- t.test(respTime$memory,respTime$control,paired = T) #p < .001
#looks like stress group has a bit slower RT; together with higher consistency could suggest a higher boundary separation (DDM)

### --- ### Bayesian analysis ### --- ###
# preparations
subjTrialIdx <- c()
stressTrialIdx <- c()
for (t in 1:nrow(decision)){
  subjTrialIdx[t] <- which(subject_ID==decision$SubjectID[t]) #counts trial-wise subjects from 1 to n
  stressTrialIdx[t] <- stressGroup[decision$SubjectID[t]==subject_ID] #counts trial-wise group assignment
}

#data definition
N1 <- sum(decision$eligTrials[stressTrialIdx==1]) #n of trials, stress group
N2 <- sum(decision$eligTrials[stressTrialIdx==0]) #n of trials, no stress group
S1 <- subjTrialIdx[decision$eligTrials&(stressTrialIdx==1)] #subjects of stress group
S2 <- subjTrialIdx[decision$eligTrials&(stressTrialIdx==0)]-max(S1) #subjects of no stress group (starting from 1)
M1 <- decision$IsControlBlock[decision$eligTrials&(stressTrialIdx==1)]==0 #memory trials, stress group
M2 <- decision$IsControlBlock[decision$eligTrials&(stressTrialIdx==0)]==0 #memory trials, no stress group
C1 <- log(decision$RespTime[decision$eligTrials&(stressTrialIdx==1)]) #RT, stress group
C2 <- log(decision$RespTime[decision$eligTrials&(stressTrialIdx==0)]) #RT, no stress group
nsubj_Stress <- length(unique(S1))
nsubj_noStress <- length(unique(S2))

# get data and initial values together and specify the model
rtData <- list('N1'=N1,'N2'=N2,'S1'=S1,'S2'=S2,'M1'=M1,'M2'=M2,'C1'=C1,'C2'=C2,
                   'priorSD'=priorSD,'nsubj_Stress'=nsubj_Stress,'nsubj_noStress'=nsubj_noStress)

rtJAGS <- jags.parallel(rtData,
                        parameters.to.save = c('mu_Baseline_Stress','sigma_Baseline_Stress','mu_Memory_Stress','sigma_Memory_Stress',
                                               'mu_Baseline_noStress','sigma_Baseline_noStress','mu_Memory_noStress','sigma_Memory_noStress',
                                               's_Baseline_Stress','s_Memory_Stress','s_Baseline_noStress','s_Memory_noStress',
                                               'g1_Baseline_Stress','g2_Baseline_Stress','g1_Memory_Stress','g2_Memory_Stress',
                                               'g1_Baseline_noStress','g2_Baseline_noStress','g1_Memory_noStress','g2_Memory_noStress',
                                               't_Baseline_Stress','t_Memory_Stress','t_Baseline_noStress','t_Memory_noStress'),
                        model.file = "BayesModel_RT.txt",working.directory = 'BayesModels',
                        n.chains = 8, n.iter = 100000, n.burnin = 60000,n.thin = 20, n.cluster= 8, DIC = TRUE)

# check convergence
rtRhats <- rtJAGS$BUGSoutput$summary[,8] #check with max(Rhats), which should ideally be < 1.01 (1.05 would also be okay)
max(rtRhats)

# plot performance in all 4 conditions (ChoiceType by Stress)
rt11 <- exp(rtJAGS$BUGSoutput$sims.list$mu_Baseline_Stress+rtJAGS$BUGSoutput$sims.list$mu_Memory_Stress)
rt12 <- exp(rtJAGS$BUGSoutput$sims.list$mu_Baseline_noStress+rtJAGS$BUGSoutput$sims.list$mu_Memory_noStress)
rt21 <- exp(rtJAGS$BUGSoutput$sims.list$mu_Baseline_Stress)
rt22 <- exp(rtJAGS$BUGSoutput$sims.list$mu_Baseline_noStress)

png('Figures/RTs',width = 600, height = 700)
histGran <- (0:300)/100
par(fig=c(0,1,0,1))
hist(rt21,histGran,xlim=c(1,2.3),ylim=c(0,2000),xaxt='n',yaxt='n',
     col=rgb(1,0,.5,1/2),xlab = 'Mean RT in s',main = 'Response times')
hist(rt11,histGran,add=T,col=rgb(1,0,0,1/2))
hist(rt12,histGran,add=T,col=rgb(0,0,1,1/2))
hist(rt22,histGran,add=T,col=rgb(0,.5,1,1/2))
axis(1,at = (0:6)/2)
axis(2,at = seq(0,3000,1000))
legend(1.75,2000,legend = c('no Stress/Memory','Stress/Memory','no Stress/Control','Stress/Control'),
       fill = c(rgb(0,0,1,1/2),rgb(1,0,0,1/2),rgb(0,.5,1,1/2),rgb(1,0,.5,1/2)),bty='n')
dev.off()

#plot effects of memory, stress and memory x stress
rtMemory <- rtJAGS$BUGSoutput$sims.list$mu_Memory_Stress + rtJAGS$BUGSoutput$sims.list$mu_Memory_noStress
rtStressBaseline <- rtJAGS$BUGSoutput$sims.list$mu_Baseline_Stress - rtJAGS$BUGSoutput$sims.list$mu_Baseline_noStress
rtStressMemory <- rtJAGS$BUGSoutput$sims.list$mu_Memory_Stress - rtJAGS$BUGSoutput$sims.list$mu_Memory_noStress

png('Figures/RTsEffects',width = 600, height = 700)
histGran <- ((0:200)-100)/50
par(fig=c(0,1,.6,1))
h1 <- hist(rtMemory,histGran,xlim = c(-1,1),main='Memory effect',xlab = '')
lines(c(0,0),c(0,max(h1$counts)),col='red',lwd=2)
HDI <- round(quantile(rtMemory,c(.025,.975)),2)
text(-.75,max(h1$counts)/2,paste('95% HDI = [',HDI[1],',',HDI[2],']'))
par(fig=c(0,1,.3,.7),new=T)
h1 <- hist(rtStressBaseline,histGran,xlim = c(-1,1),main='Stress effect on RT (in general)',xlab = '')
lines(c(0,0),c(0,max(h1$counts)),col='red',lwd=2)
HDI <- round(quantile(rtStressBaseline,c(.025,.975)),2)
text(-.75,max(h1$counts)/2,paste('95% HDI = [',HDI[1],',',HDI[2],']'))
par(fig=c(0,1,0,.4),new=T)
h1 <- hist(rtStressMemory,histGran,xlim = c(-1,1),main='Stress effect specific to memory-based decisions',xlab = '')
lines(c(0,0),c(0,max(h1$counts)),col='red',lwd=2)
HDI <- round(quantile(rtStressMemory,c(.025,.975)),2)
text(-.75,max(h1$counts)/2,paste('95% HDI = [',HDI[1],',',HDI[2],']'))
dev.off()

### --- ### Exploratory / Post-hoc analyses ### --- ###

#within no-stress group: higher T3-T1 change in blood pressure associated with higher RT: r = .221, p = .056
cs <- match(subject_ID[(nsubj_Stress+1):nsubj],Qd$VpNr)
rt_and_Blood <- cor.test(Sys[cs,3]+Dia[cs,3]-Sys[cs,1]-Dia[cs,1],rtJAGS$BUGSoutput$mean$s_Baseline_noStress)
#within no-stress group: higher T3-T1 change in STAI associated with lower RT: r = -.281, p = .015
rt_and_STAI_base <- cor.test(STAI[cs,3]-STAI[cs,1],rtJAGS$BUGSoutput$mean$s_Baseline_noStress)
#within no-stress group: higher T3-T1 change in STAI associated with memory-specific RT reduction: r = -.357, p = .002
rt_and_STAI_mem <- cor.test(STAI[cs,3]-STAI[cs,1],rtJAGS$BUGSoutput$mean$s_Memory_noStress)

#within stress group: higher T4-T1 change in Cortisol associated with memory-specific RT INCREASE: r = .246, p = .036
cs <- match(subject_ID[1:nsubj_Stress],cortMat$sID)
rt_and_Cortisol1 <- cor.test(cortMat[cs,4]-cortMat[cs,1],rtJAGS$BUGSoutput$mean$s_Memory_Stress)
#within no-stress group: higher T4-T1 change in Cortisol associated with memory-specific RT REDUCTION: r = .220, p = .057
cs <- match(subject_ID[(nsubj_Stress+1):(nsubj_Stress+nsubj_noStress)],cortMat$sID)
rt_and_Cortisol2 <- cor.test(cortMat[cs,4]-cortMat[cs,1],rtJAGS$BUGSoutput$mean$s_Memory_noStress)
#within no-stress group: higher T4-T1 change in Cortisol trend-associated with lower RT: r = -.213, p = .067
rt_and_Cortisol3 <- cor.test(cortMat[cs,4]-cortMat[cs,1],rtJAGS$BUGSoutput$mean$s_Baseline_noStress)


### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### ---
### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### ---
### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- 
### --- DDM
### --- ### --- ### --- ###
### --- ### --- ### --- ###

#data definition
decision$zValDiff <- decision$zValLeft-decision$zValRight #compute value difference of left vs. right
VD1 <- decision$zValDiff[decision$eligTrials&(stressTrialIdx==1)]
VD2 <- decision$zValDiff[decision$eligTrials&(stressTrialIdx==0)]
RT1 <- decision$RespTime[decision$eligTrials&(stressTrialIdx==1)]
RT2 <- decision$RespTime[decision$eligTrials&(stressTrialIdx==0)]
RT1[decision$Choice[decision$eligTrials&(stressTrialIdx==1)]==2] = -RT1[decision$Choice[decision$eligTrials&(stressTrialIdx==1)]==2] #for dwiener function in JAGS
RT2[decision$Choice[decision$eligTrials&(stressTrialIdx==0)]==2] = -RT2[decision$Choice[decision$eligTrials&(stressTrialIdx==0)]==2] #for dwiener function in JAGS

# get data and initial values together and specify the model
ddmData <- list('N1'=N1,'N2'=N2,'S1'=S1,'S2'=S2,'M1'=M1,'M2'=M2,'VD1'=VD1,'VD2'=VD2,'RT1'=RT1,'RT2'=RT2,
               'nsubj_Stress'=nsubj_Stress,'nsubj_noStress'=nsubj_noStress,'priorSD'=priorSD)

#this time, initial values are required
nChains <- 8
ddmInits <- vector("list", nChains)
for (c in 1:nChains){
  ddmInits[[c]] <- list(
    'mu_alphaC_Stress'=rnorm(1,3),'sigma_alphaC_Stress'=abs(rnorm(1)),'mu_alphaM_Stress'=rnorm(1),'sigma_alphaM_Stress'=abs(rnorm(1)),
    'mu_tauC_Stress'=rnorm(1,-1),'sigma_tauC_Stress'=abs(rnorm(1)),'mu_tauM_Stress'=rnorm(1),'sigma_tauM_Stress'=abs(rnorm(1)),
    'mu_deltaC_Stress'=rnorm(1),'sigma_deltaC_Stress'=abs(rnorm(1)),'mu_deltaM_Stress'=rnorm(1),'sigma_deltaM_Stress'=abs(rnorm(1)),
    'mu_alphaC_noStress'=rnorm(1,3),'sigma_alphaC_noStress'=abs(rnorm(1)),'mu_alphaM_noStress'=rnorm(1),'sigma_alphaM_noStress'=abs(rnorm(1)),
    'mu_tauC_noStress'=rnorm(1,-1),'sigma_tauC_noStress'=abs(rnorm(1)),'mu_tauM_noStress'=rnorm(1),'sigma_tauM_noStress'=abs(rnorm(1)),
    'mu_deltaC_noStress'=rnorm(1),'sigma_deltaC_noStress'=abs(rnorm(1)),'mu_deltaM_noStress'=rnorm(1),'sigma_deltaM_noStress'=abs(rnorm(1)),
    'alphaC_Stress'=rnorm(nsubj_Stress,3,1),'tauC_Stress'=runif(nsubj_Stress,-10,-2),'deltaC_Stress'=rnorm(nsubj_Stress),
    'alphaM_Stress'=rnorm(nsubj_Stress),'tauM_Stress'=rnorm(nsubj_Stress,0,1/6),'deltaM_Stress'=rnorm(nsubj_Stress),
    'alphaC_noStress'=rnorm(nsubj_noStress,3,1),'tauC_noStress'=runif(nsubj_noStress,-10,-2),'deltaC_noStress'=rnorm(nsubj_noStress),    
    'alphaM_noStress'=rnorm(nsubj_noStress),'tauM_noStress'=rnorm(nsubj_noStress,0,1/6),'deltaM_noStress'=rnorm(nsubj_noStress))
}

#parallel computing apparently not possible
T1<-Sys.time()
ddmJAGS <- jags(ddmData,inits = ddmInits,
#ddmJAGS <- jags.parallel(ddmData,inits = NULL,jags.seed = sample(100:999,1),
                    parameters.to.save = c('mu_alphaC_Stress','sigma_alphaC_Stress','mu_alphaM_Stress','sigma_alphaM_Stress',
                                           'mu_tauC_Stress','sigma_tauC_Stress ','mu_tauM_Stress','sigma_tauM_Stress',
                                           'mu_deltaC_Stress','sigma_deltaC_Stress','mu_deltaM_Stress','sigma_deltaM_Stress',
                                           'mu_alphaC_noStress','sigma_alphaC_noStress','mu_alphaM_noStress','sigma_alphaM_noStress',
                                           'mu_tauC_noStress','sigma_tauC_noStress ','mu_tauM_noStress','sigma_tauM_noStress',
                                           'mu_deltaC_noStress','sigma_deltaC_noStress','mu_deltaM_noStress','sigma_deltaM_noStress',
                                           'alphaC_Stress','tauC_Stress','deltaC_Stress','alphaM_Stress','tauM_Stress','deltaM_Stress',
                                           'alphaC_noStress','tauC_noStress','deltaC_noStress','alphaM_noStress','tauM_noStress','deltaM_noStress'),
                    model.file = "BayesModel_DDM.txt",working.directory = 'BayesModels',
                    n.chains = nChains, n.iter = 40000, n.burnin = 20000,n.thin = 10, DIC = TRUE,jags.module = c("glm","dic","wiener"))
T2<-Sys.time()
#T2-T1

# check convergence
ddmRhats <- ddmJAGS$BUGSoutput$summary[,8] #check with max(Rhats), which should ideally be < 1.01 (1.05 would also be okay)
max(ddmRhats)

#get group parameter estimates in transformed form (C = control, M = memory, S = Stress, N = no Stress, ME = memory effect)
threshCS <- log(1+exp(ddmJAGS$BUGSoutput$sims.list$mu_alphaC_Stress))
threshMS <- log(1+exp(ddmJAGS$BUGSoutput$sims.list$mu_alphaC_Stress+ddmJAGS$BUGSoutput$sims.list$mu_alphaM_Stress))
threshMES <- log(1+exp(ddmJAGS$BUGSoutput$sims.list$mu_alphaM_Stress))-log(2)
threshCN <- log(1+exp(ddmJAGS$BUGSoutput$sims.list$mu_alphaC_noStress))
threshMN <- log(1+exp(ddmJAGS$BUGSoutput$sims.list$mu_alphaC_noStress+ddmJAGS$BUGSoutput$sims.list$mu_alphaM_noStress))
threshMEN <- log(1+exp(ddmJAGS$BUGSoutput$sims.list$mu_alphaM_noStress))-log(2)
ndtCS <- log(1+exp(ddmJAGS$BUGSoutput$sims.list$mu_tauC_Stress))
ndtMS <- log(1+exp(ddmJAGS$BUGSoutput$sims.list$mu_tauC_Stress+ddmJAGS$BUGSoutput$sims.list$mu_tauM_Stress))
ndtMES <- log(1+exp(ddmJAGS$BUGSoutput$sims.list$mu_tauM_Stress))-log(2)
ndtCN <- log(1+exp(ddmJAGS$BUGSoutput$sims.list$mu_tauC_noStress))
ndtMN <- log(1+exp(ddmJAGS$BUGSoutput$sims.list$mu_tauC_noStress+ddmJAGS$BUGSoutput$sims.list$mu_tauM_noStress))
ndtMEN <- log(1+exp(ddmJAGS$BUGSoutput$sims.list$mu_tauM_noStress))-log(2)
driftCS <- log(1+exp(ddmJAGS$BUGSoutput$sims.list$mu_deltaC_Stress))
driftMS <- log(1+exp(ddmJAGS$BUGSoutput$sims.list$mu_deltaC_Stress+ddmJAGS$BUGSoutput$sims.list$mu_deltaM_Stress))
driftMES <- log(1+exp(ddmJAGS$BUGSoutput$sims.list$mu_deltaM_Stress))-log(2)
driftCN <- log(1+exp(ddmJAGS$BUGSoutput$sims.list$mu_deltaC_noStress))
driftMN <- log(1+exp(ddmJAGS$BUGSoutput$sims.list$mu_deltaC_noStress+ddmJAGS$BUGSoutput$sims.list$mu_deltaM_noStress))
driftMEN <- log(1+exp(ddmJAGS$BUGSoutput$sims.list$mu_deltaM_noStress))-log(2)

##plot each parameter with memory and stress effects
#threshold / boundary separation
png('Figures/DDMthresh',width = 600, height = 700)
histGran <- (0:400)/100
par(fig=c(0,1,.4,1))
hist(threshCS,histGran,xlim=c(1.8,3.2),ylim=c(0,1500),xaxt='n',yaxt='n',
     col=rgb(1,0,.5,1/2),xlab = '',main = 'Threshold')
hist(threshMS,histGran,add=T,col=rgb(1,0,0,1/2))
hist(threshMN,histGran,add=T,col=rgb(0,0,1,1/2))
hist(threshCN,histGran,add=T,col=rgb(0,.5,1,1/2))
axis(1,at = (1:6)/2)
axis(2,at = seq(0,1500,500))
legend(2.75,1500,legend = c('no Stress/Memory','Stress/Memory','no Stress/Control','Stress/Control'),
       fill = c(rgb(0,0,1,1/2),rgb(1,0,0,1/2),rgb(0,.5,1,1/2),rgb(1,0,.5,1/2)),bty='n')
memEffect <- threshMS+threshMN-threshCS-threshCN
stressEffect <- threshCS-threshCN+threshMS-threshMN
par(fig=c(0,.5,0,.6),new=T)
histGran <- (-13:50)/25
h1 <- hist(memEffect,histGran,xlim = c(-0.5,1.5),main='',xlab = 'Memory effect')
lines(c(0,0),c(0,max(h1$counts)),col='red',lwd=2)
HDI <- round(quantile(memEffect,c(.025,.975)),2)
text(.1,max(h1$counts)*2/3,paste('95% HDI = [',HDI[1],',',HDI[2],']'),cex=0.7)
par(fig=c(.5,1,0,.6),new=T)
histGran <- (-13:38)/20
h1 <- hist(stressEffect,histGran,xlim = c(-0.5,1.5),main='',xlab = 'Stress effect')
lines(c(0,0),c(0,max(h1$counts)),col='red',lwd=2)
HDI <- round(quantile(stressEffect,c(.025,.975)),2)
text(1.1,max(h1$counts)*2/3,paste('95% HDI = [',HDI[1],',',HDI[2],']'),cex=0.7)
dev.off()

#non-decision time
png('Figures/DDMndt',width = 600, height = 700)
histGran <- (40:140)/200
par(fig=c(0,1,.4,1))
hist(ndtCS,histGran,xlim=c(.2,.7),ylim=c(0,2000),xaxt='n',yaxt='n',
     col=rgb(1,0,.5,1/2),xlab = '',main = 'Non-decision time')
hist(ndtMS,histGran,add=T,col=rgb(1,0,0,1/2))
hist(ndtMN,histGran,add=T,col=rgb(0,0,1,1/2))
hist(ndtCN,histGran,add=T,col=rgb(0,.5,1,1/2))
axis(1,at = (2:7)/10)
axis(2,at = seq(0,2000,500))
legend(.25,2000,legend = c('no Stress/Memory','Stress/Memory','no Stress/Control','Stress/Control'),
       fill = c(rgb(0,0,1,1/2),rgb(1,0,0,1/2),rgb(0,.5,1,1/2),rgb(1,0,.5,1/2)),bty='n')
memEffect <- ndtMS+ndtMN-ndtCS-ndtCN
stressEffect <- ndtCS-ndtCN+ndtMS-ndtMN
par(fig=c(0,.5,0,.6),new=T)
histGran <- (-35:15)/50
h1 <- hist(memEffect,histGran,xlim = c(-.7,.3),main='',xlab = 'Memory effect')
lines(c(0,0),c(0,max(h1$counts)),col='red',lwd=2)
HDI <- round(quantile(memEffect,c(.025,.975)),2)
text(.05,max(h1$counts)*2/3,paste('95% HDI = [',HDI[1],',',HDI[2],']'),cex=0.7)
par(fig=c(.5,1,0,.6),new=T)
histGran <- (-25:25)/50
h1 <- hist(stressEffect,histGran,xlim = c(-0.5,0.5),main='',xlab = 'Stress effect')
lines(c(0,0),c(0,max(h1$counts)),col='red',lwd=2)
HDI <- round(quantile(stressEffect,c(.025,.975)),2)
text(-0.28,max(h1$counts)*2/3,paste('95% HDI = [',HDI[1],',',HDI[2],']'),cex=0.7)
dev.off()

#drift rate
png('Figures/DDMdrift',width = 600, height = 700)
histGran <- (0:150)/150
par(fig=c(0,1,.4,1))
hist(driftCS,histGran,xlim=c(0,1),ylim=c(0,4000),xaxt='n',yaxt='n',
     col=rgb(1,0,.5,1/2),xlab = '',main = 'Drift rate')
hist(driftMS,histGran,add=T,col=rgb(1,0,0,1/2))
hist(driftMN,histGran,add=T,col=rgb(0,0,1,1/2))
hist(driftCN,histGran,add=T,col=rgb(0,.5,1,1/2))
axis(1,at = (0:5)/5)
axis(2,at = seq(0,4000,1000))
legend(.7,4000,legend = c('no Stress/Memory','Stress/Memory','no Stress/Control','Stress/Control'),
       fill = c(rgb(0,0,1,1/2),rgb(1,0,0,1/2),rgb(0,.5,1,1/2),rgb(1,0,.5,1/2)),bty='n')
memEffect <- driftMS+driftMN-driftCS-driftCN
stressEffect <- driftCS-driftCN+driftMS-driftMN
par(fig=c(0,.5,0,.6),new=T)
histGran <- (-150:-50)/100
h1 <- hist(memEffect,histGran,xlim = c(-1.5,-.5),main='',xlab = 'Memory effect')
lines(c(0,0),c(0,max(h1$counts)),col='red',lwd=2)
HDI <- round(quantile(memEffect,c(.025,.975)),2)
text(-.7,max(h1$counts)*2/3,paste('95% HDI = [',HDI[1],',',HDI[2],']'),cex=0.7)
par(fig=c(.5,1,0,.6),new=T)
histGran <- (-70:60)/100
h1 <- hist(stressEffect,histGran,xlim = c(-0.7,0.6),main='',xlab = 'Stress effect')
lines(c(0,0),c(0,max(h1$counts)),col='red',lwd=2)
HDI <- round(quantile(stressEffect,c(.025,.975)),2)
text(-0.36,max(h1$counts)*2/3,paste('95% HDI = [',HDI[1],',',HDI[2],']'),cex=0.7)
dev.off()

### --- ### --- ### --- ### --- ###
### with starting point
### --- ### --- ### --- ### --- ###
nChains <- 8
ddmInits <- vector("list", nChains)
for (c in 1:nChains){
  ddmInits[[c]] <- list(
    'mu_alphaC_Stress'=rnorm(1,3),'sigma_alphaC_Stress'=abs(rnorm(1)),'mu_alphaM_Stress'=rnorm(1),'sigma_alphaM_Stress'=abs(rnorm(1)),
    'mu_tauC_Stress'=rnorm(1,-1),'sigma_tauC_Stress'=abs(rnorm(1)),'mu_tauM_Stress'=rnorm(1),'sigma_tauM_Stress'=abs(rnorm(1)),
    'mu_betaC_Stress'=rnorm(1),'sigma_betaC_Stress'=abs(rnorm(1)),'mu_betaM_Stress'=rnorm(1),'sigma_betaM_Stress'=abs(rnorm(1)),  
    'mu_deltaC_Stress'=rnorm(1),'sigma_deltaC_Stress'=abs(rnorm(1)),'mu_deltaM_Stress'=rnorm(1),'sigma_deltaM_Stress'=abs(rnorm(1)),
    'mu_alphaC_noStress'=rnorm(1,3),'sigma_alphaC_noStress'=abs(rnorm(1)),'mu_alphaM_noStress'=rnorm(1),'sigma_alphaM_noStress'=abs(rnorm(1)),
    'mu_tauC_noStress'=rnorm(1,-1),'sigma_tauC_noStress'=abs(rnorm(1)),'mu_tauM_noStress'=rnorm(1),'sigma_tauM_noStress'=abs(rnorm(1)),
    'mu_betaC_noStress'=rnorm(1),'sigma_betaC_noStress'=abs(rnorm(1)),'mu_betaM_noStress'=rnorm(1),'sigma_betaM_noStress'=abs(rnorm(1)),
    'mu_deltaC_noStress'=rnorm(1),'sigma_deltaC_noStress'=abs(rnorm(1)),'mu_deltaM_noStress'=rnorm(1),'sigma_deltaM_noStress'=abs(rnorm(1)),
    'alphaC_Stress'=rnorm(nsubj_Stress,3,1),'tauC_Stress'=runif(nsubj_Stress,-10,-2),'betaC_Stress'=rnorm(nsubj_Stress),'deltaC_Stress'=rnorm(nsubj_Stress),
    'alphaM_Stress'=rnorm(nsubj_Stress),'tauM_Stress'=rnorm(nsubj_Stress,0,1/6),'betaM_Stress'=rnorm(nsubj_Stress),'deltaM_Stress'=rnorm(nsubj_Stress),
    'alphaC_noStress'=rnorm(nsubj_noStress,3,1),'tauC_noStress'=runif(nsubj_noStress,-10,-2),'betaC_noStress'=rnorm(nsubj_noStress),'deltaC_noStress'=rnorm(nsubj_noStress),    
    'alphaM_noStress'=rnorm(nsubj_noStress),'tauM_noStress'=rnorm(nsubj_noStress,0,1/6),'betaM_noStress'=rnorm(nsubj_noStress),'deltaM_noStress'=rnorm(nsubj_noStress))
}

#parallel computing apparently not possible
T1<-Sys.time()
ddm2JAGS <- jags(ddmData,inits = ddmInits,
                #ddmJAGS <- jags.parallel(ddmData,inits = NULL,jags.seed = sample(100:999,1),
                parameters.to.save = c('mu_alphaC_Stress','sigma_alphaC_Stress','mu_alphaM_Stress','sigma_alphaM_Stress',
                                       'mu_tauC_Stress','sigma_tauC_Stress ','mu_tauM_Stress','sigma_tauM_Stress',
                                       'mu_betaC_Stress','sigma_betaC_Stress','mu_betaM_Stress','sigma_betaM_Stress',
                                       'mu_deltaC_Stress','sigma_deltaC_Stress','mu_deltaM_Stress','sigma_deltaM_Stress',
                                       'mu_alphaC_noStress','sigma_alphaC_noStress','mu_alphaM_noStress','sigma_alphaM_noStress',
                                       'mu_tauC_noStress','sigma_tauC_noStress ','mu_tauM_noStress','sigma_tauM_noStress',
                                       'mu_betaC_noStress','sigma_betaC_noStress','mu_betaM_noStress','sigma_betaM_noStress',
                                       'mu_deltaC_noStress','sigma_deltaC_noStress','mu_deltaM_noStress','sigma_deltaM_noStress',
                                       'alphaC_Stress','tauC_Stress','betaC_Stress','deltaC_Stress','alphaM_Stress','tauM_Stress','betaM_Stress','deltaM_Stress',
                                       'alphaC_noStress','tauC_noStress','betaC_noStress','deltaC_noStress','alphaM_noStress','tauM_noStress','betaM_noStress','deltaM_noStress'),
                model.file = "BayesModel_DDM_sp.txt",working.directory = 'BayesModels',
                n.chains = nChains, n.iter = 40000, n.burnin = 20000,n.thin = 10, DIC = TRUE,jags.module = c("glm","dic","wiener"))
T2<-Sys.time()
#T2-T1

# check convergence
ddmRhats <- ddm2JAGS$BUGSoutput$summary[,8] #check with max(Rhats), which should ideally be < 1.01 (1.05 would also be okay)
max(ddmRhats)

#get group parameter estimates in transformed form (C = control, M = memory, S = Stress, N = no Stress, ME = memory effect)
threshCS <- log(1+exp(ddm2JAGS$BUGSoutput$sims.list$mu_alphaC_Stress))
threshMS <- log(1+exp(ddm2JAGS$BUGSoutput$sims.list$mu_alphaC_Stress+ddm2JAGS$BUGSoutput$sims.list$mu_alphaM_Stress))
threshCN <- log(1+exp(ddm2JAGS$BUGSoutput$sims.list$mu_alphaC_noStress))
threshMN <- log(1+exp(ddm2JAGS$BUGSoutput$sims.list$mu_alphaC_noStress+ddm2JAGS$BUGSoutput$sims.list$mu_alphaM_noStress))
ndtCS <- log(1+exp(ddm2JAGS$BUGSoutput$sims.list$mu_tauC_Stress))
ndtMS <- log(1+exp(ddm2JAGS$BUGSoutput$sims.list$mu_tauC_Stress+ddm2JAGS$BUGSoutput$sims.list$mu_tauM_Stress))
ndtCN <- log(1+exp(ddm2JAGS$BUGSoutput$sims.list$mu_tauC_noStress))
ndtMN <- log(1+exp(ddm2JAGS$BUGSoutput$sims.list$mu_tauC_noStress+ddm2JAGS$BUGSoutput$sims.list$mu_tauM_noStress))
spCS <- pnorm(ddm2JAGS$BUGSoutput$sims.list$mu_betaC_Stress)
spMS <- pnorm(ddm2JAGS$BUGSoutput$sims.list$mu_betaC_Stress+ddm2JAGS$BUGSoutput$sims.list$mu_betaM_Stress)
spCN <- pnorm(ddm2JAGS$BUGSoutput$sims.list$mu_betaC_noStress)
spMN <- pnorm(ddm2JAGS$BUGSoutput$sims.list$mu_betaC_noStress+ddm2JAGS$BUGSoutput$sims.list$mu_betaM_noStress)
driftCS <- log(1+exp(ddm2JAGS$BUGSoutput$sims.list$mu_deltaC_Stress))
driftMS <- log(1+exp(ddm2JAGS$BUGSoutput$sims.list$mu_deltaC_Stress+ddm2JAGS$BUGSoutput$sims.list$mu_deltaM_Stress))
driftCN <- log(1+exp(ddm2JAGS$BUGSoutput$sims.list$mu_deltaC_noStress))
driftMN <- log(1+exp(ddm2JAGS$BUGSoutput$sims.list$mu_deltaC_noStress+ddm2JAGS$BUGSoutput$sims.list$mu_deltaM_noStress))

##plot each parameter with memory and stress effects
#threshold / boundary separation
png('Figures/DDM2thresh',width = 600, height = 700)
histGran <- (0:400)/100
par(fig=c(0,1,.4,1))
hist(threshCS,histGran,xlim=c(1.8,3.2),ylim=c(0,1500),xaxt='n',yaxt='n',
     col=rgb(1,0,.5,1/2),xlab = '',main = 'Threshold')
hist(threshMS,histGran,add=T,col=rgb(1,0,0,1/2))
hist(threshMN,histGran,add=T,col=rgb(0,0,1,1/2))
hist(threshCN,histGran,add=T,col=rgb(0,.5,1,1/2))
axis(1,at = (1:6)/2)
axis(2,at = seq(0,1500,500))
legend(2.75,1500,legend = c('no Stress/Memory','Stress/Memory','no Stress/Control','Stress/Control'),
       fill = c(rgb(0,0,1,1/2),rgb(1,0,0,1/2),rgb(0,.5,1,1/2),rgb(1,0,.5,1/2)),bty='n')
memEffect <- threshMS+threshMN-threshCS-threshCN
stressEffect <- threshCS-threshCN+threshMS-threshMN
par(fig=c(0,.5,0,.6),new=T)
histGran <- (-13:50)/25
h1 <- hist(memEffect,histGran,xlim = c(-0.5,1.5),main='',xlab = 'Memory effect')
lines(c(0,0),c(0,max(h1$counts)),col='red',lwd=2)
HDI <- round(quantile(memEffect,c(.025,.975)),2)
text(.1,max(h1$counts)*2/3,paste('95% HDI = [',HDI[1],',',HDI[2],']'),cex=0.7)
par(fig=c(.5,1,0,.6),new=T)
histGran <- (-13:38)/25
h1 <- hist(stressEffect,histGran,xlim = c(-0.5,1.5),main='',xlab = 'Stress effect')
lines(c(0,0),c(0,max(h1$counts)),col='red',lwd=2)
HDI <- round(quantile(stressEffect,c(.025,.975)),2)
text(1.1,max(h1$counts)*2/3,paste('95% HDI = [',HDI[1],',',HDI[2],']'),cex=0.7)
dev.off()

#non-decision time
png('Figures/DDM2ndt',width = 600, height = 700)
histGran <- (40:140)/200
par(fig=c(0,1,.4,1))
hist(ndtCS,histGran,xlim=c(.2,.7),ylim=c(0,2000),xaxt='n',yaxt='n',
     col=rgb(1,0,.5,1/2),xlab = '',main = 'Non-decision time')
hist(ndtMS,histGran,add=T,col=rgb(1,0,0,1/2))
hist(ndtMN,histGran,add=T,col=rgb(0,0,1,1/2))
hist(ndtCN,histGran,add=T,col=rgb(0,.5,1,1/2))
axis(1,at = (2:7)/10)
axis(2,at = seq(0,2000,500))
legend(.25,2000,legend = c('no Stress/Memory','Stress/Memory','no Stress/Control','Stress/Control'),
       fill = c(rgb(0,0,1,1/2),rgb(1,0,0,1/2),rgb(0,.5,1,1/2),rgb(1,0,.5,1/2)),bty='n')
memEffect <- ndtMS+ndtMN-ndtCS-ndtCN
stressEffect <- ndtCS-ndtCN+ndtMS-ndtMN
par(fig=c(0,.5,0,.6),new=T)
histGran <- (-35:15)/50
h1 <- hist(memEffect,histGran,xlim = c(-.7,.3),main='',xlab = 'Memory effect')
lines(c(0,0),c(0,max(h1$counts)),col='red',lwd=2)
HDI <- round(quantile(memEffect,c(.025,.975)),2)
text(.05,max(h1$counts)*2/3,paste('95% HDI = [',HDI[1],',',HDI[2],']'),cex=0.7)
par(fig=c(.5,1,0,.6),new=T)
histGran <- (-25:25)/50
h1 <- hist(stressEffect,histGran,xlim = c(-0.5,0.5),main='',xlab = 'Stress effect')
lines(c(0,0),c(0,max(h1$counts)),col='red',lwd=2)
HDI <- round(quantile(stressEffect,c(.025,.975)),2)
text(-0.28,max(h1$counts)*2/3,paste('95% HDI = [',HDI[1],',',HDI[2],']'),cex=0.7)
dev.off()

#starting point
png('Figures/DDM2sp',width = 600, height = 700)
histGran <- (0:500)/500
par(fig=c(0,1,.4,1))
hist(spCS,histGran,xlim=c(.4,.6),ylim=c(0,2000),xaxt='n',yaxt='n',
     col=rgb(1,0,.5,1/2),xlab = '',main = 'Starting Point')
hist(spMS,histGran,add=T,col=rgb(1,0,0,1/2))
hist(spMN,histGran,add=T,col=rgb(0,0,1,1/2))
hist(spCN,histGran,add=T,col=rgb(0,.5,1,1/2))
axis(1,at = (0:20)/20)
axis(2,at = seq(0,2000,1000))
legend(.52,2000,legend = c('no Stress/Memory','Stress/Memory','no Stress/Control','Stress/Control'),
       fill = c(rgb(0,0,1,1/2),rgb(1,0,0,1/2),rgb(0,.5,1,1/2),rgb(1,0,.5,1/2)),bty='n')
memEffect <- spMS+spMN-spCS-spCN
stressEffect <- spCS-spCN+spMS-spMN
par(fig=c(0,.5,0,.6),new=T)
histGran <- (-250:250)/500
h1 <- hist(memEffect,histGran,xlim = c(-.1,.1),main='',xlab = 'Memory effect')
lines(c(0,0),c(0,max(h1$counts)),col='red',lwd=2)
HDI <- round(quantile(memEffect,c(.025,.975)),2)
text(-.7,max(h1$counts)*2/3,paste('95% HDI = [',HDI[1],',',HDI[2],']'),cex=0.7)
par(fig=c(.5,1,0,.6),new=T)
histGran <- (-250:250)/500
h1 <- hist(stressEffect,histGran,xlim = c(-0.1,0.1),main='',xlab = 'Stress effect')
lines(c(0,0),c(0,max(h1$counts)),col='red',lwd=2)
HDI <- round(quantile(stressEffect,c(.025,.975)),2)
text(-0.36,max(h1$counts)*2/3,paste('95% HDI = [',HDI[1],',',HDI[2],']'),cex=0.7)
dev.off()

#drift rate
png('Figures/DDM2drift',width = 600, height = 700)
histGran <- (0:150)/150
par(fig=c(0,1,.4,1))
hist(driftCS,histGran,xlim=c(0,1),ylim=c(0,4000),xaxt='n',yaxt='n',
     col=rgb(1,0,.5,1/2),xlab = '',main = 'Drift rate')
hist(driftMS,histGran,add=T,col=rgb(1,0,0,1/2))
hist(driftMN,histGran,add=T,col=rgb(0,0,1,1/2))
hist(driftCN,histGran,add=T,col=rgb(0,.5,1,1/2))
axis(1,at = (0:5)/5)
axis(2,at = seq(0,4000,1000))
legend(.7,4000,legend = c('no Stress/Memory','Stress/Memory','no Stress/Control','Stress/Control'),
       fill = c(rgb(0,0,1,1/2),rgb(1,0,0,1/2),rgb(0,.5,1,1/2),rgb(1,0,.5,1/2)),bty='n')
memEffect <- driftMS+driftMN-driftCS-driftCN
stressEffect <- driftCS-driftCN+driftMS-driftMN
par(fig=c(0,.5,0,.6),new=T)
histGran <- (-150:-50)/100
h1 <- hist(memEffect,histGran,xlim = c(-1.5,-.5),main='',xlab = 'Memory effect')
lines(c(0,0),c(0,max(h1$counts)),col='red',lwd=2)
HDI <- round(quantile(memEffect,c(.025,.975)),2)
text(-.7,max(h1$counts)*2/3,paste('95% HDI = [',HDI[1],',',HDI[2],']'),cex=0.7)
par(fig=c(.5,1,0,.6),new=T)
histGran <- (-70:60)/100
h1 <- hist(stressEffect,histGran,xlim = c(-0.7,0.6),main='',xlab = 'Stress effect')
lines(c(0,0),c(0,max(h1$counts)),col='red',lwd=2)
HDI <- round(quantile(stressEffect,c(.025,.975)),2)
text(-0.36,max(h1$counts)*2/3,paste('95% HDI = [',HDI[1],',',HDI[2],']'),cex=0.7)
dev.off()


### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### ---
### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- 
### --- DDM with memory state and memory bias
### --- ### --- ### --- ###
### --- ### --- ### --- ###

#data definition
VL1 <- decision$zValLeft[decision$eligTrials&(stressTrialIdx==1)] #left item value; stress group
VL2 <- decision$zValLeft[decision$eligTrials&(stressTrialIdx==0)]
VR1 <- decision$zValRight[decision$eligTrials&(stressTrialIdx==1)]
VR2 <- decision$zValRight[decision$eligTrials&(stressTrialIdx==0)]
ML1 <- decision$LeftMem[decision$eligTrials&(stressTrialIdx==1)] #left item remembered; stress group
ML2 <- decision$LeftMem[decision$eligTrials&(stressTrialIdx==0)]
MR1 <- decision$RightMem[decision$eligTrials&(stressTrialIdx==1)] 
MR2 <- decision$RightMem[decision$eligTrials&(stressTrialIdx==0)]
RT1 <- decision$RespTime[decision$eligTrials&(stressTrialIdx==1)]
RT2 <- decision$RespTime[decision$eligTrials&(stressTrialIdx==0)]
RT1[decision$Choice[decision$eligTrials&(stressTrialIdx==1)]==2] = -RT1[decision$Choice[decision$eligTrials&(stressTrialIdx==1)]==2] #for dwiener function in JAGS
RT2[decision$Choice[decision$eligTrials&(stressTrialIdx==0)]==2] = -RT2[decision$Choice[decision$eligTrials&(stressTrialIdx==0)]==2] #for dwiener function in JAGS
M1 <- decision$IsControlBlock[decision$eligTrials&(stressTrialIdx==1)]==0 #memory trials, stress group
M2 <- decision$IsControlBlock[decision$eligTrials&(stressTrialIdx==0)]==0 #memory trials, no stress group

# get data and initial values together and specify the model
ddmData <- list('N1'=N1,'N2'=N2,'S1'=S1,'S2'=S2,'M1'=M1,'M2'=M2,
                'VL1'=VL1,'VL2'=VL2,'VR1'=VR1,'VR2'=VR2,'ML1'=ML1,'ML2'=ML2,'MR1'=MR1,'MR2'=MR2,'RT1'=RT1,'RT2'=RT2,
                'nsubj_Stress'=nsubj_Stress,'nsubj_noStress'=nsubj_noStress,'priorSD'=priorSD)

#this time, initial values are required
nChains <- 8
ddmInits <- vector("list", nChains)
for (c in 1:nChains){
  ddmInits[[c]] <- list(
    'mu_alphaC_Stress'=rnorm(1,3),'sigma_alphaC_Stress'=abs(rnorm(1)),'mu_alphaM_Stress'=rnorm(1),'sigma_alphaM_Stress'=abs(rnorm(1)),
    'mu_tauC_Stress'=rnorm(1,-1),'sigma_tauC_Stress'=abs(rnorm(1)),'mu_tauM_Stress'=rnorm(1),'sigma_tauM_Stress'=abs(rnorm(1)),
    'mu_betaC_Stress'=rnorm(1),'sigma_betaC_Stress'=abs(rnorm(1)),'mu_betaM_Stress'=rnorm(1),'sigma_betaM_Stress'=abs(rnorm(1)),  
    'mu_deltaC_Stress'=rnorm(1),'sigma_deltaC_Stress'=abs(rnorm(1)),'mu_deltaM_Stress'=rnorm(1),'sigma_deltaM_Stress'=abs(rnorm(1)),
    'mu_MB_Stress'=rnorm(1),'sigma_MB_Stress'=abs(rnorm(1)),
    'mu_alphaC_noStress'=rnorm(1,3),'sigma_alphaC_noStress'=abs(rnorm(1)),'mu_alphaM_noStress'=rnorm(1),'sigma_alphaM_noStress'=abs(rnorm(1)),
    'mu_tauC_noStress'=rnorm(1,-1),'sigma_tauC_noStress'=abs(rnorm(1)),'mu_tauM_noStress'=rnorm(1),'sigma_tauM_noStress'=abs(rnorm(1)),
    'mu_betaC_noStress'=rnorm(1),'sigma_betaC_noStress'=abs(rnorm(1)),'mu_betaM_noStress'=rnorm(1),'sigma_betaM_noStress'=abs(rnorm(1)),
    'mu_deltaC_noStress'=rnorm(1),'sigma_deltaC_noStress'=abs(rnorm(1)),'mu_deltaM_noStress'=rnorm(1),'sigma_deltaM_noStress'=abs(rnorm(1)),
    'mu_MB_noStress'=rnorm(1),'sigma_MB_noStress'=abs(rnorm(1)),
    'alphaC_Stress'=rnorm(nsubj_Stress,3,1),'tauC_Stress'=runif(nsubj_Stress,-10,-2),'betaC_Stress'=rnorm(nsubj_Stress),'deltaC_Stress'=rnorm(nsubj_Stress),
    'alphaM_Stress'=rnorm(nsubj_Stress),'tauM_Stress'=rnorm(nsubj_Stress,0,1/6),'betaM_Stress'=rnorm(nsubj_Stress),'deltaM_Stress'=rnorm(nsubj_Stress),
    'MBS'=rnorm(nsubj_Stress),
    'alphaC_noStress'=rnorm(nsubj_noStress,3,1),'tauC_noStress'=runif(nsubj_noStress,-10,-2),'betaC_noStress'=rnorm(nsubj_noStress),'deltaC_noStress'=rnorm(nsubj_noStress),    
    'alphaM_noStress'=rnorm(nsubj_noStress),'tauM_noStress'=rnorm(nsubj_noStress,0,1/6),'betaM_noStress'=rnorm(nsubj_noStress),'deltaM_noStress'=rnorm(nsubj_noStress),
    'MBN'=rnorm(nsubj_noStress))
}

#parallel computing apparently not possible
T1<-Sys.time()
ddm_MB_JAGS <- jags(ddmData,inits = ddmInits,
                #ddmJAGS <- jags.parallel(ddmData,inits = NULL,jags.seed = sample(100:999,1),
                parameters.to.save = c('mu_alphaC_Stress','sigma_alphaC_Stress','mu_alphaM_Stress','sigma_alphaM_Stress',
                                       'mu_tauC_Stress','sigma_tauC_Stress ','mu_tauM_Stress','sigma_tauM_Stress',
                                       'mu_betaC_Stress','sigma_betaC_Stress','mu_betaM_Stress','sigma_betaM_Stress',
                                       'mu_deltaC_Stress','sigma_deltaC_Stress','mu_deltaM_Stress','sigma_deltaM_Stress',
                                       'mu_MB_Stress','sigma_MB_Stress',
                                       'mu_alphaC_noStress','sigma_alphaC_noStress','mu_alphaM_noStress','sigma_alphaM_noStress',
                                       'mu_tauC_noStress','sigma_tauC_noStress ','mu_tauM_noStress','sigma_tauM_noStress',
                                       'mu_betaC_noStress','sigma_betaC_noStress','mu_betaM_noStress','sigma_betaM_noStress',
                                       'mu_deltaC_noStress','sigma_deltaC_noStress','mu_deltaM_noStress','sigma_deltaM_noStress',
                                       'mu_MB_noStress','sigma_MB_noStress',
                                       'alphaC_Stress','tauC_Stress','betaC_Stress','deltaC_Stress','alphaM_Stress','tauM_Stress','betaM_Stress','deltaM_Stress',
                                       'MBS',
                                       'alphaC_noStress','tauC_noStress','betaC_noStress','deltaC_noStress','alphaM_noStress','tauM_noStress','betaM_noStress','deltaM_noStress',
                                       'MBN'),
                model.file = "BayesModel_DDM_MB.txt",working.directory = 'BayesModels',
                n.chains = nChains, n.iter = 40000, n.burnin = 20000,n.thin = 10, DIC = TRUE,jags.module = c("glm","dic","wiener"))
T2<-Sys.time()
#T2-T1

ddmRhats <- ddm_MB_JAGS$BUGSoutput$summary[,8] #check with max(Rhats), which should ideally be < 1.01 (1.05 would also be okay)
max(ddmRhats)

#get group parameter estimates in transformed form (C = control, M = memory, S = Stress, N = no Stress, ME = memory effect)
threshCS <- log(1+exp(ddm_MB_JAGS$BUGSoutput$sims.list$mu_alphaC_Stress))
threshMS <- log(1+exp(ddm_MB_JAGS$BUGSoutput$sims.list$mu_alphaC_Stress+ddm_MB_JAGS$BUGSoutput$sims.list$mu_alphaM_Stress))
threshCN <- log(1+exp(ddm_MB_JAGS$BUGSoutput$sims.list$mu_alphaC_noStress))
threshMN <- log(1+exp(ddm_MB_JAGS$BUGSoutput$sims.list$mu_alphaC_noStress+ddm_MB_JAGS$BUGSoutput$sims.list$mu_alphaM_noStress))
ndtCS <- log(1+exp(ddm_MB_JAGS$BUGSoutput$sims.list$mu_tauC_Stress))
ndtMS <- log(1+exp(ddm_MB_JAGS$BUGSoutput$sims.list$mu_tauC_Stress+ddm_MB_JAGS$BUGSoutput$sims.list$mu_tauM_Stress))
ndtCN <- log(1+exp(ddm_MB_JAGS$BUGSoutput$sims.list$mu_tauC_noStress))
ndtMN <- log(1+exp(ddm_MB_JAGS$BUGSoutput$sims.list$mu_tauC_noStress+ddm_MB_JAGS$BUGSoutput$sims.list$mu_tauM_noStress))
spCS <- pnorm(ddm_MB_JAGS$BUGSoutput$sims.list$mu_betaC_Stress)
spMS <- pnorm(ddm_MB_JAGS$BUGSoutput$sims.list$mu_betaC_Stress+ddm_MB_JAGS$BUGSoutput$sims.list$mu_betaM_Stress)
spCN <- pnorm(ddm_MB_JAGS$BUGSoutput$sims.list$mu_betaC_noStress)
spMN <- pnorm(ddm_MB_JAGS$BUGSoutput$sims.list$mu_betaC_noStress+ddm_MB_JAGS$BUGSoutput$sims.list$mu_betaM_noStress)
driftCS <- log(1+exp(ddm_MB_JAGS$BUGSoutput$sims.list$mu_deltaC_Stress))
driftMS <- log(1+exp(ddm_MB_JAGS$BUGSoutput$sims.list$mu_deltaC_Stress+ddm_MB_JAGS$BUGSoutput$sims.list$mu_deltaM_Stress))
driftCN <- log(1+exp(ddm_MB_JAGS$BUGSoutput$sims.list$mu_deltaC_noStress))
driftMN <- log(1+exp(ddm_MB_JAGS$BUGSoutput$sims.list$mu_deltaC_noStress+ddm_MB_JAGS$BUGSoutput$sims.list$mu_deltaM_noStress))
membiasS <- ddm_MB_JAGS$BUGSoutput$sims.list$mu_MB_Stress
membiasN <- ddm_MB_JAGS$BUGSoutput$sims.list$mu_MB_noStress

##plot each parameter with memory and stress effects
#threshold / boundary separation
png('Figures/DDM_MB_thresh',width = 600, height = 700)
histGran <- (0:400)/100
par(fig=c(0,1,.4,1))
hist(threshCS,histGran,xlim=c(1.8,3.2),ylim=c(0,1500),xaxt='n',yaxt='n',
     col=rgb(1,0,.5,1/2),xlab = '',main = 'Boundary separation')
hist(threshMS,histGran,add=T,col=rgb(1,0,0,1/2))
hist(threshMN,histGran,add=T,col=rgb(0,0,1,1/2))
hist(threshCN,histGran,add=T,col=rgb(0,.5,1,1/2))
axis(1,at = (1:6)/2)
axis(2,at = seq(0,1500,500))
legend(2.75,1500,legend = c('no Stress/Memory','Stress/Memory','no Stress/Control','Stress/Control'),
       fill = c(rgb(0,0,1,1/2),rgb(1,0,0,1/2),rgb(0,.5,1,1/2),rgb(1,0,.5,1/2)),bty='n')
memEffect <- threshMS+threshMN-threshCS-threshCN
stressEffect <- threshCS-threshCN+threshMS-threshMN
par(fig=c(0,.5,0,.6),new=T)
histGran <- (-13:50)/25
h1 <- hist(memEffect,histGran,xlim = c(-0.5,1.5),main='',xlab = 'Memory effect')
lines(c(0,0),c(0,max(h1$counts)),col='red',lwd=2)
HDI <- round(quantile(memEffect,c(.025,.975)),2)
text(.1,max(h1$counts)*2/3,paste('95% HDI = [',HDI[1],',',HDI[2],']'),cex=0.7)
par(fig=c(.5,1,0,.6),new=T)
histGran <- (-20:40)/25
h1 <- hist(stressEffect,histGran,xlim = c(-0.5,1.5),main='',xlab = 'Stress effect')
lines(c(0,0),c(0,max(h1$counts)),col='red',lwd=2)
HDI <- round(quantile(stressEffect,c(.025,.975)),2)
text(1.1,max(h1$counts)*2/3,paste('95% HDI = [',HDI[1],',',HDI[2],']'),cex=0.7)
dev.off()

#non-decision time
png('Figures/DDM_MB_ndt',width = 600, height = 700)
histGran <- (40:140)/200
par(fig=c(0,1,.4,1))
hist(ndtCS,histGran,xlim=c(.2,.7),ylim=c(0,2000),xaxt='n',yaxt='n',
     col=rgb(1,0,.5,1/2),xlab = '',main = 'Non-decision time')
hist(ndtMS,histGran,add=T,col=rgb(1,0,0,1/2))
hist(ndtMN,histGran,add=T,col=rgb(0,0,1,1/2))
hist(ndtCN,histGran,add=T,col=rgb(0,.5,1,1/2))
axis(1,at = (2:7)/10)
axis(2,at = seq(0,2000,500))
legend(.25,2000,legend = c('no Stress/Memory','Stress/Memory','no Stress/Control','Stress/Control'),
       fill = c(rgb(0,0,1,1/2),rgb(1,0,0,1/2),rgb(0,.5,1,1/2),rgb(1,0,.5,1/2)),bty='n')
memEffect <- ndtMS+ndtMN-ndtCS-ndtCN
stressEffect <- ndtCS-ndtCN+ndtMS-ndtMN
par(fig=c(0,.5,0,.6),new=T)
histGran <- (-35:15)/50
h1 <- hist(memEffect,histGran,xlim = c(-.7,.3),main='',xlab = 'Memory effect')
lines(c(0,0),c(0,max(h1$counts)),col='red',lwd=2)
HDI <- round(quantile(memEffect,c(.025,.975)),2)
text(.05,max(h1$counts)*2/3,paste('95% HDI = [',HDI[1],',',HDI[2],']'),cex=0.7)
par(fig=c(.5,1,0,.6),new=T)
histGran <- (-25:25)/50
h1 <- hist(stressEffect,histGran,xlim = c(-0.5,0.5),main='',xlab = 'Stress effect')
lines(c(0,0),c(0,max(h1$counts)),col='red',lwd=2)
HDI <- round(quantile(stressEffect,c(.025,.975)),2)
text(-0.28,max(h1$counts)*2/3,paste('95% HDI = [',HDI[1],',',HDI[2],']'),cex=0.7)
dev.off()

#starting point
png('Figures/DDM_MB_sp',width = 600, height = 700)
histGran <- (0:500)/500
par(fig=c(0,1,.4,1))
hist(spCS,histGran,xlim=c(.4,.6),ylim=c(0,2000),xaxt='n',yaxt='n',
     col=rgb(1,0,.5,1/2),xlab = '',main = 'Starting Point')
hist(spMS,histGran,add=T,col=rgb(1,0,0,1/2))
hist(spMN,histGran,add=T,col=rgb(0,0,1,1/2))
hist(spCN,histGran,add=T,col=rgb(0,.5,1,1/2))
axis(1,at = (0:20)/20)
axis(2,at = seq(0,2000,1000))
legend(.52,2000,legend = c('no Stress/Memory','Stress/Memory','no Stress/Control','Stress/Control'),
       fill = c(rgb(0,0,1,1/2),rgb(1,0,0,1/2),rgb(0,.5,1,1/2),rgb(1,0,.5,1/2)),bty='n')
memEffect <- spMS+spMN-spCS-spCN
stressEffect <- spCS-spCN+spMS-spMN
par(fig=c(0,.5,0,.6),new=T)
histGran <- (-250:250)/500
h1 <- hist(memEffect,histGran,xlim = c(-.1,.1),main='',xlab = 'Memory effect')
lines(c(0,0),c(0,max(h1$counts)),col='red',lwd=2)
HDI <- round(quantile(memEffect,c(.025,.975)),2)
text(-.7,max(h1$counts)*2/3,paste('95% HDI = [',HDI[1],',',HDI[2],']'),cex=0.7)
par(fig=c(.5,1,0,.6),new=T)
histGran <- (-250:250)/500
h1 <- hist(stressEffect,histGran,xlim = c(-0.1,0.1),main='',xlab = 'Stress effect')
lines(c(0,0),c(0,max(h1$counts)),col='red',lwd=2)
HDI <- round(quantile(stressEffect,c(.025,.975)),2)
text(-0.36,max(h1$counts)*2/3,paste('95% HDI = [',HDI[1],',',HDI[2],']'),cex=0.7)
dev.off()

#drift rate
png('Figures/DDM_MB_drift',width = 600, height = 700)
histGran <- (0:150)/150
par(fig=c(0,1,.4,1))
hist(driftCS,histGran,xlim=c(0,1),ylim=c(0,4000),xaxt='n',yaxt='n',
     col=rgb(1,0,.5,1/2),xlab = '',main = 'Drift slope')
hist(driftMS,histGran,add=T,col=rgb(1,0,0,1/2))
hist(driftMN,histGran,add=T,col=rgb(0,0,1,1/2))
hist(driftCN,histGran,add=T,col=rgb(0,.5,1,1/2))
axis(1,at = (0:5)/5)
axis(2,at = seq(0,4000,1000))
legend(.7,4000,legend = c('no Stress/Memory','Stress/Memory','no Stress/Control','Stress/Control'),
       fill = c(rgb(0,0,1,1/2),rgb(1,0,0,1/2),rgb(0,.5,1,1/2),rgb(1,0,.5,1/2)),bty='n')
memEffect <- driftMS+driftMN-driftCS-driftCN
stressEffect <- driftCS-driftCN+driftMS-driftMN
par(fig=c(0,.5,0,.6),new=T)
histGran <- (-150:-50)/100
h1 <- hist(memEffect,histGran,xlim = c(-1.5,-.5),main='',xlab = 'Memory effect')
lines(c(0,0),c(0,max(h1$counts)),col='red',lwd=2)
HDI <- round(quantile(memEffect,c(.025,.975)),2)
text(-.7,max(h1$counts)*2/3,paste('95% HDI = [',HDI[1],',',HDI[2],']'),cex=0.7)
par(fig=c(.5,1,0,.6),new=T)
histGran <- (-70:60)/100
h1 <- hist(stressEffect,histGran,xlim = c(-0.7,0.6),main='',xlab = 'Stress effect')
lines(c(0,0),c(0,max(h1$counts)),col='red',lwd=2)
HDI <- round(quantile(stressEffect,c(.025,.975)),2)
text(-0.36,max(h1$counts)*2/3,paste('95% HDI = [',HDI[1],',',HDI[2],']'),cex=0.7)
dev.off()

#memory bias
png('Figures/DDM_MB_membias',width = 600, height = 700)
histGran <- (-120:120)/100
par(fig=c(0,1,.4,1))
hist(membiasS,histGran,xlim=c(-1,1),ylim=c(0,800),xaxt='n',yaxt='n',
     col=rgb(1,0,0,1/2),xlab = '',main = 'Memory Bias')
hist(membiasN,histGran,add=T,col=rgb(0,0,1,1/2))
axis(1,at = (-2:2)/2)
axis(2,at = seq(0,1000,200))
legend(-1,800,legend = c('no Stress','Stress'),
       fill = c(rgb(0,0,1,1/2),rgb(1,0,0,1/2)),bty='n')
stressEffect <- membiasS-membiasN
par(fig=c(0,.5,0,.6),new=T)
histGran <-  (-60:60)/50
h1 <- hist(stressEffect,histGran,xlim = c(-1.2,1.2),main='',xlab = 'Stress effect')
lines(c(0,0),c(0,max(h1$counts)),col='red',lwd=2)
HDI <- round(quantile(stressEffect,c(.025,.975)),2)
text(-.7,max(h1$counts)*2/3,paste('95% HDI = [',HDI[1],',',HDI[2],']'),cex=0.7)
dev.off()



###posterior predictives
#function for simulation (this is just copy-pasted from below with minor changes; so, rtdists package is used to simulated the DDM)
sim_ddm_MB  <- function(parameters,vC,vI,mC,mI){

  #extract free parameters; note that order of parameters follows rtdists, not dwiener
  alpha <- parameters[1]
  beta <- parameters[2]*alpha #don't forget: sp is NOT relative in rtdists (but it is in dwiener; gosh!)
  tau <- parameters[3]
  delta <- parameters[4]
  memBias <- parameters[5]

  #compute trial-wise value difference
  VD <- delta*(mC*(vC-memBias)-mI*(vI-memBias))

  #Simulate the DDM and store choices and RTs
  nTrials <- length(VD)
  choices <- rep(NA, nTrials)
  rts <- rep(NA, nTrials)
  for (i in 1:nTrials) {
    result <- rdiffusion(1,alpha,VD[i],tau,beta)

    choices[i] <- ifelse(result$response == "upper", 2, 1)
    rts[i] <- result$rt
  }
  simdata<-cbind(choices, rts)

    # Return the choices and RTs
  return(simdata)
}
 
nPostPred <- 40
ppSample <- sample(ddm_MB_JAGS$BUGSoutput$n.sims,nPostPred) #draw samples from posterior (without replacement)
ppDat_M_Stress_choices <- matrix(NA,nrow = sum(M1==TRUE),ncol = nPostPred)
ppDat_M_Stress_rts <- matrix(NA,nrow = sum(M1==TRUE),ncol = nPostPred)
ppDat_C_Stress_choices <- matrix(NA,nrow = sum(M1==FALSE),ncol = nPostPred)
ppDat_C_Stress_rts <- matrix(NA,nrow = sum(M1==FALSE),ncol = nPostPred)
ppDat_M_noStress_choices <- matrix(NA,nrow = sum(M2==TRUE),ncol = nPostPred)
ppDat_M_noStress_rts <- matrix(NA,nrow = sum(M2==TRUE),ncol = nPostPred)
ppDat_C_noStress_choices <- matrix(NA,nrow = sum(M2==FALSE),ncol = nPostPred)
ppDat_C_noStress_rts <- matrix(NA,nrow = sum(M2==FALSE),ncol = nPostPred)
for (p in 1:nPostPred){
  x = ppSample[p] #current sample
  simDat_M_Stress<-data.frame(choices = numeric(), rts = numeric())
  simDat_C_Stress<-data.frame(choices = numeric(), rts = numeric())
  simDat_M_noStress<-data.frame(choices = numeric(), rts = numeric())
  simDat_C_noStress<-data.frame(choices = numeric(), rts = numeric())
  #start with Stress group
  for  (s in unique(S1)){
    for (t in 1:2){ #loop over choice types (memory, control)
      if (t == 1){ #memory trials
        ctrials = (S1==s)&(M1==TRUE)
        parameters <- c(log(1+exp(ddm_MB_JAGS$BUGSoutput$sims.list$alphaC_Stress[x,s]+ddm_MB_JAGS$BUGSoutput$sims.list$alphaM_Stress[x,s])),
                        pnorm(ddm_MB_JAGS$BUGSoutput$sims.list$betaC_Stress[x,s]+ddm_MB_JAGS$BUGSoutput$sims.list$betaM_Stress[x,s]),
                        log(1+exp(ddm_MB_JAGS$BUGSoutput$sims.list$tauC_Stress[x,s]+ddm_MB_JAGS$BUGSoutput$sims.list$tauM_Stress[x,s])),
                        log(1+exp(ddm_MB_JAGS$BUGSoutput$sims.list$deltaM_Stress[x,s])),
                        ddm_MB_JAGS$BUGSoutput$sims.list$MBS[x,s])
      } else {
        ctrials = (S1==s)&(M1==FALSE)
        parameters <- c(log(1+exp(ddm_MB_JAGS$BUGSoutput$sims.list$alphaC_Stress[x,s])),
                        pnorm(ddm_MB_JAGS$BUGSoutput$sims.list$betaC_Stress[x,s]),
                        log(1+exp(ddm_MB_JAGS$BUGSoutput$sims.list$tauC_Stress[x,s])),
                        log(1+exp(ddm_MB_JAGS$BUGSoutput$sims.list$deltaC_Stress[x,s])),
                        0) #no memory bias in conventional decisions
      }

      #arrange input data
      vL <- VL1[ctrials]
      vR <- VR1[ctrials]
      mL <- ML1[ctrials]
      mR <- MR1[ctrials]
      if (t==2){
        mL[mL==FALSE] = TRUE #no forgetting in control decisions
        mR[mR==FALSE] = TRUE #no forgetting in control decisions
      }
      #recode everything into consistent vs. inconsistent choices
      correctL <- (vL>vR) #whether left option is better
      vC <- vL*(correctL==TRUE)+vR*(correctL==FALSE)
      vI <- vL*(correctL==FALSE)+vR*(correctL==TRUE)
      mC <- mL*(correctL==TRUE)+mR*(correctL==FALSE)
      mI <- mL*(correctL==FALSE)+mR*(correctL==TRUE)

      #simulate
      simDat <- sim_ddm_MB(parameters,vC,vI,mC,mI)
      if (t==1){
        simDat_M_Stress <- rbind(simDat_M_Stress, simDat)
      } else {
        simDat_C_Stress <- rbind(simDat_C_Stress, simDat)
      }
    }
  }
  
  #go on with noStress group
  for  (s in unique(S2)){
    for (t in 1:2){ #loop over choice types (memory, control)
      if (t == 1){ #memory trials
        ctrials = (S2==s)&(M2==TRUE)
        parameters <- c(log(1+exp(ddm_MB_JAGS$BUGSoutput$sims.list$alphaC_noStress[x,s]+ddm_MB_JAGS$BUGSoutput$sims.list$alphaM_noStress[x,s])),
                        pnorm(ddm_MB_JAGS$BUGSoutput$sims.list$betaC_noStress[x,s]+ddm_MB_JAGS$BUGSoutput$sims.list$betaM_noStress[x,s]),
                        log(1+exp(ddm_MB_JAGS$BUGSoutput$sims.list$tauC_noStress[x,s]+ddm_MB_JAGS$BUGSoutput$sims.list$tauM_noStress[x,s])),
                        log(1+exp(ddm_MB_JAGS$BUGSoutput$sims.list$deltaM_noStress[x,s])),
                        ddm_MB_JAGS$BUGSoutput$sims.list$MBN[x,s])
      } else {
        ctrials = (S2==s)&(M2==FALSE)
        parameters <- c(log(1+exp(ddm_MB_JAGS$BUGSoutput$sims.list$alphaC_noStress[x,s])),
                        pnorm(ddm_MB_JAGS$BUGSoutput$sims.list$betaC_noStress[x,s]),
                        log(1+exp(ddm_MB_JAGS$BUGSoutput$sims.list$tauC_noStress[x,s])),
                        log(1+exp(ddm_MB_JAGS$BUGSoutput$sims.list$deltaC_noStress[x,s])),
                        0) #no memory bias in conventional decisions
      }
      
      #arrange input data
      vL <- VL2[ctrials]
      vR <- VR2[ctrials]
      mL <- ML2[ctrials]
      mR <- MR2[ctrials]
      if (t==2){
        mL[mL==FALSE] = TRUE #no forgetting in control decisions
        mR[mR==FALSE] = TRUE #no forgetting in control decisions
      }      
      #recode everything into consistent vs. inconsistent choices
      correctL <- (vL>vR) #whether left option is better
      vC <- vL*(correctL==TRUE)+vR*(correctL==FALSE)
      vI <- vL*(correctL==FALSE)+vR*(correctL==TRUE)
      mC <- mL*(correctL==TRUE)+mR*(correctL==FALSE)
      mI <- mL*(correctL==FALSE)+mR*(correctL==TRUE)
      
      #simulate
      simDat <- sim_ddm_MB(parameters,vC,vI,mC,mI)
      if (t==1){
        simDat_M_noStress <- rbind(simDat_M_noStress, simDat)
      } else {
        simDat_C_noStress <- rbind(simDat_C_noStress, simDat)
      }
    }
  }  
  
  #save posterior predictives
  ppDat_M_Stress_choices[,p] = simDat_M_Stress$choices
  ppDat_M_Stress_rts[,p] = simDat_M_Stress$rts
  ppDat_C_Stress_choices[,p] = simDat_C_Stress$choices
  ppDat_C_Stress_rts[,p] = simDat_C_Stress$rts  
  ppDat_M_noStress_choices[,p] = simDat_M_noStress$choices
  ppDat_M_noStress_rts[,p] = simDat_M_noStress$rts
  ppDat_C_noStress_choices[,p] = simDat_C_noStress$choices
  ppDat_C_noStress_rts[,p] = simDat_C_noStress$rts  

  flush.console()
  msg = sprintf('Done with posterior predictive sample: %d',p)
  print(msg)
}


#model: quantiles
DDMquants <- seq(.1,.9,.2)
quants_c_M_S <- matrix(NA,nrow = nPostPred,ncol = 5)
quants_i_M_S <- matrix(NA,nrow = nPostPred,ncol = 5)
quants_c_M_N <- matrix(NA,nrow = nPostPred,ncol = 5)
quants_i_M_N <- matrix(NA,nrow = nPostPred,ncol = 5)
quants_c_C_S <- matrix(NA,nrow = nPostPred,ncol = 5)
quants_i_C_S <- matrix(NA,nrow = nPostPred,ncol = 5)
quants_c_C_N <- matrix(NA,nrow = nPostPred,ncol = 5)
quants_i_C_N <- matrix(NA,nrow = nPostPred,ncol = 5)
for (p in 1:nPostPred){
  quants_c_M_S[p,] <- quantile(ppDat_M_Stress_rts[ppDat_M_Stress_choices[,p]==2,p],DDMquants,na.rm = T)
  quants_i_M_S[p,] <- quantile(ppDat_M_Stress_rts[ppDat_M_Stress_choices[,p]==1,p],DDMquants,na.rm = T)
  quants_c_M_N[p,] <- quantile(ppDat_M_noStress_rts[ppDat_M_noStress_choices[,p]==2,p],DDMquants,na.rm = T)
  quants_i_M_N[p,] <- quantile(ppDat_M_noStress_rts[ppDat_M_noStress_choices[,p]==1,p],DDMquants,na.rm = T)  
  quants_c_C_S[p,] <- quantile(ppDat_C_Stress_rts[ppDat_C_Stress_choices[,p]==2,p],DDMquants,na.rm = T)
  quants_i_C_S[p,] <- quantile(ppDat_C_Stress_rts[ppDat_C_Stress_choices[,p]==1,p],DDMquants,na.rm = T)
  quants_c_C_N[p,] <- quantile(ppDat_C_noStress_rts[ppDat_C_noStress_choices[,p]==2,p],DDMquants,na.rm = T)
  quants_i_C_N[p,] <- quantile(ppDat_C_noStress_rts[ppDat_C_noStress_choices[,p]==1,p],DDMquants,na.rm = T)  
}
#model: overall choice probability
prob_c_M_S <- colMeans(ppDat_M_Stress_choices==2,na.rm = T)
prob_i_M_S <- colMeans(ppDat_M_Stress_choices==1,na.rm = T)
prob_c_M_N <- colMeans(ppDat_M_noStress_choices==2,na.rm = T)
prob_i_M_N <- colMeans(ppDat_M_noStress_choices==1,na.rm = T)
prob_c_C_S <- colMeans(ppDat_C_Stress_choices==2,na.rm = T)
prob_i_C_S <- colMeans(ppDat_C_Stress_choices==1,na.rm = T)
prob_c_C_N <- colMeans(ppDat_C_noStress_choices==2,na.rm = T)
prob_i_C_N <- colMeans(ppDat_C_noStress_choices==1,na.rm = T)

#get 95% HDIs (highest density interval)
HDI_prob_c_M_S <- quantile(prob_c_M_S,c(.025,.975))
HDI_prob_i_M_S <- quantile(prob_i_M_S,c(.025,.975))
HDI_prob_c_M_N <- quantile(prob_c_M_N,c(.025,.975))
HDI_prob_i_M_N <- quantile(prob_i_M_N,c(.025,.975))
HDI_quants_c_M_S <- matrix(nrow = length(DDMquants),ncol = 2)
HDI_quants_i_M_S <- matrix(nrow = length(DDMquants),ncol = 2)
HDI_quants_c_M_N <- matrix(nrow = length(DDMquants),ncol = 2)
HDI_quants_i_M_N <- matrix(nrow = length(DDMquants),ncol = 2)
HDI_prob_c_C_S <- quantile(prob_c_C_S,c(.025,.975))
HDI_prob_i_C_S <- quantile(prob_i_C_S,c(.025,.975))
HDI_prob_c_C_N <- quantile(prob_c_C_N,c(.025,.975))
HDI_prob_i_C_N <- quantile(prob_i_C_N,c(.025,.975))
HDI_quants_c_C_S <- matrix(nrow = length(DDMquants),ncol = 2)
HDI_quants_i_C_S <- matrix(nrow = length(DDMquants),ncol = 2)
HDI_quants_c_C_N <- matrix(nrow = length(DDMquants),ncol = 2)
HDI_quants_i_C_N <- matrix(nrow = length(DDMquants),ncol = 2)
for (q in 1:length(DDMquants)){
  HDI_quants_c_M_S[q,] <- quantile(quants_c_M_S[,q],c(.025,.975))
  HDI_quants_i_M_S[q,] <- quantile(quants_i_M_S[,q],c(.025,.975))
  HDI_quants_c_M_N[q,] <- quantile(quants_c_M_N[,q],c(.025,.975))
  HDI_quants_i_M_N[q,] <- quantile(quants_i_M_N[,q],c(.025,.975))  
  HDI_quants_c_C_S[q,] <- quantile(quants_c_C_S[,q],c(.025,.975))
  HDI_quants_i_C_S[q,] <- quantile(quants_i_C_S[,q],c(.025,.975))
  HDI_quants_c_C_N[q,] <- quantile(quants_c_C_N[,q],c(.025,.975))
  HDI_quants_i_C_N[q,] <- quantile(quants_i_C_N[,q],c(.025,.975))    
}

#data
cL <- (VL1>VR1) #whether left option is better
cRT1 <- abs(RT1)*(cL==TRUE)+abs(RT1)*(cL==TRUE)
vC <- vL*(correctL==TRUE)+vR*(correctL==FALSE)
vI <- vL*(correctL==FALSE)+vR*(correctL==TRUE)
mC <- mL*(correctL==TRUE)+mR*(correctL==FALSE)
mI <- mL*(correctL==FALSE)+mR*(correctL==TRUE)

cDat1 <- ((VL1>VR1)&(RT1>0))|((VL1<VR1)&(RT1<0))
cDat2 <- ((VL2>VR2)&(RT2>0))|((VL2<VR2)&(RT2<0))
dC_c_M_S <- mean(cDat1[M1==TRUE])
dQuants_c_M_S <- quantile(abs(RT1[(M1==TRUE)&(cDat1==TRUE)]),DDMquants)
dC_i_M_S <- mean(cDat1[M1==TRUE]==FALSE)
dQuants_i_M_S <- quantile(abs(RT1[(M1==TRUE)&(cDat1==FALSE)]),DDMquants)
dC_c_M_N <- mean(cDat2[M2==TRUE])
dQuants_c_M_N <- quantile(abs(RT2[(M2==TRUE)&(cDat2==TRUE)]),DDMquants)
dC_i_M_N <- mean(cDat2[M2==TRUE]==FALSE)
dQuants_i_M_N <- quantile(abs(RT2[(M2==TRUE)&(cDat2==FALSE)]),DDMquants)
dC_c_C_S <- mean(cDat1[M1==FALSE])
dQuants_c_C_S <- quantile(abs(RT1[(M1==FALSE)&(cDat1==TRUE)]),DDMquants)
dC_i_C_S <- mean(cDat1[M1==FALSE]==FALSE)
dQuants_i_C_S <- quantile(abs(RT1[(M1==FALSE)&(cDat1==FALSE)]),DDMquants)
dC_c_C_N <- mean(cDat2[M2==FALSE])
dQuants_c_C_N <- quantile(abs(RT2[(M2==FALSE)&(cDat2==TRUE)]),DDMquants)
dC_i_C_N <- mean(cDat2[M2==FALSE]==FALSE)
dQuants_i_C_N <- quantile(abs(RT2[(M2==FALSE)&(cDat2==FALSE)]),DDMquants)

#plot it (memory-based decisions)
png(filename = "PosteriorPredictives_QuantilePlot_M.png", width = 800, height = 600)
par(fig = c(0, 0.8, 0, 1), mar = c(5, 6, 4, 4), oma = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")

#draw (thin) lines for the "midpoints" of the HDIs
plot(rowMeans(HDI_quants_c_M_S), mean(HDI_prob_c_M_S) * DDMquants, type = 'l', col = "red",
     ylim = c(0, 1), xlim = c(0, 6),
     xlab =  substitute(paste(bold('RT (sec)'))), ylab = substitute(paste(bold('Cumulative Probability'))), lwd = 5, pch = 2, bty = 'L',
     main = 'Memory-based choices', cex.axis = 2, cex.lab = 2, cex.main = 2.5, axes = FALSE)
lines(rowMeans(HDI_quants_i_M_S), mean(HDI_prob_i_M_S) * DDMquants, type = 'l', col = "red", lwd = 5, pch = 2, lty = 2)
lines(rowMeans(HDI_quants_c_M_N), mean(HDI_prob_c_M_N) * DDMquants, type = 'l', col = "blue", lwd = 5, pch = 2, lty = 1)
lines(rowMeans(HDI_quants_i_M_N), mean(HDI_prob_i_M_N) * DDMquants, type = 'l', col = "blue", lwd = 5, pch = 2, lty = 2)

#draw HDIs as errorbars (horizontal for RT, vertical for choice probability)
for (q in 1:length(DDMquants)){
  lines(HDI_quants_c_M_S[q,],rep(mean(HDI_prob_c_M_S) * DDMquants[q],2),type = 'l', col = "red",lwd=5)
  lines(rep(mean(HDI_quants_c_M_S[q,]),2),HDI_prob_c_M_S * DDMquants[q],type = 'l', col = "red",lwd=5)
  lines(HDI_quants_i_M_S[q,],rep(mean(HDI_prob_i_M_S) * DDMquants[q],2),type = 'l', col = "red",lwd=5)
  lines(rep(mean(HDI_quants_i_M_S[q,]),2),HDI_prob_i_M_S * DDMquants[q],type = 'l', col = "red",lwd=5)
  
  lines(HDI_quants_c_M_N[q,],rep(mean(HDI_prob_c_M_N) * DDMquants[q],2),type = 'l', col = "blue",lwd=5)
  lines(rep(mean(HDI_quants_c_M_N[q,]),2),HDI_prob_c_M_N * DDMquants[q],type = 'l', col = "blue",lwd=5)
  lines(HDI_quants_i_M_N[q,],rep(mean(HDI_prob_i_M_N) * DDMquants[q],2),type = 'l', col = "blue",lwd=5)
  lines(rep(mean(HDI_quants_i_M_N[q,]),2),HDI_prob_i_M_N * DDMquants[q],type = 'l', col = "blue",lwd=5)
}

#draw the datapoints as black lines
lines(dQuants_c_M_S, dC_c_M_S * DDMquants, type = 'b', col = "black", lwd = 2, pch = 1, lty = 1)
lines(dQuants_i_M_S, dC_i_M_S * DDMquants, type = 'b', col = "black", lwd = 2, pch = 1, lty = 2)
lines(dQuants_c_M_N, dC_c_M_N * DDMquants, type = 'b', col = "black", lwd = 2, pch = 1, lty = 1)
lines(dQuants_i_M_N, dC_i_M_N * DDMquants, type = 'b', col = "black", lwd = 2, pch = 1, lty = 2)

legend("topright", legend = c('Stress/consistent', 'No stress/consistent', 'Stress/consistent', 'No stress/consistent','Data/consistent','Data/inconsistent'),
       col = c("red", "blue", "red", "blue","black","black"),
       lty = c(1, 1, 2, 2,1,2), lwd = c(4, 4, 4, 4,4,4), bty = 'n', cex = 2)

# Add y-axis with upright labels
axis(1, lwd=3, cex.axis=1.5)
axis(2, lwd=3, cex.axis=1.5, las=1)

dev.off()


#conventional choices
png(filename = "PosteriorPredictives_QuantilePlot_C.png", width = 800, height = 600)
par(fig = c(0, 0.8, 0, 1), mar = c(5, 6, 4, 4), oma = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")

#draw (thin) lines for the "midpoints" of the HDIs
plot(rowMeans(HDI_quants_c_C_S), mean(HDI_prob_c_C_S) * DDMquants, type = 'l', col = "magenta",
     ylim = c(0, 1), xlim = c(0, 6),
     xlab =  substitute(paste(bold('RT (sec)'))), ylab = substitute(paste(bold('Cumulative Probability'))), lwd = 5, pch = 2, bty = 'L',
     main = 'Conventional choices', cex.axis = 2, cex.lab = 2, cex.main = 2.5, axes = FALSE)
lines(rowMeans(HDI_quants_i_C_S), mean(HDI_prob_i_C_S) * DDMquants, type = 'l', col = "magenta", lwd = 5, pch = 2, lty = 2)
lines(rowMeans(HDI_quants_c_C_N), mean(HDI_prob_c_C_N) * DDMquants, type = 'l', col = "cyan", lwd = 5, pch = 2, lty = 1)
lines(rowMeans(HDI_quants_i_C_N), mean(HDI_prob_i_C_N) * DDMquants, type = 'l', col = "cyan", lwd = 5, pch = 2, lty = 2)

#draw HDIs as errorbars (horizontal for RT, vertical for choice probability)
for (q in 1:length(DDMquants)){
  lines(HDI_quants_c_C_S[q,],rep(mean(HDI_prob_c_C_S) * DDMquants[q],2),type = 'l', col = "magenta",lwd=5)
  lines(rep(mean(HDI_quants_c_C_S[q,]),2),HDI_prob_c_C_S * DDMquants[q],type = 'l', col = "magenta",lwd=5)
  lines(HDI_quants_i_C_S[q,],rep(mean(HDI_prob_i_C_S) * DDMquants[q],2),type = 'l', col = "magenta",lwd=5)
  lines(rep(mean(HDI_quants_i_C_S[q,]),2),HDI_prob_i_C_S * DDMquants[q],type = 'l', col = "magenta",lwd=5)
  
  lines(HDI_quants_c_C_N[q,],rep(mean(HDI_prob_c_C_N) * DDMquants[q],2),type = 'l', col = "cyan",lwd=5)
  lines(rep(mean(HDI_quants_c_C_N[q,]),2),HDI_prob_c_C_N * DDMquants[q],type = 'l', col = "cyan",lwd=5)
  lines(HDI_quants_i_C_N[q,],rep(mean(HDI_prob_i_C_N) * DDMquants[q],2),type = 'l', col = "cyan",lwd=5)
  lines(rep(mean(HDI_quants_i_C_N[q,]),2),HDI_prob_i_C_N * DDMquants[q],type = 'l', col = "cyan",lwd=5)  
}

#draw the datapoints as black lines
lines(dQuants_c_C_S, dC_c_C_S * DDMquants, type = 'b', col = "black", lwd = 2, pch = 1, lty = 1)
lines(dQuants_i_C_S, dC_i_C_S * DDMquants, type = 'b', col = "black", lwd = 2, pch = 1, lty = 2)
lines(dQuants_c_C_N, dC_c_C_N * DDMquants, type = 'b', col = "black", lwd = 2, pch = 1, lty = 1)
lines(dQuants_i_C_N, dC_i_C_N * DDMquants, type = 'b', col = "black", lwd = 2, pch = 1, lty = 2)


legend("topright", legend = c('Stress/consistent', 'No stress/consistent', 'Stress/consistent', 'No stress/consistent','Data/consistent','Data/inconsistent'),
       col = c("magenta", "cyan", "magenta", "cyan","black","black"),
       lty = c(1, 1, 2, 2,1,2), lwd = c(4, 4, 4, 4,4,4), bty = 'n', cex = 2)

# Add y-axis with upright labels
axis(1, lwd=3, cex.axis=1.5)
axis(2, lwd=3, cex.axis=1.5, las=1)

dev.off()