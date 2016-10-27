#This is my attempt to pull together a more organized code document specifically related to cell type index analyses using DLPFC Affymetrix data from the Pritzker foundation
##Megan Hagenauer, October 27 2016


#**********************************************************************************************************

#Characterizing subjects and looking for relationships between subject variables:


SubjectFactorVariables<-cbind(DiagnosisNoNA3, GenderNoNA3, SuicideNoNA3)

SubjectContinuousVariables<-cbind(BrainPHNoNA3, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, HoursFinalNoNA3, AgeNoNA3, cosTODNoNA3, sinTODNoNA3)


for (i in 1:length(SubjectContinuousVariables[1,])){
png(paste(paste("Histogram of", colnames(SubjectContinuousVariables)[i], sep="  "), "png", sep="."))	
hist(SubjectContinuousVariables[, i], col=i+1)
dev.off()		
}





#Using a scatterplot with best fit line to visually examine the relationships between the continuous subject variables:
for (i in 1:length(SubjectContinuousVariables[1,])){
for(j in 1:length(SubjectContinuousVariables[1,])){
png(paste("14", paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectContinuousVariables)[j], sep="  "), "png", sep="."))	
plot(SubjectContinuousVariables[,i]~SubjectContinuousVariables[,j], main=paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectContinuousVariables)[j], sep="  "), xlab=colnames(SubjectContinuousVariables)[j], ylab=colnames(SubjectContinuousVariables)[i])
RegressionLine<-lm(SubjectContinuousVariables[,i]~SubjectContinuousVariables[,j])
abline(RegressionLine, col=2)
mtext(paste("p-value=", round(summary.lm(RegressionLine)$coefficients[8], digits=4)))
dev.off()		
}		
}



#Using boxplots to visually examine the relationships between the continuous subject variables and categorical subject variables:
for (i in 1:length(SubjectContinuousVariables[1,])){
for(j in 1:length(SubjectFactorVariables[1,])){
png(paste("14", paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "), "png", sep="."))	
boxplot(SubjectContinuousVariables[,i]~SubjectFactorVariables[,j], main=paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "), xlab=colnames(SubjectFactorVariables)[j], ylab=colnames(SubjectContinuousVariables)[i])
mtext(paste("p-value=", round(summary(aov(SubjectContinuousVariables[,i]~SubjectFactorVariables[,j]))[[1]][["Pr(>F)"]][1], digits=4)))
dev.off()		
}		
}



#Creating a text file of contingency tables to visually examine the relationships between categorical subject variables:
CrossTabsIV<-file("14 Cross Tabs Between Subject Factors.txt")
out<-c(
capture.output(

summary(DiagnosisNoNA3),
summary(GenderNoNA3),

for (i in 1:length(SubjectFactorVariables[1,])){
for(j in 1:length(SubjectFactorVariables[1,])){
ContingencyTable<-table(SubjectFactorVariables[,i],SubjectFactorVariables[,j])
print(paste(colnames(SubjectFactorVariables)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "))
print(ContingencyTable)
print(paste("p-value=", chisq.test(ContingencyTable)$p.value))	
}		
}
)
)
cat(out, file="14 Cross Tabs Between Subject Factors.txt", sep="\n", append=TRUE)
close(CrossTabsIV)
rm(out)





StatisticalRelationshipsIV<-file("14 Statistical Relationships between Subject Variables.txt")
out<-c(

capture.output(
#Calculating the variance inflation factor (vif) to determine which subject variables are highly related to other subject variables in the data set. Most important, of course, is whether any of the subject variables strongly correlate with Diagnosis. 
vif(lm(SignalSortedNoNA3[1,]~BrainPHNoNA3 + AgonalFactorNoNA3 +HoursFinalCorrectedNoNA3+DiagnosisNoNA3+GenderNoNA3+AgeNoNA3))

),

#Using linear regression to examine the statistical relationships between the continuous subject variables:

capture.output(
for (i in 1:length(SubjectContinuousVariables[1,])){
for(j in 1:length(SubjectContinuousVariables[1,])){
print(paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectContinuousVariables)[j], sep="  "))
print(summary.lm(lm(SubjectContinuousVariables[,i]~SubjectContinuousVariables[,j])))
}		
}
),
#Using anova to examine the statistical relationships between the continuous subject variables and categorical subject variables:

capture.output(
for (i in 1:length(SubjectContinuousVariables[1,])){
for(j in 1:length(SubjectFactorVariables[1,])){
print(paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "))
print(summary(aov(SubjectContinuousVariables[,i]~SubjectFactorVariables[,j])))		
}		
}
),

#Using chi-square to examine the statistical relationships between the categorical subject variables:

capture.output(
for (i in 1:length(SubjectFactorVariables[1,])){
for(j in 1:length(SubjectFactorVariables[1,])){
print(paste(colnames(SubjectFactorVariables)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "))
print(chisq.test(ContingencyTable))		
}		
}
)

)
cat(out, file="14 Statistical Relationships between Subject Variables.txt", sep="\n", append=TRUE)
close(StatisticalRelationshipsIV)
rm(out)


#Flagging variables that are collinear with other subject variables:
FlaggedRelationshipsBetweenIV<-file("14 Flagged Relationships Between Subject Variables.txt")
out<-c(

#Using linear regression to examine the statistical relationships between the continuous subject variables:
capture.output(
for (i in 1:length(SubjectContinuousVariables[1,])){
for(j in 1:length(SubjectContinuousVariables[1,])){
if(summary.lm(lm(SubjectContinuousVariables[,i]~SubjectContinuousVariables[,j]))$coefficient[8]<0.05){
print(paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectContinuousVariables)[j], "p-value=", summary.lm(lm(SubjectContinuousVariables[,i]~SubjectContinuousVariables[,j]))$coefficient[8], sep="  "))}else{}
}		
}
),

#Using anova to examine the statistical relationships between the continuous subject variables and categorical subject variables:
capture.output(
for (i in 1:length(SubjectContinuousVariables[1,])){
for(j in 1:length(SubjectFactorVariables[1,])){
if(summary(aov(SubjectContinuousVariables[,i]~SubjectFactorVariables[,j]))[[1]][["Pr(>F)"]][1]<0.05){
print(paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectFactorVariables)[j], "p-value=", summary(aov(SubjectContinuousVariables[,i]~SubjectFactorVariables[,j]))[[1]][["Pr(>F)"]][1], sep="  "))	
}else{}		
}		
}
),

#Using chi-square to examine the statistical relationships between the categorical subject variables:
capture.output(
for (i in 1:length(SubjectFactorVariables[1,])){
for(j in 1:length(SubjectFactorVariables[1,])){
ContingencyTable<-table(SubjectFactorVariables[,i], SubjectFactorVariables[,j])
if(chisq.test(ContingencyTable)$p.value<0.05){
print(paste(colnames(SubjectFactorVariables)[i], "vs", colnames(SubjectFactorVariables)[j], "p-value=", chisq.test(ContingencyTable)$p.value, sep="  "))	
}else{}		
}		
}
)
)
cat(out, file="14 Flagged Relationships Between Subject Variables.txt", sep="\n", append=TRUE)
close(FlaggedRelationshipsBetweenIV)
rm(out)


