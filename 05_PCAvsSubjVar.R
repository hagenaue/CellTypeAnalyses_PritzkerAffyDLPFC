#This is my attempt to pull together a more organized code document specifically related to cell type index analyses using DLPFC Affymetrix data from the Pritzker foundation
##Megan Hagenauer, October 27 2016

#**********************************************

#Looking at the relationships between independent variables and principal components:

##1. To examine both visual and statistical relationships between the independent variables (with the outlier subjects removed, if there were any!) and the principal components of variation in the samples, just highlight this code and click Ctrl+R:


#***CODE TO RUN FOR STEP 15****

SubjectPCA<-cbind(PC1noOutliers, PC2noOutliers, PC3noOutliers, PC4noOutliers)

PCAoutput<-cbind(SubjectFactorVariables, SubjectContinuousVariables, SubjectPCA)
write.csv(PCAoutput, "PCAoutput.csv")

#Using a scatterplot with best fit line to visually examine the relationships between the continuous subject variables and SubjectPCA:
for (i in 1:length(SubjectPCA[1,])){
for(j in 1:length(SubjectContinuousVariables[1,])){
png(paste("15", paste(colnames(SubjectPCA)[i], "vs", colnames(SubjectContinuousVariables)[j], sep="  "), "png", sep="."))	
plot(SubjectPCA[,i]~SubjectContinuousVariables[,j], main=paste(colnames(SubjectPCA)[i], "vs", colnames(SubjectContinuousVariables)[j], sep="  "), xlab=colnames(SubjectContinuousVariables)[j], ylab=colnames(SubjectPCA)[i])
RegressionLine<-lm(SubjectPCA[,i]~SubjectContinuousVariables[,j])
abline(RegressionLine, col=2)
mtext(paste("p-value=", round(summary.lm(RegressionLine)$coefficients[8], digits=4)))
dev.off()		
}		
}



#Using boxplots to visually examine the relationships between the PCA and categorical subject variables:
for (i in 1:length(SubjectPCA[1,])){
for(j in 1:length(SubjectFactorVariables[1,])){
png(paste("15", paste(colnames(SubjectPCA)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "), "png", sep="."))	
boxplot(SubjectPCA[,i]~SubjectFactorVariables[,j], main=paste(colnames(SubjectPCA)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "), xlab=colnames(SubjectFactorVariables)[j], ylab=colnames(SubjectPCA)[i])
mtext(paste("p-value=", round(summary(aov(SubjectPCA[,i]~SubjectFactorVariables[,j]))[[1]][["Pr(>F)"]][1], digits=4)))
dev.off()		
}		
}




#Outputting a text file containing the statistical relationships between all of the subject variables and PCA:

StatisticalRelationshipsIVvsPCA<-file("15 Statistical Relationships between Subject Variables and PCA.txt")
out<-c(

capture.output(
#Calculating the variance inflation factor (vif) to determine which subject variables are highly related to other subject variables in the data set. Most important, of course, is whether any of the subject variables strongly correlate with Diagnosis. Note that "Location on Chip" has been removed as a variable because it is partially redundant with gender. 
summary.lm(lm(PC1noOutliers~BrainPHNoNA3 + AgonalFactorNoNA3 +HoursFinalCorrectedNoNA3+DiagnosisNoNA3+GenderNoNA3+AgeNoNA3+SuicideNoNA3))
),

capture.output(
summary.lm(lm(PC2noOutliers~BrainPHNoNA3 + AgonalFactorNoNA3 +HoursFinalCorrectedNoNA3+DiagnosisNoNA3+GenderNoNA3+AgeNoNA3+SuicideNoNA3))
),

capture.output(
summary.lm(lm(PC3noOutliers~BrainPHNoNA3 + AgonalFactorNoNA3 +HoursFinalCorrectedNoNA3+DiagnosisNoNA3+GenderNoNA3+AgeNoNA3+SuicideNoNA3))
),

capture.output(
summary.lm(lm(PC4noOutliers~BrainPHNoNA3 + AgonalFactorNoNA3 +HoursFinalCorrectedNoNA3+DiagnosisNoNA3+GenderNoNA3+AgeNoNA3+SuicideNoNA3))
),


#Using linear regression to examine the statistical relationships between PCA and the continuous subject variables:

capture.output(
for (i in 1:length(SubjectPCA[1,])){
for(j in 1:length(SubjectContinuousVariables[1,])){
print(paste(colnames(SubjectPCA)[i], "vs", colnames(SubjectContinuousVariables)[j], sep="  "))
print(summary.lm(lm(SubjectPCA[,i]~SubjectContinuousVariables[,j])))
}		
}
),

#Using anova to examine the statistical relationships between PCA and categorical subject variables:

capture.output(
for (i in 1:length(SubjectPCA[1,])){
for(j in 1:length(SubjectFactorVariables[1,])){
print(paste(colnames(SubjectPCA)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "))
print(summary(aov(SubjectPCA[,i]~SubjectFactorVariables[,j])))		
}		
}
)

)
cat(out, file="15 Statistical Relationships between Subject Variables and PCA.txt", sep="\n", append=TRUE)
close(StatisticalRelationshipsIVvsPCA)
rm(out)



#Flagging variables that are collinear with other subject variables:
FlaggedRelationshipsBetweenIVandPCA<-file("15 Flagged Relationships Between Subject Variables and PCA.txt")
out<-c(

#Using linear regression to examine the statistical relationships between the continuous subject variables:
capture.output(
for (i in 1:length(SubjectPCA[1,])){
for(j in 1:length(SubjectContinuousVariables[1,])){
if(summary.lm(lm(SubjectPCA[,i]~SubjectContinuousVariables[,j]))$coefficient[8]<0.05){
print(paste(colnames(SubjectPCA)[i], "vs", colnames(SubjectContinuousVariables)[j], "p-value=", summary.lm(lm(SubjectPCA[,i]~SubjectContinuousVariables[,j]))$coefficient[8], sep="  "))}else{}
}		
}
),

#Using anova to examine the statistical relationships between the continuous subject variables and categorical subject variables:
capture.output(
for (i in 1:length(SubjectPCA[1,])){
for(j in 1:length(SubjectFactorVariables[1,])){
if(summary(aov(SubjectPCA[,i]~SubjectFactorVariables[,j]))[[1]][["Pr(>F)"]][1]<0.05){
print(paste(colnames(SubjectPCA)[i], "vs", colnames(SubjectFactorVariables)[j], "p-value=", summary(aov(SubjectPCA[,i]~SubjectFactorVariables[,j]))[[1]][["Pr(>F)"]][1], sep="  "))	
}else{}		
}		
}
)
)
cat(out, file="15 Flagged Relationships Between Subject Variables and PCA.txt", sep="\n", append=TRUE)
close(FlaggedRelationshipsBetweenIVandPCA)
rm(out)
