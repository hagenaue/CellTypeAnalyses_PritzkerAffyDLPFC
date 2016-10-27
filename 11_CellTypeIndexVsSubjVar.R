#This is my attempt to pull together a more organized code document specifically related to cell type index analyses using DLPFC Affymetrix data from the Pritzker foundation
##Megan Hagenauer, October 27 2016

#**********************************************



#Using a scatterplot with best fit line to visually examine the relationships between the continuous subject variables and cell type indices:

for (i in 1:length(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[,1])){
for(j in 1:length(SubjectContinuousVariables[1,])){
png(paste(paste(row.names(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean)[i], "vs", colnames(SubjectContinuousVariables)[j], sep="  "), "png", sep="."))	
plot(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[i,]~SubjectContinuousVariables[,j], main=paste(row.names(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean)[i], "vs", colnames(SubjectContinuousVariables)[j], sep="  "), xlab=colnames(SubjectContinuousVariables)[j], ylab=row.names(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean)[i])
RegressionLine<-lm(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[i,]~SubjectContinuousVariables[,j])
abline(RegressionLine, col=2)
mtext(paste("p-value=", round(summary.lm(RegressionLine)$coefficients[8], digits=4), ", r-squared=",  round(summary.lm(lm(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[i,]~SubjectContinuousVariables[,j]))[8][[1]], digits=3)))
dev.off()		
}		
}



#Using boxplots to visually examine the relationships between the cell type indices and categorical subject variables:
for (i in 1:length(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[,1])){
for(j in 1:length(SubjectFactorVariables[1,])){
png(paste(paste(row.names(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "), "png", sep="."))	
boxplot(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[i,]~SubjectFactorVariables[,j], main=paste(row.names(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "), xlab=colnames(SubjectFactorVariables)[j], ylab=row.names(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean)[i])
mtext(paste("p-value=", round(summary(aov(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[i,]~SubjectFactorVariables[,j]))[[1]][["Pr(>F)"]][1], digits=4), ", r-squared=",  round(summary.lm(lm(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[i,]~SubjectFactorVariables[,j]))[8][[1]], digits=3)))
dev.off()		
}		
}

#Flagging variables that are collinear with other subject variables:
FlaggedRelationshipsBetweenIV<-file("14 Flagged Relationships Between Subject Variables.txt")
out<-c(

#Using linear regression to examine the statistical relationships between the cell type indices and continuous subject variables:
capture.output(
for (i in 1:length(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[,1])){
for(j in 1:length(SubjectContinuousVariables[1,])){
if(summary.lm(lm(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[i,]~SubjectContinuousVariables[,j]))$coefficient[8]<0.05){
print(paste(row.names(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean)[i], "vs", colnames(SubjectContinuousVariables)[j], "p-value=", summary.lm(lm(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[i,]~SubjectContinuousVariables[,j]))$coefficient[8], ", r-squared=",  round(summary.lm(lm(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[i,]~SubjectContinuousVariables[,j]))[8][[1]], digits=3), sep="  "))}else{}
}		
}
),

#Using anova to examine the statistical relationships between the cell type indices and categorical subject variables:
capture.output(
for (i in 1:length(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[,1])){
for(j in 1:length(SubjectFactorVariables[1,])){
if(summary(aov(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[i,]~SubjectFactorVariables[,j]))[[1]][["Pr(>F)"]][1]<0.05){
print(paste(row.names(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean)[i], "vs", colnames(SubjectFactorVariables)[j], "p-value=", summary(aov(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[i,]~SubjectFactorVariables[,j]))[[1]][["Pr(>F)"]][1], ", r-squared=",  round(summary.lm(lm(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[i,]~SubjectFactorVariables[,j]))[8][[1]], digits=3), sep="  "))	
}else{}		
}		
}
)


)
cat(out, file="14 Flagged Relationships Between Subject Variables.txt", sep="\n", append=TRUE)
close(FlaggedRelationshipsBetweenIV)
rm(out)

