#This is my attempt to pull together a more organized code document specifically related to cell type index analyses using DLPFC Affymetrix data from the Pritzker foundation
##Megan Hagenauer, October 27 2016

#**********************************************


#Using a scatterplot with best fit line to visually examine the relationships between the the publication-specific cell type indices and principal components of variation (SubjectPCA):

for (i in 1:length(SubjectPCA[1,])){
for(j in 1:length(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[,1])){
png(paste(paste(colnames(SubjectPCA)[i], "vs", row.names(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean)[j], sep="  "), "png", sep="."))	
plot(SubjectPCA[,i]~CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[j,], main=paste(colnames(SubjectPCA)[i], "vs", row.names(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean)[j], sep="  "), xlab=row.names(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean)[j], ylab=colnames(SubjectPCA)[i])
RegressionLine<-lm(SubjectPCA[,i]~CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[j,])
abline(RegressionLine, col=2)
mtext(paste("p-value=", round(summary.lm(RegressionLine)$coefficients[8], digits=8), ", r-squared = ", round(summary.lm(RegressionLine)[8][[1]], digits=2)), line=0)
dev.off()		
}		
}


FlaggedRelationshipsBetweenIVandPCA<-file("15 Flagged Relationships Between Subject Variables and PCA.txt")
out<-c(

#Using linear regression to examine the statistical relationships between the cell type indices and PCA:
capture.output(
for (i in 1:length(SubjectPCA[1,])){
for(j in 1:length(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[,1])){
if(summary.lm(lm(SubjectPCA[,i]~CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[j,]))$coefficient[8]<0.05){
print(paste(colnames(SubjectPCA)[i], "vs", row.names(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean)[j], "p-value=", summary.lm(lm(SubjectPCA[,i]~CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[j,]))$coefficient[8], "r-squared=",  round(summary.lm(lm(SubjectPCA[,i]~CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[j,]))[8][[1]], digits=3), sep="  "))}else{}
}		
}
)

)
cat(out, file="15 Flagged Relationships Between Subject Variables and PCA.txt", sep="\n", append=TRUE)
close(FlaggedRelationshipsBetweenIVandPCA)
rm(out)


