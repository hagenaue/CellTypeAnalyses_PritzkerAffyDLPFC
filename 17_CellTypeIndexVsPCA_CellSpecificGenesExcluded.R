#This is my attempt to pull together a more organized code document specifically related to cell type index analyses using DLPFC Affymetrix data from the Pritzker foundation
##Megan Hagenauer, October 27 2016

#**********************************************

#Re-Running the same PCA code to examine the relationships between principal components of variation and the consolidated cell type indices, but using a PCA analysis that excludes all cell type specific genes:

dim(SignalSortedNoNA3DF)
SignalSortedNoNA3DF[c(1:3), c(1:3)] #2nd column is gene symbol
dim(CellTypeSpecificGenes_Master3)
colnames(CellTypeSpecificGenes_Master3)#Human Gene symbol is column 4

sum(SignalSortedNoNA3DF[,2]%in%CellTypeSpecificGenes_Master3[,4])
#[1] 1683

sum((SignalSortedNoNA3DF[,2]%in%CellTypeSpecificGenes_Master3[,4])==F)
[1] 10296

SignalSortedNoNA3NoCellTypeGenes_SD<-apply(SignalSortedNoNA3DF[(SignalSortedNoNA3DF[,2]%in%CellTypeSpecificGenes_Master3[,4])==F,-c(1:2)], 1, sd)
mean(SignalSortedNoNA3NoCellTypeGenes_SD)
[1] 0.1784711

SignalSortedNoNA3CellTypeGenes_SD<-apply(SignalSortedNoNA3DF[(SignalSortedNoNA3DF[,2]%in%CellTypeSpecificGenes_Master3[,4])==T,-c(1:2)], 1, sd)
mean(SignalSortedNoNA3CellTypeGenes_SD)
[1] 0.2102225

png("StDev_CellTypeGenesVsOther.png")
boxplot(c(SignalSortedNoNA3CellTypeGenes_SD, SignalSortedNoNA3NoCellTypeGenes_SD)~c(rep("Cell Type Specific Genes", length(SignalSortedNoNA3CellTypeGenes_SD)), rep("Other", length(SignalSortedNoNA3NoCellTypeGenes_SD))), col=3, ylab="Standard Deviation of Normalized Probe Signal", cex.lab=1.4, cex=1.4)
mtext(paste("P-value:", format(t.test(c(SignalSortedNoNA3CellTypeGenes_SD, SignalSortedNoNA3NoCellTypeGenes_SD)~c(rep("Cell Type Specific Genes", length(SignalSortedNoNA3CellTypeGenes_SD)), rep("Other", length(SignalSortedNoNA3NoCellTypeGenes_SD))))$p.value, digits=3, scientific=T), sep=" "), cex=1.4)
dev.off()


SignalSortedNoNA3NoCellTypeGenes<-as.matrix(SignalSortedNoNA3DF[(SignalSortedNoNA3DF[,2]%in%CellTypeSpecificGenes_Master3[,4])==F,-c(1:2)])
dim(SignalSortedNoNA3NoCellTypeGenes)
[1] 10296   157
is.numeric(SignalSortedNoNA3NoCellTypeGenes)

# #Run principal components analysis (PCA) to determine which major gradients of sample-sample correlations exist in the data (i.e. who is similar to whom):
pcaNormFilterednoOutliers<-prcomp(t(SignalSortedNoNA3NoCellTypeGenes))
tmp<-pcaNormFilterednoOutliers$x[,1:4]
write.table(tmp, "PCA_1_4.txt", sep="\t")


PCeigenvectors<-pcaNormFilterednoOutliers$rotation[ ,c(1:4)]
PCeigenvectors2<-cbind(PCeigenvectors, GeneNames)
write.csv(PCeigenvectors2, "PCeigenvectors.csv")

PC1noOutliers<-pcaNormFilterednoOutliers$x[,1]
PC2noOutliers<-pcaNormFilterednoOutliers$x[,2]

PC3noOutliers<-pcaNormFilterednoOutliers$x[,3]
PC4noOutliers<-pcaNormFilterednoOutliers$x[,4]

# #Output a scree plot for the PCA (no outliers):
png("10 PCA Scree Plot1.png")
plot(summary(pcaNormFilterednoOutliers)$importance[2,]~(c(1:length(summary(pcaNormFilterednoOutliers)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC #", ylab="Proportion of Variance Explained", col=2)
dev.off()

png("10 PCA Scree Plot2.png")
plot(summary(pcaNormFilterednoOutliers)$importance[3,]~(c(1:length(summary(pcaNormFilterednoOutliers)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC #", ylab="Cumulative Proportion of Variance Explained", col=3)
dev.off()

# #Output a scatterplot illustrating the relationship between Principal components 1 & 2 (PC1 & PC2):
png("10 PC1 vs PC2.png")
plot(PC1noOutliers~PC2noOutliers, main="Principal Components Analysis of Normalized Filtered Data No Outliers")
points(PC1noOutliers[DiagnosisNoNA3=="Control"]~PC2noOutliers[DiagnosisNoNA3=="Control"], col=3)
points(PC1noOutliers[DiagnosisNoNA3=="MD"]~PC2noOutliers[DiagnosisNoNA3=="MD"], col=2)
points(PC1noOutliers[DiagnosisNoNA3=="BP"]~PC2noOutliers[DiagnosisNoNA3=="BP"], col=4)
points(PC1noOutliers[DiagnosisNoNA3=="Schiz"]~PC2noOutliers[DiagnosisNoNA3=="Schiz"], col=5)
legend(min(PC1noOutliers), max(PC2noOutliers)+10, c("Control", "MDD", "BP", "Schiz"), text.col=c(3, 2, 4, 5), pch=19, col=c(3, 2, 4, 5))
dev.off()

# #Output a scatterplot illustrating the relationship between Principal components 3 & 4 (PC3 & PC4):
png("10 PC3 vs PC4.png")
plot(PC3noOutliers~PC4noOutliers, main="Principal Components Analysis of Normalized Filtered Data No Outliers")
points(PC3noOutliers[DiagnosisNoNA3=="Control"]~PC4noOutliers[DiagnosisNoNA3=="Control"], col=3)
points(PC3noOutliers[DiagnosisNoNA3=="MD"]~PC4noOutliers[DiagnosisNoNA3=="MD"], col=2)
points(PC3noOutliers[DiagnosisNoNA3=="BP"]~PC4noOutliers[DiagnosisNoNA3=="BP"], col=4)
points(PC3noOutliers[DiagnosisNoNA3=="Schiz"]~PC4noOutliers[DiagnosisNoNA3=="Schiz"], col=5)
legend(min(PC3noOutliers), max(PC4noOutliers)+10, c("Control", "MDD", "BP", "Schiz"), text.col=c(3, 2, 4, 5), pch=19, col=c(3, 2, 4, 5))
dev.off()

#****Looking at the relationships between the consolidated cell type indices and principal components:


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

#Making Nicer Figures:

setwd("~/Documents/Affy/NoPC1correct/DLPFC circadian/MoreFigsForCellTypePaper/PCA wo CellTypeGenes")

temp<-(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[1,]+CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[2,])/2


pdf("PC1vsAstrocyteEndothelial.pdf", width=4, height=4, pointsize=10)
plot(SubjectPCA[,1]~temp, xlab="Combined Astrocyte & Endothelial Indices", ylab="PC1", col=1, font.lab=2, lwd=1, cex.lab=1.3, cex.axis=1)
RegressionLine<-lm(SubjectPCA[,1]~temp)
abline(RegressionLine, col="red", lwd=4)
mtext(paste("p-value = ", print(formatC(signif(summary.lm(RegressionLine)$coefficients[8],digits=3), digits=3, flag="#")), sep=""), line=0.5, cex=1.5, font=2)
mtext(paste("r-squared = ", round(summary.lm(RegressionLine)[8][[1]], digits=2), sep=""), line=2, cex=1.5, font=2)
dev.off()


pdf("PC1vsNeuron_All.pdf", width=4, height=4, pointsize=10)
plot(SubjectPCA[,1]~CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[5,], xlab="Neuron (All) Index", ylab="PC1", col=1, font.lab=2, lwd=1, cex.lab=1.3, cex.axis=1)
RegressionLine<-lm(SubjectPCA[,1]~CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[5,])
abline(RegressionLine, col="dodgerblue3", lwd=4)
mtext(paste("p-value = ", print(formatC(signif(summary.lm(RegressionLine)$coefficients[8],digits=3), digits=3, flag="#")), sep=""), line=0.5, cex=1.5, font=2)
mtext(paste("r-squared = ", round(summary.lm(RegressionLine)[8][[1]], digits=2), sep=""), line=2, cex=1.5, font=2)
dev.off()


pdf("PC2vsNeuron_Projection.pdf", width=4, height=4, pointsize=10)
plot(SubjectPCA[,2]~CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[7,], xlab="Projection Neuron Index", ylab="PC2", col=1, font.lab=2, lwd=1, cex.lab=1.3, cex.axis=1)
RegressionLine<-lm(SubjectPCA[,2]~CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[7,])
abline(RegressionLine, col="dodgerblue3", lwd=4)
mtext(paste("p-value = ", print(formatC(signif(summary.lm(RegressionLine)$coefficients[8],digits=3), digits=3, flag="#")), sep=""), line=0.5, cex=1.5, font=2)
mtext(paste("r-squared = ", round(summary.lm(RegressionLine)[8][[1]], digits=2), sep=""), line=2, cex=1.5, font=2)
dev.off()

pdf("PC2vsOligodendrocyte.pdf", width=4, height=4, pointsize=10)
plot(SubjectPCA[,2]~CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[8,], xlab="Oligodendrocyte Index", ylab="PC2", col=1, font.lab=2, lwd=1, cex.lab=1.3, cex.axis=1)
RegressionLine<-lm(SubjectPCA[,2]~CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[8,])
abline(RegressionLine, col="red", lwd=4)
mtext(paste("p-value = ", print(formatC(signif(summary.lm(RegressionLine)$coefficients[8],digits=3), digits=3, flag="#")), sep=""), line=0.5, cex=1.5, font=2)
mtext(paste("r-squared = ", round(summary.lm(RegressionLine)[8][[1]], digits=2), sep=""), line=2, cex=1.5, font=2)
dev.off()

