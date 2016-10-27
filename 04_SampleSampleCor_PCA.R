#This is my attempt to pull together a more organized code document specifically related to cell type index analyses using DLPFC Affymetrix data from the Pritzker foundation
##Megan Hagenauer, October 27 2016

#**********************************************

#Run principal components analysis (PCA) to determine which major gradients of sample-sample correlations exist in the data (i.e. who is similar to whom):

# # #Visualize the sample-sample correlations using a heatmap:
png("10 Sample Sample Correlations Heatmap.png")
image(cor(SignalSortedNoNA3), main="Visualizing the correlations between entire samples (by index#)", xlab="Red=Less correlated, Light yellow=Highly correlated")
dev.off()
# #Note that the heatmap can be tailored to focus on a certain level of correlation by using the command zlim=c(lower limit, upper limit)

# #Visualize the sample-sample correlations using a boxplot:
png("10 Boxplot Sample Sample Correlations.png")
boxplot(data.frame(cor(SignalSortedNoNA3)), cex=0.25, las=3, par(cex.axis=0.75))
Median10thQuantile<-median(apply((cor(SignalSortedNoNA3)), 1, quantile, 0.1))
abline(a=Median10thQuantile, b=0, col=2)
dev.off()
#Are these samples just really really uncorrelated due to the signal data being normalized?  

plot(SignalSortedNoNA3[,1]~SignalSortedNoNA3[,2])
#Yep, that seems to be the case - the normalized signal data has made it so that the sample correlations are sometimes even negative!


# #Run principal components analysis (PCA) to determine which major gradients of sample-sample correlations exist in the data (i.e. who is similar to whom):
pcaNormFilterednoOutliers<-prcomp(t(SignalSortedNoNA3))
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
