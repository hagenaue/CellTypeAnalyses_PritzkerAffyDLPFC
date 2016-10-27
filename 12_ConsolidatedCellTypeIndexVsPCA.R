#This is my attempt to pull together a more organized code document specifically related to cell type index analyses using DLPFC Affymetrix data from the Pritzker foundation
##Megan Hagenauer, October 27 2016

#**********************************************

#Running the same code to examine the relationships between principal components of variation and cell type balance, but for primary category indices with no overlapping genes between primary categories:

#Using a scatterplot with best fit line to visually examine the relationships between the consolidated cell type indices and SubjectPCA:
for (i in 1:length(SubjectPCA[1,])){
for(j in 1:length(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[,1])){
png(paste(paste(colnames(SubjectPCA)[i], "vs", row.names(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean)[j], sep="  "), "png", sep="."))	
plot(SubjectPCA[,i]~CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[j,], main=paste(colnames(SubjectPCA)[i], "vs", row.names(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean)[j], sep="  "), xlab=row.names(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean)[j], ylab=colnames(SubjectPCA)[i])
RegressionLine<-lm(SubjectPCA[,i]~CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[j,])
abline(RegressionLine, col=2)
mtext(paste("p-value=", round(summary.lm(RegressionLine)$coefficients[8], digits=8), ", r-squared = ", round(summary.lm(RegressionLine)[8][[1]], digits=2)), line=0)
dev.off()		
}		
}


FlaggedRelationshipsBetweenIVandPCA<-file("15 Flagged Relationships Between Subject Variables and PCA.txt")
out<-c(

#Using linear regression to examine the statistical relationships between the the consolidated cell type indices and PCA:
capture.output(
for (i in 1:length(SubjectPCA[1,])){
for(j in 1:length(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[,1])){
if(summary.lm(lm(SubjectPCA[,i]~CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[j,]))$coefficient[8]<0.05){
print(paste(colnames(SubjectPCA)[i], "vs", row.names(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean)[j], "p-value=", summary.lm(lm(SubjectPCA[,i]~CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[j,]))$coefficient[8], "r-squared=",  round(summary.lm(lm(SubjectPCA[,i]~CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[j,]))[8][[1]], digits=3), sep="  "))}else{}
}		
}
)

)
cat(out, file="15 Flagged Relationships Between Subject Variables and PCA.txt", sep="\n", append=TRUE)
close(FlaggedRelationshipsBetweenIVandPCA)
rm(out)


#PC1 is most strongly related to both astrocytes and endothelial cells. What if I combine the measure?

row.names(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean)
 [1] "Astrocyte"                "Endothelial"              "Microglia"                "Mural"                   
 [5] "Neuron_All"               "Neuron_Interneuron"       "Neuron_Projection"        "Oligodendrocyte"         
 [9] "Oligodendrocyte_Immature" "RBC" 


summary.lm(lm(SubjectPCA[,1]~CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[1,]+CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[2,]))
#Multiple R-squared:  0.9094, same as a linear combination of the two variables.

summary.lm(lm(SubjectPCA[,1]~CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[1,]+CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[2,]+CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[5,]))
#Adding neurons into the model really doesn't make a difference

summary.lm(lm(SubjectPCA[,1]~CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[1,]+CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[2,]+AgonalFactorNoNA3))
#And agonal factor no longer matters for PC1

summary.lm(lm(SubjectPCA[,1]~CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[1,]+CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[2,]+BrainPHNoNA3))
#And brain pH no longer matters for PC1 either


temp<-CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[1,]+CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[2,]

summary.lm(lm(SubjectPCA[,1]~temp))

Call:
lm(formula = SubjectPCA[, 1] ~ temp)

Residuals:
    Min      1Q  Median      3Q     Max 
-7.2320 -2.4853 -0.2628  2.0007  8.6132 

Coefficients:
             Estimate Std. Error t value Pr(>|t|)    
(Intercept) 2.414e-16  2.619e-01    0.00        1    
temp        2.076e+01  5.291e-01   39.23   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 3.281 on 155 degrees of freedom
Multiple R-squared:  0.9085,	Adjusted R-squared:  0.9079 
F-statistic:  1539 on 1 and 155 DF,  p-value: < 2.2e-16

#So basically 90% of PC1 can be explained by a combined measure of Astrocyte/Endothelial Indices


png("CombinedAstrocyteEndothelialvsPC1.png")
plot(SubjectPCA[,1]~temp, xlab="Combined Astrocyte & Endothelial Indices", ylab="PC1", main="91% of PC1 Explained by Support Cells")
dev.off()

#Fancier graphs:
png("CombinedAstrocyteEndothelialvsPC1.png")
plot(SubjectPCA[,1]~temp, xlab="Combined Astrocyte & Endothelial Indices", ylab="PC1", cex.lab=1.4)
BestFitLine<-lm(SubjectPCA[,1]~temp)
abline(BestFitLine, col=2, lwd=2)
mtext(paste("P-value: ", format(summary.lm(BestFitLine)$coefficients[2,4], digits=3, scientific=T), ", R-squared: ", format(summary.lm(BestFitLine)$r.squared, digits=3), sep=""), cex=1.4)
dev.off()

png("NeuronAllvsPC1.png")
plot(SubjectPCA[,1]~CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[5,], xlab="Neuron (All) Index", ylab="PC1", cex.lab=1.4)
BestFitLine<-lm(SubjectPCA[,1]~CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[5,])
abline(BestFitLine, col=4, lwd=2)
mtext(paste("P-value: ", format(summary.lm(BestFitLine)$coefficients[2,4], digits=3, scientific=T), ", R-squared: ", format(summary.lm(BestFitLine)$r.squared, digits=3), sep=""), cex=1.4)
dev.off()


png("OligodendrocytevsPC2.png")
plot(SubjectPCA[,2]~CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[8,], xlab="Oligodendrocyte Index", ylab="PC2", cex.lab=1.4)
BestFitLine<-lm(SubjectPCA[,2]~CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[8,])
abline(BestFitLine, col=2, lwd=2)
mtext(paste("P-value: ", format(summary.lm(BestFitLine)$coefficients[2,4], digits=3, scientific=T), ", R-squared: ", format(summary.lm(BestFitLine)$r.squared, digits=3), sep=""), cex=1.4)
dev.off()

png("ProjectionNeuronvsPC2.png")
plot(SubjectPCA[,2]~CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[7,], xlab="Neuron (Projection) Index", ylab="PC2", cex.lab=1.4)
BestFitLine<-lm(SubjectPCA[,2]~CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[7,])
abline(BestFitLine, col=4, lwd=2)
mtext(paste("P-value: ", format(summary.lm(BestFitLine)$coefficients[2,4], digits=3, scientific=T), ", R-squared: ", format(summary.lm(BestFitLine)$r.squared, digits=3), sep=""), cex=1.4)
dev.off()



#Outputting higher resolution versions for figures:

row.names(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean)

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





