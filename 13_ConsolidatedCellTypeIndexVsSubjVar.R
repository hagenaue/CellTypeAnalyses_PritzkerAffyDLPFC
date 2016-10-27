#This is my attempt to pull together a more organized code document specifically related to cell type index analyses using DLPFC Affymetrix data from the Pritzker foundation
##Megan Hagenauer, October 27 2016

#**********************************************

#Running the same code to examine the relationships between subject variables and cell type balance, but for primary category indices with no overlapping genes between primary categories:

#Using a scatterplot with best fit line to visually examine the relationships between the continuous subject variables and consolidated cell type indices:
for (i in 1:length(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[,1])){
for(j in 1:length(SubjectContinuousVariables[1,])){
png(paste(paste(row.names(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean)[i], "vs", colnames(SubjectContinuousVariables)[j], sep="  "), "png", sep="."))	
plot(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[i,]~SubjectContinuousVariables[,j], main=paste(row.names(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean)[i], "vs", colnames(SubjectContinuousVariables)[j], sep="  "), xlab=colnames(SubjectContinuousVariables)[j], ylab=row.names(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean)[i])
RegressionLine<-lm(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[i,]~SubjectContinuousVariables[,j])
abline(RegressionLine, col=2)
mtext(paste("p-value=", round(summary.lm(RegressionLine)$coefficients[8], digits=4), ", r-squared=",  round(summary.lm(lm(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[i,]~SubjectContinuousVariables[,j]))[8][[1]], digits=3)))
dev.off()		
}		
}



#Using boxplots to visually examine the relationships between consolidated cell type indices and categorical subject variables:
for (i in 1:length(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[,1])){
for(j in 1:length(SubjectFactorVariables[1,])){
png(paste(paste(row.names(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "), "png", sep="."))	
boxplot(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[i,]~SubjectFactorVariables[,j], main=paste(row.names(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "), xlab=colnames(SubjectFactorVariables)[j], ylab=row.names(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean)[i])
mtext(paste("p-value=", round(summary(aov(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[i,]~SubjectFactorVariables[,j]))[[1]][["Pr(>F)"]][1], digits=4), ", r-squared=",  round(summary.lm(lm(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[i,]~SubjectFactorVariables[,j]))[8][[1]], digits=3)))
dev.off()		
}		
}

#Flagging variables that are collinear with other subject variables:
FlaggedRelationshipsBetweenIV<-file("14 Flagged Relationships Between Subject Variables.txt")
out<-c(

#Using linear regression to examine the statistical relationships between consolidated cell type indices and the continuous subject variables:
capture.output(
for (i in 1:length(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[,1])){
for(j in 1:length(SubjectContinuousVariables[1,])){
if(summary.lm(lm(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[i,]~SubjectContinuousVariables[,j]))$coefficient[8]<0.05){
print(paste(row.names(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean)[i], "vs", colnames(SubjectContinuousVariables)[j], "p-value=", summary.lm(lm(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[i,]~SubjectContinuousVariables[,j]))$coefficient[8], ", r-squared=",  round(summary.lm(lm(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[i,]~SubjectContinuousVariables[,j]))[8][[1]], digits=3), sep="  "))}else{}
}		
}
),

#Using anova to examine the statistical relationships between consolidated cell type indices and categorical subject variables:
capture.output(
for (i in 1:length(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[,1])){
for(j in 1:length(SubjectFactorVariables[1,])){
if(summary(aov(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[i,]~SubjectFactorVariables[,j]))[[1]][["Pr(>F)"]][1]<0.05){
print(paste(row.names(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean)[i], "vs", colnames(SubjectFactorVariables)[j], "p-value=", summary(aov(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[i,]~SubjectFactorVariables[,j]))[[1]][["Pr(>F)"]][1], ", r-squared=",  round(summary.lm(lm(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[i,]~SubjectFactorVariables[,j]))[8][[1]], digits=3), sep="  "))	
}else{}		
}		
}
)


)
cat(out, file="14 Flagged Relationships Between Subject Variables.txt", sep="\n", append=TRUE)
close(FlaggedRelationshipsBetweenIV)
rm(out)


#************************************************************

#Making some prettier plots:

png("CellTypeSpecificGenes_Master2_ExsanguinatedBoxPlot.png", width=300, height=480)
boxplot(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[56,]~ExsanguinatedNoNA3, xlab="Exsanguinated Or Not", ylab="Red Blood Cell Index", col="red")
mtext("Exsanguination Decreases", line=3)
mtext("RBC-Specific Expression", line=2)
mtext(paste("p-value = ", round(summary.lm(lm(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[56,]~ExsanguinatedNoNA3))$coefficients[2,4], digits=4), sep=""))
dev.off()


Call:
lm(formula = CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[56, 
    ] ~ ExsanguinatedNoNA3)

Residuals:
     Min       1Q   Median       3Q      Max 
-0.87790 -0.29376 -0.00281  0.32672  1.27069 

Coefficients:
                   Estimate Std. Error t value Pr(>|t|)    
(Intercept)         0.03576    0.03314   1.079 0.282243    
ExsanguinatedNoNA3 -0.40107    0.11099  -3.613 0.000408 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.3964 on 155 degrees of freedom
Multiple R-squared:  0.0777,	Adjusted R-squared:  0.07175 
F-statistic: 13.06 on 1 and 155 DF,  p-value: 0.0004078




png("CellTypeSpecificGenes_Master2_ExsanguinatedBoxPlot.png", width=300, height=480)
boxplot(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[10,]~ExsanguinatedNoNA3, xlab="Exsanguinated Or Not", ylab="Red Blood Cell Index", col="red")
mtext("Exsanguination Decreases", line=2)
mtext("RBC-Specific Expression", line=1)
dev.off()


png("AgonalFactorVsEndothelialBoxPlot.png", width=300, height=480)
boxplot(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[2,]~AgonalFactorNoNA3, xlab="Agonal Factor", ylab="Endothelial Index", col="blue")
mtext("Hypoxia Increases", line=2)
mtext("Endothelial-Specific Expression", line=1)
dev.off()

summary.lm(lm(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[2,]~AgonalFactorNoNA3+BrainPHNoNA3))
#R-squared: 0.36

png("AgonalFactorVsAstrocyteBoxPlot.png", width=300, height=480)
boxplot(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[1,]~AgonalFactorNoNA3, xlab="Agonal Factor", ylab="Astrocyte Index", col="blue")
mtext("Hypoxia Increases", line=2)
mtext("Astrocyte-Specific Expression", line=1)
dev.off()

summary.lm(lm(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[1,]~AgonalFactorNoNA3+BrainPHNoNA3))
# R-squared: 0.3025

png("AgonalFactorVsNeuronBoxPlot.png", width=300, height=480)
boxplot(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[5,]~AgonalFactorNoNA3, xlab="Agonal Factor", ylab="Neuron_All Index", col="blue")
mtext("Hypoxia Decreases", line=2)
mtext("Neuron-Specific Expression", line=1)
dev.off()


png("CellTypeSpecificGenes_Master2_ExsanguinatedBoxPlot.png", width=270, height=640)
boxplot(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[10,]~ExsanguinatedNoNA3, xlab="Exsanguinated Or Not", ylab="Red Blood Cell Index", col="blue", cex.lab=1.5, cex.axis=1.5, lwd=2.5)
mtext("Exsanguination Decreases", line=2.5, cex=1.5, font=2)
mtext("RBC-Specific Expression", line=1, cex=1.5, font=2)
dev.off()

png("AgonalFactorVsEndothelialBoxPlot.png", width=350, height=640)
boxplot(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[2,]~AgonalFactorNoNA3, xlab="Agonal Factor", ylab="Endothelial Index", col="red", cex.lab=1.5, cex.axis=1.5, lwd=2.5)
mtext("Hypoxia Increases", line=2.5, cex=1.5, font=2)
mtext("Endothelial-Specific Expression", line=1, cex=1.5, font=2)
dev.off()

png("AgonalFactorVsAstrocyteBoxPlot.png", width=350, height=640)
boxplot(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[1,]~AgonalFactorNoNA3, xlab="Agonal Factor", ylab="Astrocyte Index", col="red", cex.lab=1.5, cex.axis=1.5, lwd=2.5)
mtext("Hypoxia Increases", line=2.5, cex=1.5, font=2)
mtext("Astrocyte-Specific Expression", line=1, cex=1.5, font=2)
dev.off()

png("AgonalFactorVsNeuronBoxPlot.png", width=350, height=640)
boxplot(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[5,]~AgonalFactorNoNA3, xlab="Agonal Factor", ylab="Neuron (All) Index", col="blue", cex.lab=1.5, cex.axis=1.5, lwd=2.5)
mtext("Hypoxia Decreases", line=2.5, cex=1.5, font=2)
mtext("Neuron-Specific Expression", line=1, cex=1.5, font=2)
dev.off()


summary.lm(lm(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[5,]~AgonalFactorNoNA3+BrainPHNoNA3))
#R-squared: 0.3715

png("DiagnosisVsAstrocyteBoxPlot.png", width=350, height=640)
boxplot(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[1,]~DiagnosisNoNA3, xlab="Diagnosis", ylab="Astrocyte Index", col="blue", cex.lab=1.5, cex.axis=1.5, lwd=2.5)
mtext("Major Depressive Disorder Decreases", line=2.5, cex=1.5, font=2)
mtext("Astrocyte-Specific Expression", line=1, cex=1.5, font=2)
dev.off()

png("AgeVsInterneuronPlot.png", width=350, height=640)
plot(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[6,]~AgeNoNA3, xlab="Age", ylab="Neuron (Interneuron) Index", cex.lab=1.5, cex.axis=1.5, lwd=2.5)
BestFitLine<-lm(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[6,]~AgeNoNA3)
abline(BestFitLine, col=4, lwd=4)
mtext("Age Decreases", line=2.5, cex=1.5, font=2)
mtext("Interneuron-Specific Expression", line=1, cex=1.5, font=2)
dev.off()


#I'm going to try making a higher resolution verison:

pdf(file="AgeVsInterneuronPlot_hres.pdf", width=3, height=5.5, pointsize=10)
plot(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[6,]~AgeNoNA3, xlab="Age", ylab="Neuron (Interneuron) Index", cex.lab=1.3, cex.axis=1.3)
BestFitLine<-lm(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[6,]~AgeNoNA3)
abline(BestFitLine, col="dodgerblue3", lwd=4)
mtext("Age Decreases", line=3, font=2, cex=1.3)
mtext("Interneuron-Specific", line=2, font=2, cex=1.3)
mtext("Expression", line=1, font=2, cex=1.3)
dev.off()

pdf(file="DiagnosisVsAstrocyteBoxPlot.pdf", width=3.5, height=5.5, pointsize=10)
boxplot(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[1,]~DiagnosisNoNA3, xlab="Diagnosis", ylab="Astrocyte Index", col="dodgerblue3", font.lab=2, lwd=2, cex.lab=1.3, cex.axis=1.3)
mtext("Major Depressive Disorder", line=3,  font=2, cex=1.3)
mtext("Decreases Astrocyte-Specific", line=2,  font=2, cex=1.3)
mtext("Expression", line=1,  font=2, cex=1.3)
dev.off()

pdf(file="AgonalFactorVsNeuronBoxPlot.pdf", width=3, height=5.5, pointsize=10)
boxplot(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[5,]~AgonalFactorNoNA3, xlab="Agonal Factor", ylab="Neuron (All) Index", col="dodgerblue3", font.lab=2, lwd=2, cex.lab=1.3, cex.axis=1.3)
mtext("Hypoxia Decreases", line=3,  font=2, cex=1.3)
mtext("Neuron-Specific", line=2,  font=2, cex=1.3)
mtext("Expression", line=1,  font=2, cex=1.3)
dev.off()

pdf("AgonalFactorVsEndothelialBoxPlot.pdf", width=3, height=5.5, pointsize=10)
boxplot(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[2,]~AgonalFactorNoNA3, xlab="Agonal Factor", ylab="Endothelial Index", col="red", font.lab=2, lwd=2, cex.lab=1.3, cex.axis=1.3)
mtext("Hypoxia Increases", line=3, font=2, cex=1.3)
mtext("Endothelial-Specific", line=2, font=2, cex=1.3)
mtext("Expression", line=1, font=2, cex=1.3)
dev.off()

pdf("ExsanguinationVSRBC.pdf", width=3, height=3, pointsize=10)
boxplot(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[10,]~ExsanguinatedNoNA3, xlab="Exsanguinated Or Not", ylab="Red Blood Cell Index", col="dodgerblue3", font.lab=2, lwd=2, cex.lab=1.3, cex.axis=1.3)
mtext("Exsanguination Decreases", line=3, font=2, cex=1.3)
mtext("RBC-Specific", line=2, font=2, cex=1.3)
mtext("Expression", line=1, font=2, cex=1.3)
dev.off()


#**********************************************************


#Examining the relationship between the consolidated cell type indices and subject variables while controlling for other subject variables:


StatModelsCellIndicesVsSubjVar<-file("14 StatisticalModels_CellIndicesvsSubjVar.txt")
out<-c(

capture.output(
for(i in c(1:length(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[,1]))){
	
print(row.names(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean)[i])

print(summary.lm(lm(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[i,] ~ BrainPHNoNA3 + AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+AgeNoNA3+ GenderNoNA3 + DiagnosisNoNA3)))
    
  }
)
)
cat(out, file="14 StatisticalModels_CellIndicesvsSubjVar.txt", sep="\n", append=TRUE)
close(StatModelsCellIndicesVsSubjVar)
rm(out)


StatModelsCellIndicesVsSubjVar2Coefficients<-matrix(0, length(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[,1]), 9)

StatModelsCellIndicesVsSubjVar2<-file("14 StatisticalModels_CellIndicesvsSubjVar2.txt")
out<-c(

capture.output(

for(i in c(1:length(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[,1]))){

print(row.names(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean)[i])
	
print(summary.lm(lm(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[i,] ~ BrainPHNoNA3 + AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+AgeNoNA3+ GenderNoNA3 + DiagnosisNoNA3+ExsanguinatedNoNA3)))
   
  }
)
)
cat(out, file="14 StatisticalModels_CellIndicesvsSubjVar2.txt", sep="\n", append=TRUE)
close(StatModelsCellIndicesVsSubjVar2)
rm(out)



StatModelsCellIndicesVsSubjVar3<-file("14 StatisticalModels_CellIndicesvsSubjVar3.txt")
out<-c(

capture.output(

for(i in c(1:length(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[,1]))){

print(row.names(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean)[i])
	
print(summary.lm(lm(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[i,] ~ BrainPHNoNA3 + AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+AgeNoNA3+ GenderNoNA3 + DiagnosisNoNA3+ExsanguinatedNoNA3+SuicideNoNA3)))
   
  }
)
)
cat(out, file="14 StatisticalModels_CellIndicesvsSubjVar3.txt", sep="\n", append=TRUE)
close(StatModelsCellIndicesVsSubjVar3)
rm(out)



StatModelsCellIndicesVsSubjVar4<-file("14 StatisticalModels_CellIndicesvsSubjVar4.txt")
out<-c(

capture.output(

for(i in c(1:length(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[,1]))){

print(row.names(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean)[i])
	
print(summary.lm(lm(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[i,] ~ BrainPHNoNA3 + AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+AgeNoNA3+ GenderNoNA3 + DiagnosisNoNA3+ExsanguinatedNoNA3+sinTODNoNA3+cosTODNoNA3)))
   
  }
)
)
cat(out, file="14 StatisticalModels_CellIndicesvsSubjVar4.txt", sep="\n", append=TRUE)
close(StatModelsCellIndicesVsSubjVar4)
rm(out)


StatModelsCellIndicesVsSubjVar5_Psych<-file("14 StatisticalModels_CellIndicesvsSubjVar5_Psych.txt")
out<-c(
  
  capture.output(
    
    for(i in c(1:length(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[,1]))){
      
      print(row.names(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean)[i])
      
      print(summary.lm(lm(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[i,] ~ BrainPHNoNA3 + AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+AgeNoNA3+ GenderNoNA3 + PsychiatricNoNA3+ExsanguinatedNoNA3)))
      
    }
  )
)
cat(out, file="14 StatisticalModels_CellIndicesvsSubjVar5_Psych.txt", sep="\n", append=TRUE)
close(StatModelsCellIndicesVsSubjVar5_Psych)
rm(out)


