#This is my attempt to pull together a more organized code document specifically related to cell type index analyses using DLPFC Affymetrix data from the Pritzker foundation
##Megan Hagenauer, October 27 2016

#**********************************************

#This is just an annotated version of the basic linear model output for each cell type specific gene vs. all major confounds and diagnosis:

#I need to output results from a basic LM4 for comparison:
GeneByCellTypeSubjVar2_Pvalues<-matrix(0, length(CellTypeSpecificGenes_Master3_SignalNoNA3_Norm[,1]), 9)
GeneByCellTypeSubjVar2_Betas<-matrix(0, length(CellTypeSpecificGenes_Master3_SignalNoNA3_Norm[,1]), 9)
colnames(GeneByCellTypeSubjVar2_Pvalues)<-c("Intercept", "BrainPH", "AgonalFactor", "PMI", "Age", "Gender", "BP", "MDD", "SCHIZ")
colnames(GeneByCellTypeSubjVar2_Betas)<-c("Intercept", "BrainPH", "AgonalFactor", "PMI", "Age", "Gender", "BP", "MDD", "SCHIZ")
row.names(GeneByCellTypeSubjVar2_Pvalues)<-row.names(CellTypeSpecificGenes_Master3_SignalNoNA3_Norm)
row.names(GeneByCellTypeSubjVar2_Betas)<-row.names(CellTypeSpecificGenes_Master3_SignalNoNA3_Norm)
head(GeneByCellTypeSubjVar2_Pvalues)


for(i in c(1:length(CellTypeSpecificGenes_Master3_SignalNoNA3_Norm[,1]))){
	
	temp<-summary.lm(lm(as.numeric(CellTypeSpecificGenes_Master3_SignalNoNA3_Norm[i,c(15:171)])~BrainPHCentered +AgonalFactorNoNA3 + PMICentered+ AgeCentered+ GenderNoNA3 + DiagnosisNoNA3))

GeneByCellTypeSubjVar2_Betas[i,]<-temp$coefficients[,1]
GeneByCellTypeSubjVar2_Pvalues[i,]<-temp$coefficients[,4]

}

GeneByCellTypeSubjVar2_Pvalues2<-cbind(GeneByCellTypeSubjVar2_Pvalues, CellTypeSpecificGenes_Master3_SignalNoNA3_Norm[,c(1:14)])
GeneByCellTypeSubjVar2_Betas2<-cbind(GeneByCellTypeSubjVar2_Betas, CellTypeSpecificGenes_Master3_SignalNoNA3_Norm[,c(1:14)])


write.csv(GeneByCellTypeSubjVar2_Pvalues2, "GeneByCellTypeSubjVar2_Pvalues.csv")
write.csv(GeneByCellTypeSubjVar2_Betas2, "GeneByCellTypeSubjVar2_Betas.csv")

write.csv(CellTypeSpecificGenes_Master3_SignalNoNA3_Norm[,c(1:15)], "Annotation.csv")
