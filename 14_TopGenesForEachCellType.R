#This is my attempt to pull together a more organized code document specifically related to cell type index analyses using DLPFC Affymetrix data from the Pritzker foundation
##Megan Hagenauer, October 27 2016

#**********************************************

CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap<-as.data.frame(t(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean))

GeneByCellTypeSubjVar_Pvalues<-matrix(0, length(SignalSortedNoNA3[,1]), 18)
GeneByCellTypeSubjVar_Betas<-matrix(0, length(SignalSortedNoNA3[,1]), 18)
colnames(GeneByCellTypeSubjVar_Pvalues)<-c("Intercept", "Astrocyte", "Endothelial", "Microglia", "Mural", "Neuron_All", "Neuron_Interneuron", "Neuron_Projection", "Oligodendrocyte", "RBC", "BrainPH", "AgonalFactor", "PMI", "Age", "Gender", "BP", "MDD", "SCHIZ")
colnames(GeneByCellTypeSubjVar_Betas)<-c("Intercept", "Astrocyte", "Endothelial", "Microglia", "Mural", "Neuron_All", "Neuron_Interneuron", "Neuron_Projection", "Oligodendrocyte", "RBC", "BrainPH", "AgonalFactor", "PMI", "Age", "Gender", "BP", "MDD", "SCHIZ")
row.names(GeneByCellTypeSubjVar_Pvalues)<-row.names(SignalSortedNoNA3)
row.names(GeneByCellTypeSubjVar_Betas)<-row.names(SignalSortedNoNA3)
head(GeneByCellTypeSubjVar_Pvalues)

for(i in c(1:length(SignalSortedNoNA3[,1]))){
	
	temp<-summary.lm(lm(SignalSortedNoNA3[i,]~CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap$Astrocyte+CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap$Endothelial+CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap$Microglia+CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap$Mural+CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap$Neuron_All+CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap$Neuron_Interneuron+CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap$Neuron_Projection+CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap$Oligodendrocyte+CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap$RBC+BrainPHNoNA3 +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeNoNA3+ GenderNoNA3 + DiagnosisNoNA3))

GeneByCellTypeSubjVar_Betas[i,]<-temp$coefficients[,1]
GeneByCellTypeSubjVar_Pvalues[i,]<-temp$coefficients[,4]

}

GeneByCellTypeSubjVar_Pvalues2<-cbind(GeneByCellTypeSubjVar_Pvalues, GeneNames)
GeneByCellTypeSubjVar_Betas2<-cbind(GeneByCellTypeSubjVar_Betas, GeneNames)


write.csv(GeneByCellTypeSubjVar_Pvalues2, "GeneByCellTypeSubjVar_Pvalues.csv")
write.csv(GeneByCellTypeSubjVar_Betas2, "GeneByCellTypeSubjVar_Betas.csv")

#I should really output BH-corrected p-values and histograms for those, and automatically output the top 50 positive relationships.


for (i in c(1:length(GeneByCellTypeSubjVar_Pvalues[1,]))){
png(paste(paste("17 Histogram of Raw Pvalues for", colnames(GeneByCellTypeSubjVar_Pvalues)[i], sep="  "), "png", sep="."))	
hist(GeneByCellTypeSubjVar_Pvalues[,i], breaks=100, col=i, main=paste("Raw P-values for", colnames(GeneByCellTypeSubjVar_Pvalues)[i], sep="  "), xlab="Unadjusted p-value", ylab="Count")
abline(a=(length(GeneByCellTypeSubjVar_Pvalues[,1])/100), b=0)
dev.off()		
}	



GeneByCellTypeSubjVar_PvaluesAdj<-matrix(0, length(GeneByCellTypeSubjVar_Pvalues[,1]), length(GeneByCellTypeSubjVar_Pvalues[1,]))
colnames(GeneByCellTypeSubjVar_PvaluesAdj)<-colnames(GeneByCellTypeSubjVar_Pvalues)
row.names(GeneByCellTypeSubjVar_PvaluesAdj)<-row.names(GeneByCellTypeSubjVar_Pvalues)


for (i in c(1:length(GeneByCellTypeSubjVar_Pvalues[1,]))){

#Applying two different types of multiple-comparison corrections to the raw p-values (Benjamini-Hochberg and Benjamini-Yekutieli):
TempPvalAdj<-mt.rawp2adjp(GeneByCellTypeSubjVar_Pvalues[,i], proc=c("BH"))
GeneByCellTypeSubjVar_PvaluesAdj[,i]<-TempPvalAdj$adjp[order(TempPvalAdj$index),2]

}

GeneByCellTypeSubjVar_PvaluesAdj2<-cbind(GeneByCellTypeSubjVar_PvaluesAdj, GeneNames)
write.csv(GeneByCellTypeSubjVar_PvaluesAdj2, "GeneByCellTypeSubjVar_PvaluesAdj.csv")

GeneByCellTypeSubjVarDF<-as.data.frame(cbind(GeneByCellTypeSubjVar_Betas, GeneByCellTypeSubjVar_Pvalues, GeneByCellTypeSubjVar_PvaluesAdj))

write.csv(GeneByCellTypeSubjVarDF, "GeneByCellTypeSubjVarDF.csv")


for (i in c(1:length(GeneByCellTypeSubjVar_Pvalues[1,]))){
	temp<-data.frame(GeneNames,GeneByCellTypeSubjVarDF[,i], GeneByCellTypeSubjVarDF[,i+18], GeneByCellTypeSubjVarDF[,i+36])
	tempUp<-temp[temp[,3]>0,]
    tempUpSorted<-tempUp[order(tempUp[,4]),]
    tempUpSorted2<-data.frame( tempUpSorted, as.numeric(tempUpSorted[,2]%in%CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[,4]))
    colnames(tempUpSorted2)<-c("Probe", "Gene.Symbol", "Beta", "Raw.Pvalue", "Adj.Pvalue","CellTypeSpecific")
    tempOutput<-tempUpSorted2[c(1:100),]
    tempFileName<-paste("Top100Probes_", colnames(GeneByCellTypeSubjVarDF)[i], ".csv")
    write.csv(tempOutput, tempFileName)    
    }

