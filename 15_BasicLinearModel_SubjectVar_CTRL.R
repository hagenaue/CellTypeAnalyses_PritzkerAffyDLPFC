#This is my attempt to pull together a more organized code document specifically related to cell type index analyses using DLPFC Affymetrix data from the Pritzker foundation
##Megan Hagenauer, October 27 2016

#**********************************************

#I need to output results from a basic linear model for just controls (for David Lyons):
  
 #I need to output results from a basic LM4 for comparison:
  GeneByCellTypeSubjVar2_Pvalues<-matrix(0, length(SignalSortedNoNA3[,1]), 6)
GeneByCellTypeSubjVar2_Betas<-matrix(0, length(SignalSortedNoNA3[,1]), 6)
colnames(GeneByCellTypeSubjVar2_Pvalues)<-c("Intercept", "BrainPH", "AgonalFactor", "PMI", "Age", "Gender")
colnames(GeneByCellTypeSubjVar2_Betas)<-c("Intercept", "BrainPH", "AgonalFactor", "PMI", "Age", "Gender")
row.names(GeneByCellTypeSubjVar2_Pvalues)<-row.names(SignalSortedNoNA3)
row.names(GeneByCellTypeSubjVar2_Betas)<-row.names(SignalSortedNoNA3)
head(GeneByCellTypeSubjVar2_Pvalues)


for(i in c(1:length(SignalSortedNoNA3[,1]))){
  
  temp<-summary.lm(lm(SignalSortedNoNA3[i,DiagnosisNoNA3=="Control"]~BrainPHCentered[DiagnosisNoNA3=="Control"]+AgonalFactorNoNA3[DiagnosisNoNA3=="Control"]+PMICentered[DiagnosisNoNA3=="Control"]+ AgeCentered[DiagnosisNoNA3=="Control"]+ GenderNoNA3[DiagnosisNoNA3=="Control"]))
  
  GeneByCellTypeSubjVar2_Betas[i,]<-temp$coefficients[,1]
  GeneByCellTypeSubjVar2_Pvalues[i,]<-temp$coefficients[,4]
  
}

GeneByCellTypeSubjVar2_Pvalues2<-cbind(GeneByCellTypeSubjVar2_Pvalues, GeneNames)
GeneByCellTypeSubjVar2_Betas2<-cbind(GeneByCellTypeSubjVar2_Betas, GeneNames)


write.csv(GeneByCellTypeSubjVar2_Pvalues2, "GeneByCellTypeSubjVar2_Pvalues.csv")
write.csv(GeneByCellTypeSubjVar2_Betas2, "GeneByCellTypeSubjVar2_Betas.csv")


for (i in c(1:length(GeneByCellTypeSubjVar2_Pvalues[1,]))){
  png(paste(paste("17 Histogram of Raw Pvalues for", colnames(GeneByCellTypeSubjVar2_Pvalues)[i], sep="  "), "png", sep="."))	
  hist(GeneByCellTypeSubjVar2_Pvalues[,i], breaks=100, col=i, main=paste("Raw P-values for", colnames(GeneByCellTypeSubjVar2_Pvalues)[i], sep="  "), xlab="Unadjusted p-value", ylab="Count")
  abline(a=(length(GeneByCellTypeSubjVar2_Pvalues[,1])/100), b=0)
  dev.off()		
}	


GeneByCellTypeSubjVar2_PvaluesAdj<-matrix(0, length(GeneByCellTypeSubjVar2_Pvalues[,1]), length(GeneByCellTypeSubjVar2_Pvalues[1,]))
colnames(GeneByCellTypeSubjVar2_PvaluesAdj)<-colnames(GeneByCellTypeSubjVar2_Pvalues)
row.names(GeneByCellTypeSubjVar2_PvaluesAdj)<-row.names(GeneByCellTypeSubjVar2_Pvalues)



for (i in c(1:length(GeneByCellTypeSubjVar2_Pvalues[1,]))){
  
  #Applying two different types of multiple-comparison corrections to the raw p-values (Benjamini-Hochberg and Benjamini-Yekutieli):
  TempPvalAdj<-mt.rawp2adjp(GeneByCellTypeSubjVar2_Pvalues[,i], proc=c("BH"))
  GeneByCellTypeSubjVar2_PvaluesAdj[,i]<-TempPvalAdj$adjp[order(TempPvalAdj$index),2]
  
}

GeneByCellTypeSubjVar2_PvaluesAdj2<-cbind(GeneByCellTypeSubjVar2_PvaluesAdj, GeneNames)
write.csv(GeneByCellTypeSubjVar2_PvaluesAdj2, "GeneByCellTypeSubjVar2_PvaluesAdj.csv")

GeneByCellTypeSubjVar2DF<-as.data.frame(cbind(GeneByCellTypeSubjVar2_Betas, GeneByCellTypeSubjVar2_Pvalues, GeneByCellTypeSubjVar2_PvaluesAdj))

temp<-cbind(GeneNames, GeneByCellTypeSubjVar2DF)
write.csv(temp, "GeneByCellTypeSubjVar2DF.csv" )
  
ControlDemographics<-cbind(BrainPHNoNA3[DiagnosisNoNA3=="Control"],AgonalFactorNoNA3[DiagnosisNoNA3=="Control"],HoursFinalCorrectedNoNA3[DiagnosisNoNA3=="Control"],AgeNoNA3[DiagnosisNoNA3=="Control"],GenderNoNA3[DiagnosisNoNA3=="Control"])

colnames(ControlDemographics)<-c("BrainPH", "AgonalFactor", "PMI", "Age", "Gender")

write.csv(ControlDemographics, "ControlDemographics.csv")
