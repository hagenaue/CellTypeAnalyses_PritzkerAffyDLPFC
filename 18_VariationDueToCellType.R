#This is my attempt to pull together a more organized code document specifically related to cell type index analyses using DLPFC Affymetrix data from the Pritzker foundation
##Megan Hagenauer, October 27 2016

#**********************************************

#Estimating how much variation tends to be accounted for by cell type:

GeneByCellTypeSubjVar_RSquared<-matrix(0, length(SignalSortedNoNA3[,1]), 2)
colnames(GeneByCellTypeSubjVar_RSquared)<-c("R-squared", "AdjRsquared")
row.names(GeneByCellTypeSubjVar_RSquared)<-row.names(SignalSortedNoNA3)

head(GeneByCellTypeSubjVar_RSquared)

for(i in c(1:length(SignalSortedNoNA3[,1]))){
	
	temp<-summary.lm(lm(SignalSortedNoNA3[i,]~CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap$Astrocyte+CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap$Endothelial+CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap$Microglia+CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap$Mural+CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap$Neuron_All+CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap$Neuron_Interneuron+CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap$Neuron_Projection+CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap$Oligodendrocyte+CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap$RBC+BrainPHNoNA3 +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeNoNA3+ GenderNoNA3 + DiagnosisNoNA3))

GeneByCellTypeSubjVar_RSquared[i,1]<-temp$r.squared
GeneByCellTypeSubjVar_RSquared[i,2]<-temp$adj.r.squared

}

GeneByCellTypeSubjVar_RSquared2<-cbind(GeneByCellTypeSubjVar_RSquared, GeneNames)

write.csv(GeneByCellTypeSubjVar_RSquared2, "GeneByCellTypeSubjVar_RSquared2_AllPrimaryCellTypesSubjVar.csv")

png("Histogram_Rsquared_AllPrimaryCellTypesSubjVar.png")
hist(GeneByCellTypeSubjVar_RSquared[,1], col=2, xlab="R-squared", main="Linear Model with Primary Cell Types and Subject Variables")
dev.off()

mean(GeneByCellTypeSubjVar_RSquared[,1])
#[1] 0.3535082
median(GeneByCellTypeSubjVar_RSquared[,1])
#[1] 0.322509
sd(GeneByCellTypeSubjVar_RSquared[,1])
#[1] 0.1702802

mean(GeneByCellTypeSubjVar_RSquared[,2])
#[1] 0.2744408
median(GeneByCellTypeSubjVar_RSquared[,2])
#[1] 0.2396503
sd(GeneByCellTypeSubjVar_RSquared[,2])
#[1] 0.1911059

png("Histogram_AdjRsquared_AllPrimaryCellTypesSubjVar.png")
hist(GeneByCellTypeSubjVar_RSquared[,2], col=3, xlab="Adjusted R-squared", main="Linear Model with Primary Cell Types and Subject Variables")
dev.off()


#Just Cell Type:
for(i in c(1:length(SignalSortedNoNA3[,1]))){
	
	temp<-summary.lm(lm(SignalSortedNoNA3[i,]~CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap$Astrocyte+CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap$Endothelial+CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap$Microglia+CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap$Mural+CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap$Neuron_All+CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap$Neuron_Interneuron+CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap$Neuron_Projection+CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap$Oligodendrocyte+CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap$RBC))

GeneByCellTypeSubjVar_RSquared[i,1]<-temp$r.squared
GeneByCellTypeSubjVar_RSquared[i,2]<-temp$adj.r.squared

}

GeneByCellTypeSubjVar_RSquared2<-cbind(GeneByCellTypeSubjVar_RSquared, GeneByCellTypeSubjVar_RSquared2)

#write.csv(GeneByCellTypeSubjVar_RSquared2, "GeneByCellTypeSubjVar_RSquared2_AllPrimaryCellTypesSubjVar.csv")


mean(GeneByCellTypeSubjVar_RSquared[,1])
#[1] 0.3050965
median(GeneByCellTypeSubjVar_RSquared[,1])
#[1] 0.2698552
sd(GeneByCellTypeSubjVar_RSquared[,1])
#[1] 0.1737719
sum(GeneByCellTypeSubjVar_RSquared[,1]>0.5)/length(GeneByCellTypeSubjVar_RSquared[,1])

mean(GeneByCellTypeSubjVar_RSquared[,2])
#[1] 0.2625514
median(GeneByCellTypeSubjVar_RSquared[,2])
#[1] 0.2251525
sd(GeneByCellTypeSubjVar_RSquared[,2])
#[1] 0.184411

length(GeneByCellTypeSubjVar_RSquared[,2])

png("Histogram_Rsquared_AllPrimaryCellTypes.png")
hist(GeneByCellTypeSubjVar_RSquared[,1], col=2, xlab="R-squared", xlim=c(0,1), main="Linear Model with All 10 Cell Types")
dev.off()

png("Histogram_AdjRsquared_AllPrimaryCellTypes.png")
hist(GeneByCellTypeSubjVar_RSquared[,2], col=3, xlab="Adjusted R-squared", xlim=c(0,1), main="Linear Model with All 10 Cell Types")
dev.off()

#Just the most prevalent (non-multicollinear) Cell Types:
for(i in c(1:length(SignalSortedNoNA3[,1]))){
  
  temp<-summary.lm(lm(SignalSortedNoNA3[i,]~CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap$Astrocyte+CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap$Microglia+CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap$Neuron_Interneuron+CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap$Neuron_Projection+CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap$Oligodendrocyte))
  
  GeneByCellTypeSubjVar_RSquared[i,1]<-temp$r.squared
  GeneByCellTypeSubjVar_RSquared[i,2]<-temp$adj.r.squared
  
}

GeneByCellTypeSubjVar_RSquared2<-cbind(GeneByCellTypeSubjVar_RSquared, GeneNames)
                                       
write.csv(GeneByCellTypeSubjVar_RSquared2, "GeneByCellTypeSubjVar_RSquared2_Mostprevalent.csv.csv")

png("Histogram_Rsquared_5PrevalentCellTypes.png")
hist(GeneByCellTypeSubjVar_RSquared[,1], col=2, xlab="R-squared", main="Linear Model with 5 Most Prevalent Cell Types")
dev.off()

mean(GeneByCellTypeSubjVar_RSquared[,2])

#Just 3 Cell Types:
for(i in c(1:length(SignalSortedNoNA3[,1]))){
	
	temp<-summary.lm(lm(SignalSortedNoNA3[i,]~CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap$Astrocyte++CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap$Neuron_All+CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap$Oligodendrocyte))

GeneByCellTypeSubjVar_RSquared[i,1]<-temp$r.squared
GeneByCellTypeSubjVar_RSquared[i,2]<-temp$adj.r.squared

}

GeneByCellTypeSubjVar_RSquared2<-cbind(GeneByCellTypeSubjVar_RSquared, GeneByCellTypeSubjVar_RSquared2)

#write.csv(GeneByCellTypeSubjVar_RSquared2, "GeneByCellTypeSubjVar_RSquared2_AllPrimaryCellTypesSubjVar.csv")



mean(GeneByCellTypeSubjVar_RSquared[,1])
#[1] 0.1714855
median(GeneByCellTypeSubjVar_RSquared[,1])
[1] 0.1157748
sd(GeneByCellTypeSubjVar_RSquared[,1])
[1] 0.1619384

mean(GeneByCellTypeSubjVar_RSquared[,2])
[1] 0.1552402
median(GeneByCellTypeSubjVar_RSquared[,2])
[1] 0.09843703
sd(GeneByCellTypeSubjVar_RSquared[,2])
[1] 0.1651137


png("Histogram_Rsquared_3CellTypes_OligoAstroNeuron.png")
hist(GeneByCellTypeSubjVar_RSquared[,1], col=2, xlab="R-squared", main="Linear Model with 3 Primary Cell Types: Oligodendrocytes, Astrocytes, Neurons")
dev.off()

png("Histogram_AdjRsquared_3CellTypes_OligoAstroNeuron.png")
hist(GeneByCellTypeSubjVar_RSquared[,2], col=3, xlab="Adjusted R-squared", main="Linear Model with 3 Primary Cell Types: Oligodendrocytes, Astrocytes, Neurons")
dev.off()


for(i in c(1:length(SignalSortedNoNA3[,1]))){
  
  temp<-summary.lm(lm(SignalSortedNoNA3[i,]~CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap$Astrocyte+CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap$Neuron_Projection))
  
  GeneByCellTypeSubjVar_RSquared[i,1]<-temp$r.squared
  GeneByCellTypeSubjVar_RSquared[i,2]<-temp$adj.r.squared
  
}

GeneByCellTypeSubjVar_RSquared2<-cbind(GeneByCellTypeSubjVar_RSquared, GeneByCellTypeSubjVar_RSquared2)

#Just the Astrocyte Index alone:
mean(GeneByCellTypeSubjVar_RSquared[,1])
[1] 0.09938388
median(GeneByCellTypeSubjVar_RSquared[,1])
[1] 0.04256976

mean(GeneByCellTypeSubjVar_RSquared[,2])
[1] 0.09357346
median(GeneByCellTypeSubjVar_RSquared[,2])
[1] 0.03639279

#Just the Astrocyte Index & Oligodendrocyte index alone:
mean(GeneByCellTypeSubjVar_RSquared[,1])
[1] 0.1431834
median(GeneByCellTypeSubjVar_RSquared[,1])
[1] 0.0857138

mean(GeneByCellTypeSubjVar_RSquared[,2])
[1] 0.1320559

#Just the Astrocyte Index & Projection Neuron indices alone:
mean(GeneByCellTypeSubjVar_RSquared[,1])
[1] 0.1708456
median(GeneByCellTypeSubjVar_RSquared[,1])
[1] 0.1252154

mean(GeneByCellTypeSubjVar_RSquared[,2])
[1] 0.1600773



#Just 4 Cell Types:
for(i in c(1:length(SignalSortedNoNA3[,1]))){
	
	temp<-summary.lm(lm(SignalSortedNoNA3[i,]~CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap$Astrocyte++CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap$Neuron_All+CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap$Oligodendrocyte+CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap$Neuron_Projection))

GeneByCellTypeSubjVar_RSquared[i,1]<-temp$r.squared
GeneByCellTypeSubjVar_RSquared[i,2]<-temp$adj.r.squared

}

mean(GeneByCellTypeSubjVar_RSquared[,1])
[1] 0.2246758
median(GeneByCellTypeSubjVar_RSquared[,1])
[1] 0.1763614
sd(GeneByCellTypeSubjVar_RSquared[,1])
[1] 0.1733974

 
mean(GeneByCellTypeSubjVar_RSquared[,2])
[1] 0.2042725
median(GeneByCellTypeSubjVar_RSquared[,2])
[1] 0.1546867
> sd(GeneByCellTypeSubjVar_RSquared[,2])
[1] 0.1779605

png("Histogram_Rsquared_3CellTypes_OligoAstroNeuronProjection.png")
hist(GeneByCellTypeSubjVar_RSquared[,1], col=2, xlab="R-squared", main="Linear Model with 4 Primary Cell Types: Oligodendrocytes, Astrocytes, Neurons, Projection Neurons")
dev.off()

png("Histogram_AdjRsquared_3CellTypes_OligoAstroNeuronProjection.png")
hist(GeneByCellTypeSubjVar_RSquared[,2], col=3, xlab="Adjusted R-squared", main="Linear Model with 4 Primary Cell Types: Oligodendrocytes, Astrocytes, Neurons, Projection Neurons")
dev.off()

with microglia:
mean(GeneByCellTypeSubjVar_RSquared[,1])
[1] 0.2522603

with endothelial:
mean(GeneByCellTypeSubjVar_RSquared[,1])
[1] 0.2489281

with interneurons:
[1] 0.2446624

with mural cells:
[1] 0.2445007

with RBC:
[1] 0.240354


#Just subject variables:

for(i in c(1:length(SignalSortedNoNA3[,1]))){
	
	temp<-summary.lm(lm(SignalSortedNoNA3[i,]~BrainPHNoNA3 +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeNoNA3+ GenderNoNA3 + DiagnosisNoNA3))

GeneByCellTypeSubjVar_RSquared[i,1]<-temp$r.squared
GeneByCellTypeSubjVar_RSquared[i,2]<-temp$adj.r.squared

}

mean(GeneByCellTypeSubjVar_RSquared[,1])
[1] 0.1191567
> median(GeneByCellTypeSubjVar_RSquared[,1])
[1] 0.0977273
> sd(GeneByCellTypeSubjVar_RSquared[,1])
[1] 0.08108838


> mean(GeneByCellTypeSubjVar_RSquared[,2])
[1] 0.07154358
> median(GeneByCellTypeSubjVar_RSquared[,2])
[1] 0.0489558
> sd(GeneByCellTypeSubjVar_RSquared[,2])
[1] 0.08547153

png("Histogram_Rsquared_SubjVar.png")
hist(GeneByCellTypeSubjVar_RSquared[,1], col=2, xlab="R-squared", xlim=c(0,1), main="Linear Model with Just Subject Variables")
dev.off()

png("Histogram_AdjRsquared_SubjVar.png")
hist(GeneByCellTypeSubjVar_RSquared[,2], col=3, xlab="Adjusted R-squared", xlim=c(0,1), main="Linear Model with Just Subject Variables")
dev.off()

#Subj Var & 3 cell types:
for(i in c(1:length(SignalSortedNoNA3[,1]))){
	
	temp<-summary.lm(lm(SignalSortedNoNA3[i,]~CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap$Astrocyte+CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap$Neuron_All+CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap$Oligodendrocyte+BrainPHNoNA3 +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeNoNA3+ GenderNoNA3 + DiagnosisNoNA3))

GeneByCellTypeSubjVar_RSquared[i,1]<-temp$r.squared
GeneByCellTypeSubjVar_RSquared[i,2]<-temp$adj.r.squared

}

#Not as good as if I just added projection neurons & microglia to the model:
> mean(GeneByCellTypeSubjVar_RSquared[,1])
[1] 0.2340306
> median(GeneByCellTypeSubjVar_RSquared[,1])
[1] 0.1862999
> sd(GeneByCellTypeSubjVar_RSquared[,1])
[1] 0.1610512

> 
> mean(GeneByCellTypeSubjVar_RSquared[,2])
[1] 0.1759226
> median(GeneByCellTypeSubjVar_RSquared[,2])
[1] 0.1245709
> sd(GeneByCellTypeSubjVar_RSquared[,2])
[1] 0.1732689

png("Histogram_Rsquared_SubjVar3celltypes.png")
hist(GeneByCellTypeSubjVar_RSquared[,1], col=2, xlab="R-squared", main="Linear Model with just subject variables & 3 cell types")
dev.off()

png("Histogram_AdjRsquared_SubjVar3celltypes.png")
hist(GeneByCellTypeSubjVar_RSquared[,2], col=3, xlab="Adjusted R-squared", main="Linear Model with just subject variables & 3 cell types")

