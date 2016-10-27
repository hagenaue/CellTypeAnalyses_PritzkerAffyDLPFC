#This is my attempt to pull together a more organized code document specifically related to cell type index analyses using DLPFC Affymetrix data from the Pritzker foundation
##Megan Hagenauer, October 27 2016

#**********************************************


#Do the cell type indices from different publications still correlate if we remove overlap?

#First, calculating overlap:
CellTypeSpecificGenes_Master2Best_Overlap<-matrix(0, length(table(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[,13])), length(table(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[,13])) )
colnames(CellTypeSpecificGenes_Master2Best_Overlap)<-names(table(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[,13]))
row.names(CellTypeSpecificGenes_Master2Best_Overlap)<-names(table(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[,13]))

for(i in 1: length(table(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[,13]))){
	for(j in 1: length(table(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[,13]))){

CellTypeSpecificGenes_Master2Best_Overlap[i,j]<-sum(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[,13]==names(table(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[,13])[i]), 4]%in%CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[,13]==names(table(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[,13])[j]), 4])/length(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[,13]==names(table(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[,13])[i]), 4])

}
}

write.csv(CellTypeSpecificGenes_Master2Best_Overlap, "CellTypeSpecificGenes_Master2Best_Overlap.csv")

#RemovingOverlap:

CellTypeSpecificGenes_Master2Best_NoOverlap<-matrix(0, length(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[,1]), length(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[1,]))
dim(CellTypeSpecificGenes_Master2Best_NoOverlap)

colnames(CellTypeSpecificGenes_Master2Best_NoOverlap)<-colnames(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest)
row.names(CellTypeSpecificGenes_Master2Best_NoOverlap)<-row.names(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest)

CellTypeSpecificGenes_Master2Best_OverlapbyGene<-matrix(0, length(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[,1]), 1)
row.names(CellTypeSpecificGenes_Master2Best_OverlapbyGene)<-row.names(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest)

for(i in 1:length(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[,1])){
CellTypeSpecificGenes_Master2Best_OverlapbyGene[i,1]<-sum(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[,4]%in%CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[i,4])	
}


plot.new()

png("CellTypeSpecificGenes_Master2Best_OverlapbyGene.png")
plot(sort(CellTypeSpecificGenes_Master2Best_OverlapbyGene))
dev.off()

#Removing overlap:


CellTypeSpecificGenes_Master2Best_NoOverlap<-CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[CellTypeSpecificGenes_Master2Best_OverlapbyGene[,1]==1,]
dim(CellTypeSpecificGenes_Master2Best_NoOverlap)
[1] 1393  172

CellTypeSpecificGenes_Master2Best_NoOverlapTable<-table(CellTypeSpecificGenes_Master2Best_NoOverlap[,13])
write.csv(CellTypeSpecificGenes_Master2Best_NoOverlapTable, "CellTypeSpecificGenes_Master2Best_NoOverlapTable.csv")


CellTypeSpecificGenes_Master2Best_NoOverlap_Mean<-matrix(0, nrow=length(table(CellTypeSpecificGenes_Master2Best_NoOverlap[,13])), ncol=length(CellTypeSpecificGenes_Master2Best_NoOverlap[1, c(1:172)]))

for(i in c(1:172)){
CellTypeSpecificGenes_Master2Best_NoOverlap_Mean[,i]<-tapply(CellTypeSpecificGenes_Master2Best_NoOverlap[,i], CellTypeSpecificGenes_Master2Best_NoOverlap[,13], mean)
}

head(CellTypeSpecificGenes_Master2Best_NoOverlap_Mean)
row.names(CellTypeSpecificGenes_Master2Best_NoOverlap_Mean)<-names(table(CellTypeSpecificGenes_Master2Best_NoOverlap[,13]))
colnames(CellTypeSpecificGenes_Master2Best_NoOverlap_Mean)<-colnames(CellTypeSpecificGenes_Master2Best_NoOverlap)

CellTypeSpecificGenes_Master2Best_NoOverlap_Mean<-CellTypeSpecificGenes_Master2Best_NoOverlap_Mean[,c(15:171)]


CellTypeSpecificGenes_Master2_BestIndexNoOverlapCorrMatrix<-cor(t(CellTypeSpecificGenes_Master2Best_NoOverlap_Mean))
write.csv(CellTypeSpecificGenes_Master2_BestIndexNoOverlapCorrMatrix, "CellTypeSpecificGenes_Master2_BestIndexNoOverlapCorrMatrix.csv")


