#This is my attempt to pull together a more organized code document specifically related to cell type index analyses using DLPFC Affymetrix data from the Pritzker foundation
##Megan Hagenauer, October 27 2016

#**********************************************

#Quantifying the overlap between the cell type specific gene lists from different primary categories of cells:

CellTypeSpecificGenes_Master2_Overlap<-matrix(0, length(table(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[,172])), length(table(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[,172])) )
colnames(CellTypeSpecificGenes_Master2_Overlap)<-names(table(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[,172]))
row.names(CellTypeSpecificGenes_Master2_Overlap)<-names(table(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[,172]))

for(i in 1: length(table(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[,172]))){
	for(j in 1: length(table(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[,172]))){

CellTypeSpecificGenes_Master2_Overlap[i,j]<-sum(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[,172]==names(table(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[,172])[i]), 4]%in%CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[,172]==names(table(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[,172])[j]), 4])/length(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[,172]==names(table(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[,172])[i]), 4])

}
}

write.csv(CellTypeSpecificGenes_Master2_Overlap, "CellTypeSpecificGenes_Master2_Overlap.csv")


#What happens if we eliminate overlap between primary categories and then make master indices:

dim(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest)
[1] 2066  172

CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap<-matrix(0, 1, 172)
colnames(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap)<-colnames(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest)

for(i in 1: length(table(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[,172]))){

#Choosing all data for a particular primary cell type:
TempCurrentIndexAllInfo<-CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[,172]==names(table(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[,172])[i]), ] 

#Choosing all of the gene symbols listed as specific to the current primary cell type:
TempCurrentIndex<-CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[,172]==names(table(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[,172])[i]), 4] 

#Choosing all of the gene symbols listed as specific to all other primary cell types:
TempAllOtherIndices<-CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[,172]%in%names(table(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[,172])[-i]), 4]

#Only grabs rows of data w/ gene symbols not found to be specific to other primary cell types:
CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap<-rbind(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap, TempCurrentIndexAllInfo[(TempCurrentIndex%in%TempAllOtherIndices)==F,])

}

dim(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap)

#removing that one dummy row:
CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap<-CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap[-1,]

dim(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap)


write.csv(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap, "CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap.csv")

CellTypeSpecificGenes_Master2_PrimaryTable_NoPrimaryOverlap<-table(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap[,172])

write.csv(CellTypeSpecificGenes_Master2_PrimaryTable_NoPrimaryOverlap, "CellTypeSpecificGenes_Master2_PrimaryTable_NoPrimaryOverlap.csv")


#Creating Consolidated Primary Cell Type Indices:
CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean<-matrix(0, nrow=length(table(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap[,172])), ncol=length(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap[1, c(1:172)]))

row.names(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean)<-names(table(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap[,172]))

colnames(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean)<-colnames(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap)

for(i in c(1:172)){
CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[,i]<-tapply(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap[,i], CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap[,172], mean)}
}

CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean<-CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean[, -c(c(1:14), 172) ]

is.numeric(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean)

write.csv(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean, "CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean.csv")

#Outputting a hierarchically-clustered heatmap illustrating the correlation between those consolidated indices:

heatmap(cor(t(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean))

#Outputting a correlation matrix for those consolidated indices:

CellTypeIndexNoNA3_NormBest_NoPrimaryOverlap_CorrMatrix<-cor(t(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_NoPrimaryOverlap_Mean))

write.csv(CellTypeIndexNoNA3_NormBest_NoPrimaryOverlap_CorrMatrix, "CellTypeIndexNoNA3_NormBest_NoPrimaryOverlap_CorrMatrix.csv")


#**************************************
