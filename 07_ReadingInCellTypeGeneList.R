
#This is my attempt to pull together a more organized code document specifically related to cell type index analyses using DLPFC Affymetrix data from the Pritzker foundation
##Megan Hagenauer, October 27 2016

#**********************************************
#This code uses an earlier version of the cell type specific gene database (Master2) that included some less useful cell type specific gene lists (e.g., gene lists from fetal cells, ependymal cells)
#Later code (and figure preparation) removes these less-useful cell types

#**********************************************
#Reading in the master database of cell type specific genes and running some basic characterization:

#This is an earlier version of the database that included some less useful cell type specific gene lists (e.g., gene lists from fetal cells, ependymal cells):
CellTypeSpecificGenes_Master2<-read.csv("CellTypeSpecificGenes_Master2.csv", header=T,  na.strings = "#N/A")
dim(CellTypeSpecificGenes_Master2)
#[1] 4504   13

colnames(CellTypeSpecificGenes_Master2)
 # [1] "Umbrella.Cell.Type"    "Specific.Cell.Type"    "Brain.Region"          "Gene.Symbol..Human."   "Gene.Symbol..Mouse."  
 # [6] "Species"               "Age"                   "Statistical.Criterion" "Specificity"           "Comparison"           
# [11] "Platform"              "Citation"              "Tag"    

sum(is.na(CellTypeSpecificGenes_Master2[,4]))
#[1] 464
sum(is.na(CellTypeSpecificGenes_Master2[,5]))
#[1] 19

head(CellTypeSpecificGenes_Master2)

table(CellTypeSpecificGenes_Master2$Umbrella.Cell.Type, CellTypeSpecificGenes_Master2$Specific.Cell.Type)

#Remove genes that have no human equivalent:

CellTypeSpecificGenes_Master2_NoNA<-CellTypeSpecificGenes_Master2[is.na(CellTypeSpecificGenes_Master2[,4])==F,]
dim(CellTypeSpecificGenes_Master2_NoNA)
#[1] 4040   13

str(CellTypeSpecificGenes_Master2_NoNA)

CellTypeSpecificGenes_Master2_NoNA[,4]<-as.character(CellTypeSpecificGenes_Master2_NoNA[,4])
CellTypeSpecificGenes_Master2_NoNA[,5]<-as.character(CellTypeSpecificGenes_Master2_NoNA[,5])
str(CellTypeSpecificGenes_Master2_NoNA)

library(plyr)

SignalSortedNoNA3DF<-data.frame(GeneNames, SignalSortedNoNA3 )
str(SignalSortedNoNA3DF)
SignalSortedNoNA3DF[,2]<-as.character(SignalSortedNoNA3DF[,2])

colnames(SignalSortedNoNA3DF)[2]<-"Gene.Symbol..Human."

CellTypeSpecificGenes_Master2_SignalNoNA3<-join(CellTypeSpecificGenes_Master2_NoNA, SignalSortedNoNA3DF, by="Gene.Symbol..Human.", type="inner", match="all")

dim(CellTypeSpecificGenes_Master2_SignalNoNA3)
#[1] 2678  171

colnames(CellTypeSpecificGenes_Master2_SignalNoNA3)
str(CellTypeSpecificGenes_Master2_SignalNoNA3)

table(CellTypeSpecificGenes_Master2_SignalNoNA3[,1])

      # Astrocyte     Endothelial       Ependymal    FetalGroup10     FetalGroup9       Microglia           Mural          Neuron 
            # 243             313             215              17              17             334             158             927 
# Oligodendrocyte             RBC   UnknownGroup2   UnknownGroup6 
            # 417              11               9              17 
            
CellTypeSpecificGenes_Master2_table<-table(CellTypeSpecificGenes_Master2_SignalNoNA3[,1], CellTypeSpecificGenes_Master2_SignalNoNA3[,2])
write.csv(CellTypeSpecificGenes_Master2_table, "CellTypeSpecificGenes_Master2_table.csv")

#This converts the signal data into z-scores:
temp<-t(scale(t(CellTypeSpecificGenes_Master2_SignalNoNA3[,c(15:171)]), center=T, scale=T))

CellTypeSpecificGenes_Master2_SignalNoNA3_Norm<-data.frame(CellTypeSpecificGenes_Master2_SignalNoNA3[,c(1:14)], temp)

plot(sort(apply(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm[,c(15:171)], 1, mean)))
plot(sort(apply(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm[,c(15:171)], 1, sd)))
#Looks like the normalization worked fine. :)


#This averages all of the normalized signal values for each publication-specific list of cell type specific genes to create teh cell type indices:
CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean<-matrix(0, nrow=length(table(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm[,13])), ncol=length(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm[1, c(15:171)]))
dim(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean)
#[1]  58 157

CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean<-matrix(0, nrow=length(table(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm[,13])), ncol=length(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm[1, c(1:171)]))

for(i in c(1:171)){
CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[,i]<-tapply(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm[,i], CellTypeSpecificGenes_Master2_SignalNoNA3_Norm[,13], mean)}
}

head(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean)
row.names(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean)<-names(table(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm[,13]))
colnames(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean)<-colnames(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm)

str(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean)


CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean<-CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[,c(15:171)]


#This creates a correlation matrix comparing all of the publication-specific cell type indices:

CellTypeSpecificGenes_Master2_IndexCorrelationMatrix<-cor(t(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean))
write.csv(CellTypeSpecificGenes_Master2_IndexCorrelationMatrix, "CellTypeSpecificGenes_Master2_IndexCorrelationMatrix.csv")

#This outputs a hierarchically-clustered heatmap illustrating the correlation matrix:

png("CellTypeSpecificGenes_Master2_IndexHeatmap.png", width=1000, height=1000)
heatmap(CellTypeSpecificGenes_Master2_IndexCorrelationMatrix)
dev.off()

row.names(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean)


#This outputs a consensus-clustering of the cell type indices:

library(ConsensusClusterPlus)

ClusterMasterCellTypeIndices<-ConsensusClusterPlus(t(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean), maxK=10, reps=50, clusterAlg="km", title="kMeans_ConsCluster_MasterCellTypeIndices", distance="euclidean", plot="png", writeTable=TRUE)

#What if I try re-scaling it again (since the averaging made some of the indices have larger ranges than others):
CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_MeanScaled<-t(scale(t(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean), center=T, scale=T))

ClusterMasterCellTypeIndicesScaled<-ConsensusClusterPlus(t(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_MeanScaled), maxK=15, reps=50, clusterAlg="km", title="kMeans_ConsCluster_MasterCellTypeIndicesScaled", distance="euclidean", plot="png", writeTable=TRUE)

#Clustering by gene doesn't seem to be working:
ClusterMasterCellTypeGenes<-ConsensusClusterPlus(t(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm[,c(15:171)]), maxK=10, reps=50, clusterAlg="km", title="kMeans_ConsCluster_MasterCellTypeGenes", distance="euclidean", plot="png", writeTable=TRUE)

colnames(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm)
  # [1] "Umbrella.Cell.Type"    "Specific.Cell.Type"    "Brain.Region"          "Gene.Symbol..Human."   "Gene.Symbol..Mouse."  
  # [6] "Species"               "Age"                   "Statistical.Criterion" "Specificity"           "Comparison"           
 # [11] "Platform"              "Citation"              "Tag"                   "Probe"                 "X1834"  ...
 
 
 CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_MainCat<-CellTypeSpecificGenes_Master2_SignalNoNA3_Norm[CellTypeSpecificGenes_Master2_SignalNoNA3_Norm[,2]=="All",]
 dim(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_MainCat)
#[1] 1660  171
#Hmm... that misses some of the mural cells, as well as some of the nice large broad categories of neurons.



#This is the code for throwing out all of the less-useful cell type indices (e.g., fetal cell indices, ependymal indices)
CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest<-join(as.data.frame(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm), as.data.frame(ListOfBestCellTypeIndices), by="Tag", type="inner")

library(plyr)



#*********************************************************
#Throwing out all of the less useful cell type specific gene lists (e.g., gene lists from fetal cells, ependymal cells):


ListOfBestCellTypeIndices<-read.csv("ListOfBestCellTypeIndices_wCategory.csv", header=T)
str(ListOfBestCellTypeIndices)

CellTypeGenesMaster_TableTag<-table(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm[,13])
write.csv(CellTypeGenesMaster_TableTag, "CellTypeGenesMaster_TableTag.csv")


CellTypeSpecificGenes_Master3<-join(CellTypeSpecificGenes_Master2, as.data.frame(ListOfBestCellTypeIndices), by="Tag", type="inner")

write.csv(CellTypeSpecificGenes_Master3, "CellTypeSpecificGenes_Master3.csv")


#This code characterizes the distribution of cell type specific genes that were included in the trimmed-down database, or that are included in our actual signal data:

CellTypeSpecificGenes_Master3_SignalNoNA3_Norm<-CellTypeSpecificGenes_Master2_SignalNoNA3_Norm[CellTypeSpecificGenes_Master2_SignalNoNA3_Norm[,13]%in%CellTypeSpecificGenes_Master3[,13],]

write.csv(CellTypeSpecificGenes_Master3_SignalNoNA3_Norm, "CellTypeSpecificGenes_Master3_SignalNoNA3_Norm.csv")

CellTypeSpecificGenes_Master3_PrimaryTable<-table(CellTypeSpecificGenes_Master3[,14])
write.csv(CellTypeSpecificGenes_Master3_PrimaryTable, "CellTypeSpecificGenes_Master3_PrimaryTable.csv")

CellTypeSpecificGenes_Master3_IndicesTable<-table(CellTypeSpecificGenes_Master3[,13])
write.csv(CellTypeSpecificGenes_Master3_IndicesTable, "CellTypeSpecificGenes_Master3_IndicesTable.csv")

CellTypeSpecificGenes_Master3_SourceTable<-table(CellTypeSpecificGenes_Master3[,14], CellTypeSpecificGenes_Master3[,12])
write.csv(CellTypeSpecificGenes_Master3_SourceTable, "CellTypeSpecificGenes_Master3_SourceTable.csv")


dim(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest)
#[1] 2066  172

CellTypeSpecificGenes_Master2_PrimaryTable<-table(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[,172])

               # Astrocyte              Endothelial                Microglia                    Mural               Neuron_All 
                     # 243                      313                      334                      158                       90 
      # Neuron_Interneuron        Neuron_Projection          Oligodendrocyte Oligodendrocyte_Immature                      RBC 
                     # 260                      240                      359                       58                       11 

write.csv(CellTypeSpecificGenes_Master2_PrimaryTable, "CellTypeSpecificGenes_Master2_PrimaryTable.csv")


#This code recreates the cell type indices just for the trimmed-down database:

CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_Mean<-matrix(0, nrow=length(table(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[,13])), ncol=length(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[1, c(1:172)]))

for(i in c(1:172)){
CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_Mean[,i]<-tapply(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[,i], CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[,13], mean)}
}


head(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_Mean)
row.names(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_Mean)<-names(table(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[,13]))
colnames(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_Mean)<-colnames(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest)

str(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_Mean)

#This code recreates table of cell type specific genes for each primary cell type for each publication using just the trimmed down table:

CellTypeSpecificGenes_Master2_IndexBestTable<-table(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest[,13])
write.csv(CellTypeSpecificGenes_Master2_IndexBestTable, "CellTypeSpecificGenes_Master2_IndexBestTable.csv")


head(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_Mean)

CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_Mean<-CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_Mean[,c(15:171)]

#This code recreates the cell type indice correlation matrix using just for the trimmed-down database:

CellTypeSpecificGenes_Master2_BestIndexCorrelationMatrix<-cor(t(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_Mean))
write.csv(CellTypeSpecificGenes_Master2_BestIndexCorrelationMatrix, "CellTypeSpecificGenes_Master2_BestIndexCorrelationMatrix.csv")
#All of the tags that don't exist anymore still show up in the correlation matrix as NA :(


sum(is.na(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_Mean[,1])==T)
CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_MeanNoNA<-CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_Mean[is.na(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_Mean[,1])==F,]

#This code recreates the consensus clustering of the cell type indices using just for the trimmed-down database:

ClusterMasterCellTypeIndicesBest<-ConsensusClusterPlus(t(CellTypeSpecificGenes_Master2_SignalNoNA3_NormBest_MeanNoNA), maxK=10, reps=50, clusterAlg="km", title="kMeans_ConsCluster_MasterCellTypeIndicesBest", distance="euclidean", plot="png", writeTable=TRUE)

#****************************************************************
#April 20, 2016:
#Outputting some additional figures for the paper

row.names(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean)


png("AstrocyteIndices_CahoyVsDamanis.png")
plot(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[2,]~ CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[1,], xlab="Astrocyte_All_Cahoy_JNeuro_2008", ylab="Astrocyte_All_Darmanis_PNAS_2015", cex.lab=1.4)
BestFitLine<-lm(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[2,]~ CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[1,])
abline(BestFitLine, col=2, lwd=2)
mtext(paste("P-value: ", format(summary.lm(BestFitLine)$coefficients[2,4], digits=3, scientific=T), ", R-squared: ", format(summary.lm(BestFitLine)$r.squared, digits=3), sep=""), cex=1.4)
dev.off()

png("OligodendrocyteIndices_ZeiselVsDamanis.png")
plot(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[50,]~ CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[49,], xlab="Oligodendrocyte_All_Zeisel_Science_2015", ylab="Oligodendrocyte_Mature_Darmanis_PNAS_2015", cex.lab=1.4)
BestFitLine<-lm(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[50,]~ CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[49,])
abline(BestFitLine, col=2, lwd=2)
mtext(paste("P-value: ", format(summary.lm(BestFitLine)$coefficients[2,4], digits=3, scientific=T), ", R-squared: ", format(summary.lm(BestFitLine)$r.squared, digits=3), sep=""), cex=1.4)
dev.off()

png("EndothelialIndices_ZeiselVsDamarnis.png")
plot(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[8,]~ CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[7,], xlab="Endothelial_All_Zeisel_Science_2015", ylab="Endothelial_All_Darmanis_PNAS_2015", cex.lab=1.4)
BestFitLine<-lm(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[8,]~ CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[7,])
abline(BestFitLine, col=2, lwd=2)
mtext(paste("P-value: ", format(summary.lm(BestFitLine)$coefficients[2,4], digits=3, scientific=T), ", R-squared: ", format(summary.lm(BestFitLine)$r.squared, digits=3), sep=""), cex=1.4)
dev.off()

png("MicrogliaIndices_ZeiselVsDamarnis.png")
plot(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[14,]~ CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[13,], xlab="Microglia_All_Zeisel_Science_2015", ylab="Microglia_All_Darmanis_PNAS_2015", cex.lab=1.4)
BestFitLine<-lm(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[14,]~ CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[13,])
abline(BestFitLine, col=2, lwd=2)
mtext(paste("P-value: ", format(summary.lm(BestFitLine)$coefficients[2,4], digits=3, scientific=T), ", R-squared: ", format(summary.lm(BestFitLine)$r.squared, digits=3), sep=""), cex=1.4)
dev.off()

png("NeuronIndices_CahoyVsDarmanis.png")
plot(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[27,]~ CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[26,], xlab="Neuron_All_Cahoy_JNeuro_2008", ylab="Neuron_All_Darmanis_PNAS_2015", cex.lab=1.4)
BestFitLine<-lm(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[27,]~ CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[26,])
abline(BestFitLine, col=2, lwd=2)
mtext(paste("P-value: ", format(summary.lm(BestFitLine)$coefficients[2,4], digits=3, scientific=T), ", R-squared: ", format(summary.lm(BestFitLine)$r.squared, digits=3), sep=""), cex=1.4)
dev.off()


png("NeuronvsAstrocyteIndices_DarmanisVsZeisel.png")
plot(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[27,]~ CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[4,], xlab="Astrocyte_All_Zeisel_Science_2015", ylab="Neuron_All_Darmanis_PNAS_2015", cex.lab=1.4)
BestFitLine<-lm(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[27,]~ CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[4,])
abline(BestFitLine, col=4, lwd=2)
mtext(paste("P-value: ", format(summary.lm(BestFitLine)$coefficients[2,4], digits=3, scientific=T), ", R-squared: ", format(summary.lm(BestFitLine)$r.squared, digits=3), sep=""), cex=1.4)
dev.off()

png("NeuronvsEndothelialIndices_CahoyVsZeisel.png")
plot(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[26,]~ CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[8,], xlab="Endothelial_All_Zeisel_Science_2015", ylab="Neuron_All_Cahoy_JNeuro_2008", cex.lab=1.4)
BestFitLine<-lm(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[26,]~ CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[8,])
abline(BestFitLine, col=4, lwd=2)
mtext(paste("P-value: ", format(summary.lm(BestFitLine)$coefficients[2,4], digits=3, scientific=T), ", R-squared: ", format(summary.lm(BestFitLine)$r.squared, digits=3), sep=""), cex=1.4)
dev.off()


#Making these figures higher resolution:


pdf("AstrocyteIndices_CahoyVsDamanis.pdf", width=4, height=4.5, pointsize=10)
plot(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[2,]~ CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[1,], xlab="Astrocyte_All_Cahoy_JNeuro_2008", ylab="Astrocyte_All_Darmanis_PNAS_2015", col=1, font.lab=2, lwd=1, cex.lab=1.1, cex.axis=1)
RegressionLine<-lm(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[2,]~ CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[1,])
abline(RegressionLine, col="red", lwd=4)
mtext(paste("p-value = ", print(formatC(signif(summary.lm(RegressionLine)$coefficients[8],digits=3), digits=3, flag="#")), sep=""), line=0.5, cex=1.5, font=2)
mtext(paste("r-squared = ", round(summary.lm(RegressionLine)[8][[1]], digits=2), sep=""), line=2, cex=1.5, font=2)
dev.off()

pdf("OligodendrocyteIndices_ZeiselVsDamanis.pdf", width=4, height=4.5, pointsize=10)
plot(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[50,]~ CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[49,], xlab="Oligodendrocyte_All_Zeisel_Science_2015", ylab="Oligodendrocyte_Mature_Darmanis_PNAS_2015", col=1, font.lab=2, lwd=1, cex.lab=1.1, cex.axis=1)
RegressionLine<-lm(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[50,]~ CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[49,])
abline(RegressionLine, col="red", lwd=4)
mtext(paste("p-value = ", print(formatC(signif(summary.lm(RegressionLine)$coefficients[8],digits=3), digits=3, flag="#")), sep=""), line=0.5, cex=1.5, font=2)
mtext(paste("r-squared = ", round(summary.lm(RegressionLine)[8][[1]], digits=2), sep=""), line=2, cex=1.5, font=2)
dev.off()

pdf("NeuronvsAstrocyteIndices_DarmanisVsZeisel.pdf", width=4, height=4.5, pointsize=10)
plot(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[27,]~ CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[4,], xlab="Astrocyte_All_Zeisel_Science_2015", ylab="Neuron_All_Darmanis_PNAS_2015", col=1, font.lab=2, lwd=1, cex.lab=1.1, cex.axis=1)
RegressionLine<-lm(CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[27,]~ CellTypeSpecificGenes_Master2_SignalNoNA3_Norm_Mean[4,])
abline(RegressionLine, col="dodgerblue3", lwd=4)
mtext(paste("p-value = ", print(formatC(signif(summary.lm(RegressionLine)$coefficients[8],digits=3), digits=3, flag="#")), sep=""), line=0.5, cex=1.5, font=2)
mtext(paste("r-squared = ", round(summary.lm(RegressionLine)[8][[1]], digits=2), sep=""), line=2, cex=1.5, font=2)
dev.off()






