library(fgsea)

setwd("~/Documents/Affy/NoPC1correct/DLPFC circadian/fGSEA/GMTfiles_Human")


temp<-gmtPathways("CombinedGMT_C2_C5_JustTraditional_CellType_20170928.gmt.txt")
names(temp)
temp[[1]]


setwd("~/Documents/Affy/NoPC1correct/DLPFC circadian/fGSEA")


Age_betas_forGSEA<-read.csv("Age_betas_forGSEAmoredecimals.csv", header=T, stringsAsFactors = F)
colnames(Age_betas_forGSEA)<-c()

Age_betas_forGSEA_asVector<-Age_betas_forGSEA[,2]
names(Age_betas_forGSEA_asVector)<-Age_betas_forGSEA[,1]

temp1<-fgsea(temp, Age_betas_forGSEA_asVector, nperm=2000, minSize = 15, maxSize = 500)
str(temp1)
#Worked.

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "Age_betas_GSEA_Results.csv")



AgonalFactor_betas_forGSEA<-read.csv("AgonalFactor_betas_forGSEA.csv", header=T, stringsAsFactors = F)
colnames(AgonalFactor_betas_forGSEA)<-c()

AgonalFactor_betas_forGSEA_asVector<-AgonalFactor_betas_forGSEA[,2]
names(AgonalFactor_betas_forGSEA_asVector)<-AgonalFactor_betas_forGSEA[,1]

temp1<-fgsea(temp, AgonalFactor_betas_forGSEA_asVector, nperm=2000, minSize = 15, maxSize = 500)
str(temp1)
#Worked.

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "AgonalFactor_betas_GSEA_Results.csv")


BP_betas_forGSEA<-read.csv("BP_betas_forGSEA.csv", header=T, stringsAsFactors = F)
colnames(BP_betas_forGSEA)<-c()

BP_betas_forGSEA_asVector<-BP_betas_forGSEA[,2]
names(BP_betas_forGSEA_asVector)<-BP_betas_forGSEA[,1]

temp1<-fgsea(temp, BP_betas_forGSEA_asVector, nperm=2000, minSize = 15, maxSize = 500)
str(temp1)
#Worked.

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "BP_betas_GSEA_Results.csv")


BrainPH_betas_forGSEA<-read.csv("BrainPH_betas_forGSEA.csv", header=T, stringsAsFactors = F)
colnames(BrainPH_betas_forGSEA)<-c()

BrainPH_betas_forGSEA_asVector<-BrainPH_betas_forGSEA[,2]
names(BrainPH_betas_forGSEA_asVector)<-BrainPH_betas_forGSEA[,1]

temp1<-fgsea(temp, BrainPH_betas_forGSEA_asVector, nperm=2000, minSize = 15, maxSize = 500)
str(temp1)
#Worked.

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "BrainPH_betas_GSEA_Results.csv")


Gender_betas_forGSEA<-read.csv("Gender_betas_forGSEA.csv", header=T, stringsAsFactors = F)
colnames(Gender_betas_forGSEA)<-c()

Gender_betas_forGSEA_asVector<-Gender_betas_forGSEA[,2]
names(Gender_betas_forGSEA_asVector)<-Gender_betas_forGSEA[,1]

temp1<-fgsea(temp, Gender_betas_forGSEA_asVector, nperm=2000, minSize = 15, maxSize = 500)
str(temp1)
#Worked.

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "Gender_betas_GSEA_Results.csv")


MDD_betas_forGSEA<-read.csv("MDD_betas_forGSEA.csv", header=T, stringsAsFactors = F)
colnames(MDD_betas_forGSEA)<-c()

MDD_betas_forGSEA_asVector<-MDD_betas_forGSEA[,2]
names(MDD_betas_forGSEA_asVector)<-MDD_betas_forGSEA[,1]

temp1<-fgsea(temp, MDD_betas_forGSEA_asVector, nperm=2000, minSize = 15, maxSize = 500)
str(temp1)
#Worked.

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "MDD_betas_GSEA_Results.csv")


PMI_betas_forGSEA<-read.csv("PMI_betas_forGSEA.csv", header=T, stringsAsFactors = F)
colnames(PMI_betas_forGSEA)<-c()

PMI_betas_forGSEA_asVector<-PMI_betas_forGSEA[,2]
names(PMI_betas_forGSEA_asVector)<-PMI_betas_forGSEA[,1]

temp1<-fgsea(temp, PMI_betas_forGSEA_asVector, nperm=2000, minSize = 15, maxSize = 500)
str(temp1)
#Worked.

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "PMI_betas_GSEA_Results.csv")


Schiz_betas_forGSEA<-read.csv("Schiz_betas_forGSEA.csv", header=T, stringsAsFactors = F)
colnames(Schiz_betas_forGSEA)<-c()

Schiz_betas_forGSEA_asVector<-Schiz_betas_forGSEA[,2]
names(Schiz_betas_forGSEA_asVector)<-Schiz_betas_forGSEA[,1]

temp1<-fgsea(temp, Schiz_betas_forGSEA_asVector, nperm=2000, minSize = 15, maxSize = 500)
str(temp1)
#Worked.

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "Schiz_betas_GSEA_Results.csv")