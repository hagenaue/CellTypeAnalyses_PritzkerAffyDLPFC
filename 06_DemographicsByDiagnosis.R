
#This is my attempt to pull together a more organized code document specifically related to cell type index analyses using DLPFC Affymetrix data from the Pritzker foundation
##Megan Hagenauer, October 27 2016

#**********************************************

#Outputting demographics by diagnosis:

DemographicsByDiagnosis<-matrix(0, 4, 15)
colnames(DemographicsByDiagnosis)<-c("Age", "Age SD", "Age Min", "Age Max" , "Gender Percent F", "Percent Suicide", "pH", "pH SD", "pH min", "pH Max", "PMI", "PMI sd", "PMI min", "PMI max", "Agonal Factor Percent 0")
row.names(DemographicsByDiagnosis)<-c("Control", "MDD", "BP", "SCHIZ")

DemographicsByDiagnosis[,1]<-c(median(AgeNoNA3[DiagnosisNoNA3=="Control"]), median(AgeNoNA3[DiagnosisNoNA3=="MD"]), median(AgeNoNA3[DiagnosisNoNA3=="BP"]), median(AgeNoNA3[DiagnosisNoNA3=="Schiz"]))

DemographicsByDiagnosis[,2]<-c(sd(AgeNoNA3[DiagnosisNoNA3=="Control"]), sd(AgeNoNA3[DiagnosisNoNA3=="MD"]), sd(AgeNoNA3[DiagnosisNoNA3=="BP"]), sd(AgeNoNA3[DiagnosisNoNA3=="Schiz"]))

DemographicsByDiagnosis[,3]<-c(min(AgeNoNA3[DiagnosisNoNA3=="Control"]), min(AgeNoNA3[DiagnosisNoNA3=="MD"]), min(AgeNoNA3[DiagnosisNoNA3=="BP"]), min(AgeNoNA3[DiagnosisNoNA3=="Schiz"]))

DemographicsByDiagnosis[,4]<-c(max(AgeNoNA3[DiagnosisNoNA3=="Control"]), max(AgeNoNA3[DiagnosisNoNA3=="MD"]), max(AgeNoNA3[DiagnosisNoNA3=="BP"]), max(AgeNoNA3[DiagnosisNoNA3=="Schiz"]))

DemographicsByDiagnosis[,5]<-c(sum(GenderNoNA3[DiagnosisNoNA3=="Control"]=="F")/length(GenderNoNA3[DiagnosisNoNA3=="Control"]), sum(GenderNoNA3[DiagnosisNoNA3=="MD"]=="F")/length(GenderNoNA3[DiagnosisNoNA3=="MD"]), sum(GenderNoNA3[DiagnosisNoNA3=="BP"]=="F")/length(GenderNoNA3[DiagnosisNoNA3=="BP"]), sum(GenderNoNA3[DiagnosisNoNA3=="Schiz"]=="F")/length(GenderNoNA3[DiagnosisNoNA3=="Schiz"]))

DemographicsByDiagnosis[,6]<-c(sum(SuicideNoNA3[DiagnosisNoNA3=="Control"]==1)/length(SuicideNoNA3[DiagnosisNoNA3=="Control"]), sum(SuicideNoNA3[DiagnosisNoNA3=="MD"]==1)/length(SuicideNoNA3[DiagnosisNoNA3=="MD"]), sum(SuicideNoNA3[DiagnosisNoNA3=="BP"]==1)/length(SuicideNoNA3[DiagnosisNoNA3=="BP"]), sum(SuicideNoNA3[DiagnosisNoNA3=="Schiz"]==1)/length(SuicideNoNA3[DiagnosisNoNA3=="Schiz"]))

DemographicsByDiagnosis[,7]<-c(median(BrainPHNoNA3[DiagnosisNoNA3=="Control"]), median(BrainPHNoNA3[DiagnosisNoNA3=="MD"]), median(BrainPHNoNA3[DiagnosisNoNA3=="BP"]), median(BrainPHNoNA3[DiagnosisNoNA3=="Schiz"]))

DemographicsByDiagnosis[,8]<-c(sd(BrainPHNoNA3[DiagnosisNoNA3=="Control"]), sd(BrainPHNoNA3[DiagnosisNoNA3=="MD"]), sd(BrainPHNoNA3[DiagnosisNoNA3=="BP"]), sd(BrainPHNoNA3[DiagnosisNoNA3=="Schiz"]))

DemographicsByDiagnosis[,9]<-c(min(BrainPHNoNA3[DiagnosisNoNA3=="Control"]), min(BrainPHNoNA3[DiagnosisNoNA3=="MD"]), min(BrainPHNoNA3[DiagnosisNoNA3=="BP"]), min(BrainPHNoNA3[DiagnosisNoNA3=="Schiz"]))

DemographicsByDiagnosis[,10]<-c(max(BrainPHNoNA3[DiagnosisNoNA3=="Control"]), max(BrainPHNoNA3[DiagnosisNoNA3=="MD"]), max(BrainPHNoNA3[DiagnosisNoNA3=="BP"]), max(BrainPHNoNA3[DiagnosisNoNA3=="Schiz"]))

DemographicsByDiagnosis[,11]<-c(median(HoursFinalCorrectedNoNA3[DiagnosisNoNA3=="Control"]), median(HoursFinalCorrectedNoNA3[DiagnosisNoNA3=="MD"]), median(HoursFinalCorrectedNoNA3[DiagnosisNoNA3=="BP"]), median(HoursFinalCorrectedNoNA3[DiagnosisNoNA3=="Schiz"]))

DemographicsByDiagnosis[,12]<-c(sd(HoursFinalCorrectedNoNA3[DiagnosisNoNA3=="Control"]), sd(HoursFinalCorrectedNoNA3[DiagnosisNoNA3=="MD"]), sd(HoursFinalCorrectedNoNA3[DiagnosisNoNA3=="BP"]), sd(HoursFinalCorrectedNoNA3[DiagnosisNoNA3=="Schiz"]))

DemographicsByDiagnosis[,13]<-c(min(HoursFinalCorrectedNoNA3[DiagnosisNoNA3=="Control"]), min(HoursFinalCorrectedNoNA3[DiagnosisNoNA3=="MD"]), min(HoursFinalCorrectedNoNA3[DiagnosisNoNA3=="BP"]), min(HoursFinalCorrectedNoNA3[DiagnosisNoNA3=="Schiz"]))

DemographicsByDiagnosis[,14]<-c(max(HoursFinalCorrectedNoNA3[DiagnosisNoNA3=="Control"]), max(HoursFinalCorrectedNoNA3[DiagnosisNoNA3=="MD"]), max(HoursFinalCorrectedNoNA3[DiagnosisNoNA3=="BP"]), max(HoursFinalCorrectedNoNA3[DiagnosisNoNA3=="Schiz"]))

DemographicsByDiagnosis[,15]<-c(sum(AgonalFactorNoNA3[DiagnosisNoNA3=="Control"]==0)/length(AgonalFactorNoNA3[DiagnosisNoNA3=="Control"]), sum(AgonalFactorNoNA3[DiagnosisNoNA3=="MD"]==0)/length(AgonalFactorNoNA3[DiagnosisNoNA3=="MD"]), sum(AgonalFactorNoNA3[DiagnosisNoNA3=="BP"]==0)/length(AgonalFactorNoNA3[DiagnosisNoNA3=="BP"]), sum(AgonalFactorNoNA3[DiagnosisNoNA3=="Schiz"]==0)/length(AgonalFactorNoNA3[DiagnosisNoNA3=="Schiz"]))

write.csv(DemographicsByDiagnosis, "DemographicsByDiagnosis.csv")

