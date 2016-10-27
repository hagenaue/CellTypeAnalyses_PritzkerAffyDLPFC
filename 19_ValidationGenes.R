#This is my attempt to pull together a more organized code document specifically related to cell type index analyses using DLPFC Affymetrix data from the Pritzker foundation
##Megan Hagenauer, October 27 2016

#**********************************************

#This code is a pretty serious mess - there are multiple times here where I should have functionalized but instead just ran the code over and over again in different ways using different input/ouput... and then just commented out the code detailing the earlier input/output.  That means that everything is documented here, but not in a format where it can be simply highlighted and rerun - it needs to be picked through.

#Testing out various linear models using a list of validation genes:

GenesStrongGoodDiagnosis2<-read.csv("GenesStrongGoodDiagnosis.csv", stringsAsFactors =F)
dim(GenesStrongGoodDiagnosis2)
#[1] 83  5
colnames(GenesStrongGoodDiagnosis2)
# [1] "Quality.of.Evidence"                     
# [2] "Validation.genes"                        
# [3] "Subject.Variable.with.known.relationship"
# [4] "Cell.type.in.which.the.effect.occurs"    
# [5] "Nature.of.Effect.s." 

colnames(GenesStrongGoodDiagnosis2)[2]<-"Gene.Symbol"

temp<-data.frame(GeneNames,SignalSortedNoNA3)
ValidationGeneSignalNoNA3wAnnotation<-join(GenesStrongGoodDiagnosis2, temp, by="Gene.Symbol", type="inner")

dim(ValidationGeneSignalNoNA3wAnnotation)
#[1]  62 163
colnames(ValidationGeneSignalNoNA3wAnnotation)
# [1] "Quality.of.Evidence"                     
# [2] "Gene.Symbol"                             
# [3] "Subject.Variable.with.known.relationship"
# [4] "Cell.type.in.which.the.effect.occurs"    
# [5] "Nature.of.Effect.s."                     
# [6] "Probe"                                   
# [7] "X1834"                                   
# [8] "X1856" 
# ...

ValidationGeneSignalNoNA3<-ValidationGeneSignalNoNA3wAnnotation[,c(7:163)]
str(ValidationGeneSignalNoNA3)
#'data.frame':	62 obs. of  157 variables:
#   $ X1834: num  0.3666 0.2339 1.0212 0.0925 0.7847 ...
# $ X1856: num  -0.02841 -0.64541 -1.45368 0.21383 -0.00463 ...

row.names(ValidationGeneSignalNoNA3)<-ValidationGeneSignalNoNA3wAnnotation[,2]
head(ValidationGeneSignalNoNA3)
str(ValidationGeneSignalNoNA3)
ValidationGeneSignalNoNA3<-as.matrix(ValidationGeneSignalNoNA3)
str(ValidationGeneSignalNoNA3)
head(ValidationGeneSignalNoNA3)


#*********First attempt: Stepwise Regression for model selection****************


#These are notes about all of the different predictor lists that I explored, followed by the actual stepwise code:

###Note: For any models with interaction terms, it is essential that the input data for continuous variables be centered  - I'm not sure if I already included that code in this set of code documents, so I'm going to add it here again:

AgeCentered<-scale(AgeNoNA3, center=T, scale=F)
BrainPHCentered<-scale(BrainPHNoNA3, center=T, scale=F)



#I did two fwd & backward runs using these predictors, and entitled the output: ValidationGenes_Stepwise_Bwd_WDiagnosis & ValidationGenes_Stepwise_Fwd_WDiagnosis
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 + PsychiatricNoNA3+DiagnosisNoNA3+SuicideNoNA3+Astrocyte + Endothelial + Microglia + Mural + Neuron_All + Neuron_Interneuron + Neuron_Projection + Oligodendrocyte + Oligodendrocyte_Immature + RBC + BrainPHCentered*Astrocyte + BrainPHCentered*Endothelial + BrainPHCentered*Microglia + BrainPHCentered*Mural + BrainPHCentered*Neuron_All + BrainPHCentered*Neuron_Interneuron + BrainPHCentered*Neuron_Projection + BrainPHCentered*Oligodendrocyte + BrainPHCentered*Oligodendrocyte_Immature + BrainPHCentered*RBC + AgeCentered*Astrocyte + AgeCentered*Endothelial + AgeCentered*Microglia + AgeCentered*Mural + AgeCentered*Neuron_All + AgeCentered*Neuron_Interneuron + AgeCentered*Neuron_Projection + AgeCentered*Oligodendrocyte + AgeCentered*Oligodendrocyte_Immature + AgeCentered*RBC + PsychiatricNoNA3*Astrocyte + PsychiatricNoNA3*Endothelial + PsychiatricNoNA3*Microglia + PsychiatricNoNA3*Mural + PsychiatricNoNA3*Neuron_All + PsychiatricNoNA3*Neuron_Interneuron + PsychiatricNoNA3*Neuron_Projection + PsychiatricNoNA3*Oligodendrocyte + PsychiatricNoNA3*Oligodendrocyte_Immature + PsychiatricNoNA3*RBC + DiagnosisNoNA3*Astrocyte + DiagnosisNoNA3*Endothelial + DiagnosisNoNA3*Microglia + DiagnosisNoNA3*Mural + DiagnosisNoNA3*Neuron_All + DiagnosisNoNA3*Neuron_Interneuron + DiagnosisNoNA3*Neuron_Projection + DiagnosisNoNA3*Oligodendrocyte + DiagnosisNoNA3*Oligodendrocyte_Immature + DiagnosisNoNA3*RBC + SuicideNoNA3*Astrocyte + SuicideNoNA3*Endothelial + SuicideNoNA3*Microglia + SuicideNoNA3*Mural + SuicideNoNA3*Neuron_All + SuicideNoNA3*Neuron_Interneuron + SuicideNoNA3*Neuron_Projection + SuicideNoNA3*Oligodendrocyte + SuicideNoNA3*Oligodendrocyte_Immature + SuicideNoNA3*RBC, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 + PsychiatricNoNA3+DiagnosisNoNA3+SuicideNoNA3+Astrocyte + Endothelial + Microglia + Mural + Neuron_All + Neuron_Interneuron + Neuron_Projection + Oligodendrocyte + Oligodendrocyte_Immature + RBC + BrainPHCentered*Astrocyte + BrainPHCentered*Endothelial + BrainPHCentered*Microglia + BrainPHCentered*Mural + BrainPHCentered*Neuron_All + BrainPHCentered*Neuron_Interneuron + BrainPHCentered*Neuron_Projection + BrainPHCentered*Oligodendrocyte + BrainPHCentered*Oligodendrocyte_Immature + BrainPHCentered*RBC + AgeCentered*Astrocyte + AgeCentered*Endothelial + AgeCentered*Microglia + AgeCentered*Mural + AgeCentered*Neuron_All + AgeCentered*Neuron_Interneuron + AgeCentered*Neuron_Projection + AgeCentered*Oligodendrocyte + AgeCentered*Oligodendrocyte_Immature + AgeCentered*RBC + PsychiatricNoNA3*Astrocyte + PsychiatricNoNA3*Endothelial + PsychiatricNoNA3*Microglia + PsychiatricNoNA3*Mural + PsychiatricNoNA3*Neuron_All + PsychiatricNoNA3*Neuron_Interneuron + PsychiatricNoNA3*Neuron_Projection + PsychiatricNoNA3*Oligodendrocyte + PsychiatricNoNA3*Oligodendrocyte_Immature + PsychiatricNoNA3*RBC + DiagnosisNoNA3*Astrocyte + DiagnosisNoNA3*Endothelial + DiagnosisNoNA3*Microglia + DiagnosisNoNA3*Mural + DiagnosisNoNA3*Neuron_All + DiagnosisNoNA3*Neuron_Interneuron + DiagnosisNoNA3*Neuron_Projection + DiagnosisNoNA3*Oligodendrocyte + DiagnosisNoNA3*Oligodendrocyte_Immature + DiagnosisNoNA3*RBC + SuicideNoNA3*Astrocyte + SuicideNoNA3*Endothelial + SuicideNoNA3*Microglia + SuicideNoNA3*Mural + SuicideNoNA3*Neuron_All + SuicideNoNA3*Neuron_Interneuron + SuicideNoNA3*Neuron_Projection + SuicideNoNA3*Oligodendrocyte + SuicideNoNA3*Oligodendrocyte_Immature + SuicideNoNA3*RBC, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))

#Out of curiousity, I then decided to see if the results I got would be different if I didn't include the interaction terms. My expectation was that the models would differ for backwards but not forwards stepwise regresssion (since all but one gene lacked interaction terms in the output from forward stepwise regression). This expectation was met. The output files are called: "ValidationGenes_Stepwise_Bwd_WDiagnosisNoInteract" and "ValidationGenes_Stepwise_Fwd_WDiagnosisNoInteract".
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 + PsychiatricNoNA3+DiagnosisNoNA3+SuicideNoNA3+Astrocyte + Endothelial + Microglia + Mural + Neuron_All + Neuron_Interneuron + Neuron_Projection + Oligodendrocyte + Oligodendrocyte_Immature + RBC, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 + PsychiatricNoNA3+DiagnosisNoNA3+SuicideNoNA3+Astrocyte + Endothelial + Microglia + Mural + Neuron_All + Neuron_Interneuron + Neuron_Projection + Oligodendrocyte + Oligodendrocyte_Immature + RBC, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))

#Next, I wondered about the performance in models in which we lacked cell type variables. As expected, many more genes were related to pH and agonal factor when cell type wasn't included in the model.  Output called: "ValidationGenes_Stepwise_Bwd_WDiagnosisJustConfouds.csv"

AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 + PsychiatricNoNA3+DiagnosisNoNA3+SuicideNoNA3, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 + PsychiatricNoNA3+DiagnosisNoNA3, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))

#Next, I wondered about the performance in models in which we had cell type variables but no confounds - with psychosis or mood as options as well.  Output called: "ValidationGenes_Stepwise_Bwd_WDiagnosisJustConfoudsPsychosisMood.csv" Psychosis and mood don't ever seem to matter."ValidationGenes_Stepwise_Fwd_WDiagnosisJustConfoudsPsychosisMood.csv" - in the forward version, psychosis matters for one gene.
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 + PsychiatricNoNA3+DiagnosisNoNA3+SuicideNoNA3+PsychosisNoNA3+MoodNoNA3, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 + PsychiatricNoNA3+DiagnosisNoNA3+SuicideNoNA3+PsychosisNoNA3+MoodNoNA3, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))


#Let's double check that in the version that includes controlling for cell type: In the forward version, psychosis and mood are included in the model more frequently than diagnosis, although the results are still pretty underwhelming. In the backwards version, the model primarily includes diagnosis and suicide. Interesting. I'm guessing that maybe the forward version has problems with multi-level factor variables and the backwards version has trouble with multicollinearity.
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 + PsychiatricNoNA3+DiagnosisNoNA3+SuicideNoNA3+PsychosisNoNA3+MoodNoNA3+Astrocyte + Endothelial + Microglia + Mural + Neuron_All + Neuron_Interneuron + Neuron_Projection + Oligodendrocyte + Oligodendrocyte_Immature + RBC, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 + PsychiatricNoNA3+DiagnosisNoNA3+SuicideNoNA3+PsychosisNoNA3+MoodNoNA3+Astrocyte + Endothelial + Microglia + Mural + Neuron_All + Neuron_Interneuron + Neuron_Projection + Oligodendrocyte + Oligodendrocyte_Immature + RBC, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))


#Here is the actual stepwise code:

StepwiseBetas<-matrix(0, nrow(ValidationGeneSignalNoNA3), length(AllPredictorsNames))
StepwisePvalues<-matrix(0, nrow(ValidationGeneSignalNoNA3), length(AllPredictorsNames))
StepwiseSE<-matrix(0, nrow(ValidationGeneSignalNoNA3), length(AllPredictorsNames))
StepwiseT<-matrix(0, nrow(ValidationGeneSignalNoNA3), length(AllPredictorsNames))

colnames(StepwiseBetas)<-AllPredictorsNames
colnames(StepwisePvalues)<-AllPredictorsNames
colnames(StepwiseSE)<-AllPredictorsNames
colnames(StepwiseT)<-AllPredictorsNames

system.time(
  for(i in 1:nrow(ValidationGeneSignalNoNA3)){
    StepwiseOutput<-stepwise(ModeltoUseForStepwise, direction="backward/forward")
    
    StepwiseOutputTemp<-data.frame(row.names(summary.lm(lm(StepwiseOutput$model))$coefficients), summary.lm(lm(StepwiseOutput$model))$coefficients)
    colnames(StepwiseOutputTemp)[1]<-"AllPredictorsNames"
    StepwiseOutputTemp2<-join(as.data.frame(AllPredictorsNames), StepwiseOutputTemp, by="AllPredictorsNames")
    
    StepwiseBetas[i,]<-StepwiseOutputTemp2[,2]
    StepwiseSE[i,]<-StepwiseOutputTemp2[,3]
    StepwiseT[i,]<-StepwiseOutputTemp2[,4]
    StepwisePvalues[i,]<-StepwiseOutputTemp2[,5]
  }
)


NumberOfTimesInModels<-apply(StepwisePvalues, 2, function(y) length(y)-sum(is.na(y)))
names(NumberOfTimesInModels)<-colnames(StepwisePvalues)

NumberOfTimesSig<-apply(StepwisePvalues, 2, function(y) sum(y<0.05, na.rm=T))
names(NumberOfTimesSig)<-colnames(StepwisePvalues)

temp<-data.frame(ValidationGeneSignalNoNA3wAnnotation[,c(1:6)], StepwiseBetas, StepwisePvalues, StepwiseT, StepwiseSE)

#And here are the various outputs: (currently commented out)

#write.csv(temp, "ValidationGenes_Stepwise_Fwd_WDiagnosis.csv")
#write.csv(temp, "ValidationGenes_Stepwise_Bwd_WDiagnosis.csv")
#write.csv(temp, "ValidationGenes_Stepwise_Bwd_WDiagnosisNoInteract.csv")
#write.csv(temp, "ValidationGenes_Stepwise_Fwd_WDiagnosisNoInteract.csv")
#write.csv(temp, "ValidationGenes_Stepwise_Bwd_WDiagnosisJustConfouds.csv")
#write.csv(temp, "ValidationGenes_Stepwise_Bwd_WDiagnosisJustConfoudsPsychosisMood.csv")
#write.csv(temp, "ValidationGenes_Stepwise_Fwd_WDiagnosisJustConfoudsPsychosisMood.csv")
write.csv(temp, "ValidationGenes_Stepwise_Fwd_WDiagnosisPsychosisMood.csv")

#write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_Stepwise_Fwd_WDiagnosis_CountInModelAndSig.csv")
#write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_Stepwise_Bwd_WDiagnosis_CountInModelAndSig.csv")
#write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_Stepwise_Bwd_WDiagnosisNoInteract_CountInModelAndSig.csv")
#write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_Stepwise_Fwd_WDiagnosisNoInteract_CountInModelAndSig.csv")
#write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_Stepwise_Bwd_WDiagnosisJustConfounds_CountInModelAndSig.csv")
#write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_Stepwise_Bwd_WDiagnosisJustConfoundsPsychosisMood_CountInModelAndSig.csv")
#write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_Stepwise_Fwd_WDiagnosisJustConfoundsPsychosisMood_CountInModelAndSig.csv")
#write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_Stepwise_Fwd_WDiagnosisPsychosisMood_CountInModelAndSig.csv")
write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_Stepwise_Bwd_WDiagnosisPsychosisMood_CountInModelAndSig.csv")

#O.k., so that was pretty interesting, but stepwise methods are typically frowned upon without including cross-validation (a train/test version). That unfortunately means using half the sample size to test effects. :(

#****************************

#Applying the simplest version of stepwise to the full dataset:

#Full dataset version:
AllPredictorsNames<-names(lm(SignalSortedNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 + PsychiatricNoNA3+DiagnosisNoNA3+SuicideNoNA3+Astrocyte + Endothelial + Microglia + Mural + Neuron_All + Neuron_Interneuron + Neuron_Projection + Oligodendrocyte + Oligodendrocyte_Immature + RBC, data=data.frame(SignalSortedNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-lm(SignalSortedNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 + PsychiatricNoNA3+DiagnosisNoNA3+SuicideNoNA3+Astrocyte + Endothelial + Microglia + Mural + Neuron_All + Neuron_Interneuron + Neuron_Projection + Oligodendrocyte + Oligodendrocyte_Immature + RBC, data=data.frame(SignalSortedNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))


StepwiseBetas<-matrix(0, nrow(SignalSortedNoNA3), length(AllPredictorsNames))
StepwisePvalues<-matrix(0, nrow(SignalSortedNoNA3), length(AllPredictorsNames))
StepwiseSE<-matrix(0, nrow(SignalSortedNoNA3), length(AllPredictorsNames))
StepwiseT<-matrix(0, nrow(SignalSortedNoNA3), length(AllPredictorsNames))

colnames(StepwiseBetas)<-AllPredictorsNames
colnames(StepwisePvalues)<-AllPredictorsNames
colnames(StepwiseSE)<-AllPredictorsNames
colnames(StepwiseT)<-AllPredictorsNames

system.time(
  for(i in 1:nrow(SignalSortedNoNA3)){
    StepwiseOutput<-stepwise(ModeltoUseForStepwise, direction="forward/backward")
    
    StepwiseOutputTemp<-data.frame(row.names(summary.lm(lm(StepwiseOutput$model))$coefficients), summary.lm(lm(StepwiseOutput$model))$coefficients)
    colnames(StepwiseOutputTemp)[1]<-"AllPredictorsNames"
    StepwiseOutputTemp2<-join(as.data.frame(AllPredictorsNames), StepwiseOutputTemp, by="AllPredictorsNames")
    
    StepwiseBetas[i,]<-StepwiseOutputTemp2[,2]
    StepwiseSE[i,]<-StepwiseOutputTemp2[,3]
    StepwiseT[i,]<-StepwiseOutputTemp2[,4]
    StepwisePvalues[i,]<-StepwiseOutputTemp2[,5]
  }
)


NumberOfTimesInModels<-apply(StepwisePvalues, 2, function(y) length(y)-sum(is.na(y)))
names(NumberOfTimesInModels)<-colnames(StepwisePvalues)

NumberOfTimesSig<-apply(StepwisePvalues, 2, function(y) sum(y<0.05, na.rm=T))
names(NumberOfTimesSig)<-colnames(StepwisePvalues)

temp<-data.frame(GeneNames, StepwiseBetas, StepwisePvalues, StepwiseT, StepwiseSE)

write.csv(temp, "AllGenes_Stepwise_Fwd_WDiagnosisNoInteract.csv")
write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "AllGenes_Stepwise_Fwd_WDiagnosisNoInteract_CountInModelAndSig.csv")



#*****************************************************************************************

#Moving on to theory-based models (not stepwise):

#Please note - this code is pretty confusing because I was lazy and recycled some of my stepwise code from above... but these are not actually stepwise models, just regular linear regression models.  So the code says "stepwise" even though it is not.

#Also note that there are more models here than discussed in the paper - I chose some representative ones for the paper.


#Different models explored:


#LMWDiagnosisSuicideCelltypeConfounds
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +DiagnosisNoNA3+SuicideNoNA3+Astrocyte + Endothelial + Microglia + Mural + Neuron_All + Neuron_Interneuron + Neuron_Projection + Oligodendrocyte + Oligodendrocyte_Immature + RBC, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-function(i){lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +DiagnosisNoNA3+SuicideNoNA3+Astrocyte + Endothelial + Microglia + Mural + Neuron_All + Neuron_Interneuron + Neuron_Projection + Oligodendrocyte + Oligodendrocyte_Immature + RBC, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))}

#LMWDiagnosisCelltypeConfounds:
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +DiagnosisNoNA3+Astrocyte + Endothelial + Microglia + Mural + Neuron_All + Neuron_Interneuron + Neuron_Projection + Oligodendrocyte + Oligodendrocyte_Immature + RBC, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-function(i){lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +DiagnosisNoNA3+Astrocyte + Endothelial + Microglia + Mural + Neuron_All + Neuron_Interneuron + Neuron_Projection + Oligodendrocyte + Oligodendrocyte_Immature + RBC, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))}

#LM4:
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +DiagnosisNoNA3, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-function(i){lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +DiagnosisNoNA3, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))}

#LM4wSuicide:
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +DiagnosisNoNA3+SuicideNoNA3, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-function(i){lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +DiagnosisNoNA3+SuicideNoNA3, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))}

#LMWDiagnosisPrevelantCelltypesConfounds:
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +DiagnosisNoNA3+Astrocyte + Microglia +  Neuron_Interneuron + Neuron_Projection + Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-function(i){lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +DiagnosisNoNA3+ Astrocyte + Microglia +  Neuron_Interneuron + Neuron_Projection + Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))}

#LMWDiagnosis2CelltypesConfounds:
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +DiagnosisNoNA3+Astrocyte + Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-function(i){lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +DiagnosisNoNA3+ Astrocyte + Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))}

#LMWDiagnosisGliaCelltypesConfounds:
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +DiagnosisNoNA3+Astrocyte + Oligodendrocyte+ Microglia, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-function(i){lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +DiagnosisNoNA3+ Astrocyte + Oligodendrocyte+ Microglia, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))}


#LMWDiagnosis3CelltypesConfounds:
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +DiagnosisNoNA3+Astrocyte + Oligodendrocyte+ Neuron_All, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-function(i){lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +DiagnosisNoNA3+ Astrocyte + Oligodendrocyte+ Neuron_All, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))}

#Side note: I ran some variance inflation factor calculations to try to determine which of the cell type indices were the most multicollinear
library(car)
vif(lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +DiagnosisNoNA3+Astrocyte + Endothelial + Microglia + Mural + Neuron_All + Neuron_Interneuron + Neuron_Projection + Oligodendrocyte + Oligodendrocyte_Immature + RBC, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap)))

# GVIF Df GVIF^(1/(2*Df))
# BrainPHCentered          1.655267  1        1.286572
# AgonalFactorNoNA3        1.864435  1        1.365443
# HoursFinalCorrectedNoNA3 1.511193  1        1.229306
# AgeCentered              1.366539  1        1.168991
# GenderNoNA3              1.213165  1        1.101438
# DiagnosisNoNA3           1.516908  3        1.071914
# Astrocyte                3.805356  1        1.950732
# Endothelial              7.518195  1        2.741933####### Redundant
# Microglia                2.101712  1        1.449728
# Mural                    3.035078  1        1.742147
# Neuron_All               8.154045  1        2.855529####### Redundant
# Neuron_Interneuron       4.310703  1        2.076223
# Neuron_Projection        6.692459  1        2.586979
# Oligodendrocyte          2.600468  1        1.612597
# Oligodendrocyte_Immature 2.297947  1        1.515898##Probably not relevant - cut?
# RBC                      1.889213  1        1.374486

#Just cell types:
vif(lm(ValidationGeneSignalNoNA3[i,]~Astrocyte + Endothelial + Microglia + Mural + Neuron_All + Neuron_Interneuron + Neuron_Projection + Oligodendrocyte + Oligodendrocyte_Immature + RBC, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap)))

# Astrocyte              Endothelial                Microglia                    Mural               Neuron_All 
# 3.024576                 7.120087                 1.954662                 2.815714                 6.599726 
# Neuron_Interneuron        Neuron_Projection          Oligodendrocyte Oligodendrocyte_Immature                      RBC 
# 3.509660                 5.792388                 2.171777                 2.064843                 1.704132 

#Same conclusion - cut Neuron_All and Endothelial, and oligodendrocyte_immature since it may not be valid:

vif(lm(ValidationGeneSignalNoNA3[i,]~Astrocyte + Microglia + Mural + Neuron_Interneuron + Neuron_Projection + Oligodendrocyte + RBC, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap)))

# Astrocyte          Microglia              Mural Neuron_Interneuron  Neuron_Projection    Oligodendrocyte                RBC 
# 2.351019           1.397098           1.864927           3.315728           3.698774           2.026621           1.667946 
#Much better!

vif(lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +DiagnosisNoNA3+Astrocyte + Microglia + Mural + Neuron_Interneuron + Neuron_Projection + Oligodendrocyte + RBC, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap)))

# GVIF Df GVIF^(1/(2*Df))
# BrainPHCentered          1.538839  1        1.240499
# AgonalFactorNoNA3        1.598681  1        1.264390
# HoursFinalCorrectedNoNA3 1.377544  1        1.173688
# AgeCentered              1.326010  1        1.151525
# GenderNoNA3              1.201177  1        1.095982
# DiagnosisNoNA3           1.390568  3        1.056490
# Astrocyte                3.127530  1        1.768483
# Microglia                1.583904  1        1.258533
# Mural                    1.948214  1        1.395784
# Neuron_Interneuron       4.005375  1        2.001343
# Neuron_Projection        4.303278  1        2.074434
# Oligodendrocyte          2.241298  1        1.497097
# RBC                      1.841856  1        1.357150
#Yep- looking good. Let's try that one.

#LMWDiagnosisVIFCelltypesConfounds:
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +DiagnosisNoNA3+Astrocyte + Microglia + Mural + Neuron_Interneuron + Neuron_Projection + Oligodendrocyte + RBC, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-function(i){lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +DiagnosisNoNA3+ Astrocyte + Microglia + Mural + Neuron_Interneuron + Neuron_Projection + Oligodendrocyte + RBC, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))}


#JustDiagnosis:
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~ DiagnosisNoNA3, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-function(i){lm(ValidationGeneSignalNoNA3[i,]~ DiagnosisNoNA3, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))}

#JustAge:
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~ AgeCentered, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-function(i){lm(ValidationGeneSignalNoNA3[i,]~ AgeCentered, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))}

#JustHypoxia:
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~ BrainPHCentered+AgonalFactorNoNA3, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-function(i){lm(ValidationGeneSignalNoNA3[i,]~ BrainPHCentered+AgonalFactorNoNA3, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))}

#DiagnosisPrevelantCellTypesNoConfounds:
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~ DiagnosisNoNA3+Astrocyte + Microglia +  Neuron_Interneuron + Neuron_Projection + Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-function(i){lm(ValidationGeneSignalNoNA3[i,]~ DiagnosisNoNA3+Astrocyte + Microglia +  Neuron_Interneuron + Neuron_Projection + Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))}

#AgePrevelantCellTypesNoConfounds:
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~ AgeCentered+Astrocyte + Microglia +  Neuron_Interneuron + Neuron_Projection + Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-function(i){lm(ValidationGeneSignalNoNA3[i,]~ AgeCentered+Astrocyte + Microglia +  Neuron_Interneuron + Neuron_Projection + Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))}

#HypoxiaPrevelantCellTypesNoConfounds:
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~ BrainPHCentered+AgonalFactorNoNA3+Astrocyte + Microglia +  Neuron_Interneuron + Neuron_Projection + Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-function(i){lm(ValidationGeneSignalNoNA3[i,]~ BrainPHCentered+AgonalFactorNoNA3+Astrocyte + Microglia +  Neuron_Interneuron + Neuron_Projection + Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))}

#DiagnosisThreeCellTypesNoConfounds:
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~ DiagnosisNoNA3+ Astrocyte +  Neuron_All + Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-function(i){lm(ValidationGeneSignalNoNA3[i,]~ DiagnosisNoNA3+ Astrocyte +  Neuron_All + Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))}

#DiagnosisThreeCellTypesConfounds:
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~ BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +DiagnosisNoNA3+ Astrocyte +  Neuron_All + Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-function(i){lm(ValidationGeneSignalNoNA3[i,]~ BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +DiagnosisNoNA3+ Astrocyte +  Neuron_All + Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))}

#DiagnosisThreeCellTypesInteraction:
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~ DiagnosisNoNA3+ Astrocyte +  Neuron_All + Oligodendrocyte+DiagnosisNoNA3* Astrocyte+DiagnosisNoNA3*Neuron_All+DiagnosisNoNA3*Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-function(i){lm(ValidationGeneSignalNoNA3[i,]~ DiagnosisNoNA3+ Astrocyte +  Neuron_All + Oligodendrocyte+DiagnosisNoNA3* Astrocyte+DiagnosisNoNA3*Neuron_All+DiagnosisNoNA3*Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))}


#DiagnosisNeuronInteraction:
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~ DiagnosisNoNA3+ Astrocyte +  Neuron_All+ Oligodendrocyte +DiagnosisNoNA3*Neuron_All, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-function(i){lm(ValidationGeneSignalNoNA3[i,]~ DiagnosisNoNA3+ Astrocyte +  Neuron_All + Oligodendrocyte+DiagnosisNoNA3*Neuron_All, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))}

#DiagnosisInterNeuronInteraction:
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~ DiagnosisNoNA3+ Astrocyte +  Neuron_Interneuron+ Oligodendrocyte +DiagnosisNoNA3*Neuron_Interneuron, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-function(i){lm(ValidationGeneSignalNoNA3[i,]~ DiagnosisNoNA3+ Astrocyte + Neuron_Interneuron + Oligodendrocyte+DiagnosisNoNA3*Neuron_Interneuron, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))}

#DiagnosisInterNeuronProjection:
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~ DiagnosisNoNA3+ Astrocyte +  Neuron_Projection+ Oligodendrocyte +DiagnosisNoNA3*Neuron_Projection, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-function(i){lm(ValidationGeneSignalNoNA3[i,]~ DiagnosisNoNA3+ Astrocyte + Neuron_Projection + Oligodendrocyte+DiagnosisNoNA3*Neuron_Projection, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))}

#DiagnosisThreeCellTypesConfoundsInteraction:
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~ BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +DiagnosisNoNA3+ Astrocyte +  Neuron_All + Oligodendrocyte+DiagnosisNoNA3* Astrocyte+DiagnosisNoNA3*Neuron_All+DiagnosisNoNA3*Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-function(i){lm(ValidationGeneSignalNoNA3[i,]~ BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +DiagnosisNoNA3+ Astrocyte +  Neuron_All + Oligodendrocyte+DiagnosisNoNA3* Astrocyte+DiagnosisNoNA3*Neuron_All+DiagnosisNoNA3*Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))}

#DiagnosisThreeCellTypesConfoundsInteractionswConfounds:
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~ BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +DiagnosisNoNA3+ Astrocyte +  Neuron_All + Oligodendrocyte+DiagnosisNoNA3*Astrocyte+DiagnosisNoNA3*Neuron_All+DiagnosisNoNA3*Oligodendrocyte+BrainPHCentered*Astrocyte+BrainPHCentered*Neuron_All+BrainPHCentered*Oligodendrocyte+AgeCentered*Astrocyte+AgeCentered*Neuron_All+AgeCentered*Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-function(i){lm(ValidationGeneSignalNoNA3[i,]~ BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +DiagnosisNoNA3+ Astrocyte +  Neuron_All + Oligodendrocyte+DiagnosisNoNA3*Astrocyte+DiagnosisNoNA3*Neuron_All+DiagnosisNoNA3*Oligodendrocyte+BrainPHCentered*Astrocyte+BrainPHCentered*Neuron_All+BrainPHCentered*Oligodendrocyte+AgeCentered*Astrocyte+AgeCentered*Neuron_All+AgeCentered*Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))}


#DiagnosisPrevelantCellTypesNoConfoundsInteractions:
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~ DiagnosisNoNA3+Astrocyte + Microglia +  Neuron_Interneuron + Neuron_Projection + Oligodendrocyte+DiagnosisNoNA3*Astrocyte+DiagnosisNoNA3*Microglia+DiagnosisNoNA3*Neuron_Interneuron+DiagnosisNoNA3*Neuron_Projection+DiagnosisNoNA3*Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-function(i){lm(ValidationGeneSignalNoNA3[i,]~ DiagnosisNoNA3+Astrocyte + Microglia +  Neuron_Interneuron + Neuron_Projection + Oligodendrocyte+DiagnosisNoNA3*Astrocyte+DiagnosisNoNA3*Microglia+DiagnosisNoNA3*Neuron_Interneuron+DiagnosisNoNA3*Neuron_Projection+DiagnosisNoNA3*Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))}


#LMWPsychosisMoodPrevelantCelltypesConfounds:
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +PsychosisNoNA3+ MoodNoNA3+Astrocyte + Microglia +  Neuron_Interneuron + Neuron_Projection + Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-function(i){lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +PsychosisNoNA3+ MoodNoNA3+ Astrocyte + Microglia +  Neuron_Interneuron + Neuron_Projection + Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))}

#LMWPsychiatricSuicidePrevelantCelltypesConfounds:
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3+PsychiatricNoNA3+SuicideNoNA3+Astrocyte + Microglia +  Neuron_Interneuron + Neuron_Projection + Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-function(i){lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +PsychiatricNoNA3+SuicideNoNA3+ Astrocyte + Microglia +  Neuron_Interneuron + Neuron_Projection + Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))}

#LMWPsychiatricPrevelantCelltypesConfounds:
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3+PsychiatricNoNA3+Astrocyte + Microglia +  Neuron_Interneuron + Neuron_Projection + Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-function(i){lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +PsychiatricNoNA3+ Astrocyte + Microglia +  Neuron_Interneuron + Neuron_Projection + Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))}


#LMWPsychiatricPrevelantCelltypes:
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~PsychiatricNoNA3+Astrocyte + Microglia +  Neuron_Interneuron + Neuron_Projection + Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-function(i){lm(ValidationGeneSignalNoNA3[i,]~PsychiatricNoNA3+ Astrocyte + Microglia +  Neuron_Interneuron + Neuron_Projection + Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))}


#LM4WPsychiatric:
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3+PsychiatricNoNA3, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-function(i){lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +PsychiatricNoNA3, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))}

#JustPsychiatric
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~PsychiatricNoNA3, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-function(i){lm(ValidationGeneSignalNoNA3[i,]~PsychiatricNoNA3, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))}

#LMWPsychiatricAllCelltypesConfounds:
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3+PsychiatricNoNA3+Astrocyte + Endothelial + Microglia + Mural + Neuron_All + Neuron_Interneuron + Neuron_Projection + Oligodendrocyte + Oligodendrocyte_Immature + RBC, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-function(i){lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +PsychiatricNoNA3+ Astrocyte + Endothelial + Microglia + Mural + Neuron_All + Neuron_Interneuron + Neuron_Projection + Oligodendrocyte + Oligodendrocyte_Immature + RBC, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))}


#PsychiatricThreeCellTypesConfoundsInteractionswConfounds:
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~ BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +PsychiatricNoNA3+ Astrocyte +  Neuron_All + Oligodendrocyte+PsychiatricNoNA3*Astrocyte+PsychiatricNoNA3*Neuron_All+PsychiatricNoNA3*Oligodendrocyte+BrainPHCentered*Astrocyte+BrainPHCentered*Neuron_All+BrainPHCentered*Oligodendrocyte+AgeCentered*Astrocyte+AgeCentered*Neuron_All+AgeCentered*Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-function(i){lm(ValidationGeneSignalNoNA3[i,]~ BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +PsychiatricNoNA3+ Astrocyte +  Neuron_All + Oligodendrocyte+PsychiatricNoNA3*Astrocyte+PsychiatricNoNA3*Neuron_All+PsychiatricNoNA3*Oligodendrocyte+BrainPHCentered*Astrocyte+BrainPHCentered*Neuron_All+BrainPHCentered*Oligodendrocyte+AgeCentered*Astrocyte+AgeCentered*Neuron_All+AgeCentered*Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))}

#PsychiatricThreeCellTypesConfoundsInteractions:
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~ BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +PsychiatricNoNA3+ Astrocyte +  Neuron_All + Oligodendrocyte+PsychiatricNoNA3*Astrocyte+PsychiatricNoNA3*Neuron_All+PsychiatricNoNA3*Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-function(i){lm(ValidationGeneSignalNoNA3[i,]~ BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +PsychiatricNoNA3+ Astrocyte +  Neuron_All + Oligodendrocyte+PsychiatricNoNA3*Astrocyte+PsychiatricNoNA3*Neuron_All+PsychiatricNoNA3*Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))}

#PsychiatricThreeCellTypesInteractions:
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~ PsychiatricNoNA3+ Astrocyte +  Neuron_All + Oligodendrocyte+PsychiatricNoNA3*Astrocyte+PsychiatricNoNA3*Neuron_All+PsychiatricNoNA3*Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-function(i){lm(ValidationGeneSignalNoNA3[i,]~ PsychiatricNoNA3+ Astrocyte +  Neuron_All + Oligodendrocyte+PsychiatricNoNA3*Astrocyte+PsychiatricNoNA3*Neuron_All+PsychiatricNoNA3*Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))}


#PsychiatricPrevelantCellTypesNoConfoundsInteractions:
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~ PsychiatricNoNA3+Astrocyte + Microglia +  Neuron_Interneuron + Neuron_Projection + Oligodendrocyte+PsychiatricNoNA3*Astrocyte+PsychiatricNoNA3*Microglia+PsychiatricNoNA3*Neuron_Interneuron+PsychiatricNoNA3*Neuron_Projection+PsychiatricNoNA3*Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-function(i){lm(ValidationGeneSignalNoNA3[i,]~ PsychiatricNoNA3+Astrocyte + Microglia +  Neuron_Interneuron + Neuron_Projection + Oligodendrocyte+PsychiatricNoNA3*Astrocyte+PsychiatricNoNA3*Microglia+PsychiatricNoNA3*Neuron_Interneuron+PsychiatricNoNA3*Neuron_Projection+PsychiatricNoNA3*Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))}

#PsychiatricPrevelantCellTypeswConfoundsInteractions:
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~ BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +PsychiatricNoNA3+Astrocyte + Microglia +  Neuron_Interneuron + Neuron_Projection + Oligodendrocyte+PsychiatricNoNA3*Astrocyte+PsychiatricNoNA3*Microglia+PsychiatricNoNA3*Neuron_Interneuron+PsychiatricNoNA3*Neuron_Projection+PsychiatricNoNA3*Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-function(i){lm(ValidationGeneSignalNoNA3[i,]~ BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +PsychiatricNoNA3+Astrocyte + Microglia +  Neuron_Interneuron + Neuron_Projection + Oligodendrocyte+PsychiatricNoNA3*Astrocyte+PsychiatricNoNA3*Microglia+PsychiatricNoNA3*Neuron_Interneuron+PsychiatricNoNA3*Neuron_Projection+PsychiatricNoNA3*Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))}


#****Actual code for running those models:


StepwiseBetas<-matrix(0, nrow(ValidationGeneSignalNoNA3), length(AllPredictorsNames))
StepwisePvalues<-matrix(0, nrow(ValidationGeneSignalNoNA3), length(AllPredictorsNames))
StepwiseSE<-matrix(0, nrow(ValidationGeneSignalNoNA3), length(AllPredictorsNames))
StepwiseT<-matrix(0, nrow(ValidationGeneSignalNoNA3), length(AllPredictorsNames))
StepwiseRsquared<-matrix(0, nrow(ValidationGeneSignalNoNA3), 5)

colnames(StepwiseBetas)<-AllPredictorsNames
colnames(StepwisePvalues)<-AllPredictorsNames
colnames(StepwiseSE)<-AllPredictorsNames
colnames(StepwiseT)<-AllPredictorsNames
colnames(StepwiseRsquared)<-c("Rsquared", "AdjRsquared", "ResidualSE", "AIC", "BIC")

system.time(
  for(i in 1:nrow(ValidationGeneSignalNoNA3)){
    StepwiseOutput<-ModeltoUseForStepwise(i)
    StepwiseOutputTemp<-data.frame(row.names(summary.lm(StepwiseOutput)$coefficients), summary.lm(StepwiseOutput)$coefficients)
    colnames(StepwiseOutputTemp)[1]<-"AllPredictorsNames"
    StepwiseOutputTemp2<-join(as.data.frame(AllPredictorsNames), StepwiseOutputTemp, by="AllPredictorsNames")
    
    StepwiseBetas[i,]<-StepwiseOutputTemp2[,2]
    StepwiseSE[i,]<-StepwiseOutputTemp2[,3]
    StepwiseT[i,]<-StepwiseOutputTemp2[,4]
    StepwisePvalues[i,]<-StepwiseOutputTemp2[,5]
    StepwiseRsquared[i,1]<-summary.lm(StepwiseOutput)$r.squared
    StepwiseRsquared[i,2]<-summary.lm(StepwiseOutput)$adj.r.squared
    StepwiseRsquared[i,3]<-summary.lm(StepwiseOutput)[[6]]
    StepwiseRsquared[i,4]<-AIC(StepwiseOutput)
    StepwiseRsquared[i,5]<-BIC(StepwiseOutput)
  }
)


NumberOfTimesInModels<-apply(StepwisePvalues, 2, function(y) length(y)-sum(is.na(y)))
names(NumberOfTimesInModels)<-colnames(StepwisePvalues)

NumberOfTimesSig<-apply(StepwisePvalues, 2, function(y) sum(y<0.05, na.rm=T))
names(NumberOfTimesSig)<-colnames(StepwisePvalues)

temp<-data.frame(ValidationGeneSignalNoNA3wAnnotation[,c(1:6)], StepwiseBetas, StepwisePvalues, StepwiseT, StepwiseSE, StepwiseRsquared)

str(temp)

#And the various output code:

#write.csv(temp, "ValidationGenes_LMWDiagnosisSuicideCelltypeConfounds.csv")
#write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_LMWDiagnosisSuicideCelltypeConfounds_CountInModelAndSig.csv")

write.csv(temp, "ValidationGenes_LMWDiagnosisCelltypeConfounds.csv")
write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_LMWDiagnosisCelltypeConfounds_CountInModelAndSig.csv")

write.csv(temp, "ValidationGenes_LM4.csv")
write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_LM4_CountInModelAndSig.csv")

#write.csv(temp, "ValidationGenes_LM4wSuicide.csv")
#write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_LM4wSuicide_CountInModelAndSig.csv")

#LMWDiagnosisPrevelantCelltypesConfounds:
write.csv(temp, "ValidationGenes_LMWDiagnosisPrevelantCelltypeConfounds.csv")
write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_LMWDiagnosisPrevelantCelltypeConfounds_CountInModelAndSig.csv")

#write.csv(temp, "ValidationGenes_LMWDiagnosis2CelltypeConfounds.csv")
#write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_LMWDiagnosis2CelltypeConfounds_CountInModelAndSig.csv")

#write.csv(temp, "ValidationGenes_LMWDiagnosisGliaCelltypeConfounds.csv")
#write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_LMWDiagnosisGliaCelltypeConfounds_CountInModelAndSig.csv")

#write.csv(temp, "ValidationGenes_LMWDiagnosis3CelltypeConfounds.csv")
#write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_LMWDiagnosis3CelltypeConfounds_CountInModelAndSig.csv")

#write.csv(temp, "ValidationGenes_LMWDiagnosisVIFCelltypeConfounds.csv")
#write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_LMWDiagnosisVIFCelltypeConfounds_CountInModelAndSig.csv")

write.csv(temp, "ValidationGenes_JustDiagnosis.csv")
write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_JustDiagnosis_CountInModelAndSig.csv")

write.csv(temp, "ValidationGenes_JustAge.csv")
write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_JustAge_CountInModelAndSig.csv")

#write.csv(temp, "ValidationGenes_JustHypoxia.csv")
#write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_JustHypoxia_CountInModelAndSig.csv")

write.csv(temp, "ValidationGenes_JustDiagnosisPrevelantCellTypesNoConfounds.csv")
write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_JustDiagnosisPrevelantCellTypesNoConfounds_CountInModelAndSig.csv")

write.csv(temp, "ValidationGenes_JustAgePrevelantCellTypesNoConfounds.csv")
write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_JustAgePrevelantCellTypesNoConfounds_CountInModelAndSig.csv")

#write.csv(temp, "ValidationGenes_JustHypoxiaPrevelantCellTypesNoConfounds.csv")
#write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_JustHypoxiaPrevelantCellTypesNoConfounds_CountInModelAndSig.csv")

#write.csv(temp, "ValidationGenes_JustDiagnosis3CellTypesNoConfounds.csv")
#write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_JustDiagnosis3CellTypesNoConfounds_CountInModelAndSig.csv")

#write.csv(temp, "ValidationGenes_JustDiagnosis3CellTypesConfounds.csv")
#write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_JustDiagnosis3CellTypesConfounds_CountInModelAndSig.csv")

#write.csv(temp, "ValidationGenes_JustDiagnosis3CellTypesInteractions.csv")
#write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_JustDiagnosis3CellTypesInteractions_CountInModelAndSig.csv")

#write.csv(temp, "ValidationGenes_JustDiagnosisNeuronInteraction.csv")
#write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_JustDiagnosisNeuronInteraction_CountInModelAndSig.csv")

#write.csv(temp, "ValidationGenes_JustDiagnosisInterneuronInteraction.csv")
#write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_JustDiagnosisInterneuronInteraction_CountInModelAndSig.csv")

#write.csv(temp, "ValidationGenes_JustDiagnosisProjectionInteraction.csv")
#write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_JustDiagnosisProjectionInteraction_CountInModelAndSig.csv")

#write.csv(temp, "ValidationGenes_DiagnosisConfounds3CellTypesInteractions.csv")
#write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_DiagnosisConfounds3CellTypesInteractions_CountInModelAndSig.csv")

#write.csv(temp, "ValidationGenes_DiagnosisConfounds3CellTypesInteractwConfounds.csv")
#write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_DiagnosisConfounds3CellTypesInteractwConfounds_CountInModelAndSig.csv")

#write.csv(temp, "ValidationGenes_JustDiagnosisPrevelantCellTypesNoConfoundsInteract.csv")
#write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_JustDiagnosisPrevelantCellTypesNoConfoundsInteract_CountInModelAndSig.csv")

#write.csv(temp, "ValidationGenes_LMWPsychosisMoodPrevelantCelltypeConfounds.csv")
#write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_LMWPsychosisMoodPrevelantCelltypeConfounds_CountInModelAndSig.csv")

#write.csv(temp, "ValidationGenes_LMWPsychiatricSuicidePrevelantCelltypeConfounds.csv")
#write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_LMWPsychiatricSuicidePrevelantCelltypeConfounds_CountInModelAndSig.csv")

write.csv(temp, "ValidationGenes_LMWPsychiatricPrevelantCelltypeConfounds.csv")
write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_LMWPsychiatricPrevelantCelltypeConfounds_CountInModelAndSig.csv")

write.csv(temp, "ValidationGenes_LMWPsychiatricPrevelantCelltype.csv")
write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_LMWPsychiatricPrevelantCelltype_CountInModelAndSig.csv")

write.csv(temp, "ValidationGenes_LM4WPsychiatric.csv")
write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_LM4WPsychiatric_CountInModelAndSig.csv")

write.csv(temp, "ValidationGenes_JustPsychiatric.csv")
write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_JustPsychiatric_CountInModelAndSig.csv")

write.csv(temp, "ValidationGenes_LMWPsychiatricAllCelltypeConfounds.csv")
write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_LMWPsychiatricAllCelltypeConfounds_CountInModelAndSig.csv")

write.csv(temp, "ValidationGenes_PsychiatricThreeCellTypesConfoundsInteractionswConfounds.csv")
write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_PsychiatricThreeCellTypesConfoundsInteractionswConfounds_CountInModelAndSig.csv")

write.csv(temp, "ValidationGenes_PsychiatricThreeCellTypesConfoundsInteractions.csv")
write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_PsychiatricThreeCellTypesConfoundsInteractions_CountInModelAndSig.csv")

write.csv(temp, "ValidationGenes_PsychiatricThreeCellTypesInteractions.csv")
write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_PsychiatricThreeCellTypesInteractions_CountInModelAndSig.csv")

write.csv(temp, "ValidationGenes_PsychiatricPrevelantCellTypesNoConfoundsInteractions.csv")
write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_PsychiatricPrevelantCellTypesNoConfoundsInteractions_CountInModelAndSig.csv")

write.csv(temp, "ValidationGenes_PsychiatricPrevelantCellTypeswConfoundsInteractions.csv")
write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_PsychiatricPrevelantCellTypeswConfoundsInteractions_CountInModelAndSig.csv")


temp<-apply(ValidationGeneSignalNoNA3, 1, sd)

write.csv(temp, "STDEV_forValidationGenes.csv")



#************Trying out permutation-based p-values


#I'm going to try a permutation version and see if it helps:
library(lmPerm)

#Basic test run:
Output<-lmp(formula=ValidationGeneSignalNoNA3[i,]~PsychiatricNoNA3+ Astrocyte + Microglia +  Neuron_Interneuron + Neuron_Projection + Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap), perm="Exact")

summary(Output)

Output<-lm(ValidationGeneSignalNoNA3[i,]~PsychiatricNoNA3+ Astrocyte + Microglia +  Neuron_Interneuron + Neuron_Projection + Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))

summary.lm(Output)
summary(Output)


#LMWPsychiatricPrevelantCelltypes+Confounds:
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+AgeCentered+ GenderNoNA3+PsychiatricNoNA3+Astrocyte + Microglia +  Neuron_Interneuron + Neuron_Projection + Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-function(i){lmp(formula=ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3+PsychiatricNoNA3+Astrocyte + Microglia +  Neuron_Interneuron + Neuron_Projection + Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap), perm="Exact")}


#LMWDiagnosisPrevelantCelltypes+Confounds:
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+AgeCentered+ GenderNoNA3+DiagnosisNoNA3+Astrocyte + Microglia +  Neuron_Interneuron + Neuron_Projection + Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-function(i){lmp(formula=ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3+DiagnosisNoNA3+Astrocyte + Microglia +  Neuron_Interneuron + Neuron_Projection + Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap), perm="Exact")}



#I tweaked the stepwise code for permutation based p-values - it is named stepwise but I'm not using it for stepwise procedures
StepwiseBetas<-matrix(0, nrow(ValidationGeneSignalNoNA3), length(AllPredictorsNames))
StepwisePvalues<-matrix(0, nrow(ValidationGeneSignalNoNA3), length(AllPredictorsNames))
StepwiseIterations<-matrix(0, nrow(ValidationGeneSignalNoNA3), length(AllPredictorsNames))
StepwiseRsquared<-matrix(0, nrow(ValidationGeneSignalNoNA3), 5)

colnames(StepwiseBetas)<-AllPredictorsNames
colnames(StepwisePvalues)<-AllPredictorsNames
colnames(StepwiseIterations)<-AllPredictorsNames
colnames(StepwiseRsquared)<-c("Rsquared", "AdjRsquared", "ResidualSE", "AIC", "BIC")

system.time(
  for(i in 1:nrow(ValidationGeneSignalNoNA3)){
    StepwiseOutput<-ModeltoUseForStepwise(i)
    StepwiseOutputTemp2<-data.frame(row.names(summary(StepwiseOutput)$coefficients), summary(StepwiseOutput)$coefficients)
    colnames(StepwiseOutputTemp)[1]<-"AllPredictorsNames"
    #StepwiseOutputTemp2<-join(as.data.frame(AllPredictorsNames), StepwiseOutputTemp, by="AllPredictorsNames")
    
    StepwiseBetas[i,]<-summary(StepwiseOutput)$coefficients[,1]
    StepwiseIterations[i,]<-summary(StepwiseOutput)$coefficients[,2]
    StepwisePvalues[i,]<-summary(StepwiseOutput)$coefficients[,3]
    StepwiseRsquared[i,1]<-summary.lm(StepwiseOutput)$r.squared
    StepwiseRsquared[i,2]<-summary.lm(StepwiseOutput)$adj.r.squared
    StepwiseRsquared[i,3]<-summary.lm(StepwiseOutput)[[6]]
    StepwiseRsquared[i,4]<-AIC(StepwiseOutput)
    StepwiseRsquared[i,5]<-BIC(StepwiseOutput)
  }
)


NumberOfTimesInModels<-apply(StepwisePvalues, 2, function(y) length(y)-sum(is.na(y)))
names(NumberOfTimesInModels)<-colnames(StepwisePvalues)

NumberOfTimesSig<-apply(StepwisePvalues, 2, function(y) sum(y<0.05, na.rm=T))
names(NumberOfTimesSig)<-colnames(StepwisePvalues)

temp<-data.frame(ValidationGeneSignalNoNA3wAnnotation[,c(1:6)], StepwiseBetas, StepwisePvalues, StepwiseIterations, StepwiseRsquared)

str(temp)

write.csv(temp, "ValidationGenes_PsychiatricPrevelantCellTypeswConfoundsPermutation.csv")
write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_PsychiatricPrevelantCellTypeswConfoundsPermutation_CountInModelAndSig.csv")

write.csv(temp, "ValidationGenes_DiagnosisPrevelantCellTypeswConfoundsPermutation.csv")
write.csv(rbind(NumberOfTimesInModels, NumberOfTimesSig), "ValidationGenes_DiagnosisPrevelantCellTypeswConfoundsPermutation_CountInModelAndSig.csv")


#***********************************************************************

#Outputting the best (or most interesting) models for all genes:


#LMWDiagnosisPrevelantCelltypesConfounds:
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +DiagnosisNoNA3+Astrocyte + Microglia +  Neuron_Interneuron + Neuron_Projection + Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-function(i){lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +DiagnosisNoNA3+ Astrocyte + Microglia +  Neuron_Interneuron + Neuron_Projection + Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))}

#LMWPsychiatricAllCelltypesConfounds:
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3+PsychiatricNoNA3+Astrocyte + Endothelial + Microglia + Mural + Neuron_All + Neuron_Interneuron + Neuron_Projection + Oligodendrocyte + Oligodendrocyte_Immature + RBC, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-function(i){lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +PsychiatricNoNA3+ Astrocyte + Endothelial + Microglia + Mural + Neuron_All + Neuron_Interneuron + Neuron_Projection + Oligodendrocyte + Oligodendrocyte_Immature + RBC, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))}

#LM4WPsychiatric:
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3+PsychiatricNoNA3, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-function(i){lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +PsychiatricNoNA3, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))}

#LMWPsychiatricPrevelantCelltypesConfounds:
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3+PsychiatricNoNA3+Astrocyte + Microglia +  Neuron_Interneuron + Neuron_Projection + Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)

ModeltoUseForStepwise<-function(i){lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +PsychiatricNoNA3+ Astrocyte + Microglia +  Neuron_Interneuron + Neuron_Projection + Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))}

#PsychiatricPrevelantCellTypeswConfoundsInteractions:
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~ BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +PsychiatricNoNA3+Astrocyte + Microglia +  Neuron_Interneuron + Neuron_Projection + Oligodendrocyte+PsychiatricNoNA3*Astrocyte+PsychiatricNoNA3*Microglia+PsychiatricNoNA3*Neuron_Interneuron+PsychiatricNoNA3*Neuron_Projection+PsychiatricNoNA3*Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)


ModeltoUseForStepwise<-function(i){lm(ValidationGeneSignalNoNA3[i,]~ BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +PsychiatricNoNA3+Astrocyte + Microglia +  Neuron_Interneuron + Neuron_Projection + Oligodendrocyte+PsychiatricNoNA3*Astrocyte+PsychiatricNoNA3*Microglia+PsychiatricNoNA3*Neuron_Interneuron+PsychiatricNoNA3*Neuron_Projection+PsychiatricNoNA3*Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))}


#LMWDiagnosisPrevelantCelltypesConfounds:
AllPredictorsNames<-names(lm(ValidationGeneSignalNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +DiagnosisNoNA3+Astrocyte + Microglia +  Neuron_Interneuron + Neuron_Projection + Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap))$coefficients)


GeneByCellTypeSubjVar2_Pvalues<-matrix(0, length(SignalSortedNoNA3[,1]), length(AllPredictorsNames))
GeneByCellTypeSubjVar2_Betas<-matrix(0, length(SignalSortedNoNA3[,1]), length(AllPredictorsNames))
colnames(GeneByCellTypeSubjVar2_Pvalues)<-AllPredictorsNames
colnames(GeneByCellTypeSubjVar2_Betas)<-AllPredictorsNames
row.names(GeneByCellTypeSubjVar2_Pvalues)<-row.names(SignalSortedNoNA3)
row.names(GeneByCellTypeSubjVar2_Betas)<-row.names(SignalSortedNoNA3)
head(GeneByCellTypeSubjVar2_Pvalues)


for(i in c(1:length(SignalSortedNoNA3[,1]))){
  
  #temp<-summary.lm(lm(SignalSortedNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +PsychiatricNoNA3+Astrocyte + Microglia +  Neuron_Interneuron + Neuron_Projection + Oligodendrocyte+PsychiatricNoNA3*Astrocyte+PsychiatricNoNA3*Microglia+PsychiatricNoNA3*Neuron_Interneuron+PsychiatricNoNA3*Neuron_Projection+PsychiatricNoNA3*Oligodendrocyte, data=data.frame(SignalSortedNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap)))
 
  temp<-summary(lmp(SignalSortedNoNA3[i,]~BrainPHCentered +AgonalFactorNoNA3 + HoursFinalCorrectedNoNA3+ AgeCentered+ GenderNoNA3 +PsychiatricNoNA3+Astrocyte + Microglia +  Neuron_Interneuron + Neuron_Projection + Oligodendrocyte, data=data.frame(SignalSortedNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap), maxIter=9999, center=F))
   
  GeneByCellTypeSubjVar2_Betas[i,]<-temp$coefficients[,1]
  #GeneByCellTypeSubjVar2_Pvalues[i,]<-temp$coefficients[,4]
  GeneByCellTypeSubjVar2_Pvalues[i,]<-temp$coefficients[,3]
}

#specifically necessary for permutation p-values:
sum(GeneByCellTypeSubjVar2_Pvalues[,7]==0)
GeneByCellTypeSubjVar2_Pvalues[GeneByCellTypeSubjVar2_Pvalues==0]<-1/10000
sum(GeneByCellTypeSubjVar2_Pvalues[,7]==0)

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


 
library(multtest)

GeneByCellTypeSubjVar2_PvaluesAdj<-matrix(0, length(GeneByCellTypeSubjVar2_Pvalues[,1]), length(GeneByCellTypeSubjVar2_Pvalues[1,]))
colnames(GeneByCellTypeSubjVar2_PvaluesAdj)<-colnames(GeneByCellTypeSubjVar2_Pvalues)
row.names(GeneByCellTypeSubjVar2_PvaluesAdj)<-row.names(GeneByCellTypeSubjVar2_Pvalues)


for (i in c(1:length(GeneByCellTypeSubjVar2_Pvalues[1,]))){
  
  #Applying two different types of multiple-comparison corrections to the raw p-values (Benjamini-Hochberg and Benjamini-Yekutieli):
  TempPvalAdj<-mt.rawp2adjp(GeneByCellTypeSubjVar2_Pvalues[,i], proc=c("BH"))
  GeneByCellTypeSubjVar2_PvaluesAdj[,i]<-TempPvalAdj$adjp[order(TempPvalAdj$index),2]
  
}

GeneByCellTypeSubjVar2_PvaluesAdj2<-cbind(GeneByCellTypeSubjVar2_PvaluesAdj, GeneNames)
# write.csv(GeneByCellTypeSubjVar2_PvaluesAdj2, "GeneByCellTypeSubjVar2_PvaluesAdj.csv")

GeneByCellTypeSubjVar2DF<-as.data.frame(cbind(GeneByCellTypeSubjVar2_Betas, GeneByCellTypeSubjVar2_Pvalues, GeneByCellTypeSubjVar2_PvaluesAdj))

temp<-cbind(GeneNames, GeneByCellTypeSubjVar2DF)
write.csv(temp, paste("GeneByCellTypeSubjVar2DF.csv", sep=""))


#####Adding some chi-square stats to back-up improved sensitivity: I need to generate some random gene lists to determine whether to top genes in association with diagnosis in our models are more associated with diagnosis in the literature than simple lists of randomly-selected genes:

write.csv(GeneNames[sample(1:11979, 40, replace=T),], "RandomGeneGenerator2.csv")

