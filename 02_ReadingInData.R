#This is my attempt to pull together a more organized code document specifically related to cell type index analyses using DLPFC Affymetrix data from the Pritzker foundation
##Megan Hagenauer, October 27 2016

#**************************************

#Data source:

#I found the appropriate Affy data on the Pritzker server (/home/GROUPS/pritzker/pritzker_umich/Affy_reanalysis_2010). At the end of the "dlpfc_summary_new.doc" description, I found this text:

#The analysis of dlpfc data created the following R objects.  Some are not exported, but are available upon request:
##1. Raw data, dlpfc, exported as "dlpfc_5708_rma.txt";
##2. Normalized, dlpfc.norm, not exported;
##3. Median centered, dlpfc.mc, not exported;
##4. Re-normalized, dlpfc.mc2, not exported;
##5. Averaged based on dlpfc.mc2, dlpfc.avg2, 172 unique samples, exported as "dlpfc_avg2.txt";
##6. dlpfc.avg2 adjusted by PC1, dlpfc.avg2.PC1, exported as "dlpfc_avg2_PC1.txt";


#Sample information files:
##11. "dlpfc_sample.txt", for the original 367 samples;
##12. "dlpfc_sample_filter.txt", for the 337 samples after outlier removal;
##13. "dlpfc_avg_sample.txt", for 172 unique samples.
#Summary file:
##14. "dlpfc_ summary_new.docx".

#For the DLPFC, 
#I chose #5 ("dlpfc_avg2.txt") and #13 ("dlpfc_avg_sample.txt") for the cell type analyses because the data in these files should be fully preprocessed with the exception of correcting for the effects of pH and age.

#Later, after preprocessing other microarray datasets, I decided that HoursFinal and Suicide were likely to matter, and so I went back to the master Pritzker database and added this information to "dlpfc_avg_sample.txt".  The new file is called "DLPFC_SubjectInfoUpdated.csv"

#**************************************

#Reading in the SampleInfo:

SampleInfo<-as.matrix(read.csv("DLPFC_SubjectInfoUpdated_IrvineZT.csv",header=T))
row.names(SampleInfo)<-as.numeric(SampleInfo[,2])
SubjectInfoSorted<-SampleInfo[order(row.names(SampleInfo)),]

colnames(SubjectInfoSorted)

 [1] "X"                   "SUBJECT"             "COHORT"              "SITE"               
 [5] "DIAGNOSIS"           "DIAGNOSIS_SSRI"      "AGONAL_FACTOR"       "PH"                 
 [9] "GENDER"              "RACE"                "AGE"                 "AGE2"               
[13] "HOURS.FINAL"         "HoursFinalCorrected" "Suicide"             "ACI"                
[17] "SAMPLE_TYPE"         "ORIGIN"              "CHIP_TYPE"           "BLOCK"              
[21] "BLOCK_FULLNAME"      "TOD" 


rownames(SubjectInfoSorted)[c(1:4)]
[1] "1834" "1856" "1858" "1881"

is.numeric(SubjectInfoSorted)
[1] FALSE

SubjectInfoSorted<-SubjectInfoSorted[order(row.names(SubjectInfoSorted)),]
SubjectInfoSorted[c(1:3), c(1:3)]

     X                    SUBJECT COHORT          
1834 "D_133A_1834_HC_1_C" "1834"  "Schiz Cohort 1"
1856 "D_133A_1856_HC_1_C" "1856"  "Schiz Cohort 1"
1858 "D_133A_1858_HC_1_C" "1858"  "Schiz Cohort 1"


Diagnosis<-as.factor(SubjectInfoSorted[,5])
Diagnosis<-relevel(Diagnosis, ref="Control")
str(Diagnosis)

AgonalFactor<-as.numeric(SubjectInfoSorted[,7])

BrainPH<-as.numeric(SubjectInfoSorted[,8])

Gender<-as.factor(SubjectInfoSorted[,9])
Gender<-relevel(Gender, ref="M")
str(Gender)

Age<-as.numeric(SubjectInfoSorted[,11])

HoursFinal<-as.numeric(SubjectInfoSorted[,13])

HoursFinalCorrected<-as.numeric(SubjectInfoSorted[,14])

Suicide<-as.factor(SubjectInfoSorted[,15])

TOD<-as.numeric(SubjectInfoSorted[,22])

ExsanguinatedSubjects<-c("1881", "3038", "3281", "3426", "3618", "3850", "3952", "4063", "4069")
#note: 3952 is iffy - they died from blunt trauma, but presumably that was accompanied by massive internal bleeding



#*********************************

## Moving on to the (already pre-processed) Signal information:

Signal<-as.matrix(read.delim("dlpfc_avg2.txt",header=T,row.names=1,sep="\t"))
SignalSorted<-Signal[,order(colnames(Signal))]

SignalSorted[c(1:3),] 
                X1834       X1856      X1858       X1881       X1964       X2208
10000_at -0.096410058  0.07127575  0.1876630 -0.07429976 -0.04345254  0.12822757
10001_at  0.355015060 -0.05120225 -0.1396510  0.20049059 -0.25016275 -0.14930513
10002_at  0.008055338 -0.02904995  0.0783977 -0.02886939  0.08380615  0.09241613

#Note: this data is already log(2)transformed, quantile normalized, and median centered

##  hmmm... matching the subject info and signal is going to be a pain in the butt because both arent numeric - we have to figure out how to extract the numeric sample number information

colnames(SignalSorted)<-as.numeric(gsub("[^0-9]", "", colnames(SignalSorted)))
SignalSorted[c(1:3),c(1:3)]

                 1834        1856       1858
10000_at -0.096410058  0.07127575  0.1876630
10001_at  0.355015060 -0.05120225 -0.1396510
10002_at  0.008055338 -0.02904995  0.0783977

Yes!

SignalSorted<-SignalSorted[,order(colnames(SignalSorted))]
SignalSorted[c(1:3),c(1:3)]
                 1834        1856       1858
10000_at -0.096410058  0.07127575  0.1876630
10001_at  0.355015060 -0.05120225 -0.1396510
10002_at  0.008055338 -0.02904995  0.0783977

cbind(colnames(SignalSorted), row.names(SubjectInfoSorted))
plot(as.numeric(colnames(SignalSorted))~as.numeric(row.names(SubjectInfoSorted)))
#Everything matches up. Hurray!


GeneNames<-as.matrix(read.csv("AffyGeneNames.csv",header=T))
dim(GeneNames)
GeneNames[c(1:3), c(1:2)]
#I already set up this file so that it would match the order of Signal

#*********************************

#Removing data that is missing associated subject information regarding important confounds (either pH, HoursFinal, or Agonal Factor):

BrainPHNoNA<-BrainPH[is.na(BrainPH)==FALSE]
DiagnosisNoNA<-Diagnosis[is.na(BrainPH)==FALSE]
AgonalFactorNoNA<-AgonalFactor[is.na(BrainPH)==FALSE]
GenderNoNA<-Gender[is.na(BrainPH)==FALSE]
HoursFinalNoNA<-HoursFinal[is.na(BrainPH)==FALSE]
HoursFinalCorrectedNoNA<-HoursFinalCorrected[is.na(BrainPH)==FALSE]
AgeNoNA<-Age[is.na(BrainPH)==FALSE]
SuicideNoNA<-Suicide[is.na(BrainPH)==FALSE]
TODNoNA<-TOD[is.na(BrainPH)==FALSE]


SignalSortedNoNA<-SignalSorted[,is.na(BrainPH)==FALSE]
MeanSignalByProbe<-apply(SignalSortedNoNA, 1, mean)
SDSignalByProbe<-apply(SignalSortedNoNA, 1,sd)
plot(sort(MeanSignalByProbe))
plot(sort(SDSignalByProbe))
#Note: There is no way to use the Mean Signal to filter out non-detected genes because this data was already normalized to deal with batch effects.


BrainPHNoNA2<-as.numeric(BrainPHNoNA[is.na(HoursFinalNoNA)==FALSE])
DiagnosisNoNA2<-DiagnosisNoNA[is.na(HoursFinalNoNA)==FALSE]
AgonalFactorNoNA2<-AgonalFactorNoNA[is.na(HoursFinalNoNA)==FALSE]
GenderNoNA2<-GenderNoNA[is.na(HoursFinalNoNA)==FALSE]
HoursFinalCorrectedNoNA2<-as.numeric(HoursFinalCorrectedNoNA[is.na(HoursFinalNoNA)==FALSE])
HoursFinalNoNA2<-as.numeric(HoursFinalNoNA[is.na(HoursFinalNoNA)==FALSE])
AgeNoNA2<-AgeNoNA[is.na(HoursFinalNoNA)==FALSE]
SuicideNoNA2<-SuicideNoNA[is.na(HoursFinalNoNA)==FALSE]
TODNoNA2<-TODNoNA[is.na(HoursFinalNoNA)==FALSE]

SignalSortedNoNA2<-SignalSortedNoNA[,is.na(HoursFinalNoNA)==FALSE]


BrainPHNoNA3<-as.numeric(BrainPHNoNA2[is.na(AgonalFactorNoNA2)==FALSE])
DiagnosisNoNA3<-DiagnosisNoNA2[is.na(AgonalFactorNoNA2)==FALSE]
AgonalFactorNoNA3<-AgonalFactorNoNA2[is.na(AgonalFactorNoNA2)==FALSE]
GenderNoNA3<-GenderNoNA2[is.na(AgonalFactorNoNA2)==FALSE]
HoursFinalCorrectedNoNA3<-as.numeric(HoursFinalCorrectedNoNA2[is.na(AgonalFactorNoNA2)==FALSE])
HoursFinalNoNA3<-as.numeric(HoursFinalNoNA2[is.na(AgonalFactorNoNA2)==FALSE])
AgeNoNA3<-AgeNoNA2[is.na(AgonalFactorNoNA2)==FALSE]
SuicideNoNA3<-SuicideNoNA2[is.na(AgonalFactorNoNA2)==FALSE]
TODNoNA3<-TODNoNA2[is.na(AgonalFactorNoNA2)==FALSE]

SignalSortedNoNA3<-SignalSortedNoNA2[,is.na(AgonalFactorNoNA2)==FALSE]

cosTODNoNA3<-cos(TODNoNA3)
sinTODNoNA3<-sin(TODNoNA3)

PsychiatricNoNA3<-rep(0, length(DiagnosisNoNA3))

for(i in 1:length(DiagnosisNoNA3) ){
	if(DiagnosisNoNA3[i]=="Control"){
	PsychiatricNoNA3[i]<-0
	}else{
		PsychiatricNoNA3[i]<-1
	}
}


ExsanguinatedNoNA3<-as.numeric(colnames(SignalSortedNoNA3)%in%ExsanguinatedSubjects)
table(ExsanguinatedNoNA3, DiagnosisNoNA3)

ExsanguinatedNoNA3 Control BP MD Schiz
                 	0      47 	15 	33   16
                 	1       5 	 3  	1     0
                 					
                 		DiagnosisNoNA3
ExsanguinatedNoNA3 Control BP MD Schiz
                	 0      	66 	21 	39    22
                 	1      	 5 		 3  	1     0
 #A higher percentage of bipolar subjects were exsanguinated than other diagnostic categories


#This was run specifically for making figures later: 
ExsanguinatedSignalSortedNoNA3<-SignalSortedNoNA3[, colnames(SignalSortedNoNA3)%in%ExsanguinatedSubjects]
dim(ExsanguinatedSignalSortedNoNA3)
[1] 11979     9
colnames(ExsanguinatedSignalSortedNoNA3)
[1] "1881" "3038" "3281" "3426" "3618" "3850" "3952" "4063" "4069"



