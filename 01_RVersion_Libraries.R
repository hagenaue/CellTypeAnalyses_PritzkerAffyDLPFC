#This is my attempt to pull together a more organized code document specifically related to cell type index analyses using DLPFC Affymetrix data from the Pritzker foundation
##Megan Hagenauer, October 27 2016

#**************************************

#R: 

##Overall, this analysis spanned 2 years.  Therefore, the version of R used for running this code varied across the years.
##Currently I have:
### R version 3.3.0 (2016-05-03) -- "Supposedly Educational"
### Copyright (C) 2016 The R Foundation for Statistical Computing
### Platform: x86_64-apple-darwin13.4.0 (64-bit)

#I also tend to use R-studio as a GUI:
### Version 0.99.896 – © 2009-2016 RStudio, Inc.
### Mozilla/5.0 (Macintosh; Intel Mac OS X 10_10_5) AppleWebKit/601.5.17 (KHTML, like Gecko)

#****************************************************

#(Potentially) Relevant code libraries:
##I apologize for the fact that I'm not entirely sure which of these I actually used - I just tend to load them all because I use them regularly
##Each of these packages requires installation before they can be loaded.

library(gdata)
library("fields")
library(stats)
library(car)
library(affy)
library(preprocessCore)
library(multtest)
library(plyr)
library("car")
library("drc")
library("gdata")
library("gtools")
library("lattice")
library("lmtest")
library("magic")
library("MASS")
library("plotrix")
library("proto")
library("sandwich")
library("spam")

#************************
