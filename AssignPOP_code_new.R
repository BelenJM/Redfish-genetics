


library(adegenet)
library(diveRsity)
library(hierfstat)
#library(rapport)
#library(graph4lg)
#library(SNPRelate)
#library(qvalue)
#library(dartR)
#library(StAMPP)
#library(HardyWeinberg)
library(genetics)
#library(pegas)
#library(readxl)
#library(klaR)
library(dplyr)

# Load the assignPOP: 
library(assignPOP)

setwd("C:/Users/Paula Langkilde/Desktop/10 ECTS projekt/Artikler fra Belén/AssignPOP/AssignPOP_New/")
#setwd("H:/DTU/13. Rødfisk/anna-mae/")


# Loading the data (baseline) using the command form assignPOP package:
# Also naming the populations: DEEP; SHALLOW and SLOP
Baseline_total <- read.Genepop( "assignPOP_baseline_GP.gen", 
                                pop.names=c("DEEP", "SHALLOW", "SLOPE"), haploid = FALSE)
str(Baseline_total)

#______________________________________________
#### Step 1  Dimensional reduction ####
#______________________________________________

# REMOVAL of low variance loci. Still keeping all the 91 individuals)
# How should decide upon a theshold? Should we use 80% like in the data-analasis part or is 95% ok?  
YourGenepopRd <- reduce.allele( Baseline_total, p = 0.95)
#      p = 0.95 indicates the removal of loci having the frequency of an allele
#      greater than 0.95.
#     By coshing this threshold, 6 column (alleles) have been removed, and 93 SNP remains. 

#______________________________________________
####  MONTE CARLO- procedure  ####
#______________________________________________

## Using the Monte Carlo resampling method: 
accuMC <- assign.MC(YourGenepopRd, train.inds=c(0.5, 0.7), train.loci=c(0.1, 0.25, 0.5, 1),
                    loci.sample="fst", iterations=30, model="tree", dir="Result-folder1/")
# A total of 240 assignment test were completed

accuMC <- accuracy.MC( dir = "Result-folder1/" )
# See the table for the calculation
accuMC <- read.table("Result-folder1/Rate_of_240_tests_3_pops.txt", header=T)
head(accuMC)
# We see here that we have indicies below 30. This is not good for Monte Carlo resampling (I think)
# One shoud go with K-fold instad since we have baseline populations of different sizes. 16 -> 47


# Create an assignment accuracy box plot:
library(ggplot2)
accuracy.plot( accuMC, pop=c("all", "DEEP", "SHALLOW", "SLOPE")) +
  ylim(0, 1) + #Set y limit between 0 and 1
  annotate("segment",x=0.4,xend=3.6,y=0.33,yend=0.33,colour="red",size=1) +
  #Add a red horizontal line at y = 0.33 (null assignment rate for 3 populations)
  ggtitle("Monte-Carlo cross-validation using genetic loci")+
  #Add a plot title
  theme(plot.title = element_text(size=20, face="bold"))

#______________________________________________
####  K-fold - procedure  ####
#______________________________________________

#  Here I'm trying the same thing, but I'm using the K-fold re-sampling. This is 
#  probably better since I have a population in the baseline with only 16 individuals 
#  compared to 47 individual in another group. 
assign.kfold( YourGenepopRd, k.fold=c(3,4,5), train.loci=c(0.1,0.25,0.5,1),
              loci.sample="random", dir="Result-folder10/", model="lda" )
# 48 assignment test were completed with the K-fold method. 

accuKF <- accuracy.kfold(dir = "Result-folder10/")


membership.plot( dir = "Result-folder10/" )
# Showing the result of the K-fold method for the baseline: (table)
accuKF <- read.table("Result-folder10/Rate_of_48_tests_3_pops.txt", header=T)
head(accuKF)


#______________________________________________
####  Assigning unknown samples  ####
#______________________________________________

YourGenepopUnknown <- read.Genepop( "HybridLab_simu_ecotype.gen")                                                                           #  "pop_C"), haploid = FALSE)
#YourGenepopRd2 <- reduce.allele( YourGenepopUnknown, p = 0.95)

      # A-M: This command is not working. I cant figure out what the problem is. 
assign.X( x1=YourGenepopRd, x2=YourGenepopUnknown, dir="Result-folder2/", model="tree")

accuUnknown <- read.table("Result-folder3 - K-fold/AssignmentResult.txt", header=T)
head(accuUnknown)

write.table(accuUnknown, "assignPOP_result_K_fold.xlsx")







