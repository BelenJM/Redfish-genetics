#BiocManager::install(c("SNPRelate", "qvalue"))
#install.packages("BiocManager")

#_____________________________________________________________________
#### Loading Packages: #### 
#_____________________________________________________________________


library(adegenet)
library(diveRsity)
library(hierfstat)
library(rapport)
library(graph4lg)
library(SNPRelate)
library(qvalue)
library(dartR)
library(StAMPP)
library(HardyWeinberg)
library(genetics)
library(pegas)

# this is needed
library(readxl)


#setwd("C:/Users/Paula Langkilde/Desktop/10 ECTS projekt/Files for PA/3_part_RedFish/Merging part try 2/")
setwd("H:/DTU/13. Rødfisk/anna-mae/error_15032022/")

setwd("H:/DTU/13. Rødfisk/anna-mae/simulations/")


#_____________________________________________________________________
####   LOADING THE DATA  ####
#_____________________________________________________________________


#     Since we defined the working-area, we're able to just load it. Furthermore, 
#     we're using the "read.genepop" since this is the format of the Baseline-file

#     We're now defining a "object" named "GP_RedFish"

#     Furthermore, the alleles are now in a 2 "number-system". 02 or 04 etc. Therefore 
#     ncode=2
test <- read.genepop("HybridLab_baseline_ALL.gen", ncode=3)
GP_RedFish <- read.genepop("redfish_noBaseline_controlNames (2).gen",ncode = 2)


#     Showing the structure of the object "GP_RedFish" 
GP_RedFish

#_____________________________________________________________________
####                  PREPARATION OF THE DATA SET             ####
#_____________________________________________________________________


# Creation of the pop-map file with the Habitat information form GeneClass2. 
# Each individual has been allocated to one of the three habitat-types with
# a 96.6% certainty (mean) : 
pop_file <- read_excel("Pop_file_test.xlsx")
## A-M: Clones now included :) 



#  Here I'm renaming the names given by GeneClass2 for the habitats: DEEP, SHALLOW, and SLOPE 
#  RE19_GR_TB306_21= DEEP, RE19_N_TB294_86 = SHALLOW, RE19_GR_TB309_07 = SLOPE)
library(dplyr)
pop_file$Habitat <- recode_factor(pop_file$Habitat, 
                                        RE19_GR_TB306_21  = "DEEP",
                                        RE19_N_TB294_86 = "SHALLOW", 
                                        RE19_GR_TB309_07 = "SLOPE")


# Check the order of the names: 
head(pop_file)


# Storing names into a vector with the same names of individuals from my genepop
# file
names <- rownames(GP_RedFish@tab)

# We need to change the format of the vector, as it won't allow me to change
# the column name (which we need to use the merge function below)
names <- as.data.frame(names)

# Now I can change the name of the column of the "name". It should be the exact same 
# as the column I want to merge with from the other file "pop_file"
colnames(names) <- c("Individual_ID")


# Check if both files "names" and "pop_file" have the same heading
# and which individuals they have at the beginning. Should not necessarily be the same since we haven't
# merged them yet. 
head(names)
head(pop_file)
tail(names)
tail(pop_file)


# Now I can merge the names from my genepop file with their population info:
ordered_GP_pop <- merge(names, pop_file, by = "Individual_ID", sort=F)



# Check the headers to see if the merge function worked properly. It should be the 
# same order for both of them. Otherwise, it didn't work. 
head(ordered_GP_pop)
head(names)
tail(ordered_GP_pop)
tail(names)



# Now I am confident that the population information corresponds to each
# individual in the GP file

str(ordered_GP_pop$Habitat)

# Here I'm setting the pop within the GP_RedFish object to be the 3 leveled factor form 
# the habitats I just made. DEEP, SHALLOW and SLOPE
GP_RedFish@pop <- as.factor(ordered_GP_pop$Habitat)



# Here, we're setting the factors, so that we can show a table in the next command bellow: 
GP_RedFish_pop.map_fac <- GP_RedFish@pop
# Print the table to see how many individuals in each group: 
table(GP_RedFish_pop.map_fac)


# Now we can see form the table, how many individuals from each group 
# NOTE: We can see that we have a total of 943 individuals. All clones has been included
# so far. 

#         DEEP: ----------------369 
#         SHALLOW: -------------274
#         SLOPE: ---------------300

#         TOTAL: ---------------943


#_____________________________________________________________________
#### 2. PREPARATION OF THE DATA SET 2 ####
#_____________________________________________________________________

#GP_RedFish2 <- read.genepop("redfish_noBaseline_controlNames_Sampling.gen",ncode = 2)


# I don't know if it is necessary to load this one again under another name? 
## BJM: What is the difference between the two datasets? (GP_RedFish and GP_RedFish2?
## If any, then it is no problem in just overwriting GP_RedFish
## You can do this also by removing whatever is in memory now using xx
## and then loading again GP_RedFish

### A-M: Thanks! There's no difference, so i could just go on with GP_RedFish. But
### sinc I'm using GP_RedFish2 in my script from here and down, I'll just keep it. 
### But I'll keep in mind that I could have just overwritten the original one since 
### there's no difference. Thanks :) 
GP_RedFish2 <- read.genepop("redfish_noBaseline_controlNames (2).gen",ncode = 2)
GP_RedFish2@pop

# Creation of the pop-map file:
# Import the population file with information about where the individuals were sampled form
# Sampling site: Iceland, Greenland or Faroe

### A-M: I've fixed the Individual_ID so they're matching the ones used in "Pop_file_test.xlsx" :)
pop_file2 <- read_excel("Sampling_site.xlsx")
head(pop_file2)

# Check the order of the names: 
head(pop_file2)


# Storing names into a vector with the names of individuals from my genepop file
## BJM: remember to use the correct name of file,
## before you used GP_RedFish2:
#names2 <- rownames(GP_RedFish@tab)
## BJM
### A-M: Thanks, silly mistake
names2 <- rownames(GP_RedFish2@tab)


# We need to change the format of the vector, as it won't allow me to change
# the column name (which we need to use the merge function below)
names2 <- as.data.frame(names2)

# BJM: I count the number of entries of the vector, so I have an overview
# We do this by using length(): usually one needs to know the dimensions
# of the object (is it a vector? is it a matrix?) to play with counting within
# columns, rows, etc
dim(names2) # this tells me the dimensions: it has 943 rows and 1 column
# so to count the number of rows, I'd have to tell R to count the column:
length(names2[,1])
# we can also count the pop_file2
dim(pop_file2)
length(pop_file2[,1])

### A-M: Thanks, really handy. I think this method is faster than finding it in 
### global environment to the top right 


# I can change the name of the column of the "name". It should be the exact same 
# as the column I want to merge with from the other file "pop_file"
colnames(names2) <- c("Individual_ID")



# Check if both files "names" and "pop_file" have the same heading
# and which individuals they have at the beginning. Should not necessarily be the 
# same since we haven't merged them yet. 
head(names2)
head(pop_file2)
tail(names2)
tail(pop_file2)

# Now I can merge the names from my genepop file with their population info:
ordered_GP_pop2 <- merge(names2, pop_file2, by = "Individual_ID", sort=F)
length(ordered_GP_pop2[,1])

# Check the headers to see if the merge function worked properly. It should be the 
# same order for both of them. Otherwise, it didn't work. 
head(ordered_GP_pop2)
head(names2)
tail(ordered_GP_pop2)
tail(names2)


# Now I am confident that the population information corresponds to each
# individual in the GP file
GP_RedFish2@pop <- as.factor(ordered_GP_pop2$Sampling_site)


# Here I'm setting the pop within the GP_RedFish object to be the 3 level factor form 
# the habitats I just made. Greenland, Iceland and Faroe
GP_RedFish_pop.map_fac2 <- GP_RedFish2@pop

# Print the table to see how many individuals in each group: 
table(GP_RedFish_pop.map_fac2)

#         Faroe: ----------------213 
#         Greenland:-------------406
#         Iceland: --------------324

#         TOTAL: ---------------943



#__________________________________________________________________________

# I now want to merge the two ordered_GP_pop files that I just created:
# ordered_GP_pop and ordered_GP_pop2
# First I'm checking the individuals.
head(pop_file)
head(pop_file2)
tail(pop_file)
tail(pop_file2)

# I now have two pop-files I need to merge into one so that I have 
        #   Individual_ID      Habitat      Sampling site
ordered_GP_pop3 <- merge(pop_file, pop_file2, by = "Individual_ID", sort=F)


# Here I'm setting the sampling-site as a factor 
ordered_GP_pop3$Sampling_site <- as.factor(ordered_GP_pop3$Sampling_site)  
str(ordered_GP_pop3)


#  Here, I combine the two columns with 3 factores in each, so that we have 9 factores in total. 
Samp_Hab <- paste(ordered_GP_pop3$Habitat, ordered_GP_pop3$Sampling_site, 
                  sep = "_",recycle0 = ordered_GP_pop3$Individual_ID )

str(Samp_Hab)


# MERGING:
# Here I'm merging the individuals with the Samp_hab: 
ordered_GP_pop3$Samp_Hab <- paste(ordered_GP_pop3$Habitat, ordered_GP_pop3$Sampling_site, sep = "_",recycle0 = ordered_GP_pop3$Individual_ID )
Master_pop_file <- ordered_GP_pop3 %>% select(Individual_ID, Samp_Hab)



# Setting it as a factor of 9 levels instead of a character.
# We have 9 factores. DEEP, SKALLOW and SLOPE for each country (Faro, Iceland and Greenland)
Master_pop_file$Samp_Hab <- as.factor(Master_pop_file$Samp_Hab)



# Print the table to see how many individuals in each group: 
table(Master_pop_file$Samp_Hab)

#     Now we can see form the table, how many individuals from each group 
#         DEEP - Faroe. ...................67
#         DEEP - Greenland. ..............117 
#         DEEP - Iceland .................185
#         SHALLOW - Faror ................128
#         SHALLOW - Greenland ............129
#         SHALLOW - Iceland ...............17
#         SLOPE - Faroe ...................18
#         SLOPE - Greenland ..............160
#         SLOPE - Iceland ................122
#         ___________________________________
#         SUM                             943
#         ___________________________________




# Now I am confident that the population information corresponds to each
# individual in the GP file
GP_RedFish2@pop <- as.factor(ordered_GP_pop3$Samp_Hab)


# Here I'm setting the pop within the GP_RedFish object to be the 9 level factor form 
# the habitats I just made. Greenland, Iceland and Faroe
## BJM: watch out -- here you are not setting the pop within the GP_RedFish object,
## You are copying the pop file into a separate vector called ".map_fac3"
## Is this what you want?

### A-M: Thanks yes, i need it for the next part. For the: GP_hf <- genind2hierfstat(GP_RedFish2,pop=GP_RedFish_pop.map_fac3)
GP_RedFish_pop.map_fac3 <- GP_RedFish2@pop


#_____________________________________________________________________
            ####  CONTROLS  ####
#_____________________________________________________________________
# BJM: before you start with the rest of the analysis, we will first
# check if the controls have consistent genotypes. We can do this 
# by different ways, e.g. one way is to visualize their positions
# within the PCA

adeg_pca <- scaleGen(GP_RedFish2, NA.method="mean",scale=F)
pca.adeg_pca <- dudi.pca(adeg_pca, scale=F, nf = 10,scannf = F)
# this is the normal PCA, as we have been doing it:
s.class(pca.adeg_pca$li, fac=GP_RedFish2@pop, col=rep('black', 6),cpoint=1)

# Now, we want the PCA but only with the controls, just to reduce the
# number of points to visualize the plot better:
# 1)we first select the controls:
# many ways to do this, e.g we "grep" (select) the names that contain
# a certain character, in our case the controls all follow the same
# pattern, they have a "C" in their string
test_id <- pca.adeg_pca$li[grepl("_C", rownames(pca.adeg_pca$li)),]
test_names <- which(rownames(pca.adeg_pca$li) %in% rownames(test_id)) 
#fac=GP_RedFish2@pop[test_names],

# 2)we plot the labels of the controls
# We can change the size of the font with clabel
s.label(test_id, boxes = F,clabel = 0.5)

# 3)once the controls are tested,
# we proceed to remove all but once from each control
# Also different ways to do this: we could randomly take one control
# from each set, or we could take the one with least missing data
# or we can take the first one
gt <- extract.gt(vcf, return.alleles = FALSE)  # with TRUE returns the A-G-C-T instead of 1/0, 0/0,...
dataset <- as.data.frame(GP_RedFish2$tab[test_names,],stringsAsFactors = F)                   # create dataframe of genotypes 

# count the missing SNPs per individual:
snpLoci = rowSums(is.na(dataset))
hist(snpLoci, breaks=100)
# which loci have extreme values of NA (i.e. missing data? e.g. >130)
which(snpLoci>130)

# Can you fill this list with all the individuals to remove? (you keep only
# one from each set of controls:
# e.g. controls_to_remove <- c("RE19_GR112_C2", "xx", "xx")
#e.g.
controls_to_remove <- c("RE19_GR112_C2")
# We drop these individuals from the entire genind object:
GP_RedFish3 <- GP_RedFish2[!row.names(GP_RedFish2$tab) %in% controls_to_remove]


#_____________________________________________________________________
        ####  FILTERING: LOCI WITH HIGH MISSING DATA  ####
#_____________________________________________________________________
# BJM: before you start with the rest of the analysis, we will first
# check if the controls have consistent genotypes. We can do this 
# by different ways, e.g. one way is to visualize their positions
# within the PCA

# check the call rate from all loci
GP_RedFish3_gl <- gi2gl(GP_RedFish3, parallel = TRUE)
GP_RedFish3_gl <- gl.recalc.metrics(GP_RedFish3_gl) # first calculate the metrics into the @other slot for gl.filter.callrate to work...
hist(GP_RedFish3_gl$other$loc.metrics$CallRate, 
     main="Levels of call rate per marker", xlab = "Call rate", breaks=100)
# filter for loci with call rate < 0.8
length(which(GP_RedFish3_gl$other$loc.metrics$CallRate>0.8))

select_lowestNA <- as.data.frame(GP_RedFish3_gl$other$loc.metrics$CallRate)
head(select_lowestNA)
select_lowestNA[,2] <- names(GP_RedFish3_gl$other$loc.metrics$CallRate)
head(select_lowestNA)
names_split <- strsplit(select_lowestNA[,2],"_")

# select the highest value from the unique
colnames(select_lowestNA) <- c("call_rate")

# filter out markers with missing data > 0.2
markers_filt <- which(select_lowestNA$call_rate>0.8) 
length(markers_filt) # 91 markers left
index2 <- which(GP_RedFish3@loc.fac %in% rownames(select_lowestNA)[markers_filt])
GP_RedFish4 <- GP_RedFish3[,index2]



#_____________________________________________________________________
      ####  FILTERING: INDIVIDUALS WITH HIGH MISSING DATA  ####
#_____________________________________________________________________
# BJM: before you start with the rest of the analysis, we will first
# check if the controls have consistent genotypes. We can do this 
# by different ways, e.g. one way is to visualize their positions
# within the PCA
dataset <- as.data.frame(GP_RedFish4$tab,stringsAsFactors = F)                   # create dataframe of genotypes 

# Exploring the dataset:
dataset[0:3,0:3]
length(dataset[,1]) # 942 individuals
length(dataset[1,]) ## 91 loci (x2 alleles)

snpLoci2 <- matrix(ncol=2)
snpLoci2$NAs <- unlist(rowSums(is.na(dataset)))
snpLoci2$pop <- GP_RedFish4$pop
head(snpLoci2)
hist(snpLoci2$NAs, breaks=100)
thresh <- 0.2*91 # I choose a threshold of 20%
snpLoci2_filt <- snpLoci2$pop[snpLoci2$NAs<thresh]
snpLoci2_filt_indiv <- snpLoci2$NAs[snpLoci2$NAs<thresh]
indiv_keep <- names(snpLoci2_filt_indiv)
length(indiv_keep) # 838 individuals to keep

# We now remove individuals with high missing data 
# so we Keep the individuals with least missing data:
GP_RedFish5 <- GP_RedFish4[row.names(GP_RedFish4$tab) %in% indiv_keep]


#_____________________________________________________________________
            ####  CALCULATING ALLELE FREQUENCIES  ####
#_____________________________________________________________________

#       First thing is to convert it into yet another format called 
#       hierfstat
#       Convert to hierfstat
GP_hf <- genind2hierfstat(GP_RedFish5,pop=GP_RedFish5$pop)

#       Looking at the structure of the very beginning of the data by using the 
#       command "header"
#       So here, we're asking R to show us row 1-3 and column 1-3

head(GP_hf[1:3,1:3])

#       We're now ready to calculate the allele frequencies per. SNP per population
#       By using the command bellow: 

GP_stats <- basic.stats(GP_hf,diploid=TRUE,digits=4)
#       digits  = How many digits to print in the output
#       diploid = Weatehr the individuals are diploids (default) or haploids


#       We can now look at the structure of the object and the pop.freq = population frequency
str(GP_stats$pop.freq)



#_____________________________________________________________________
                  ####  HARDY-WEINBERG EQUILIBRIUM  ####
#_____________________________________________________________________
#       First we need to change the format again. This command converts a genind
#       object into a genlight object
gp_gi <- gi2gl(GP_RedFish2, verbose = NULL)

#       We need to do 2 things before we procied: 
#       1) add a metrics entry into the data-set, and convert it into a data-frame 
#          for it to work. 
genli_vcf<- gl.recalc.metrics(gp_gi)

#       2) Convert it into a data-frame for it to work. 
genli_vcf$other$loc.metrics <- as.data.frame(genli_vcf$other$loc.metrics)


#       We're now ready to run the HWE test per locus. 
gg <- gl.report.hwe(genli_vcf, p = 0.05, subset = "each",bonf = TRUE)

#       Saving it on my laptop: "HW_each_Figure_1.txt"

#_____________________________________________________________________
          ####  PRINCIPAL COMPONENT ANALYSIS (PCA) ####
#_____________________________________________________________________

#       We will visualize the population structure at an individual level

adeg_pca <- scaleGen(GP_RedFish2, NA.method="mean",scale=F)
pca.adeg_pca <- dudi.pca(adeg_pca, scale=F, nf = 10,scannf = F)
s.class(pca.adeg_pca$li, fac=GP_RedFish2@pop, col=funky(15),cpoint=1)
s.class(pca.adeg_pca$li,fac=GP_RedFish2@pop,  xax=1, yax=2,
        col=transp(funky(12),.9),
        axesel=F, cstar=0,
        cpoint=2)

#       To figure out which locus are differing the most among the populations. 
deg_genind <- gi2gl(GP_RedFish2, parallel = TRUE)

#       Remember to enter the number of axes in the console window. 
pca_genlight <- glPca(deg_genind)
#___________ A-M: How many axes do we want to investigate? How do we know? (5)


#       Loading the plot
loci_axis <- loadingplot(pca_genlight,axis = 1, lab.jitter = 1)




#_____________________________________________________________________
            ####  DEGREE OF DIFFERENTIATION (Fst) ####
#_____________________________________________________________________

#         Measuring of the genetic differentiation between the populations
#         Fst is a measure of how different a population is from another one. 


fst_test <- stamppFst(deg_genind, nboots = 10000, percent = 95)
fst_test$Fsts
fst_test$Pvalues

#         Since they're all 0 they're all significant! This support our assumption
#         of the four different populations which are genetically different. 






