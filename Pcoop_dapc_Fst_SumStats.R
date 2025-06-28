library(adegenet)
library(pegas)
library(ape)
library(vcfR)
library(SNPRelate)
library(hierfstat)
library(diveRsity)
library(dartR)

#setwd("Path/to/input/files")

########################### load in data and assign population to genind object ###########

#read in data
data<-read.genepop("pcoop.gen")
#how many populations are there?
nPop(data)
# how many individuals are there?
nInd(data)
# how many loci? 
nLoc(data)
# what are the populationS
data@pop
# name pops
PopNames <- c("OHR", "TNR")
# set pop names in genind object 
popNames(data) <- PopNames
popNames(data)
# try again 
data@pop
# summary of the data
data 

##################################################################################
######################## estimate K means clusters ###################################
##################################################################################

# For this dataset one is the most likely cluster so there is not a whole lot to be done here ass far as displaying the data 
# without as prior population/river grouping
grp1<-find.clusters(data, max.n.clust = 5)
25
1
# explore "grp1" data. See what the different columns of information are. 
names(grp1)
# Kstat, this will return the BIC value for each of the 25 possible clusters you created using find.cluster. 
grp1$Kstat
# stat, the BIC for your choosen k value
grp1$stat
# grp, which cluster (1 -4) each indivudal is assigned to 
grp1$grp
# size, number of individuals assigned to each cluster 
grp1$size


###this wont show much
dapc.fnd.clus<-dapc(data, grp1$grp)
25
1

scatter(dapc.fnd.clus, bg="white", scree.da=FALSE, scree.pca = FALSE,
        posi.pca = "topright", legend=TRUE, posi.leg = "topleft", solid=.5, cstar = 0, 
        clabel = 0, cellipse = 2.5)


##################################################################################################
################## DAPC with prior grouping assignments ###########################################
##################################################################################################

dapc.pop<-dapc(seed = 99, data, data$pop)
25
1
# the "optim.a.score" will take the information you just saved into the dapc.np.5 object and determine the 
# optimimum number of PCs to retain 
temp1<-optim.a.score(seed = 99, dapc.pop, n.sim = 20)
# the optimimum output is displayed graphifically or can be display hitting enter below
temp1$best

# rerun the DAPC using the optimial number of PCs as determined in "temp1" i.e. temp1$best
# again retain all discriminant functions, which should be 3
dapc.pop<-dapc(data, data$pop, n.pca = temp1$best)
1
# plot the results of the dapc. 
pcoop.col<-c('darkorange','purple')


set.seed(909)
scatter(dapc.pop, bg="white", scree.da=FALSE, scree.pca = FALSE,
        posi.pca = "topright", legend=TRUE, posi.leg = "topleft", solid=1, cstar = 0, 
        clabel = 0, cellipse = 2.5, col=pcoop.col)


########################## Fst #################################

# convert genind to hierstat data frame
hierdata<-genind2hierfstat(data, pop = data@pop)

# generate pairwise Fst values; this will take a few minutes with a standard RADseq SNP dataset
fst <- pairwise.WCfst(hierdata [,-2], diploid = TRUE)
fst

# bootstrap fst values for static significance; below 95CI. This might also take a few minutes to run. 
bootsfst<-hierfstat::boot.ppfst (hierdata [,-2], nboot = 1000, quant = c(0.025, 0.975), diploid = TRUE)
bootsfst


write.csv(bootsfst, file="pcoop_fst_boots.csv")

write.csv(fst, file="pcoop_fst.csv")


#########################################################################################
################## Calculate general summary statistics #################################
#########################################################################################

#### calculate most summary statistics with DartR package ###

# first convert genind to genlight object
gl <- gi2gl(data) 

# based on above analyses, we conclude there is no genetic differentiation between collection sites (rivers).
# so we need to combine the collection sites back together into a single population
gl.all<-gl.merge.pop(gl, old = c('OHR', 'TNR'), new = 'all') 

# calculate sumstats
df <- gl.report.heterozygosity(gl.all,method='pop')

#### calculate Allelic Richness with DiveRsity package
#### The Allelic Richness calculation code was modified from code published on Nathan Whelan's GitHub
#### https://github.com/nathanwhelan

##Reload in GENEPOP FILE.
test_results<-divBasic(infile="pcoop.gen",gp=2)

##Split out Allelic Richness, make into table.  NOTE: Could modify here (and later with variable name) to split out other statistics
test_results$Allelic_richness
alleicRicheness<-test_results$Allelic_richness
allelicRichness<-as.table(alleicRicheness)

alleicRicheness
dat<-alleicRicheness[-c(3528), ]
tail(dat)

dat$ar<-dat$pc10

##Data check. Note overall values.
allelicRicheness
allelicRichness["overall",]

tail(alleicRicheness)

allelicRichness.df<-as.data.frame(alleicRicheness)

dat<-allelicRichness[-c(3528), ]
tail(dat)
write.csv(dat, file = "ar.csv" )