library(adegenet)
library(pegas)
library(ape)
library(vcfR)
library(SNPRelate)
library(hierfstat)
library(diveRsity)
library(dartR)


#setwd ("/path/to/input/files")

# read in data and assign population to genind object
data<-read.genepop("pcyph.gen")
#how many populations are there?
nPop(data)
# how many individuals are there?
nInd(data)
# how many loci? 
nLoc(data)
# what are the populationS
data@pop
# name pops
PopNames <- c("MSR","TNR")
# set pop names in genind object 
popNames(data) <- PopNames
popNames(data)
# try again 
data@pop
# summary of the data
data 

##################################################################################
######################## estimate K means clusters ###############################
##################################################################################

# For this dataset one is the most likely cluster so there is not a whole lot to be done here ass far as displaying the data 
# without as prior population/river grouping
grp1<-find.clusters(data, max.n.clust = 5)
20
2
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
20
1

temp1<-optim.a.score(seed = 99, dapc.fnd.clus, n.sim = 20)
# the optimimum output is displayed graphifically or can be display hitting enter below
temp1$best
# see the percent variation explained by each PC. remember you are only keeping the first 9 though for 
# this analysis
#dapc.NoPrio.4$pca.eig/sum(dapc.NoPrio$pca.eig)*100
# rerun the DAPC using the optimial number of PCs as determined in "temp1" i.e. temp1$best
# again retain all discriminant functions, which should be 3
dapc.fnd.clus<-dapc(data, data$pop, n.pca = temp1$best)
1

pcyph.col<-c('royalblue','green')
scatter(dapc.fnd.clus, bg="white", scree.da=FALSE, scree.pca = FALSE,
        posi.pca = "topright", legend=TRUE, posi.leg = "topleft", solid=1, cstar = 0, 
        clabel = 0, cellipse = 2.5, col=pcyph.col)


##################################################################################################
################## DAPC with prior grouping assignments ##########################################
##################################################################################################

dapc.pop<-dapc(seed = 99, data, data$pop)
24
1
# the "optim.a.score" will take the information you just saved into the dapc.np.5 object and determine the 
# optimimum number of PCs to retain 
temp1<-optim.a.score(seed = 99, dapc.pop, n.sim = 20)
# the optimimum output is displayed graphifically or can be display hitting enter below
temp1$best

# rerun the DAPC using the optimial number of PCs as determined in "temp1" i.e. temp1$best
# again retain all discriminant functions, which should be 3
dapc.pop<-dapc(data, data$pop, n.pca = 2)
1
# plot the results of the dapc. 
set.seed(909)
scatter(dapc.pop, bg="white", scree.da=FALSE, scree.pca = FALSE,
        posi.pca = "topright", legend=TRUE, posi.leg = "topleft", solid=1, cstar = 0, 
        clabel = 0, cellipse = 2.5, col=pcyph.col)

#########################################################################################
#################### Calculate Fst between populations identified above #################
#########################################################################################

# convert genind to hierstat data frame
hierdata<-genind2hierfstat(data, pop = data@pop)

# generate pairwise Fst values; this will take a few minutes with a standard RADseq SNP dataset
fst <- pairwise.WCfst(hierdata [,-2], diploid = TRUE)
fst


# bootstrap fst values for static significance; below 95CI. This might also take a few minutes to run. 
bootsfst<-hierfstat::boot.ppfst (hierdata [,-2], nboot = 1000, quant = c(0.025, 0.975), diploid = TRUE)
bootsfst

#########################################################################################
################## Calculate general summary statistics #################################
#########################################################################################

# calculate most summary statistics with DartR package 

gl <- gi2gl(data) #convert to genlight object

df <- gl.report.heterozygosity(gl)
df <- gl.report.heterozygosity(gl,method='pop')

df <- gl.report.heterozygosity(gl)


#### calculate Allelic Richness with DiveRsity package
#### The Allelic Richness calculation code was modified from code published on Nathan Whelan's GitHub
#### https://github.com/nathanwhelan

##Reload in GENEPOP FILE.
test_results<-divBasic(infile="pcyph.gen",gp=2)

##Split out Allelic Richness, make into table.  NOTE: Could modify here (and later with variable name) to split out other statistics
test_results$Allelic_richness
alleicRichenes<-test_results$Allelic_richness
allelicRichness<-as.table(alleicRichenes)

##Data check. Note overall values.
allelicRichness
allelicRichness["overall",]

##seperate by population and remove loci that were not genotyped for that population (i.e., Allelic Richness of 0)
pop1_allelicR<-allelicRichness[,1]
pop1_allelicR<-as.data.frame(pop1_allelicR)
pop1_allelicR_zeroesRemoved<-pop1_allelicR[apply(pop1_allelicR!=0,1,all),]

pop2_allelicR<-allelicRichness[,2]
pop2_allelicR<-as.data.frame(pop2_allelicR)
pop2_allelicR_zeroesRemoved<-pop2_allelicR[apply(pop2_allelicR!=0,1,all),]

#Calculate Mean
pop1mean<-mean(pop1_allelicR_zeroesRemoved)
pop2mean<-mean(pop2_allelicR_zeroesRemoved)

##Calculate Sstandard Deviation
pop1sd_Ar<-sd(pop1_allelicR_zeroesRemoved)
pop2sd_Ar<-sd(pop2_allelicR_zeroesRemoved)



















