library(adegenet) 
library(gapminder)
library(gdsfmt)
library(ggplot2)
library(gmodels)
library(gridExtra)
library(lme4)
library(lmerTest)
library(poppr)
library(RColorBrewer)
library(SeqArray)
library(SNPRelate)
library(tidyverse)
library(vcfR)

############# PROTEUS SWABBING VS TISSUE ############### 
### DNA quality and quantity from swabs ###
setwd("C:/Users/hans_/Dropbox/Proteus_tissue_swab_ddRADSeq/05_data_sharing")
swabs <- read.csv("qubit_values_swabs_mod2.csv", header=TRUE)
str(swabs)

### FIGURE 2 A ###
par(mfrow=c(2,1))
# show all data
hist(swabs$DNA_swab, xaxt="n", breaks = 50)
axis(side=1, at=seq(0, 120, by=10))
abline(v = c(5), col="red", lwd=2, lty=2)
### FIGURE 2 B ###
## Relationship between size and DNA
# size
#par(mfrow=c(1,1))
summary(lm(DNA_swab~Prot_Length, data=swabs))
plot(DNA_swab~Prot_Length, pch=19, cex=1, col=ifelse(Prot_Length>0,"grey20","grey80"), data=swabs)
reg<-lm(DNA_swab~Prot_Length, data = swabs)
abline(reg, col="orange", lty=2, lwd=2)
# LM
summary(lm(DNA_swab~Prot_Length, na.action=na.omit, data=swabs))

### FIGURE 3 ###
### Nucleotide diversity  and coverage variation between sample types ###
data <- read.csv("Pi_tissue_vs_swabs_outlier.csv", header=T)
str(data)
## check if method explains variation in coverage and pi
summary(mod1<-(aov(coverage~sample_type, data=data)))
summary(mod2<-(aov(Pi~sample_type, data=data)))
## check variation in pi due to samples vs. method
mod1<-(aov(Pi~Lineage+ID+sample_type, data=data))
mod1<-(aov(Pi~ID+sample_type+coverage, data=data))
summary(mod1)
af <- anova(mod1)
afss <- af$"Sum Sq"
print(cbind(af,PctExp=afss/sum(afss)*100))
## make a boxplot
p1 <- data %>%
  ggplot(aes(x=(reorder(sample_type,Pi)), y= coverage, fill=sample_type)) +
  geom_boxplot() +
  geom_point()+ 
  geom_line(aes(group=ID),linetype = "dashed") +
  theme_bw() + 
  theme(text = element_text(size = 20)) +
  theme(panel.border = element_blank(), legend.position = "none", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p2 <- data %>%
  ggplot(aes(x=(reorder(sample_type,Pi)), y= Pi, fill=sample_type)) +
  geom_boxplot() +
  geom_point()+ 
  geom_line(aes(group=ID),linetype = "dashed") +
  ylim(0.0005,0.001)+
  theme_bw() + 
  theme(text = element_text(size = 20)) +
  theme(panel.border = element_blank(), legend.position = "none", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
grid.arrange(p1,p2,ncol=2)

### FIGURE S1 ###
## coverage vs. number of loci
setwd("C:/Users/hans_/Dropbox/Proteus_tissue_swab_ddRADSeq/05_data_sharing")
cov.loci <- read.csv("coverage_vs_loci.csv", header=TRUE)
str(cov.loci)
ggplot(cov.loci, aes(mean_cov,n_loci,colour=factor(method)))+
  geom_point(size=2) + 
  stat_smooth(method = "lm", formula = y ~ splines::bs(x, 3), col = "black") + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
summary(aov(mean_cov~method, data=cov.loci))

### TABLE S3 ###
### ERROR RATE ###
### download error rate calculation scripts from https://github.com/AliciaMstt/RAD-error-rates/blob/master/SNPs_error.R
source("C:/Users/hans_/Desktop/Projects/Proteus_methods_swabbing_vs_tissue/01_data/12_error_rate/PairsDiff.R")
source("C:/Users/hans_/Desktop/Projects/Proteus_methods_swabbing_vs_tissue/01_data/12_error_rate/SNPs_error.R") # read in the R functions, which also calls the needed packages
## DATASET 1: Unguided analysis
setwd("C:/Users/hans_/Dropbox/Proteus_tissue_swab_ddRADSeq/05_data_sharing")
proteus.error.unguide <- read.PLINK("replicates_66_p8_r1.plink.raw", parallel=FALSE)
SNP_error(proteus.error.unguide) 
str(proteus.error.unguide) # N SNPs = 149,678
## DATASET 2: tissue guided analysis
# load vcf files of both swabs and tissues into R and replace and reorder to match them
proteus.error.tissg <- read.PLINK("swab_tissue_guided_p66_sorted.plink.raw", parallel=FALSE)
SNP_error(proteus.error.tissg)
str(proteus.error.tissg)  # N SNPs = 148,780
## DATASET 3: reference genome guided analysis
proteus.error.refg <- read.PLINK("proteus_refmap_p8_r1.raw", parallel=FALSE)
SNP_error(proteus.error.refg)
str(proteus.error.refg) # N SNPs = 279,590
## compare all three
SNP_error(proteus.error.unguide)
SNP_error(proteus.error.tissg)
SNP_error(proteus.error.refg)

### TISSUE-GUIDED ANALYSIS MAPPING RESULTS ###
### FIGURE S2B ###
par(mfrow=c(1,1))
hits <- read.csv("swab_catalogue_loci_bt2_local_results.csv", header=TRUE)
str(hits)
hits_m <- hits[1:3,2:3]
row.names(hits_m)<- c("unique hits across all loci","unique hits across mapped loci", "total mapped loci")
coul <- c("#ffeda0","#feb24c","#f03b20") 
# Transform this data in %
data_percentage <- apply(hits_m, 2, function(x){x*100})
# Make a stacked barplot--> it will be in %!
barplot(data_percentage, col=coul, border="white", legend = row.names(hits_m), ylim= c(0,100), beside = TRUE)

### LOCUS ANNOTATION ###
### FIGURE 4B ###
# mapping to genome: unique
tg.annot <- read.csv("annotation_loci_blast.csv", header=TRUE)
str(tg.annot)
# calculate percentage:
gmodels::CrossTable(tg.annot$type)
tg.annot$type <- factor(tg.annot$type, levels=c("mRNA", "LTR", "LINE","DNA", "SINE"))
# Basic piechart
ggplot(tg.annot, aes(x=factor(1), fill=type))+
  geom_bar(width = 1)+
  coord_polar("y")+
  scale_fill_manual(values=c("#AA4499", "#44AA99", "#88CCEE","#DDCC77","#332288"))+ 
  theme_void() # remove background, grid, numeric labels
aggregate(overlap.length ~ type, data=tg.annot, FUN=mean)
### FIGURE 4C ###
## check genome-guided annotation
# filtered loci that were also used to calculate genotyping errors
rep.loci <- read.csv("replicate_filtered_loci_annotation.csv", header=TRUE)
str(rep.loci)
# calculate percentage:
gmodels::CrossTable(rep.loci$type)
rep.loci$type <- factor(rep.loci$type, levels=c("mRNA", "LTR", "LINE","DNA", "SINE"))
# Basic piechart
ggplot(rep.loci, aes(x=factor(1), fill=type))+
  geom_bar(width = 1)+
  coord_polar("y")+
  scale_fill_manual(values=c("#AA4499", "#44AA99", "#88CCEE","#DDCC77","#332288"))+ 
  theme_void() # remove background, grid, numeric labels
# show mean overlap between types
aggregate(overlap.length ~ type, data=rep.loci, FUN=mean)

### EXOGENOUS DNA ###
## display the average of swab and tissue samples, and all combined
### FIGURE 6A ###
## superkingdoms
setwd("C:/Users/hans_/Desktop/Projects/Proteus_methods_swabbing_vs_tissue/01_data/14_exogenous_reads/02_individual_results/blastNTindividualtaxa_tables")
suking <- read.csv("superkingdoms.csv", header=TRUE)
str(suking)
suking_m <- suking[36:37,3:6]
rownames(suking_m) <- c("Swab_Individuals","Tissue_Individuals")
colnames(suking_m) <- c("Bacteria",	"Eukaryota", "Viruses",	"Archaea")
suking_m <- t(suking_m)
coul <- brewer.pal(4, "Pastel2") 
# Transform this data in %
data_percentage <- apply(suking_m, 2, function(x){x*100/sum(x,na.rm=T)})
# Make a stacked barplot--> it will be in %!
sup <- barplot(data_percentage, col=coul , border="white", legend = TRUE, xlab="group")
### FIGURE 6B ###
## Bacteria
bacteria <- read.csv("bacteria_order.csv", header=TRUE)
str(bacteria)
bacteria_m <- bacteria[36:37,3:15]
rownames(bacteria_m) <- c("Swab_Individuals","Tissue_Individuals")
colnames(bacteria_m) <- c("Unknown", "Other", "Balneolales", "Propionibacteriales", "Micrococcales", "Flavobacteriales", "Bacillales", "Hyphomicrobiales", "Mycobacteriales", "Pseudomonadales", "Burkholderiales", "Enterobacterales", "Moraxellales")
bacteria_m <- t(bacteria_m)
coul <- c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462","#b3de69","#fccde5","#d9d9d9","#bc80bd","#ccebc5","#ffed6f") 
# Transform this data in %
data_percentage <- apply(bacteria_m, 2, function(x){x*100/sum(x,na.rm=T)})
# Make a stacked barplot--> it will be in %!
bac <- barplot(data_percentage, col=coul , border="white", legend = TRUE, xlab="group")
### FIGURE 6C ###
## Chordates
# base R stacked barplot
chordata <- read.csv("chordata_class.csv", header=TRUE)
str(chordata)
chordata_m <- chordata[36:37,3:6]
rownames(chordata_m) <- c("Swab_Individuals","Tissue_Individuals")
colnames(chordata_m) <- c("Unknown", "Birds", "Mammals", "Fishes")
str(chordata_m)
chordata_m <- t(chordata_m)
coul <- brewer.pal(6, "Pastel2") 
# Transform this data in %
data_percentage <- apply(chordata_m, 2, function(x){x*100/sum(x,na.rm=T)})
# Make a stacked barplot--> it will be in %!
chor <- barplot(data_percentage, col=coul , border="white", legend = TRUE, xlab="group")

### CONSPECIFIC DNA ###
## using vcftools heterozygosity file, filtered to include "better" loci with less genotyping errors
data<-read.csv("heterozygosity_swab_tissue_genome_ref_vcftools.csv",header=TRUE)
str(data)
mod2<-lmer(H.O.~mean_cov+type+(1 | ID), data=data)
summary(mod2)
anova(mod2)
