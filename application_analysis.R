##Analysis of metabolite data. 
#We use data of the effect sizes of each SNP on the 118 metabolites combined with the standard error of those SNP exposure associations 
#(extracted from the GWAS results avaliable at http://www.computationalmedicine.fi/data#NMR_GWAS).(1) We also use data on the SNP associations 
#with age related macular degeneration (AMD) from Fritsche et al 2016 (2). 

rm(list = ls(all=TRUE))


#functions defined for this analysis

library(remotes)
#install_github("WSpiller/MRChallenge2019")
library(data.table)
library(knitr)
library(tidyr)
library(dplyr)
library(devtools)
library(readxl)
library(MRChallenge2019)
source("app_functions.R")

dat <- Challenge_dat
dat_se <- data.frame(read.csv("data_incse.txt"))
NMRAdat <- NMRA_dat

names <- NMRAdat$Abbreviation
colnames(dat_se) <- gsub("_", ".", colnames(dat_se))

ids <- as.vector(dat_se$rsid)
row.names(dat_se) <- ids
dat_se <- dat_se[,2:(length(names)+1)]

names <- c("ldl", "hdl", "tg", names)

exp <- subset(dat, select=c(1,9,12,15,18,32:149))
pvals <- subset(dat, select=c(11,14,17,150:267))
colnames(exp) <- sub("beta_","",colnames(exp))
names(exp)[names(exp) == 'acAce'] <- 'AcAce'
colnames(pvals) <- sub("p_","",colnames(pvals))

ids <- exp$rsid
row.names(exp) <- ids
row.names(pvals) <- ids

dat_se <- data.frame(dat$se_amd, dat$se_ldl, dat$se_hdl, dat$se_tg, dat_se)
colnames(dat_se) <- gsub("dat.se_", "", colnames(dat_se))

Fstat <- data.frame()
for(x in 1:length(names)){
  for(y in 1:length(ids)){
    
    Fstat[ids[y],names[x]] <- (exp[ids[y],names[x]]/dat_se[ids[y],names[x]])^2
    
  }
}

#import and sort out correlations (NB - correltations are calculated from ALSPAC data and therefore not currently publicly avaliable)

correlations <- read_excel("correlations.xlsx")
correlations <- data.frame(correlations)
row.names(correlations) <- correlations[,1]
correlations[,1] <- NULL

#calculate the exposures with the most SNPs with an F>5 then keep all snps with individual F>5 for at least one of those exposures.

F.ind <- Fstrong(names[4:length(names)])

F.ind <- F.ind[order(-F.ind$no.snps),]

topexp <- row.names(F.ind[1:14,])


F.MR <- data.frame(Fstat[,topexp])
ex.MR <- data.frame(exp[,topexp])

maxF_row <- apply(F.MR,1,function(x) max(as.numeric(x)))
keep <- as.vector(as.numeric(maxF_row > 5))

ex.MR <- ex.MR[,1:length(topexp)]*keep
ex.MR[ex.MR == 0] <- NA

MR.all <- (summary(lm(dat$beta_amd~ -1 + ., data = ex.MR, weights = (dat$se_amd)^2)))$coefficients

b<- Fstrong(topexp)
c <- conditionalF(topexp)

Ftop <- data.frame(c, b)
colnames(Ftop)[4] <- "Ind.F.Stat" 
colnames(Ftop)[5] <- "No.snps.Ind" 


#subset by type
subset_b <- c("IDL.PL", "IDL.P", "IDL.TG", "IDL.L")
subset_c <- c("L.LDL.L", "L.LDL.P", "M.LDL.P")
subset_d <- c("S.VLDL.PL",  "S.VLDL.C",   "S.VLDL.FC")  
subset_e <- c("XS.VLDL.L", "XS.VLDL.TG", "XS.VLDL.P")  



Fstrong_b <- Fstrong(subset_b)
Fcond_b <- conditionalF(subset_b)
F.setb <- data.frame(Fcond_b,Fstrong_b)

Fstrong_c <- Fstrong(subset_c)
Fcond_c <- conditionalF(subset_c)
F.setc <- data.frame(Fcond_c,Fstrong_c)

Fstrong_d <- Fstrong(subset_d)
Fcond_d <- conditionalF(subset_d)
F.setd <- data.frame(Fcond_d,Fstrong_d)

Fstrong_e <- Fstrong(subset_e)
Fcond_e <- conditionalF(subset_e)
F.sete <- data.frame(Fcond_e,Fstrong_e)


##MR for the final set of exposures
subexp <- c("XS.VLDL.P", "S.VLDL.PL", "L.LDL.L", "IDL.TG")
subexp_se <- c("XS.VLDL.P_se", "S.VLDL.PL_se", "L.LDL.L_se", "IDL.TG_se")
subexp_f <- c("XS.VLDL.P_f", "S.VLDL.PL_f", "L.LDL.L_f", "IDL.TG_f")

F.MR <- data.frame(Fstat[,subexp])
ex.MR <- data.frame(exp[,subexp])

maxF_row <- apply(F.MR,1,function(x) max(as.numeric(x)))
keep <- as.vector(as.numeric(maxF_row > 5))

ex.MR <- ex.MR[,1:length(subexp)]*keep
ex.MR[ex.MR == 0] <- NA
MR.subset <- summary(lm(dat$beta_amd~ -1 + ., data = ex.MR, weights = (dat$se_amd)^-2))$coefficients
conditionalF(subexp)
Fstrong(subexp)

kx <- length(subexp)

#jackknife

analysis.dat_all <- data.frame(exp[,c("amd",subexp)])
analysis.dat_all <- data.frame(cbind(analysis.dat_all, data.frame(dat_se[,c("amd",subexp)]), data.frame(Fstat[,c(subexp)])))
names(analysis.dat_all) <- c("amd",subexp, "amd_se", subexp_se, subexp_f)
F.analysis <- analysis.dat_all[,c(subexp_f)]
maxF_row <- apply(F.analysis,1,function(x) max(as.numeric(x)))
keep <- as.vector(as.numeric(maxF_row > 5))

analysis.dat_all <- analysis.dat_all[,1:length(c("amd",subexp, "amd_se", subexp_se, subexp_f))]*keep
analysis.dat_all[analysis.dat_all==0] <-NA
analysis.dat_all <- na.omit(analysis.dat_all)

analysis.dat <- analysis.dat_all

MR.results <-  MRfunction_jk(subexp)

results <- MR.results


for(s in 1:68){
  analysis.dat <- analysis.dat_all
  analysis.dat[s,] <- NA
  analysis.dat <- na.omit(analysis.dat)
  
  temp <-  MRfunction_jk(subexp)
  results <- rbind(results, temp)
  
}

