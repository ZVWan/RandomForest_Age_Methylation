###########################
#     Quality Control     #
#     Pre-processing      #
###########################

##### Adjust for datasets

##### Quality control - Minfi
setwd("~/wz/ZWAge/Age_Sex/Data/Whole_Blood/GSE67444_Data/GSE67444")
untar("GSE67444_RAW.tar")
fileNames <- list.files()
sapply(fileNames, gunzip)
RGSet <- read.metharray.exp("~/wz/ZWAge/Age_Sex/Data/Whole_Blood/GSE67444_Data/GSE67444")        #Change the dataset here.
MSet <- preprocessRaw(RGSet)                                                        #A MethylSet objects contains only the methylated and unmethylated signals
RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)                            #A RatioSet object is a class designed to store Beta values and/or M values
qc <- getQC(MSet)
setwd("~/wz/ZWAge/Age_Sex/Data/Whole_Blood/GSE67444_Data")
qcReport(RGSet, pdf= "qcReport.pdf")
plotQC(qc)
controlStripPlot(RGSet, controls="BISULFITE CONVERSION II")

##### Sex prediction
GRset <- mapToGenome(RSet)
predictedSex <- getSex(GRset, cutoff = -2)$predictedSex

##### Normalization
RGset.quantile <- preprocessQuantile(RGSet, removeBadSamples = TRUE)
RGset.funnorm <- preprocessFunnorm(RGSet)
RGSet.erase <- (colnames(RGset.funnorm) %in% colnames(RGset.quantile))
RGset.erase1 <- RGset.funnorm[,RGSet.erase]

##### Obtain Probes on X&Y chromosome
setwd("~/wz/ZWAge/Age_Sex/Data/Whole_Blood/GSE67444_Data")
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
keep_X_Y <- (featureNames(RGset.erase1) %in% ann450k$Name[ann450k$chr %in%  c("chrX","chrY")])
X_Y_chromosome <- RGset.erase1[keep_X_Y,]

##### Removal of SNPdb 
snps <- getSnpInfo(X_Y_chromosome)
removed1 <- dropLociWithSnps(X_Y_chromosome, snps=c("SBE","CpG"), maf=0)
Sex_chromosome_SNP <- removed1

##### Merge 2017 probes
#For Whole blood, use:
keep_2017 <- (featureNames(RGset.funnorm) %in% ann450k$Name[ann450k$Name %in% c("cg16867657","cg12934382","cg11807280","cg02872426","cg06874016","cg08097417")])
#For buffy coat, use:,
#keep_2017 <- (featureNames(RGset.erase1) %in% ann450k$Name[ann450k$Name %in% c("cg16867657","cg08097417","cg06784991","cg07553761","cg16054275","cg18473521")])
probes_2017<- RGset.erase1[keep_2017,]

##### Saving beta value (DNA methylation level)
setwd("~/wz/ZWAge/Age_Sex/Data/Whole_Blood/GSE67444_Data")
sex_chromosome_beta <- getBeta(Sex_chromosome_SNP)
probes_beta <- getBeta(probes_2017)
remainprobes <- rbind(sex_chromosome_beta, probes_beta)
write.csv(remainprobes,file = "remainprobes67444.csv")

##### Removal of Cross-hybridize, Chen et al. 2013
probe3 <- read.csv("~/wz/ZWAge/Age_Sex/48639-non-specific-probes-Illumina450k-cg+ch.csv"); 
delete3 <- probe3[,1]
write.csv(x = delete3, file = "deleteprobe3.csv")
filter_probe3 <- read.table("deleteprobe3.csv", header = TRUE, sep = ',')
dataset_crhy <- read.table("beta_GSE67444_sex_snp.csv", header = TRUE, sep = ',')
beta_GSE67444_removed_3 <- dataset_crhy[!(dataset_crhy$X %in% filter_probe3$x),]
write.csv(x = beta_GSE67444_removed_3, file = "beta_GSE67444_removed_sex_snp_crosshybr.csv")

##### Removal of blood sub cell type specific pattern, Jaffe et al.
probe4 <- read.csv("~/wz/ZWAge/Age_Sex/Supplementary_Table_2.csv")
delete4 <- filter_at(probe4, vars(starts_with("p.value")),all_vars(.<0.05))
write.csv(x = delete4, file = "deleteprobe4.csv")
beta_GSE67444_removed_4 <- beta_GSE67444_removed_3[!(beta_GSE67444_removed_3$X %in% delete4$Name),]
write.csv(x = beta_GSE67444_removed_4, file = "beta_GSE67444_removed_sex_snp_crosshybr_blood.csv")

##### Removal of p detection value <0.01
Pdetectionvalue <- detectionP(RGSet, type = "m+u")
write.csv(x = Pdetectionvalue, file = "Pdetectionvalue_67444.csv")
filter_failure <- (rowMeans(Pdetectionvalue>0.01)>0.2)
filtered_probes <- Pdetectionvalue[filter_failure,]
dataset_5 <- read.table("beta_GSE67444_removed_sex_snp_crosshybr_blood.csv", header = TRUE, sep = ',')
filtered_probes <- as.data.frame(cbind(Row.Names = rownames(filtered_probes), filtered_probes))
beta_GSE67444_removed_5 <- dataset_5[!(dataset_5$X %in% filtered_probes$Row.Names),]
write.csv(x = beta_GSE67444_removed_5, file = "beta_GSE67444_removed_final.csv")

##### Attract X&Y probes separately
keep_Y_final <- beta_GSE67444_removed_5$X %in% ann450k$Name[ann450k$chr %in%  c("chrY")]
keep_X_final <- beta_GSE67444_removed_5$X %in% ann450k$Name[ann450k$chr %in%  c("chrX")]
Y_chromosome_probes_final <- beta_GSE67444_removed_5[keep_Y_final,]
X_chromosome_probes_final <- beta_GSE67444_removed_5[keep_X_final,]
write.csv(x = Y_chromosome_probes_final, file = "Y_chromosome_probes_final.csv")
write.csv(x = X_chromosome_probes_final, file = "X_chromosome_probes_final.csv")

