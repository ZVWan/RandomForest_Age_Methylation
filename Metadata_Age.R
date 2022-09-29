############################
#                          #
#   Metadata collection    #
#      Age     Gender      #
#                          #
############################

#Adjust for each dataset for further analysis
setwd("~/Data/Whole_Blood")
dir.create("GSE67444_Data")
setwd("~/Data/Whole_Blood/GSE67444_Data")
getGEOSuppFiles("GSE67444")
GSE67444_Meta <- getGEO("GSE67444")
save(GSE67444_Meta, file = "GSE67444.RData")
metadata <- as.data.frame(pData(GSE67444_Meta[[1]]))
names(metadata)
# Export
write.csv(metadata, "GSE67444_meta.csv")
# Age histogram
setwd("~/wz/ZWAge/Age_Sex/Data/Whole_Blood/GSE67444_Data")
GSE67444_Metadata <- as.data.frame(read.csv("GSE67444_meta.csv"),sheet=1)
colnames(GSE67444_Metadata)
hist(GSE67444_Metadata$age.ch1,las = 1)
table(GSE67444_Metadata$age.ch1)
age_count_GSE67444 <- count(GSE67444_Metadata,c("age.ch1"))
sex_count_GSE67444 <- count(GSE67444_Metadata,c("Sex.ch1"))
counts <- table(GSE67444_Metadata$Sex.ch1,GSE67444_Metadata$age.ch1)
barplot(counts, main = 'frequency by Age & Sex', xlab = "Age", legend = rownames(counts),beside = T)

