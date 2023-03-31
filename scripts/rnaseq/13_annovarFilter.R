###################################################################
#                  Annovar Variants filters
###################################################################
library(tidyverse)

#vcf files directory
path = "/srv/data/my_shared_data_folder/epilepsy/rnaseq_results/vcf_anno/"
SAMPLI_ID = list.files(path = path, pattern = ".txt")

#-----------------------------------------------------

#CADD, SIFT, PolyPhen2 and FATHMM prediction clinvar
pred_filter = list("clinvar"=c(".", "Uncertain_significance","Likely_pathogenic","Pathogenic/Likely_pathogenic",
            "Pathogenic", "risk_factor"), "silicoPred"=c(".","P","D"))
#------------------------------------------------------

for (FILES in SAMPLI_ID) {
  
#Exporting annovar annotation .txt or csv file
data = read.delim(FILES, header = T, sep = "\t", quote = NULL)

#Filtering variant that PASS quality check and read deeph
genePass = filter(data, data$Otherinfo10 == "PASS" & data$Otherinfo9 > 15)

#Will looking only for exonic and Non synonymous variation
genFun = filter(genePass, genePass$Func.refGene == "exonic" & genePass$ExonicFunc.refGene == "nonsynonymous SNV")

#SIFT, PolyPhen2 and FATHMM prediction clinvar
funPred = filter(genFun, genFun$SIFT_pred %in% pred_filter$silicoPred & genFun$Polyphen2_HDIV_pred %in% pred_filter$silicoPred & genFun$FATHMM_pred %in% pred_filter$silicoPred, genFun$CLNSIG %in% pred_filter$clinvar)

#CADD filter
caddPref = filter(funPref, genFun$CADD_phred >= 20)
    
# Filter base on allele frequencies in all population less than 1%
varAF = filter(caddPred, funPred$ExAC_ALL <= 0.01)

write.table(varAF, file = paste(unlist(strsplit(FILES, split = "\\."))[1], ".annovarFilt.txt", sep = ""), row.names = F, sep = "\t", quote = F)
}


########    Merge files     #########
samplesFilt = list.files(path = path, pattern = "annovarFilt.txt$")

#Temporaly fille that will contain merge filles
temp_data = data.frame(ModelName = character(), Object = character(), stringsAsFactors = F)

for (f in 1:length(samplesFilt)) {
  #read the file
  currentFile = read.delim(samplesFilt[f], sep = "\t", quote = NULL)
  
  #Add Accession column with SRR ID
  sample_id = data.frame("Accession"=c(rep(unlist(strsplit(samplesFilt[f], split = "\\."))[1], nrow(currentFile))))
  
  currentFile = cbind(sample_id,currentFile)
  
  #Append current file to temp_data
  temp_data = rbind(temp_data, currentFile)
  
}

# Importing final merge variants call for all samples
write.table(temp_data, file = "Epilepsy_var.txt", row.names = F, sep = "\t", quote = F)