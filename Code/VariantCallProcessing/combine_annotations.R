data_folder = "../../Data/Processed/FreeBayesCalls"
# Read variant annotation from Annovar
d = read.table(file.path(data_folder, "OVCA_TargetedSeq_BIH_FreeBayes_calls_merged_ANNOVAR.hg19_multianno.txt"), sep="\t", fill=T, header=F, skip = 1, stringsAsFactors = F, na.strings = c("NA", "."))
hh = read.table(file.path(data_folder, "OVCA_TargetedSeq_BIH_FreeBayes_calls_merged_ANNOVAR.hg19_multianno.txt"), sep="\t", fill=T, header=F, skip = 0, nrows = 1, stringsAsFactors = F)
colnames(d)[1:ncol(hh)] = hh[1,]

# Read variant annotation from SnpEff
d1 = read.table(file.path(data_folder, "OVCA_TargetedSeq_BIH_FreeBayes_calls_merged_SnpEff_table.txt"), sep="\t", stringsAsFactors = F, header=T)
rl = nchar(d1$RefAllele)
al = nchar(d1$VariantAllele)
d1$vartype = ifelse(rl==al & rl==1, "SNV", ifelse(rl<al, "INS","DEL"))

# Get sample labels
library(readxl)
samples = as.data.frame(read_excel("./BIH_to_Label.xlsx"))
rownames(samples) = make.names(samples$`BIH-ID`)

ci = which(colnames(d1) %in% rownames(samples))
colnames(d1)[ci] = apply( samples[colnames(d1)[ci],c("Label","BIH-ID")],1, paste, collapse="-")

# Merge Annovar and SnpEff annotations on same variants
d$varid = with(d, paste(Chr, Start, Ref, Alt, sep="_"))
d1$varid = with(d1, ifelse(vartype=="SNV", paste(Chromosome, Position, RefAllele, VariantAllele, sep="_"),
                    ifelse(vartype=="INS", paste(Chromosome, Position, "-", substr(VariantAllele,2,nchar(VariantAllele)), sep="_"), paste(Chromosome, Position+1, substr(RefAllele,2,nchar(RefAllele)), "-", sep="_"))))
                           
m = merge(d, d1, by="varid", all.x=T, all.y = F)

cols_to_keep = colnames(d1)
cols_to_keep = append(cols_to_keep, c("esp6500si_all","1000g2015aug_all","MetaSVM_score","MetaSVM_rankscore","MetaSVM_pred","MetaLR_score", "MetaLR_rankscore","MetaLR_pred","M-CAP_score","M-CAP_rankscore","M-CAP_pred","CADD_raw","CADD_raw_rankscore","CADD_phred","CLNALLELEID","CLNDN","CLNDISDB","CLNREVSTAT","CLNSIG", "varid"))

tab_final = m[, cols_to_keep]

# Mark variants with less than 1% minor allele frequency in general population based on 1000G and ESP6500 data
tab_final$rare_in_healthy=with(tab_final, ifelse(is.na(esp6500si_all) & is.na(tab_final$`1000g2015aug_all`), T, ifelse((!is.na(esp6500si_all) & esp6500si_all < 0.01) | (is.na(esp6500si_all) & tab_final$`1000g2015aug_all` < 0.01 ), "yes, <1% of population", "no")))
# Describe impact on protein function based on 3 predictive pipelines
tab_final$predicted_impact = apply(tab_final[, c("MetaSVM_pred", "MetaLR_pred", "M-CAP_pred")], 1, function(x) {x[is.na(x)]<-"."; paste(x,  collapse="|") })

write.table(tab_final, file=file.path(data_folder, "/OVCA_TargetedSeq_BIH_FreeBayes_calls_merged_SnpEff_and_Annovar_table.txt"), sep="\t", row.names = F, quote=F)
