bash ./SCRIPT_combine_variants
bash ./SCRIPT_annotate_Annovar
bash ./SCRIPT_annotate_SnpEff
$(which R) CMD BATCH combine_annotations.R 
