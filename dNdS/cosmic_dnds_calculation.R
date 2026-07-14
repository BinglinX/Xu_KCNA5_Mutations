library(dndscv)
library(dplyr)


cosmic_data <- read.delim("/data/cosmic_LI_AC_for_dnds.tsv")

cosmic_data_unique <-  cosmic_data %>% distinct(SAMPLE_NAME, CHROMOSOME, GENOME_START, 
                                                GENOMIC_MUT_ALLELE, .keep_all = TRUE)

coadread_data_COSMIC <- as.data.frame(cbind(
  cosmic_data_unique$SAMPLE_NAME,
  cosmic_data_unique$CHROMOSOME,
  cosmic_data_unique$GENOME_START,
  cosmic_data_unique$GENOMIC_WT_ALLELE,
  cosmic_data_unique$GENOMIC_MUT_ALLELE))

colnames(coadread_data_COSMIC) <- c('sampleID','chr','pos','ref','mut')
coadread_data_COSMIC[coadread_data_COSMIC == "None"] <- ""

coadread_data_COSMIC$chr <- as.character(coadread_data_COSMIC$chr)
coadread_data_COSMIC$pos <- as.integer(coadread_data_COSMIC$pos)

dndsout = dndscv(coadread_data_COSMIC,
                refdb = "hg38",
                max_muts_per_gene_per_sample = Inf,
                max_coding_muts_per_sample = Inf,
                outmats = T)

sel_cv = dndsout$sel_cv
ci = geneci(dndsout)

write.csv(sel_cv, file = "dnds_coad_cosmic.csv")

