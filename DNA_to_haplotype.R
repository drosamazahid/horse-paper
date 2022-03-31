#DNA to Haplotype

library(microseq)
library(plyr)
library(tidyverse)
library("writexl")
library(dplyr)
library("seqRFLP")
library(DECIPHER)
library(readxl)

setwd("D:/Osama/UK_tree")
list <- sub(".fasta", "", list.files(pattern = ".fasta"))
dir.create("aggregated")
dir.create("trimed")
dir.create("final")

#aggregate loop
for (i in list) {
  loadfile <- file.path(paste0(i, ".fasta"))
  Data <- subset(readFasta(loadfile),select = 2,)
  Data2 <- aggregate(list(numdup=rep(1,nrow(Data))), Data, length)
  Data3 <- Data2[order(-Data2$numdup),]
  Data3$name <- i
  Data3$number <- c(1:nrow(Data3))
  Data3 <- Data3 %>%
    unite("Header", name:number, remove = TRUE)
  Data3 <- Data3[, c(3,1,2)]
  Data3 <- Data3 %>% filter(numdup > 5) %>% head(1000)
  excelfile <- file.path("aggregated", paste0(i, ".xlsx"))
  write_xlsx(Data3, excelfile)
  fastafile <- file.path("aggregated", paste0(i, ".fasta"))
  Data3.fasta = dataframe2fas(Data3[,c(1:2)], file=fastafile)}


#trim loop
for (i in list) {
  loadfile <- file.path("aggregated", paste0(i, ".fasta"))
  seqs <- RemoveGaps(readDNAStringSet(loadfile, format = "fasta"))
  aligned <- AlignSeqs(seqs, iterations = 5, refinements = 5)
  adjusted <- AdjustAlignment(aligned)
  trimed <- TrimDNA(adjusted, "ACGTCT", "AACTAA", type = "sequences")
  trimedfile <- file.path("trimed", paste0(i, ".fasta"))
  writeXStringSet(trimed, file=trimedfile)}



#final loop

for (i in list) {
  loadfile <- file.path("trimed", paste0(i, ".fasta"))
  Data <- subset(readFasta(loadfile),select = 2,)
  Data2 <- aggregate(list(numdup=rep(1,nrow(Data))), Data, length)
  Data3 <- Data2[order(-Data2$numdup),]
  Data3$name <- i
  Data3$number <- c(1:nrow(Data3))
  Data3 <- Data3 %>%
    unite("Header", name:number, remove = TRUE)
  Data3 <- Data3[, c(3,1,2)]
  excelfile <- file.path("final", paste0(i, " final", ".xlsx"))
  write_xlsx(Data3, excelfile)
  Data3 <- Data3 %>% filter(numdup > 1)
  fastafile <- file.path("final", paste0(i, ".fasta"))
  Data3.fasta = dataframe2fas(Data3[,c(1:2)], file=fastafile)}

#all loop
for (i in list) {
  loadfile <- file.path("final", paste0(i, ".fasta"))
  assign(i, readDNAStringSet(loadfile))
  aligned <- AlignProfiles(Haemonchidae_unclassified, Rhabditida_unclassified)
  adjusted <- AdjustAlignment(aligned)
  adjustedfile <- file.path("final", "all.fasta")
  writeXStringSet(adjusted, file=adjustedfile)}

file.copy(adjustedfile, "K:/sargison_grp/Mike_Osama_shared/NCBI/blast-2.12.0+/bin", overwrite = TRUE)
writeLines("cd \\cmvm.datastore.ed.ac.uk\cmvm\eb\groups\sargison_grp\Mike_Osama_shared\NCBI\blast-2.12.0+\bin", "outfile.txt")
