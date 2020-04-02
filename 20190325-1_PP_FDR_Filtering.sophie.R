library("tidyverse")
library("readxl")
library(stringr)


#create files of FDR report #PUTNAMEOFTHEFILE
FDR_thresh0 <- read_excel("Data/Ryegrass_Poaceae_FDR.xlsx", sheet = "Protein FDR Summary")
                      
Protein_summary <-read_excel("Data/Ryegrass_Poaceae_FDR.xlsx", sheet = "Protein Summary")
Peptide_summary <-read_excel("Data/Ryegrass_Poaceae_FDR.xlsx", sheet = "Peptide Summary")

#Select FDR Threshold of unused
FDR_threshnumber <-  select (FDR_thresh0,c=19)
FDR5 <- FDR_threshnumber [12,]
#FDR1 <-  FDR_threshnumber [11,] #If you want 1% un lock this line

#make FDR numeric
FDR5 <- as.numeric(FDR5$c)

####Proteins
#filter FDR
Protein_list5FDR <- filter(Protein_summary, Unused >= FDR5)


#filter proteins with more than 1 peptide
Protein_over1peptide <- Protein_list5FDR %>% 
  filter_at(10, all_vars(. > 1))

#Filter take out HUMAN BOVIN PIG
Protein_filter_specie<- filter(Protein_over1peptide, !grepl
                             ('HUMAN|BOVIN|PIG', Accession))

#Filter take out Reversed Proteins
Protein_final_list <- filter(Protein_filter_specie, !grepl
                             ('RRRR', Accession))

#Export
write.csv (Protein_final_list, "Results/Protein_final_list.csv")



