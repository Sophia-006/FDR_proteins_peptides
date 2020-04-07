library("tidyverse")
library("readxl")
library(stringr)


#create files of FDR report #PUTNAMEOFTHEFILE
FDR_thresh0 <- read_excel("Data/Ryegrass multi-file pooideae__FDR.xlsx", sheet = "Protein FDR Summary")
                      
Protein_summary <-read_excel("Data/Ryegrass multi-file pooideae__FDR.xlsx", 
                             sheet = "Protein Summary")
Peptide_summary <-read_excel("Data/Ryegrass multi-file pooideae__FDR.xlsx", 
                             sheet = "Distinct Peptide Summary")
Conf_thresh_pep <- read_excel ("Data/Ryegrass multi-file pooideae__FDR.xlsx", 
                               sheet ="Distinct Peptide FDR Summary")

#Select FDR Threshold of unused
FDR_threshnumber <-  select (FDR_thresh0,c=19)
FDR5 <- FDR_threshnumber [12,]
#FDR1 <-  FDR_threshnumber [11,] #If you want 1% un lock this line

#make FDR numeric
FDR5 <- as.numeric(FDR5$c)
#FDR1 <- as.numeric(FDR1$c) #For 1 %

####Proteins

#filter FDR
Protein_list5FDR <- filter(Protein_summary, Unused >= FDR5)
#Protein_list1FDR <- filter(Protein_summary, Unused >= FDR1)

#filter proteins with more than 1 peptide
Protein_5FDR_1pep <- Protein_list5FDR %>% 
  filter_at(10, all_vars(. > 1))
#Protein_1FDR_1pep <- Protein_list1FDR %>% 
 # filter_at(10, all_vars(. > 1))

#Filter take out HUMAN BOVIN PIG
Protein_5FDR_specie<- filter(Protein_5FDR_1pep, !grepl
                             ('HUMAN|BOVIN|PIG', Accession))
#Protein_1FDR_specie<- filter(Protein_1FDR_1pep, !grepl
#                             ('HUMAN|BOVIN|PIG', Accession))

#Filter take out Reversed Proteins
Protein_5FDR_final_list <- filter(Protein_5FDR_specie, !grepl
                             ('RRRR', Accession))
#Protein_1FDR_final_list <- filter(Protein_1FDR_specie, !grepl
#                                  ('RRRR', Accession))

#Export #CHANGE_NAME_OF_FILE_HERE
write.csv (Protein_5FDR_final_list, "Results/Protein_5FDR_final_Pooideae.csv")
#write.csv (Protein_1FDR_final_list, "Results/Protein_1FDR_final_Pooideae.csv")

####Peptides
#select peptide confidence
conf_number <-  select (Conf_thresh_pep,c=19)
conf0.05 <- conf_number [6,] 

#make confidence numeric
conf5 <- conf0.05 * 100
conf5<- as.numeric(conf5$c)


#filter peptide summary by confidence
Peptide_conf <- Peptide_summary %>% filter_at(7, all_vars(. >= conf5))

#Filter peptides of proteins from other species
Peptide_conf_specie<- filter(Peptide_conf, !grepl
                             ('HUMAN|BOVIN|PIG', Accessions))

#Keep peptide locus
colnames(Peptide_conf_specie )[1] <- 'Peptide_locus' 
Peptide_conf_locus <- filter (Peptide_conf_specie, grepl ('\\.001', Peptide_locus))

#Filter tryptic peptide
Peptide_conf_tryp <- Peptide_conf_locus %>% filter(is.na(Cleavages)) 

#filter modifications
Peptide_conf_na <- Peptide_conf_tryp %>% filter(is.na(Modifications)) 


Peptide_conf_mod <-Peptide_conf_tryp %>% filter(str_detect(Modifications,
                                                           '^((Gln->pyro-Glu@N-term|Carbamidomethyl\\(C\\)@\\d*|Oxidation\\(M\\)@\\d*);?\\s?){0,5}$'))

#Finish export final data

write.csv(Peptide_list5FDR,"Results/Peptide_list5FDR.csv")
