### Libraries ----
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(survival)
library(RTCGAToolbox)
library(survminer)
library(readr)
library(tibble)
library(ggforce)
library(SOfun)
### Set Working Directory PLEASE ADJUST TO YOUR COMPUTER ----
# suggestion: create a new folder called R_Cache on desktop. It will save time
setwd("~/Desktop/R_Cache")
### Functions to download data ----
# Download Complete Sets of Data
dl.full <- function(x) {
  getFirehoseData(
    dataset = x,
    runDate = "20160128",
    gistic2_Date = "20160128",
    # gistic2 is CNA
    RNAseq_Gene = TRUE,
    # gene level expression data from RNA-seq [raw values]
    Clinic = TRUE,
    # clinical information of patient sample
    mRNA_Array = TRUE,
    # gene level expression data by array platform
    Mutation = TRUE,
    # Gene level mutation information matrix
    CNA_SNP = T, # copy number alterations in somatic cells provided by segmented sequecing
    CNA_Seq = T, # copy number alterations provided by NGS sequences
    CNA_CGH = T, # copy number alternations provided by CGH platform
    CNV_SNP = T, # copy number alterantion in germline cells
    # Methylation = T, # methylation provided by array platform
    # RPPA = T, # reverse phase protein array expression
    RNAseq2_Gene_Norm = TRUE,
    # normalized count
    fileSizeLimit = 99999,
    # getUUIDs = T,
    destdir = "FireHose Data",
    forceDownload = F
  )
}
# this S4 object can be extracted for dataframe containing information that can be plotted

# Gather mRNA in long form
dl.RNA.select <- function(y) {
  set <- getFirehoseData(
    dataset = y,
    runDate = "20160128",
    gistic2_Date = "20160128",
    # gistic2 is CNA``
    RNAseq_Gene = TRUE,
    # gene level expression data from RNA-seq [raw values]
    Clinic = F,
    # clinical information of patient sample
    mRNA_Array = F,
    # gene level expression data by array platform
    Mutation = TRUE,
    # Gene level mutation information matrix
    # CNA_SNP = T, # copy number alterations in somatic cells provided by segmented sequecing
    # CNA_Seq = T, # copy number alterations provided by NGS sequences
    # CNA_CGH = T, # copy number alternations provided by CGH platform
    # CNV_SNP = T, # copy number alterantion in germline cells
    # Methylation = T, # methylation provided by array platform
    # RPPA = T, # reverse phase protein array expression
    RNAseq2_Gene_Norm = TRUE,
    # normalized count
    fileSizeLimit = 99999,
    # getUUIDs = T,
    destdir = "FireHose Data",
    forceDownload = F
  )
  edit <- gather(rownames_to_column(as.data.frame(getData(set, type = "RNASeq2GeneNorm")), var = "Gene.Symbol"), 
                 key = "Patient.ID", value = "mRNA.Value", -Gene.Symbol) # extract RNAseq, transform row name as a colmn, then transform into long form
  edit$Patient.ID <- str_sub(edit$Patient.ID, start = 1, end = 15) # we remove the trailing barcodes so we can cross match the data later to patients
  # now remove RNA-expressions of genes we dont need
  edit2 <- edit[edit$Gene.Symbol %in% "NOX4",]
  rm(edit)
  rm(set)
  edit2
}

# Gather Clinical Data in long form
# we will use days to death for the survival
dl.surv <- function(z) {
  full.set <- getFirehoseData(
    dataset = z,
    runDate = "20160128",
    #gistic2_Date = "20160128",
    # gistic2 is CNA
    #RNAseq_Gene = TRUE,
    # gene level expression data from RNA-seq [raw values]
    Clinic = TRUE,
    # clinical information of patient sample
    # mRNA_Array = TRUE,
    # gene level expression data by array platform
    # Mutation = TRUE,
    # Gene level mutation information matrix
    # CNA_SNP = T, # copy number alterations in somatic cells provided by segmented sequecing
    # CNA_Seq = T, # copy number alterations provided by NGS sequences
    # CNA_CGH = T, # copy number alternations provided by CGH platform
    # CNV_SNP = T, # copy number alterantion in germline cells
    # Methylation = T, # methylation provided by array platform
    # RPPA = T, # reverse phase protein array expression
    # RNAseq2_Gene_Norm = TRUE,
    # normalized count
    fileSizeLimit = 99999,
    # getUUIDs = T,
    destdir = "FireHose Data",
    forceDownload = F
  )
  clin <- getData(full.set, type = "Clinical")
  # to obtain overall survival (OS) we combine three columns: vital status and days to last fol up OR days to death
  # if patient is still living (vital status = 0) we use days to last follow up
  # if patient is dead (vital status = 1) we use days to death
  clin1 <- rownames_to_column(as.data.frame(clin), var = "Patient.ID")
  clin1$Patient.ID <-
    str_replace_all(clin1$Patient.ID, "[.]", "-") # change divider to dash
  clin1$Patient.ID <- toupper(clin1$Patient.ID) # change to upper case
  
  # here we create a new column call time (time to event). 
  # since the days to last follow up and time to death given
  # are mutually exclusive, we can merge them together to get one column
  sur <- clin1 %>% 
    mutate(time = days_to_death,
           time = as.numeric(time),
           days_to_last_followup = as.numeric(days_to_last_followup)) %>% 
    select(Patient.ID, vital_status, time, everything()) # now time to event is only time to death
  # because they are mutally exclusive, we replace where NA in time to death with last follow up
  # to get full time to event
  sur$time[is.na(sur$time)] <- sur$days_to_last_followup[is.na(sur$time)]
  sur1 <- select(sur, Patient.ID, vital_status, time, pathologic_stage)
  sur1
}
### General Function to plot NOX4 exp by cancer stage in a set cancer ----

stageNOX4 <- function(cancer, pos = 0.5) {
  clin <- dl.surv(cancer)
  Nox4 <- dl.RNA.select(cancer) 
  
  # trim patient ID and inner join
  Nox4$Patient.ID <- str_sub(Nox4$Patient.ID, start = 1, end = 12)
  Nox4 <- Nox4[!duplicated(Nox4$Patient.ID),]
  Nox4 <- Nox4 %>% 
    select(
      Patient.ID, 
      NOX4 = mRNA.Value
    )
  
  # remove rows with NA as cancer stage
  FINAL <- na.omit(inner_join(Nox4, clin, by = "Patient.ID"))
  FINAL$pathologic_stage <- toupper(FINAL$pathologic_stage)
  FINAL$pathologic_stage <- str_sub(FINAL$pathologic_stage, start = 7, end = 30)
  leng <- FINAL %>% 
    group_by(pathologic_stage) %>% 
    summarise(n = length(Patient.ID))
  
  plot <- ggplot(FINAL) + aes( x = pathologic_stage, y = log(NOX4)) +
    geom_boxplot(outlier.shape = NA, lwd=1.5) +
    geom_sina(color = "royalblue", alpha = 0.8) +
    ylab("Relative Fold Change in NOX4 mRNA") +
    xlab(NULL) +
    ggtitle(cancer) +
    geom_label(size = 6, data = leng, aes(label = paste("n", n, sep = ":")), 
               y = pos, # specify location 
               label.size = 0 # specify no boarder
               ) + 
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 20, vjust=0.5),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 23),
          panel.grid.minor.y=element_blank(),
          panel.grid.minor.x=element_blank(),
          panel.grid.major.x=element_blank(),
          panel.grid.major.y=element_blank(),
          panel.background = element_rect(fill = "white", color = "black", size = 4),
          axis.ticks = element_line(size=2))
  plot
  ggsave(plot, filename = cancer, device = "pdf", dpi=900, width = 9, height = 6.5)
  
}

### Plotting the cancers by loop ----
# the follonwg code writes the pdf plots to the wd
stageNOX4("BLCA", 0.2) # bladder car
stageNOX4("BRCA") # breast car, uses sub stages, try without substsages
stageNOX4("COAD", 1) # colon adeno
stageNOX4("ESCA") # esophogeal
stageNOX4("HNSC") # head and neck sq carc
stageNOX4("KIRC", 2.8) # kidney renal clear cell car
stageNOX4("KIRP") # kidney renal pap cell car
stageNOX4("LIHC") # liver hep car
stageNOX4("LUSC", 2) # lung squ cell car
stageNOX4("LUAD") # lung adeno
stageNOX4("MESO") # mesothelioma
stageNOX4("PAAD") # pancreatic
stageNOX4("READ") # rectum, uses substages
stageNOX4("STAD", 0) # stomach adeno
stageNOX4("THCA", 1.5) # thyroid carcinoma


# # selecting the studies that have stage informations ----
# studies <- c("BLCA", # bladder car
#              "BRCA", # breast car, uses sub stages, try without substsages
#              #"CESC" # cervical # cervical only uses TNM stage markers
#              "COAD", # colon adeno
#              "ESCA", # esophogeal
#              #"GBM",  # gliblastoma multiform, uses performance score
#              "HNSC", # head and neck sq carc
#              "KIRC", # kidney renal clear cell car
#              "KIRP", # kidney renal pap cell car
#              #"LGG",# brain lower grade glioma, uses performance score q 
#              "LIHC", # liver hep car
#              "LUSC", # lung squ cell car
#              "LUAD", # lung adeno
#              "MESO",# mesothelioma
#              #"OV",   # ovarian ser cyst no stage data availiable
#              #"PRAD", # prostate adeno, TNM only
#              "PAAD", # pancreatic
#              "READ", # rectum, uses substages
#              #"SARC", # sarcoma, no tumor stage
#              "STAD", # stomach adeno
#              "THCA" # thyroid carcinoma
#              #"UCEC", # uterine corpus endometrial carcinoma, no data on tumor stage
#              #"UCS" #no tumor data
# )