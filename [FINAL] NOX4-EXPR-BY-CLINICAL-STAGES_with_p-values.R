### Libraries ----
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(FSA)
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

### Enhanced Function to plot and generate adj.p-value frames by Dunn.Test ----
# assign cancer value

NOX4.dunn.test <-
  function(cancer, pos = 1, ymax = 10) {
    clin <- dl.surv(cancer)
    Nox4 <- dl.RNA.select(cancer)
    
    # trim patient ID and inner join
    Nox4$Patient.ID <- str_sub(Nox4$Patient.ID, start = 1, end = 12)
    Nox4 <- Nox4[!duplicated(Nox4$Patient.ID), ]
    Nox4 <- Nox4 %>%
      select(Patient.ID,
             NOX4 = mRNA.Value)
    
    # remove rows with NA as cancer stage
    FINAL <- na.omit(inner_join(Nox4, clin, by = "Patient.ID"))
    FINAL$pathologic_stage <- toupper(FINAL$pathologic_stage)
    FINAL$pathologic_stage <-
      str_sub(FINAL$pathologic_stage, start = 7, end = 30)
    leng <- FINAL %>%
      group_by(pathologic_stage) %>%
      summarise(n = length(Patient.ID))
    
    plot <- ggplot(FINAL) + aes(x = pathologic_stage, y = log(NOX4)) +
      geom_boxplot(outlier.shape = NA, lwd = 1.5) +
      geom_sina(color = "royalblue", alpha = 0.6) +
      ylab("Relative Fold Change in NOX4 mRNA") +
      xlab(NULL) +
      ggtitle(cancer) +
      expand_limits(y= ymax) +
      geom_text(
        size = 6,
        data = leng,
        aes(label = paste("n", n, sep = ":")),
        y = pos,
        # specify location
        label.size = 0 # specify no boarder
      ) +
      theme(
        axis.text.x = element_text(
          angle = 0,
          hjust = 0.5,
          size = 20,
          vjust = 0.5
        ),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 23),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_rect(
          fill = "white",
          color = "black",
          size = 4
        ),
        axis.ticks = element_line(size = 2)
      )
    
    print(plot)
    ggsave(plot, filename = cancer, device = "pdf", dpi=900, width = 9, height = 6.5)
    ##Kruskal-Wallis test is a non-parametric ANOVA, basically an extended mann-whitney for >2 groups
    # first, use the kruskal test to see if any group is statistically different from another
    krus.test <- kruskal.test(x = FINAL$NOX4, g =  as.factor(FINAL$pathologic_stage))
    # p = 1.87 e-06, so now we decide which one is statistically different using the Dunn Test
    dun.results <-
      dunnTest(NOX4 ~ as.factor(pathologic_stage),
               data = FINAL,
               method = "bh") # Benjamini-Hochberg adjustment 
    print(dun.results) # these tests are special classes and difficult to export as table/excel
    print(krus.test)
  }

NOX4.dunn.test.ignore.substages <-
  function(cancer, pos = 1, ymax = 10) {
    clin <- dl.surv(cancer)
    Nox4 <- dl.RNA.select(cancer)
    
    
    # trim patient ID and inner join
    Nox4$Patient.ID <- str_sub(Nox4$Patient.ID, start = 1, end = 12)
    Nox4 <- Nox4[!duplicated(Nox4$Patient.ID), ]
    Nox4 <- Nox4 %>%
      select(Patient.ID,
             NOX4 = mRNA.Value)     
    
    # remove rows with NA as cancer stage
    FINAL <- na.omit(inner_join(Nox4, clin, by = "Patient.ID"))
    FINAL$pathologic_stage <- toupper(FINAL$pathologic_stage) # upper case
    FINAL$pathologic_stage <-
      str_sub(FINAL$pathologic_stage, start = 7, end = 30) # removes "STAGES"
    
    # remove substages
    FINAL$pathologic_stage <- 
      gsub("A", "", FINAL$pathologic_stage) # remove all A
    FINAL$pathologic_stage <- 
      gsub("B", "", FINAL$pathologic_stage) # remove all B
    FINAL$pathologic_stage <- 
      gsub("C", "", FINAL$pathologic_stage) # remove all C
    
    leng <- FINAL %>%
      group_by(pathologic_stage) %>%
      summarise(n = length(Patient.ID))
    
    plot <- ggplot(FINAL) + aes(x = pathologic_stage, y = log(NOX4)) +
      geom_boxplot(outlier.shape = NA, lwd = 1.5) +
      geom_sina(color = "royalblue", alpha = 0.6) +
      ylab("Relative Fold Change in NOX4 mRNA") +
      xlab(NULL) +
      ggtitle(cancer) +
      expand_limits(y= ymax) +
      geom_text(
        size = 6,
        data = leng,
        aes(label = paste("n", n, sep = ":")),
        y = pos,
        # specify location
        label.size = 0 # specify no boarder
      ) +
      theme(
        axis.text.x = element_text(
          angle = 0,
          hjust = 0.5,
          size = 20,
          vjust = 0.5
        ),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 23),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_rect(
          fill = "white",
          color = "black",
          size = 4
        ),
        axis.ticks = element_line(size = 2)     )
    
    print(plot)
    ggsave(plot, filename = paste(cancer, "_no_substage"), device = "pdf", dpi=900, width = 9, height = 6.5)
    ##Kruskal-Wallis test is a non-parametric ANOVA, basically an extended mann-whitney for >2 groups
    # first, use the kruskal test to see if any group is statistically different from another
    krus.test <- kruskal.test(x = FINAL$NOX4, g =  as.factor(FINAL$pathologic_stage))
    # p = 1.87 e-06, so now we decide which one is statistically different using the Dunn Test
    dun.results <-
      dunnTest(NOX4 ~ as.factor(pathologic_stage),
               data = FINAL,
               method = "bh") # Benjamini-Hochberg adjustment 
    print(dun.results) # these tests are special classes and difficult to export as table/excel
    print(krus.test)

  }

# these are the only ones with adequate stage information
# run each line at a time, and interpret the output;
  # kruskal-wallis need to have a p value <0.05 then you can look at the dun results to see which pair is different
  # KW only tells you whether or not there is a likely difference among groups (one p-values for the entire set)
  # you then follow with Dunn test to compare each group to each other and get adjusted p-values.

# significant sets
NOX4.dunn.test("BLCA", 0, ymax = 8) # signifiant
NOX4.dunn.test("THCA", 1.5, ymax = 8) # thyroid carcinoma ### significant
NOX4.dunn.test.ignore.substages("ESCA", 1, 6.5)

# barely under
NOX4.dunn.test.ignore.substages("MESO", 1) # mesothelioma

#non-significant sets
NOX4.dunn.test("READ", 1) # rectum, uses substages
NOX4.dunn.test.ignore.substages("READ", 1) # rectum, uses substages
NOX4.dunn.test("KIRC", 2.8) # kidney renal clear cell car
NOX4.dunn.test("LIHC", 1) # liver hep car
NOX4.dunn.test.ignore.substages("LIHC", 1) # liver hep car
NOX4.dunn.test("KIRP", 1) # kidney renal pap cell car
NOX4.dunn.test("PAAD", 1) # pancreatic
NOX4.dunn.test.ignore.substages("PAAD", 1) # pancreatic
NOX4.dunn.test("LUSC", 2) # lung squ cell car
NOX4.dunn.test.ignore.substages("LUSC", 2) # lung squ cell car
NOX4.dunn.test("LUAD", 1) # lung adeno
NOX4.dunn.test.ignore.substages("LUAD", 1) # lung adeno
NOX4.dunn.test("STAD", 0) # stomach adeno
NOX4.dunn.test.ignore.substages("STAD", 0) # stomach adeno
NOX4.dunn.test("COAD", 1, 6) # colon adeno
NOX4.dunn.test.ignore.substages("COAD", 1) # colon adeno
NOX4.dunn.test("HNSC", 1) # head and neck sq carc
NOX4.dunn.test("ESCA", 1, 6.5 ) # esophogeal

NOX4.dunn.test.ignore.substages("BRCA", 1, 6.5) # breast car, uses sub stages, try without substsages
