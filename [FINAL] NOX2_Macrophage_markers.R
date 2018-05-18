### Libraries ----
library(ggplot2)
library(reshape2)
library(tidyr)
library(stringr)
library(survival)
library(RTCGAToolbox)
library(survminer)
library(readr)
library(tibble)
library(ggforce)
library(superheat) 
# custom install from github, 
# install.packages("devtools")
# devtools::install_github("rlbarter/superheat")
library(dplyr)

### Set Working Directory PLEASE ADJUST TO YOUR COMPUTER ----
# suggestion: create a new folder called R_Cache on desktop. It will save time
setwd("~/Desktop/R_Cache")
### Gene Lists----
Gene.list <- c(
  #M1 markers
  "CD86",
  "CD80",
  "CD68",
  "HLA-DMA",
  "HLA-DMB",
  "HLA-DOA",
  "HLA-DOB",
  "HLA-DPA1",
  "HLA-DPB1",
  "HLA-DQA1",
  "HLA-DQB1",
  "HLA-DRA",
  "HLA-DRB1",
  "HLA-DRB3",
  "HLA-DRB4",
  "HLA-DRB5",
  "IL1R1",
  "TLR2",
  "TLR4",
  "NOS2",
  "SOCS3",
  "CD163",
  "HLA-DMA", # all the HLA are MHCII isotypes
  "HLA-DMB",
  "HLA-DOA",
  "HLA-DOB",
  "HLA-DPA1",
  "HLA-DPB1",
  "HLA-DQA1",
  "HLA-DQB1",
  "HLA-DRA",
  "HLA-DRB1",
  "HLA-DRB3",
  "HLA-DRB4",
  "HLA-DRB5",
  "MSR1", #CD204, macrophage scanvenger receptor 1
  "MRC1", # mannose receptor 
  "NOX4", "CYBB", "NOX3", "NOX5", "DUOX1", "DUOX2", "CYBA", "NOX1"
)


M1.genes <- c(  "CD86",
                "CD80",
                "CD68",
                "HLA-DMA",
                "HLA-DMB",
                "HLA-DOA",
                "HLA-DOB",
                "HLA-DPA1",
                "HLA-DPB1",
                "HLA-DQA1",
                "HLA-DQB1",
                "HLA-DRA",
                "HLA-DRB1",
                "HLA-DRB3",
                "HLA-DRB4",
                "HLA-DRB5",
                "IL1R1",
                "TLR2",
                "TLR4",
                "NOS2",
                "SOCS3",
                "NOX4", "CYBB", "NOX3", "NOX5", "DUOX1", "DUOX2", "CYBA", "NOX1"
)
M2ab.genes <-  c(
  "CD163",
  "HLA-DMA", # all the HLA are MHCII isotypes
  "HLA-DMB",
  "HLA-DOA",
  "HLA-DOB",
  "HLA-DPA1",
  "HLA-DPB1",
  "HLA-DQA1",
  "HLA-DQB1",
  "HLA-DRA",
  "HLA-DRB1",
  "HLA-DRB3",
  "HLA-DRB4",
  "HLA-DRB5",
  "MSR1", #CD204, macrophage scanvenger receptor 1
  "MRC1", # mannose receptor 
  "NOX4", "CYBB", "NOX3", "NOX5", "DUOX1", "DUOX2", "CYBA", "NOX1"
  
)

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
  edit2 <- edit[edit$Gene.Symbol %in% Gene.list,]
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
  sur1 <- select(sur, Patient.ID, vital_status, time)
  sur1
}

### Pan-Cancer mRNA of NOX4 download ----
BLCA <- dl.RNA.select("BLCA") # bladder car
BRCA <- dl.RNA.select("BRCA") # breast car
CESC <- dl.RNA.select("CESC")
COAD <- dl.RNA.select("COAD") # colon adeno
COADREAD <- dl.RNA.select("COADREAD") # colorectal adenocarcinoma
ESCA <- dl.RNA.select("ESCA")

gc() # garbage collection

GBM  <- dl.RNA.select("GBM")  # gliblastoma multiform
HNSC <- dl.RNA.select("HNSC") # head and neck sq carc
KIRC <- dl.RNA.select("KIRC") # kidney renal clear cell car
KIRP <- dl.RNA.select("KIRP") # kidney renal pap cell car

gc()

LGG  <- dl.RNA.select("LGG")
LIHC <- dl.RNA.select("LIHC") # liver hep car
LUSC <- dl.RNA.select("LUSC") # lung squ cell car
LUAD <- dl.RNA.select("LUAD") # lung adeno
MESO <- dl.RNA.select("MESO")
OV   <- dl.RNA.select("OV")   # ovarian ser cyst
PRAD <- dl.RNA.select("PRAD") # prostate adeno
PAAD <- dl.RNA.select("PAAD")

gc() 

PCPG <- dl.RNA.select("PCPG")
READ <- dl.RNA.select("READ")
SARC <- dl.RNA.select("SARC") # sarcoma
#SKCM <- dl.RNA.select("SKCM")
STAD <- dl.RNA.select("STAD") # stomach adeno
THCA <- dl.RNA.select("THCA") # thyroid carcinoma
gc()
THCA <- dl.RNA.select("THCA")
UCEC <- dl.RNA.select("UCEC") # uterine corpus endometrial carcinoma
UCS  <- dl.RNA.select("UCS")

gc()

#rbind for PanCan
Pan.RNA  <- rbind( BLCA, # bladder car
                   BRCA, # breast car
                   CESC,
                   COAD, # colon adeno
                   COADREAD, # colorectal adenocarcinoma
                   ESCA,
                   GBM,  # gliblastoma multiform
                   HNSC, # head and neck sq carc
                   KIRC, # kidney renal clear cell car
                   KIRP, # kidney renal pap cell car
                   LGG,
                   LIHC, # liver hep car
                   LUSC, # lung squ cell car
                   LUAD, # lung adeno
                   MESO,
                   OV,   # ovarian ser cyst
                   PRAD, # prostate adeno
                   PAAD,
                   PCPG,
                   READ,
                   SARC, # sarcoma
                   # SKCM,
                   STAD, # stomach adeno
                   THCA, # thyroid carcinoma
                   THCA,
                   UCEC, # uterine corpus endometrial carcinoma
                   UCS
)

Pan.RNA$Patient.ID <- str_sub(Pan.RNA$Patient.ID, start = 1, end = 15) # 13 and beyond are samples id
Pan.RNA$Gene.Symbol <- as.character(Pan.RNA$Gene.Symbol)
### Read p53 Mutation Files from cBioPortal ----
P53_exp  <-
  read.delim(file = "TP53_Expression_All.txt",
             header = T,
             na.strings = "NA")   # p53 with Expression Data
P53_exp <-
  P53_exp %>% dplyr::select(P53.Mutation = Mutation, Cancer.Study, Sample.Id)
P53_exp <-
  na.omit(P53_exp) # important to remove missing values now to avoid mixing with 'no mutation'

P53_mut <-
  read.delim("mutation_table_TP53_02022017.tsv",
             header = T,
             na.strings = "NA")
P53_mut <-
  P53_mut %>% dplyr::select(P53.Mutation = AA.change, Cancer.Study, Sample.Id = Sample.ID)
P53_mut <- na.omit(P53_mut)
### Data Clean Up p53 
#'Not Sequenced' data is now designated as NA
P53_mut$P53.Mutation[(P53_mut$P53.Mutation == "Not Sequenced")]   <-
  NA
P53_mut$P53.Mutation[(P53_mut$P53.Mutation == "[Not Available]")] <-
  NA

P53_exp$P53.Mutation[(P53_exp$P53.Mutation == "Not Sequenced")]   <-
  NA
P53_exp$P53.Mutation[(P53_exp$P53.Mutation == "[Not Available]")] <-
  NA

# Merging two p53 Data Files to get more mutant reads on all ID's
P53_com <-
  full_join(P53_mut,
            P53_exp,
            by = "Sample.Id",
            suffix = c(".mut_file" , ".exp_file"))
P53_com <- dplyr::select(P53_com, Sample.Id, everything())
# mut: 6623 NO sequence
# exp: 8075 NO sequence

# Now, if from the exp file it is NA, copy p53 mutation status from mut file
# first convert mutation from factors to characters to allow replacement operation
P53_com$P53.Mutation.mut_file <-
  as.character(P53_com$P53.Mutation.mut_file)
P53_com$P53.Mutation.exp_file <-
  as.character(P53_com$P53.Mutation.exp_file)
# replace missing mutation status
P53_com$P53.Mutation.exp_file[is.na(P53_com$P53.Mutation.exp_file)] <-
  P53_com$P53.Mutation.mut_file[is.na(P53_com$P53.Mutation.exp_file)]
# "check where exp is NA, and replace with mut where there is NA'
# mut: 2006 NO sequence
# exp: 6623 NO sequence, leaving 10013 obs at this level

# renaming NA as WT, assuming no mutation is equal to wild-type
# because we know that NA in the exp are WT, we can assume those NA in the mut are also WT
P53_com$P53.Mutation.mut_file[is.na(P53_com$P53.Mutation.mut_file)] <-
  "WT"
P53_com$P53.Mutation.exp_file[is.na(P53_com$P53.Mutation.exp_file)] <-
  "WT"

# removing non-matched mutation status rows via subsetting
P53_com.sub <-
  subset(P53_com, P53.Mutation.exp_file == P53.Mutation.mut_file)
P53.final <-
  dplyr::select(P53_com.sub,
                Patient.ID = Sample.Id,
                P53.Mutation = P53.Mutation.exp_file,
                Case.Study = Cancer.Study.mut_file) # reordering columns and mut mutation column

gc()
### NOX4-4;DUOX1-2 M1 and M2 macrophage corrl ----

# this function will extract NOX4 and a gene (fa) from the df = Pan.final and create a matrix 
# of correlation grouped by p53 mutation status
# and also create a column of
calc.rho.unclustered.no_p53 <-
  function(list, wid = 0.5, hei = 0.5, dpi = 100, low.rho = -1, high.rho = 1) {
    df <- Pan.final[Pan.final$Gene.Symbol %in% list,] 
    
    # use this  function to calculate the spearman rank corrl in a matrix
    temp <- select(df, # this object is created later after filtering for selected p53 mutations
                   -Case.Study, -P53.Mutation)
    temp.unique <- unique(temp, nmax = 2 ) # removes duplicates
    temp.wide <- spread(data = temp.unique, key = Gene.Symbol, value = mRNA.Value)
    
    rownames(temp.wide) <- temp.wide$Patient.ID
    temp.wide.IDRow <- select(temp.wide, -Patient.ID) # remove IDrow
    cor.table <- cor(temp.wide.IDRow, method = "spearman", use = "pairwise.complete.obs")
    
    noxes.unsorted <- c("NOX4", "NOX1", "CYBB", "NOX3", "NOX5", "DUOX1", "DUOX2", "CYBA") # unsorted, ordered by alphabet
    noxes.vector <- c("NOX4", "NOX1", "CYBB", "CYBA", "NOX3", "NOX5", "DUOX1", "DUOX2") # specified order
    noxes <- factor(noxes.unsorted, levels = noxes.vector) # ordered in my specific
    # just want NOX rows
    cor.select.1 <- cor.table[rownames(cor.table) %in% noxes,]
    # rm NOX from column (rho = 1's)
    cor.select.2 <- cor.select.1[, !colnames(cor.select.1) %in% noxes]
    # need to remove NOX'es from heatmap
    
    # now create a heatmap
    png(paste(deparse(substitute(list)), # deparse(substitute()) calls the name of the df as character string
              "_nogrouping.png"), width = wid, height = hei, units = 'in', res = dpi)
    
    superheat(
      cor.select.2,
      bottom.label.text.angle = 90,
      bottom.label.text.size = 3,
      scale = F,
      grid.hline.col = "white",
      grid.vline.col = "white",
      grid.hline.size = 1,
      grid.vline.size = 1,
      pretty.order.rows = TRUE,
      pretty.order.cols = TRUE,
      left.label.text.alignment = "right",
      bottom.label.text.alignment = "right",
      heat.pal = c("turquoise", "white", "violetred1"),
      heat.lim = c(-1,1))
    
    dev.off()
    
    png(paste(deparse(substitute(list)), # deparse(substitute()) calls the name of the df as character string
              "_nogrouping_withnumber.png"), width = wid, height = hei, units = 'in', res = dpi)
    
    superheat(
      cor.select.2,
      bottom.label.text.angle = 90,
      bottom.label.text.size = 3,
      scale = F,
      grid.hline.col = "white",
      grid.vline.col = "white",
      grid.hline.size = 1,
      grid.vline.size = 1,
      pretty.order.rows = TRUE,
      pretty.order.cols = TRUE,
      left.label.text.alignment = "right",
      bottom.label.text.alignment = "right",
      X.text = round(as.matrix(cor.select.2), 1),
      heat.pal = c("turquoise", "white", "violetred1"),
      heat.lim = c(-1,1))
    
    dev.off()
  }    
calc.rho.unclustered.no_p53.nos2 <-
  function(list, wid = 0.5, hei = 0.5, dpi = 100, low.rho = -1, high.rho = 1) {
    df <- Pan.final[Pan.final$Gene.Symbol %in% list,] 
    
    # use this  function to calculate the spearman rank corrl in a matrix
    temp <- select(df, # this object is created later after filtering for selected p53 mutations
                   -Case.Study, -P53.Mutation)
    temp.unique <- unique(temp, nmax = 2 ) # removes duplicates
    temp.wide <- spread(data = temp.unique, key = Gene.Symbol, value = mRNA.Value)
    
    rownames(temp.wide) <- temp.wide$Patient.ID
    temp.wide.IDRow <- select(temp.wide, -Patient.ID) # remove IDrow
    cor.table <- cor(temp.wide.IDRow, method = "spearman", use = "pairwise.complete.obs")
    
    noxes.unsorted <- c("NOX4", "NOX1", "CYBB", "NOX3", "NOX5", "DUOX1", "DUOX2", "CYBA", "NOS2") # unsorted, ordered by alphabet
    noxes.vector <- c("NOX4", "NOX1", "CYBB", "CYBA", "NOX3", "NOX5", "DUOX1", "DUOX2", "NOS2") # specified order
    noxes <- factor(noxes.unsorted, levels = noxes.vector) # ordered in my specific
    # just want NOX rows
    cor.select.1 <- cor.table[rownames(cor.table) %in% noxes,]
    # rm NOX from column (rho = 1's)
    cor.select.2 <- cor.select.1[, !colnames(cor.select.1) %in% noxes]
    # need to remove NOX'es from heatmap
    
    # now create a heatmap
    png(paste(deparse(substitute(list)), # deparse(substitute()) calls the name of the df as character string
              "_nogrouping.png"), width = wid, height = hei, units = 'in', res = dpi)
    
    superheat(
      cor.select.2,
      bottom.label.text.angle = 90,
      bottom.label.text.size = 3,
      scale = F,
      grid.hline.col = "white",
      grid.vline.col = "white",
      grid.hline.size = 1,
      grid.vline.size = 1,
      pretty.order.rows = TRUE,
      pretty.order.cols = TRUE,
      left.label.text.alignment = "right",
      bottom.label.text.alignment = "right",
      heat.pal = c("turquoise", "white", "violetred1"),
      heat.lim = c(-1,1))
    
    dev.off()
    
    png(paste(deparse(substitute(list)), # deparse(substitute()) calls the name of the df as character string
              "_nogrouping_withnumber.png"), width = wid, height = hei, units = 'in', res = dpi)
    
    superheat(
      cor.select.2,
      bottom.label.text.angle = 90,
      bottom.label.text.size = 3,
      scale = F,
      grid.hline.col = "white",
      grid.vline.col = "white",
      grid.hline.size = 1,
      grid.vline.size = 1,
      pretty.order.rows = TRUE,
      pretty.order.cols = TRUE,
      left.label.text.alignment = "right",
      bottom.label.text.alignment = "right",
      X.text = round(as.matrix(cor.select.2), 1),
      heat.pal = c("turquoise", "white", "violetred1"),
      heat.lim = c(-1,1))
    
    dev.off()
  }    


# frequency of p53 mutation in selected patients
# P53_com.sub will then be merged with pan-cancer set to obtain p53 mutation for each patient
RNA.joined.p53 <- inner_join(P53.final, Pan.RNA, by = "Patient.ID")

# Now, Name the mutations of interest 
mut.interest <- c(
  "R175H",
  "R248Q",
  "R273H",
  "R273C",
  "R248W",
  "Y220C",
  "R249S",
  "G245D",
  "R273C",
  "R248Q",
  "H179R",
  "R282W",
  "V157F",
  "H193R",
  "R158L",
  "R273L",
  # "G248S", # very few of these
  "R158L",
  # "C175F", # very few of these
  "H179R",
  "G245S",
  "WT"
)

Pan.final <- RNA.joined.p53[RNA.joined.p53$P53.Mutation %in% mut.interest,]
head(Pan.final)

calc.rho.unclustered.no_p53(list= M1.genes, wid = 7, hei = 5, dpi = 600, low.rho = -1, high.rho = 1)
calc.rho.unclustered.no_p53(list= M2ab.genes, wid = 7, hei = 5, dpi = 600, low.rho = -1, high.rho = 1)
calc.rho.unclustered.no_p53(list= Gene.list, wid = 7, hei = 5, dpi = 600, low.rho = -1, high.rho = 1)
calc.rho.unclustered.no_p53.nos2(list= Gene.list.with.NOS2, wid = 7, hei = 5, dpi = 600, low.rho = -1, high.rho = 1)
