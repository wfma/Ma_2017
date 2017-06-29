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
  edit2 <- edit[edit$Gene.Symbol %in% EMT.Genes,]
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

### Chooosing Interested EMT genes and NADPH ox ----
EMT.Genes <- c(
  "FN1",
  "VIM",
  "CDH2",
  "ZEB1",
  "ZEB2",
  "MMP9",
  "MMP2",
  "TWIST1",
  "SNAI1",
  "SNAI2",
  "ACTA2",
  "CLDN1",
  "DSP",
  "PKP1",
  "CDH1",
  "NOX4",
  "CYBB", # NOX2
  "NOX3",
  "NOX5",
  "CYBA", # p22phox
  "NOX1",
  "DUOX1",
  "DUOX2",
  ### recently added list
  "VTN", # vitronectin, valid
  "IL8", # IL8
  "COL1A1", # collagen type 1 alpha 1 chain, valid
  "ITGA5", # alpha 5 subunit of integrin, valid
  "AGT", # angiotensin, IPA
  "AKT1", # alcohol binding phosphotransferase
  "AKT2",
  "TGFB1",
  "TGFB2",
  "TGFB3",
  "MMP3",
  "WNT5A",
  "WNT5B",
  "COL1A2", # collagen type I alpha 2 chain
  "COL3A1",
  "COL5A2",
  "ILK", # integrin linked kinase
  "CXCR4",
  "CXCL12", # chemokine, also known as SDF1
  "ITGAV",
  "TGFBR1",
  "TJP1", # zonuli1
  "CRB1", # cell polarity complex crumbs1, valid
  "KRT18",
  "TLR4"# keratin 18, valid 
  
)


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
### NOX4-EMT correlations pan cancer with p53 ----

# this function will extract NOX4 and a gene (fa) from the df = Pan.final and create a matrix 
# of correlation grouped by p53 mutation status
# and also create a column of
calc.rho <- function(fa) { # use this function to calculate the rho between two genes in a long.form data frame
  gene  <- Pan.final[Pan.final$Gene.Symbol %in% fa,] # use genes in ""
  NOX4 <- Pan.final[Pan.final$Gene.Symbol %in% "NOX4",]
  co <- inner_join(x= gene, y= NOX4, by = "Patient.ID")
  co2 <- co %>% 
    select(
      Patient.ID,
      P53.Mut = P53.Mutation.x,
      gene = mRNA.Value.x,
      NOX4 = mRNA.Value.y,
      -P53.Mutation.y,
      -Gene.Symbol.x,
      -Gene.Symbol.y,
      -Case.Study.x,
      -Case.Study.y
    )
  com <- co2[!duplicated(co2),]
  corr <- com %>% 
    group_by(P53.Mut) %>% 
    summarise(
      rho = cor(NOX4, gene, method = "spearman", use = "complete.obs"))
  corr <- mutate(corr, 
                 P53.Mut = as.factor(P53.Mut))
  corr <- mutate(corr, rho = round(rho, 2))
  data.frame(corr)
}
calc.rho.n <- function(fa) { # use this function to calculate the rho between two genes in a long.form data frame
  gene  <- Pan.final[Pan.final$Gene.Symbol %in% fa,] # use genes in ""
  NOX4 <- Pan.final[Pan.final$Gene.Symbol %in% "NOX4",]
  co <- inner_join(x= gene, y= NOX4, by = "Patient.ID")
  co2 <- co %>% 
    select(
      Patient.ID,
      P53.Mut = P53.Mutation.x,
      gene = mRNA.Value.x,
      NOX4 = mRNA.Value.y,
      -P53.Mutation.y,
      -Gene.Symbol.x,
      -Gene.Symbol.y,
      -Case.Study.x,
      -Case.Study.y
    )
  com <- co2[!duplicated(co2),]
  corr <- com %>% 
    group_by(P53.Mut) %>% 
    summarise(
      rho = cor(NOX4, gene, method = "spearman", use = "complete.obs"),
      n = length(gene))
  corr <- mutate(corr, 
                 P53.Mut = as.factor(P53.Mut))
  corr <- mutate(corr, rho = round(rho, 2))
  data.frame(corr)
}

# frequency of p53 mutation in selected patients
# P53_com.sub will then be merged with pan-cancer set to obtain p53 mutation for each patient
Pan.select.p53 <- inner_join(P53.final, Pan.RNA, by = "Patient.ID")

p53.freq <- Pan.select.p53 %>% 
  group_by(P53.Mutation) %>% 
  summarise(n = length(P53.Mutation))


# Now, Name the mutations of interest 
mut.interest <- c(
  "R175H",
  "R248Q",
  "R273H",
  # "R280K",
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
  "G248S",
  "R158L",
  "C175F",
  "H179R",
  "G245S",
  # "R175G",
  "WT"
)

Pan.final <- Pan.select.p53[Pan.select.p53$P53.Mutation %in% mut.interest,]
# this is a good place to clean up the environment
#rm(list=setdiff(ls(), c("Pan.final", "calc.rho", "EMT.Genes", "Pan.RNA")))

# invidually calculating and creating object (learn how to loop better)
FN1.NOX4.n <- calc.rho.n("FN1")
FN1 <- calc.rho("FN1")
VIM <- calc.rho("VIM")
CDH2 <- calc.rho("CDH2")
ZEB1 <- calc.rho("ZEB1")
ZEB2 <- calc.rho("ZEB2")
MMP9 <- calc.rho("MMP9")
MMP2 <- calc.rho("MMP2")
TWIST1 <- calc.rho("TWIST1")
SNAI1 <- calc.rho("SNAI1")
SNAI2 <- calc.rho("SNAI2")
ACTA2 <- calc.rho("ACTA2")
CLDN1 <- calc.rho("CLDN1")
TJP1 <- calc.rho("TJP1")
CRB1 <- calc.rho("CRB1")
KRT18 <- calc.rho("KRT18")
DSP <- calc.rho("DSP")
PKP1 <- calc.rho("PKP1")
CDH1 <- calc.rho("CDH1")
NOX4 <- calc.rho("NOX4")
### recently added list

COL1A1 <- calc.rho("COL1A1") # collagen type 1 alpha 1 chain, valid
ITGA5 <- calc.rho("ITGA5") # alpha 5 subunit of integrin, valid
TGFB1 <- calc.rho("TGFB1")
TGFB2 <- calc.rho("TGFB2")
TGFB3 <- calc.rho("TGFB3")
MMP3 <- calc.rho("MMP3")
WNT5A <- calc.rho("WNT5A")
WNT5B <- calc.rho("WNT5B")
COL1A2 <- calc.rho("COL1A2") # collagen type I alpha 2 chain
COL3A1 <- calc.rho("COL3A1")
COL5A2 <- calc.rho("COL5A2")
ILK <- calc.rho("ILK") # integrin linked kinase
CXCR4 <- calc.rho("CXCR4")
CXCL12 <- calc.rho("CXCL12") # chemokine, also known as SDF1
ITGAV <- calc.rho("ITGAV")
TGFBR1 <- calc.rho("TGFBR1")
TJP1 <- calc.rho("TJP1") # zonuli1
CRB1 <- calc.rho("CRB1") # cell polarity complex crumbs1, valid
KRT18 <- calc.rho("KRT18") # keratin 18, valid 

# rename the rho in each dataframe to its respective gene_nox4
FN1 <- rename(FN1,FN1 = rho)
VIM <- rename(VIM,VIM = rho)
CDH2 <- rename(CDH2,CDH2 = rho)
ZEB1 <- rename(ZEB1,ZEB1 = rho)
ZEB2 <- rename(ZEB2,ZEB2 = rho)
MMP9 <- rename(MMP9,MMP9 = rho)
MMP2 <- rename(MMP2,MMP2 = rho)
TWIST1 <- rename(TWIST1,TWIST1 = rho)
SNAI1 <- rename(SNAI1,SNAI1 = rho)
SNAI2 <- rename(SNAI2,SNAI2 = rho)
ACTA2 <- rename(ACTA2,ACTA2 = rho)
CLDN1 <- rename(CLDN1,CLDN1 = rho)
TJP1 <- rename(TJP1,TJP1 = rho)
CRB1 <- rename(CRB1,CRB1 = rho)
KRT18 <- rename(KRT18,KRT18 = rho)
DSP <- rename(DSP,DSP = rho)
PKP1 <- rename(PKP1,PKP1 = rho)
CDH1 <- rename(CDH1,CDH1 = rho)
NOX4 <- rename(NOX4,NOX4 = rho)

# newly added

COL1A1 <- rename(COL1A1, COL1A1=rho) # collagen type 1 alpha 1 chain, valid
ITGA5 <- rename(ITGA5, ITGA5=rho) # alpha 5 subunit of integrin, valid
TGFB1 <- rename(TGFB1, TGFB1=rho)
TGFB2 <- rename(TGFB2, TGFB2=rho)
TGFB3 <- rename(TGFB3, TGFB3=rho)
MMP3 <- rename(MMP3, MMP3=rho)
WNT5A <- rename(WNT5A, WNT5A=rho)
WNT5B <- rename(WNT5B, WNT5B=rho)
COL1A2 <- rename(COL1A2, COL1A2=rho) # collagen type I alpha 2 chain
COL3A1 <- rename(COL3A1, COL3A1=rho)
COL5A2 <- rename(COL5A2, COL5A2=rho)
ILK <- rename(ILK, ILK=rho) # integrin linked kinase
CXCR4 <- rename(CXCR4, CXCR4=rho)
CXCL12 <- rename(CXCL12, CXCL12=rho) # chemokine, also known as SDF1
ITGAV <- rename(ITGAV, ITGAV=rho)
TGFBR1 <- rename(TGFBR1, TGFBR1=rho)

# combine to form one df, using inner_join to ensure mut corresponds
rho.comb <- inner_join(FN1, VIM, by= "P53.Mut") %>%
  # inner_join(., CDH2, by = "P53.Mut") %>%
  # inner_join(., ZEB1, by = "P53.Mut") %>%
  # inner_join(., ZEB2, by = "P53.Mut") %>%
  inner_join(., MMP9, by = "P53.Mut") %>%
  inner_join(., MMP2, by = "P53.Mut") %>%
  inner_join(., TWIST1, by = "P53.Mut") %>%
  inner_join(., SNAI1, by = "P53.Mut") %>%
  inner_join(., SNAI2, by = "P53.Mut") %>%
  inner_join(., ACTA2, by = "P53.Mut") %>%
  inner_join(., CLDN1, by = "P53.Mut") %>%
  inner_join(., TJP1, by = "P53.Mut") %>%
  inner_join(., CRB1, by = "P53.Mut") %>%
  inner_join(., KRT18, by = "P53.Mut") %>%
  inner_join(., DSP, by = "P53.Mut") %>%
  inner_join(., PKP1, by = "P53.Mut") %>%
  inner_join(., CDH1, by = "P53.Mut") %>% 
  ### recently added list
  
  inner_join(., COL1A1, by = "P53.Mut") %>% # collagen type 1 alpha 1 chain, valid
  inner_join(., ITGA5, by = "P53.Mut") %>% # alpha 5 subunit of integrin, valid
  inner_join(., TGFB1, by = "P53.Mut") %>%
  inner_join(., TGFB2, by = "P53.Mut") %>%
  inner_join(., TGFB3, by = "P53.Mut") %>%
  inner_join(., MMP3, by = "P53.Mut") %>%
  inner_join(., COL1A2, by = "P53.Mut") %>% # collagen type I alpha 2 chain
  inner_join(., COL3A1, by = "P53.Mut") %>%
  inner_join(., COL5A2, by = "P53.Mut") %>%
  inner_join(., ILK, by = "P53.Mut") %>% # integrin linked kinase
  inner_join(., CXCR4, by = "P53.Mut") %>%
  inner_join(., CXCL12, by = "P53.Mut") %>% # chemokine, also known as SDF1
  inner_join(., ITGAV, by = "P53.Mut") 
# ggplot geom raster wants long form  
rho.long <- gather(rho.comb, key = "Markers", value = "rho", -P53.Mut)
rho.long$Markers <- as.factor(rho.long$Markers)
# reorder EMT markers by type
# EPITHELIAL: CL1, OCC, E-CAD, PLAK, CRUMB3 CYTOKERATINS, ZO1
# MESENCH: FN1, N-CAD, MMP, TWIST1, SNAI1,2, ACTA2, VIM

rho.long$Markers <- factor(rho.long$Markers, levels = c(
  # MESENCHYMAL
  "ACTA2", 
  "COL1A1", # collagen type 1 alpha 1 chain, valid
  "COL1A2", # collagen type I alpha 2 chain
  "COL3A1",
  "COL5A2",
  "CXCR4",
  "CXCL12", # chemokine, also known as SDF1
  "FN1",
  "ILK", # integrin linked kinase
  "ITGA5", # alpha 5 subunit of integrin, valid
  "ITGAV",
  "MMP2",
  "MMP3",
  "MMP9",
  "SNAI1",
  "SNAI2",
  "TGFB1",
  "TGFB2",
  "TGFB3",
  "TWIST1",
  "VIM",
  
  # EPITHELIAL
  "CDH1",
  "CLDN1",
  "CRB1",
  "DSP",
  "KRT18", # keratin
  "PKP1",
  "TJP1"
)
)

# reorder by numeric then put WT first
# fist extract and order the mutations by gene location
p53.factors <- rho.long %>% 
  group_by(P53.Mut) %>% 
  summarize(l = length(rho))

p53.factors$position.gene <- str_sub(p53.factors$P53.Mut, start = 2, end = 4)
p53.factors <-p53.factors[order(p53.factors$position.gene),]
# extract as a vector of characters
gene.position <- as.character(p53.factors$P53.Mut)

# reorder the rho data frame with new level
rho.long$P53.Mut <- factor(rho.long$P53.Mut,
                           levels = c(gene.position))


# now we can create a heat map
NOX4.hm <- ggplot(rho.long) +
  aes(x = Markers, y = P53.Mut, fill = rho) +
  geom_raster() +
  xlab(NULL) +
  ylab(NULL) +
  scale_fill_gradient2(mid = "white", low = "royalblue4", high = "firebrick1",
                       limits = c(-1, 1)) + # limit sets the range of color to use
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 15, vjust=0.5),
        axis.text.y = element_text(size = 15)) +
  geom_text(aes(label = round(rho, 1))) 
# ggtitle("rho values of mRNA expression of various markers and NOX4")
ggsave(NOX4.hm, file = "NOX4_EMT_rho_PanCan_V1.png", dpi=900, width = 11, height = 6.5)

NOX4.hm


### NOX1-EMT correlations pan cancer with p53 ----

# this function will extract NOX1 and a gene (fa) from the df = Pan.final and create a matrix 
# of correlation grouped by p53 mutation status
# and also create a column of
calc.rho <- function(fa) { # use this function to calculate the rho between two genes in a long.form data frame
  gene  <- Pan.final[Pan.final$Gene.Symbol %in% fa,] # use genes in ""
  NOX1 <- Pan.final[Pan.final$Gene.Symbol %in% "NOX1",]
  co <- inner_join(x= gene, y= NOX1, by = "Patient.ID")
  co2 <- co %>% 
    select(
      Patient.ID,
      P53.Mut = P53.Mutation.x,
      gene = mRNA.Value.x,
      NOX1 = mRNA.Value.y,
      -P53.Mutation.y,
      -Gene.Symbol.x,
      -Gene.Symbol.y,
      -Case.Study.x,
      -Case.Study.y
    )
  com <- co2[!duplicated(co2),]
  corr <- com %>% 
    group_by(P53.Mut) %>% 
    summarise(
      rho = cor(NOX1, gene, method = "spearman", use = "complete.obs"))
  corr <- mutate(corr, 
                 P53.Mut = as.factor(P53.Mut))
  corr <- mutate(corr, rho = round(rho, 2))
  data.frame(corr)
}
calc.rho.n <- function(fa) { # use this function to calculate the rho between two genes in a long.form data frame
  gene  <- Pan.final[Pan.final$Gene.Symbol %in% fa,] # use genes in ""
  NOX1 <- Pan.final[Pan.final$Gene.Symbol %in% "NOX1",]
  co <- inner_join(x= gene, y= NOX1, by = "Patient.ID")
  co2 <- co %>% 
    select(
      Patient.ID,
      P53.Mut = P53.Mutation.x,
      gene = mRNA.Value.x,
      NOX1 = mRNA.Value.y,
      -P53.Mutation.y,
      -Gene.Symbol.x,
      -Gene.Symbol.y,
      -Case.Study.x,
      -Case.Study.y
    )
  com <- co2[!duplicated(co2),]
  corr <- com %>% 
    group_by(P53.Mut) %>% 
    summarise(
      rho = cor(NOX1, gene, method = "spearman", use = "complete.obs"),
      n = length(gene))
  corr <- mutate(corr, 
                 P53.Mut = as.factor(P53.Mut))
  corr <- mutate(corr, rho = round(rho, 2))
  data.frame(corr)
}

# frequency of p53 mutation in selected patients
# P53_com.sub will then be merged with pan-cancer set to obtain p53 mutation for each patient
Pan.select.p53 <- inner_join(P53.final, Pan.RNA, by = "Patient.ID")

p53.freq <- Pan.select.p53 %>% 
  group_by(P53.Mutation) %>% 
  summarise(n = length(P53.Mutation))


# Now, Name the mutations of interest 
mut.interest <- c(
  "R175H",
  "R248Q",
  "R273H",
  # "R280K",
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
  "G248S",
  "R158L",
  "C175F",
  "H179R",
  "G245S",
  # "R175G",
  "WT"
)

Pan.final <- Pan.select.p53[Pan.select.p53$P53.Mutation %in% mut.interest,]
# this is a good place to clean up the environment
#rm(list=setdiff(ls(), c("Pan.final", "calc.rho", "EMT.Genes", "Pan.RNA")))

# invidually calculating and creating object (learn how to loop better)
FN1.NOX1.n <- calc.rho.n("FN1")
FN1 <- calc.rho("FN1")
VIM <- calc.rho("VIM")
CDH2 <- calc.rho("CDH2")
ZEB1 <- calc.rho("ZEB1")
ZEB2 <- calc.rho("ZEB2")
MMP9 <- calc.rho("MMP9")
MMP2 <- calc.rho("MMP2")
TWIST1 <- calc.rho("TWIST1")
SNAI1 <- calc.rho("SNAI1")
SNAI2 <- calc.rho("SNAI2")
ACTA2 <- calc.rho("ACTA2")
CLDN1 <- calc.rho("CLDN1")
TJP1 <- calc.rho("TJP1")
CRB1 <- calc.rho("CRB1")
KRT18 <- calc.rho("KRT18")
DSP <- calc.rho("DSP")
PKP1 <- calc.rho("PKP1")
CDH1 <- calc.rho("CDH1")
NOX1 <- calc.rho("NOX1")
### recently added list

COL1A1 <- calc.rho("COL1A1") # collagen type 1 alpha 1 chain, valid
ITGA5 <- calc.rho("ITGA5") # alpha 5 subunit of integrin, valid
TGFB1 <- calc.rho("TGFB1")
TGFB2 <- calc.rho("TGFB2")
TGFB3 <- calc.rho("TGFB3")
MMP3 <- calc.rho("MMP3")
WNT5A <- calc.rho("WNT5A")
WNT5B <- calc.rho("WNT5B")
COL1A2 <- calc.rho("COL1A2") # collagen type I alpha 2 chain
COL3A1 <- calc.rho("COL3A1")
COL5A2 <- calc.rho("COL5A2")
ILK <- calc.rho("ILK") # integrin linked kinase
CXCR4 <- calc.rho("CXCR4")
CXCL12 <- calc.rho("CXCL12") # chemokine, also known as SDF1
ITGAV <- calc.rho("ITGAV")
TGFBR1 <- calc.rho("TGFBR1")
TJP1 <- calc.rho("TJP1") # zonuli1
CRB1 <- calc.rho("CRB1") # cell polarity complex crumbs1, valid
KRT18 <- calc.rho("KRT18") # keratin 18, valid 

# rename the rho in each dataframe to its respective gene_NOX1
FN1 <- rename(FN1,FN1 = rho)
VIM <- rename(VIM,VIM = rho)
CDH2 <- rename(CDH2,CDH2 = rho)
ZEB1 <- rename(ZEB1,ZEB1 = rho)
ZEB2 <- rename(ZEB2,ZEB2 = rho)
MMP9 <- rename(MMP9,MMP9 = rho)
MMP2 <- rename(MMP2,MMP2 = rho)
TWIST1 <- rename(TWIST1,TWIST1 = rho)
SNAI1 <- rename(SNAI1,SNAI1 = rho)
SNAI2 <- rename(SNAI2,SNAI2 = rho)
ACTA2 <- rename(ACTA2,ACTA2 = rho)
CLDN1 <- rename(CLDN1,CLDN1 = rho)
TJP1 <- rename(TJP1,TJP1 = rho)
CRB1 <- rename(CRB1,CRB1 = rho)
KRT18 <- rename(KRT18,KRT18 = rho)
DSP <- rename(DSP,DSP = rho)
PKP1 <- rename(PKP1,PKP1 = rho)
CDH1 <- rename(CDH1,CDH1 = rho)
NOX1 <- rename(NOX1,NOX1 = rho)

# newly added

COL1A1 <- rename(COL1A1, COL1A1=rho) # collagen type 1 alpha 1 chain, valid
ITGA5 <- rename(ITGA5, ITGA5=rho) # alpha 5 subunit of integrin, valid
TGFB1 <- rename(TGFB1, TGFB1=rho)
TGFB2 <- rename(TGFB2, TGFB2=rho)
TGFB3 <- rename(TGFB3, TGFB3=rho)
MMP3 <- rename(MMP3, MMP3=rho)
WNT5A <- rename(WNT5A, WNT5A=rho)
WNT5B <- rename(WNT5B, WNT5B=rho)
COL1A2 <- rename(COL1A2, COL1A2=rho) # collagen type I alpha 2 chain
COL3A1 <- rename(COL3A1, COL3A1=rho)
COL5A2 <- rename(COL5A2, COL5A2=rho)
ILK <- rename(ILK, ILK=rho) # integrin linked kinase
CXCR4 <- rename(CXCR4, CXCR4=rho)
CXCL12 <- rename(CXCL12, CXCL12=rho) # chemokine, also known as SDF1
ITGAV <- rename(ITGAV, ITGAV=rho)
TGFBR1 <- rename(TGFBR1, TGFBR1=rho)

# combine to form one df, using inner_join to ensure mut corresponds
rho.comb <- inner_join(FN1, VIM, by= "P53.Mut") %>%
  # inner_join(., CDH2, by = "P53.Mut") %>%
  # inner_join(., ZEB1, by = "P53.Mut") %>%
  # inner_join(., ZEB2, by = "P53.Mut") %>%
  inner_join(., MMP9, by = "P53.Mut") %>%
  inner_join(., MMP2, by = "P53.Mut") %>%
  inner_join(., TWIST1, by = "P53.Mut") %>%
  inner_join(., SNAI1, by = "P53.Mut") %>%
  inner_join(., SNAI2, by = "P53.Mut") %>%
  inner_join(., ACTA2, by = "P53.Mut") %>%
  inner_join(., CLDN1, by = "P53.Mut") %>%
  inner_join(., TJP1, by = "P53.Mut") %>%
  inner_join(., CRB1, by = "P53.Mut") %>%
  inner_join(., KRT18, by = "P53.Mut") %>%
  inner_join(., DSP, by = "P53.Mut") %>%
  inner_join(., PKP1, by = "P53.Mut") %>%
  inner_join(., CDH1, by = "P53.Mut") %>% 
  ### recently added list
  
  inner_join(., COL1A1, by = "P53.Mut") %>% # collagen type 1 alpha 1 chain, valid
  inner_join(., ITGA5, by = "P53.Mut") %>% # alpha 5 subunit of integrin, valid
  inner_join(., TGFB1, by = "P53.Mut") %>%
  inner_join(., TGFB2, by = "P53.Mut") %>%
  inner_join(., TGFB3, by = "P53.Mut") %>%
  inner_join(., MMP3, by = "P53.Mut") %>%
  inner_join(., COL1A2, by = "P53.Mut") %>% # collagen type I alpha 2 chain
  inner_join(., COL3A1, by = "P53.Mut") %>%
  inner_join(., COL5A2, by = "P53.Mut") %>%
  inner_join(., ILK, by = "P53.Mut") %>% # integrin linked kinase
  inner_join(., CXCR4, by = "P53.Mut") %>%
  inner_join(., CXCL12, by = "P53.Mut") %>% # chemokine, also known as SDF1
  inner_join(., ITGAV, by = "P53.Mut") 
# ggplot geom raster wants long form  
rho.long <- gather(rho.comb, key = "Markers", value = "rho", -P53.Mut)
rho.long$Markers <- as.factor(rho.long$Markers)
# reorder EMT markers by type
# EPITHELIAL: CL1, OCC, E-CAD, PLAK, CRUMB3 CYTOKERATINS, ZO1
# MESENCH: FN1, N-CAD, MMP, TWIST1, SNAI1,2, ACTA2, VIM
rho.long$Markers <- factor(rho.long$Markers, levels = c(
  # MESENCHYMAL
  "ACTA2", 
  "COL1A1", # collagen type 1 alpha 1 chain, valid
  "COL1A2", # collagen type I alpha 2 chain
  "COL3A1",
  "COL5A2",
  "CXCR4",
  "CXCL12", # chemokine, also known as SDF1
  "FN1",
  "ILK", # integrin linked kinase
  "ITGA5", # alpha 5 subunit of integrin, valid
  "ITGAV",
  "MMP2",
  "MMP3",
  "MMP9",
  "SNAI1",
  "SNAI2",
  "TGFB1",
  "TGFB2",
  "TGFB3",
  "TWIST1",
  "VIM",
  
  # EPITHELIAL
  "CDH1",
  "CLDN1",
  "CRB1",
  "DSP",
  "KRT18", # keratin
  "PKP1",
  "TJP1"
)
)
# reorder by numeric then put WT first
# fist extract and order the mutations by gene location
p53.factors <- rho.long %>% 
  group_by(P53.Mut) %>% 
  summarize(l = length(rho))

p53.factors$position.gene <- str_sub(p53.factors$P53.Mut, start = 2, end = 4)
p53.factors <-p53.factors[order(p53.factors$position.gene),]
# extract as a vector of characters
gene.position <- as.character(p53.factors$P53.Mut)

# reorder the rho data frame with new level
rho.long$P53.Mut <- factor(rho.long$P53.Mut,
                           levels = c(gene.position))


# now we can create a heat map
NOX1.hm <- ggplot(rho.long) +
  aes(x = Markers, y = P53.Mut, fill = rho) +
  geom_raster() +
  xlab(NULL) +
  ylab(NULL) +
  scale_fill_gradient2(mid = "white", low = "royalblue4", high = "firebrick1",
                       limits = c(-1, 1)) + # limit sets the range of color to use
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 15, vjust=0.5),
        axis.text.y = element_text(size = 15)) +
  geom_text(aes(label = round(rho, 1))) 
# ggtitle("rho values of mRNA expression of various markers and NOX1")
ggsave(NOX1.hm, file = "NOX1_EMT_rho_PanCan_V1.png", dpi=900, width = 11, height = 6.5)
NOX1.hm


### DUOX1-EMT correlations pan cancer with p53 ----

# this function will extract DUOX1 and a gene (fa) from the df = Pan.final and create a matrix 
# of correlation grouped by p53 mutation status
# and also create a column of
calc.rho <- function(fa) { # use this function to calculate the rho between two genes in a long.form data frame
  gene  <- Pan.final[Pan.final$Gene.Symbol %in% fa,] # use genes in ""
  DUOX1 <- Pan.final[Pan.final$Gene.Symbol %in% "DUOX1",]
  co <- inner_join(x= gene, y= DUOX1, by = "Patient.ID")
  co2 <- co %>% 
    select(
      Patient.ID,
      P53.Mut = P53.Mutation.x,
      gene = mRNA.Value.x,
      DUOX1 = mRNA.Value.y,
      -P53.Mutation.y,
      -Gene.Symbol.x,
      -Gene.Symbol.y,
      -Case.Study.x,
      -Case.Study.y
    )
  com <- co2[!duplicated(co2),]
  corr <- com %>% 
    group_by(P53.Mut) %>% 
    summarise(
      rho = cor(DUOX1, gene, method = "spearman", use = "complete.obs"))
  corr <- mutate(corr, 
                 P53.Mut = as.factor(P53.Mut))
  corr <- mutate(corr, rho = round(rho, 2))
  data.frame(corr)
}
calc.rho.n <- function(fa) { # use this function to calculate the rho between two genes in a long.form data frame
  gene  <- Pan.final[Pan.final$Gene.Symbol %in% fa,] # use genes in ""
  DUOX1 <- Pan.final[Pan.final$Gene.Symbol %in% "DUOX1",]
  co <- inner_join(x= gene, y= DUOX1, by = "Patient.ID")
  co2 <- co %>% 
    select(
      Patient.ID,
      P53.Mut = P53.Mutation.x,
      gene = mRNA.Value.x,
      DUOX1 = mRNA.Value.y,
      -P53.Mutation.y,
      -Gene.Symbol.x,
      -Gene.Symbol.y,
      -Case.Study.x,
      -Case.Study.y
    )
  com <- co2[!duplicated(co2),]
  corr <- com %>% 
    group_by(P53.Mut) %>% 
    summarise(
      rho = cor(DUOX1, gene, method = "spearman", use = "complete.obs"),
      n = length(gene))
  corr <- mutate(corr, 
                 P53.Mut = as.factor(P53.Mut))
  corr <- mutate(corr, rho = round(rho, 2))
  data.frame(corr)
}

# frequency of p53 mutation in selected patients
# P53_com.sub will then be merged with pan-cancer set to obtain p53 mutation for each patient
Pan.select.p53 <- inner_join(P53.final, Pan.RNA, by = "Patient.ID")

p53.freq <- Pan.select.p53 %>% 
  group_by(P53.Mutation) %>% 
  summarise(n = length(P53.Mutation))


# Now, Name the mutations of interest 
mut.interest <- c(
  "R175H",
  "R248Q",
  "R273H",
  # "R280K",
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
  "G248S",
  "R158L",
  "C175F",
  "H179R",
  "G245S",
  # "R175G",
  "WT"
)

Pan.final <- Pan.select.p53[Pan.select.p53$P53.Mutation %in% mut.interest,]
# this is a good place to clean up the environment
#rm(list=setdiff(ls(), c("Pan.final", "calc.rho", "EMT.Genes", "Pan.RNA")))

# invidually calculating and creating object (learn how to loop better)
FN1.DUOX1.n <- calc.rho.n("FN1")
FN1 <- calc.rho("FN1")
VIM <- calc.rho("VIM")
CDH2 <- calc.rho("CDH2")
ZEB1 <- calc.rho("ZEB1")
ZEB2 <- calc.rho("ZEB2")
MMP9 <- calc.rho("MMP9")
MMP2 <- calc.rho("MMP2")
TWIST1 <- calc.rho("TWIST1")
SNAI1 <- calc.rho("SNAI1")
SNAI2 <- calc.rho("SNAI2")
ACTA2 <- calc.rho("ACTA2")
CLDN1 <- calc.rho("CLDN1")
TJP1 <- calc.rho("TJP1")
CRB1 <- calc.rho("CRB1")
KRT18 <- calc.rho("KRT18")
DSP <- calc.rho("DSP")
PKP1 <- calc.rho("PKP1")
CDH1 <- calc.rho("CDH1")
DUOX1 <- calc.rho("DUOX1")
### recently added list

COL1A1 <- calc.rho("COL1A1") # collagen type 1 alpha 1 chain, valid
ITGA5 <- calc.rho("ITGA5") # alpha 5 subunit of integrin, valid
TGFB1 <- calc.rho("TGFB1")
TGFB2 <- calc.rho("TGFB2")
TGFB3 <- calc.rho("TGFB3")
MMP3 <- calc.rho("MMP3")
WNT5A <- calc.rho("WNT5A")
WNT5B <- calc.rho("WNT5B")
COL1A2 <- calc.rho("COL1A2") # collagen type I alpha 2 chain
COL3A1 <- calc.rho("COL3A1")
COL5A2 <- calc.rho("COL5A2")
ILK <- calc.rho("ILK") # integrin linked kinase
CXCR4 <- calc.rho("CXCR4")
CXCL12 <- calc.rho("CXCL12") # chemokine, also known as SDF1
ITGAV <- calc.rho("ITGAV")
TGFBR1 <- calc.rho("TGFBR1")
TJP1 <- calc.rho("TJP1") # zonuli1
CRB1 <- calc.rho("CRB1") # cell polarity complex crumbs1, valid
KRT18 <- calc.rho("KRT18") # keratin 18, valid 

# rename the rho in each dataframe to its respective gene_DUOX1
FN1 <- rename(FN1,FN1 = rho)
VIM <- rename(VIM,VIM = rho)
CDH2 <- rename(CDH2,CDH2 = rho)
ZEB1 <- rename(ZEB1,ZEB1 = rho)
ZEB2 <- rename(ZEB2,ZEB2 = rho)
MMP9 <- rename(MMP9,MMP9 = rho)
MMP2 <- rename(MMP2,MMP2 = rho)
TWIST1 <- rename(TWIST1,TWIST1 = rho)
SNAI1 <- rename(SNAI1,SNAI1 = rho)
SNAI2 <- rename(SNAI2,SNAI2 = rho)
ACTA2 <- rename(ACTA2,ACTA2 = rho)
CLDN1 <- rename(CLDN1,CLDN1 = rho)
TJP1 <- rename(TJP1,TJP1 = rho)
CRB1 <- rename(CRB1,CRB1 = rho)
KRT18 <- rename(KRT18,KRT18 = rho)
DSP <- rename(DSP,DSP = rho)
PKP1 <- rename(PKP1,PKP1 = rho)
CDH1 <- rename(CDH1,CDH1 = rho)
DUOX1 <- rename(DUOX1,DUOX1 = rho)

# newly added

COL1A1 <- rename(COL1A1, COL1A1=rho) # collagen type 1 alpha 1 chain, valid
ITGA5 <- rename(ITGA5, ITGA5=rho) # alpha 5 subunit of integrin, valid
TGFB1 <- rename(TGFB1, TGFB1=rho)
TGFB2 <- rename(TGFB2, TGFB2=rho)
TGFB3 <- rename(TGFB3, TGFB3=rho)
MMP3 <- rename(MMP3, MMP3=rho)
WNT5A <- rename(WNT5A, WNT5A=rho)
WNT5B <- rename(WNT5B, WNT5B=rho)
COL1A2 <- rename(COL1A2, COL1A2=rho) # collagen type I alpha 2 chain
COL3A1 <- rename(COL3A1, COL3A1=rho)
COL5A2 <- rename(COL5A2, COL5A2=rho)
ILK <- rename(ILK, ILK=rho) # integrin linked kinase
CXCR4 <- rename(CXCR4, CXCR4=rho)
CXCL12 <- rename(CXCL12, CXCL12=rho) # chemokine, also known as SDF1
ITGAV <- rename(ITGAV, ITGAV=rho)
TGFBR1 <- rename(TGFBR1, TGFBR1=rho)

# combine to form one df, using inner_join to ensure mut corresponds
rho.comb <- inner_join(FN1, VIM, by= "P53.Mut") %>%
  # inner_join(., CDH2, by = "P53.Mut") %>%
  # inner_join(., ZEB1, by = "P53.Mut") %>%
  # inner_join(., ZEB2, by = "P53.Mut") %>%
  inner_join(., MMP9, by = "P53.Mut") %>%
  inner_join(., MMP2, by = "P53.Mut") %>%
  inner_join(., TWIST1, by = "P53.Mut") %>%
  inner_join(., SNAI1, by = "P53.Mut") %>%
  inner_join(., SNAI2, by = "P53.Mut") %>%
  inner_join(., ACTA2, by = "P53.Mut") %>%
  inner_join(., CLDN1, by = "P53.Mut") %>%
  inner_join(., TJP1, by = "P53.Mut") %>%
  inner_join(., CRB1, by = "P53.Mut") %>%
  inner_join(., KRT18, by = "P53.Mut") %>%
  inner_join(., DSP, by = "P53.Mut") %>%
  inner_join(., PKP1, by = "P53.Mut") %>%
  inner_join(., CDH1, by = "P53.Mut") %>% 
  ### recently added list
  
  inner_join(., COL1A1, by = "P53.Mut") %>% # collagen type 1 alpha 1 chain, valid
  inner_join(., ITGA5, by = "P53.Mut") %>% # alpha 5 subunit of integrin, valid
  inner_join(., TGFB1, by = "P53.Mut") %>%
  inner_join(., TGFB2, by = "P53.Mut") %>%
  inner_join(., TGFB3, by = "P53.Mut") %>%
  inner_join(., MMP3, by = "P53.Mut") %>%
  inner_join(., COL1A2, by = "P53.Mut") %>% # collagen type I alpha 2 chain
  inner_join(., COL3A1, by = "P53.Mut") %>%
  inner_join(., COL5A2, by = "P53.Mut") %>%
  inner_join(., ILK, by = "P53.Mut") %>% # integrin linked kinase
  inner_join(., CXCR4, by = "P53.Mut") %>%
  inner_join(., CXCL12, by = "P53.Mut") %>% # chemokine, also known as SDF1
  inner_join(., ITGAV, by = "P53.Mut") 
# ggplot geom raster wants long form  
rho.long <- gather(rho.comb, key = "Markers", value = "rho", -P53.Mut)
rho.long$Markers <- as.factor(rho.long$Markers)
# reorder EMT markers by type
# EPITHELIAL: CL1, OCC, E-CAD, PLAK, CRUMB3 CYTOKERATINS, ZO1
# MESENCH: FN1, N-CAD, MMP, TWIST1, SNAI1,2, ACTA2, VIM

rho.long$Markers <- factor(rho.long$Markers, levels = c(
  # MESENCHYMAL
  "ACTA2", 
  "COL1A1", # collagen type 1 alpha 1 chain, valid
  "COL1A2", # collagen type I alpha 2 chain
  "COL3A1",
  "COL5A2",
  "CXCR4",
  "CXCL12", # chemokine, also known as SDF1
  "FN1",
  "ILK", # integrin linked kinase
  "ITGA5", # alpha 5 subunit of integrin, valid
  "ITGAV",
  "MMP2",
  "MMP3",
  "MMP9",
  "SNAI1",
  "SNAI2",
  "TGFB1",
  "TGFB2",
  "TGFB3",
  "TWIST1",
  "VIM",
  
  # EPITHELIAL
  "CDH1",
  "CLDN1",
  "CRB1",
  "DSP",
  "KRT18", # keratin
  "PKP1",
  "TJP1"
)
)

# reorder by numeric then put WT first
# fist extract and order the mutations by gene location
p53.factors <- rho.long %>% 
  group_by(P53.Mut) %>% 
  summarize(l = length(rho))

p53.factors$position.gene <- str_sub(p53.factors$P53.Mut, start = 2, end = 4)
p53.factors <-p53.factors[order(p53.factors$position.gene),]
# extract as a vector of characters
gene.position <- as.character(p53.factors$P53.Mut)

# reorder the rho data frame with new level
rho.long$P53.Mut <- factor(rho.long$P53.Mut,
                           levels = c(gene.position))


# now we can create a heat map
DUOX1.hm <- ggplot(rho.long) +
  aes(x = Markers, y = P53.Mut, fill = rho) +
  geom_raster() +
  xlab(NULL) +
  ylab(NULL) +
  scale_fill_gradient2(mid = "white", low = "royalblue4", high = "firebrick1",
                       limits = c(-1, 1)) + # limit sets the range of color to use
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 15, vjust=0.5),
        axis.text.y = element_text(size = 15)) +
  geom_text(aes(label = round(rho, 1))) 
# ggtitle("rho values of mRNA expression of various markers and DUOX1")
ggsave(DUOX1.hm, file = "DUOX1_EMT_rho_PanCan_V1.png", dpi=900, width = 11, height = 6.5)
DUOX1.hm


### DUOX2-EMT correlations pan cancer with p53 ----

# this function will extract DUOX2 and a gene (fa) from the df = Pan.final and create a matrix 
# of correlation grouped by p53 mutation status
# and also create a column of
calc.rho <- function(fa) { # use this function to calculate the rho between two genes in a long.form data frame
  gene  <- Pan.final[Pan.final$Gene.Symbol %in% fa,] # use genes in ""
  DUOX2 <- Pan.final[Pan.final$Gene.Symbol %in% "DUOX2",]
  co <- inner_join(x= gene, y= DUOX2, by = "Patient.ID")
  co2 <- co %>% 
    select(
      Patient.ID,
      P53.Mut = P53.Mutation.x,
      gene = mRNA.Value.x,
      DUOX2 = mRNA.Value.y,
      -P53.Mutation.y,
      -Gene.Symbol.x,
      -Gene.Symbol.y,
      -Case.Study.x,
      -Case.Study.y
    )
  com <- co2[!duplicated(co2),]
  corr <- com %>% 
    group_by(P53.Mut) %>% 
    summarise(
      rho = cor(DUOX2, gene, method = "spearman", use = "complete.obs"))
  corr <- mutate(corr, 
                 P53.Mut = as.factor(P53.Mut))
  corr <- mutate(corr, rho = round(rho, 2))
  data.frame(corr)
}
calc.rho.n <- function(fa) { # use this function to calculate the rho between two genes in a long.form data frame
  gene  <- Pan.final[Pan.final$Gene.Symbol %in% fa,] # use genes in ""
  DUOX2 <- Pan.final[Pan.final$Gene.Symbol %in% "DUOX2",]
  co <- inner_join(x= gene, y= DUOX2, by = "Patient.ID")
  co2 <- co %>% 
    select(
      Patient.ID,
      P53.Mut = P53.Mutation.x,
      gene = mRNA.Value.x,
      DUOX2 = mRNA.Value.y,
      -P53.Mutation.y,
      -Gene.Symbol.x,
      -Gene.Symbol.y,
      -Case.Study.x,
      -Case.Study.y
    )
  com <- co2[!duplicated(co2),]
  corr <- com %>% 
    group_by(P53.Mut) %>% 
    summarise(
      rho = cor(DUOX2, gene, method = "spearman", use = "complete.obs"),
      n = length(gene))
  corr <- mutate(corr, 
                 P53.Mut = as.factor(P53.Mut))
  corr <- mutate(corr, rho = round(rho, 2))
  data.frame(corr)
}

# frequency of p53 mutation in selected patients
# P53_com.sub will then be merged with pan-cancer set to obtain p53 mutation for each patient
Pan.select.p53 <- inner_join(P53.final, Pan.RNA, by = "Patient.ID")

p53.freq <- Pan.select.p53 %>% 
  group_by(P53.Mutation) %>% 
  summarise(n = length(P53.Mutation))


# Now, Name the mutations of interest 
mut.interest <- c(
  "R175H",
  "R248Q",
  "R273H",
  # "R280K",
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
  "G248S",
  "R158L",
  "C175F",
  "H179R",
  "G245S",
  # "R175G",
  "WT"
)

Pan.final <- Pan.select.p53[Pan.select.p53$P53.Mutation %in% mut.interest,]
# this is a good place to clean up the environment
#rm(list=setdiff(ls(), c("Pan.final", "calc.rho", "EMT.Genes", "Pan.RNA")))

# invidually calculating and creating object (learn how to loop better)
FN1.DUOX2.n <- calc.rho.n("FN1")
FN1 <- calc.rho("FN1")
VIM <- calc.rho("VIM")
CDH2 <- calc.rho("CDH2")
ZEB1 <- calc.rho("ZEB1")
ZEB2 <- calc.rho("ZEB2")
MMP9 <- calc.rho("MMP9")
MMP2 <- calc.rho("MMP2")
TWIST1 <- calc.rho("TWIST1")
SNAI1 <- calc.rho("SNAI1")
SNAI2 <- calc.rho("SNAI2")
ACTA2 <- calc.rho("ACTA2")
CLDN1 <- calc.rho("CLDN1")
TJP1 <- calc.rho("TJP1")
CRB1 <- calc.rho("CRB1")
KRT18 <- calc.rho("KRT18")
DSP <- calc.rho("DSP")
PKP1 <- calc.rho("PKP1")
CDH1 <- calc.rho("CDH1")
DUOX2 <- calc.rho("DUOX2")
### recently added list

COL1A1 <- calc.rho("COL1A1") # collagen type 1 alpha 1 chain, valid
ITGA5 <- calc.rho("ITGA5") # alpha 5 subunit of integrin, valid
TGFB1 <- calc.rho("TGFB1")
TGFB2 <- calc.rho("TGFB2")
TGFB3 <- calc.rho("TGFB3")
MMP3 <- calc.rho("MMP3")
WNT5A <- calc.rho("WNT5A")
WNT5B <- calc.rho("WNT5B")
COL1A2 <- calc.rho("COL1A2") # collagen type I alpha 2 chain
COL3A1 <- calc.rho("COL3A1")
COL5A2 <- calc.rho("COL5A2")
ILK <- calc.rho("ILK") # integrin linked kinase
CXCR4 <- calc.rho("CXCR4")
CXCL12 <- calc.rho("CXCL12") # chemokine, also known as SDF1
ITGAV <- calc.rho("ITGAV")
TGFBR1 <- calc.rho("TGFBR1")
TJP1 <- calc.rho("TJP1") # zonuli1
CRB1 <- calc.rho("CRB1") # cell polarity complex crumbs1, valid
KRT18 <- calc.rho("KRT18") # keratin 18, valid 

# rename the rho in each dataframe to its respective gene_DUOX2
FN1 <- rename(FN1,FN1 = rho)
VIM <- rename(VIM,VIM = rho)
CDH2 <- rename(CDH2,CDH2 = rho)
ZEB1 <- rename(ZEB1,ZEB1 = rho)
ZEB2 <- rename(ZEB2,ZEB2 = rho)
MMP9 <- rename(MMP9,MMP9 = rho)
MMP2 <- rename(MMP2,MMP2 = rho)
TWIST1 <- rename(TWIST1,TWIST1 = rho)
SNAI1 <- rename(SNAI1,SNAI1 = rho)
SNAI2 <- rename(SNAI2,SNAI2 = rho)
ACTA2 <- rename(ACTA2,ACTA2 = rho)
CLDN1 <- rename(CLDN1,CLDN1 = rho)
TJP1 <- rename(TJP1,TJP1 = rho)
CRB1 <- rename(CRB1,CRB1 = rho)
KRT18 <- rename(KRT18,KRT18 = rho)
DSP <- rename(DSP,DSP = rho)
PKP1 <- rename(PKP1,PKP1 = rho)
CDH1 <- rename(CDH1,CDH1 = rho)
DUOX2 <- rename(DUOX2,DUOX2 = rho)

# newly added

COL1A1 <- rename(COL1A1, COL1A1=rho) # collagen type 1 alpha 1 chain, valid
ITGA5 <- rename(ITGA5, ITGA5=rho) # alpha 5 subunit of integrin, valid
TGFB1 <- rename(TGFB1, TGFB1=rho)
TGFB2 <- rename(TGFB2, TGFB2=rho)
TGFB3 <- rename(TGFB3, TGFB3=rho)
MMP3 <- rename(MMP3, MMP3=rho)
WNT5A <- rename(WNT5A, WNT5A=rho)
WNT5B <- rename(WNT5B, WNT5B=rho)
COL1A2 <- rename(COL1A2, COL1A2=rho) # collagen type I alpha 2 chain
COL3A1 <- rename(COL3A1, COL3A1=rho)
COL5A2 <- rename(COL5A2, COL5A2=rho)
ILK <- rename(ILK, ILK=rho) # integrin linked kinase
CXCR4 <- rename(CXCR4, CXCR4=rho)
CXCL12 <- rename(CXCL12, CXCL12=rho) # chemokine, also known as SDF1
ITGAV <- rename(ITGAV, ITGAV=rho)
TGFBR1 <- rename(TGFBR1, TGFBR1=rho)

# combine to form one df, using inner_join to ensure mut corresponds
rho.comb <- inner_join(FN1, VIM, by= "P53.Mut") %>%
  # inner_join(., CDH2, by = "P53.Mut") %>%
  # inner_join(., ZEB1, by = "P53.Mut") %>%
  # inner_join(., ZEB2, by = "P53.Mut") %>%
  inner_join(., MMP9, by = "P53.Mut") %>%
  inner_join(., MMP2, by = "P53.Mut") %>%
  inner_join(., TWIST1, by = "P53.Mut") %>%
  inner_join(., SNAI1, by = "P53.Mut") %>%
  inner_join(., SNAI2, by = "P53.Mut") %>%
  inner_join(., ACTA2, by = "P53.Mut") %>%
  inner_join(., CLDN1, by = "P53.Mut") %>%
  inner_join(., TJP1, by = "P53.Mut") %>%
  inner_join(., CRB1, by = "P53.Mut") %>%
  inner_join(., KRT18, by = "P53.Mut") %>%
  inner_join(., DSP, by = "P53.Mut") %>%
  inner_join(., PKP1, by = "P53.Mut") %>%
  inner_join(., CDH1, by = "P53.Mut") %>% 
  ### recently added list
  
  inner_join(., COL1A1, by = "P53.Mut") %>% # collagen type 1 alpha 1 chain, valid
  inner_join(., ITGA5, by = "P53.Mut") %>% # alpha 5 subunit of integrin, valid
  inner_join(., TGFB1, by = "P53.Mut") %>%
  inner_join(., TGFB2, by = "P53.Mut") %>%
  inner_join(., TGFB3, by = "P53.Mut") %>%
  inner_join(., MMP3, by = "P53.Mut") %>%
  inner_join(., COL1A2, by = "P53.Mut") %>% # collagen type I alpha 2 chain
  inner_join(., COL3A1, by = "P53.Mut") %>%
  inner_join(., COL5A2, by = "P53.Mut") %>%
  inner_join(., ILK, by = "P53.Mut") %>% # integrin linked kinase
  inner_join(., CXCR4, by = "P53.Mut") %>%
  inner_join(., CXCL12, by = "P53.Mut") %>% # chemokine, also known as SDF1
  inner_join(., ITGAV, by = "P53.Mut") 
# ggplot geom raster wants long form  
rho.long <- gather(rho.comb, key = "Markers", value = "rho", -P53.Mut)
rho.long$Markers <- as.factor(rho.long$Markers)
# reorder EMT markers by type
# EPITHELIAL: CL1, OCC, E-CAD, PLAK, CRUMB3 CYTOKERATINS, ZO1
# MESENCH: FN1, N-CAD, MMP, TWIST1, SNAI1,2, ACTA2, VIM

rho.long$Markers <- factor(rho.long$Markers, levels = c(
  # MESENCHYMAL
  "ACTA2", 
  "COL1A1", # collagen type 1 alpha 1 chain, valid
  "COL1A2", # collagen type I alpha 2 chain
  "COL3A1",
  "COL5A2",
  "CXCR4",
  "CXCL12", # chemokine, also known as SDF1
  "FN1",
  "ILK", # integrin linked kinase
  "ITGA5", # alpha 5 subunit of integrin, valid
  "ITGAV",
  "MMP2",
  "MMP3",
  "MMP9",
  "SNAI1",
  "SNAI2",
  "TGFB1",
  "TGFB2",
  "TGFB3",
  "TWIST1",
  "VIM",
  
  # EPITHELIAL
  "CDH1",
  "CLDN1",
  "CRB1",
  "DSP",
  "KRT18", # keratin
  "PKP1",
  "TJP1"
)
)

# reorder by numeric then put WT first
# fist extract and order the mutations by gene location
p53.factors <- rho.long %>% 
  group_by(P53.Mut) %>% 
  summarize(l = length(rho))

p53.factors$position.gene <- str_sub(p53.factors$P53.Mut, start = 2, end = 4)
p53.factors <-p53.factors[order(p53.factors$position.gene),]
# extract as a vector of characters
gene.position <- as.character(p53.factors$P53.Mut)

# reorder the rho data frame with new level
rho.long$P53.Mut <- factor(rho.long$P53.Mut,
                           levels = c(gene.position))


# now we can create a heat map
DUOX2.hm <- ggplot(rho.long) +
  aes(x = Markers, y = P53.Mut, fill = rho) +
  geom_raster() +
  xlab(NULL) +
  ylab(NULL) +
  scale_fill_gradient2(mid = "white", low = "royalblue4", high = "firebrick1",
                       limits = c(-1, 1)) + # limit sets the range of color to use
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 15, vjust=0.5),
        axis.text.y = element_text(size = 15)) +
  geom_text(aes(label = round(rho, 1))) 
# ggtitle("rho values of mRNA expression of various markers and DUOX2")
ggsave(DUOX2.hm, file = "DUOX2_EMT_rho_PanCan_V1.png", dpi=900, width = 11, height = 6.5)
DUOX2.hm



### NOX3-EMT correlations pan cancer with p53 ----

# this function will extract NOX3 and a gene (fa) from the df = Pan.final and create a matrix 
# of correlation grouped by p53 mutation status
# and also create a column of
calc.rho <- function(fa) { # use this function to calculate the rho between two genes in a long.form data frame
  gene  <- Pan.final[Pan.final$Gene.Symbol %in% fa,] # use genes in ""
  NOX3 <- Pan.final[Pan.final$Gene.Symbol %in% "NOX3",]
  co <- inner_join(x= gene, y= NOX3, by = "Patient.ID")
  co2 <- co %>% 
    select(
      Patient.ID,
      P53.Mut = P53.Mutation.x,
      gene = mRNA.Value.x,
      NOX3 = mRNA.Value.y,
      -P53.Mutation.y,
      -Gene.Symbol.x,
      -Gene.Symbol.y,
      -Case.Study.x,
      -Case.Study.y
    )
  com <- co2[!duplicated(co2),]
  corr <- com %>% 
    group_by(P53.Mut) %>% 
    summarise(
      rho = cor(NOX3, gene, method = "spearman", use = "complete.obs"))
  corr <- mutate(corr, 
                 P53.Mut = as.factor(P53.Mut))
  corr <- mutate(corr, rho = round(rho, 2))
  data.frame(corr)
}
calc.rho.n <- function(fa) { # use this function to calculate the rho between two genes in a long.form data frame
  gene  <- Pan.final[Pan.final$Gene.Symbol %in% fa,] # use genes in ""
  NOX3 <- Pan.final[Pan.final$Gene.Symbol %in% "NOX3",]
  co <- inner_join(x= gene, y= NOX3, by = "Patient.ID")
  co2 <- co %>% 
    select(
      Patient.ID,
      P53.Mut = P53.Mutation.x,
      gene = mRNA.Value.x,
      NOX3 = mRNA.Value.y,
      -P53.Mutation.y,
      -Gene.Symbol.x,
      -Gene.Symbol.y,
      -Case.Study.x,
      -Case.Study.y
    )
  com <- co2[!duplicated(co2),]
  corr <- com %>% 
    group_by(P53.Mut) %>% 
    summarise(
      rho = cor(NOX3, gene, method = "spearman", use = "complete.obs"),
      n = length(gene))
  corr <- mutate(corr, 
                 P53.Mut = as.factor(P53.Mut))
  corr <- mutate(corr, rho = round(rho, 2))
  data.frame(corr)
}

# frequency of p53 mutation in selected patients
# P53_com.sub will then be merged with pan-cancer set to obtain p53 mutation for each patient
Pan.select.p53 <- inner_join(P53.final, Pan.RNA, by = "Patient.ID")

p53.freq <- Pan.select.p53 %>% 
  group_by(P53.Mutation) %>% 
  summarise(n = length(P53.Mutation))


# Now, Name the mutations of interest 
mut.interest <- c(
  "R175H",
  "R248Q",
  "R273H",
  # "R280K",
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
  "G248S",
  "R158L",
  "C175F",
  "H179R",
  "G245S",
  # "R175G",
  "WT"
)

Pan.final <- Pan.select.p53[Pan.select.p53$P53.Mutation %in% mut.interest,]
# this is a good place to clean up the environment
#rm(list=setdiff(ls(), c("Pan.final", "calc.rho", "EMT.Genes", "Pan.RNA")))

# invidually calculating and creating object (learn how to loop better)
FN1.NOX3.n <- calc.rho.n("FN1")
FN1 <- calc.rho("FN1")
VIM <- calc.rho("VIM")
CDH2 <- calc.rho("CDH2")
ZEB1 <- calc.rho("ZEB1")
ZEB2 <- calc.rho("ZEB2")
MMP9 <- calc.rho("MMP9")
MMP2 <- calc.rho("MMP2")
TWIST1 <- calc.rho("TWIST1")
SNAI1 <- calc.rho("SNAI1")
SNAI2 <- calc.rho("SNAI2")
ACTA2 <- calc.rho("ACTA2")
CLDN1 <- calc.rho("CLDN1")
TJP1 <- calc.rho("TJP1")
CRB1 <- calc.rho("CRB1")
KRT18 <- calc.rho("KRT18")
DSP <- calc.rho("DSP")
PKP1 <- calc.rho("PKP1")
CDH1 <- calc.rho("CDH1")
NOX3 <- calc.rho("NOX3")
### recently added list

COL1A1 <- calc.rho("COL1A1") # collagen type 1 alpha 1 chain, valid
ITGA5 <- calc.rho("ITGA5") # alpha 5 subunit of integrin, valid
TGFB1 <- calc.rho("TGFB1")
TGFB2 <- calc.rho("TGFB2")
TGFB3 <- calc.rho("TGFB3")
MMP3 <- calc.rho("MMP3")
WNT5A <- calc.rho("WNT5A")
WNT5B <- calc.rho("WNT5B")
COL1A2 <- calc.rho("COL1A2") # collagen type I alpha 2 chain
COL3A1 <- calc.rho("COL3A1")
COL5A2 <- calc.rho("COL5A2")
ILK <- calc.rho("ILK") # integrin linked kinase
CXCR4 <- calc.rho("CXCR4")
CXCL12 <- calc.rho("CXCL12") # chemokine, also known as SDF1
ITGAV <- calc.rho("ITGAV")
TGFBR1 <- calc.rho("TGFBR1")
TJP1 <- calc.rho("TJP1") # zonuli1
CRB1 <- calc.rho("CRB1") # cell polarity complex crumbs1, valid
KRT18 <- calc.rho("KRT18") # keratin 18, valid 

# rename the rho in each dataframe to its respective gene_NOX3
FN1 <- rename(FN1,FN1 = rho)
VIM <- rename(VIM,VIM = rho)
CDH2 <- rename(CDH2,CDH2 = rho)
ZEB1 <- rename(ZEB1,ZEB1 = rho)
ZEB2 <- rename(ZEB2,ZEB2 = rho)
MMP9 <- rename(MMP9,MMP9 = rho)
MMP2 <- rename(MMP2,MMP2 = rho)
TWIST1 <- rename(TWIST1,TWIST1 = rho)
SNAI1 <- rename(SNAI1,SNAI1 = rho)
SNAI2 <- rename(SNAI2,SNAI2 = rho)
ACTA2 <- rename(ACTA2,ACTA2 = rho)
CLDN1 <- rename(CLDN1,CLDN1 = rho)
TJP1 <- rename(TJP1,TJP1 = rho)
CRB1 <- rename(CRB1,CRB1 = rho)
KRT18 <- rename(KRT18,KRT18 = rho)
DSP <- rename(DSP,DSP = rho)
PKP1 <- rename(PKP1,PKP1 = rho)
CDH1 <- rename(CDH1,CDH1 = rho)
NOX3 <- rename(NOX3,NOX3 = rho)

# newly added

COL1A1 <- rename(COL1A1, COL1A1=rho) # collagen type 1 alpha 1 chain, valid
ITGA5 <- rename(ITGA5, ITGA5=rho) # alpha 5 subunit of integrin, valid
TGFB1 <- rename(TGFB1, TGFB1=rho)
TGFB2 <- rename(TGFB2, TGFB2=rho)
TGFB3 <- rename(TGFB3, TGFB3=rho)
MMP3 <- rename(MMP3, MMP3=rho)
WNT5A <- rename(WNT5A, WNT5A=rho)
WNT5B <- rename(WNT5B, WNT5B=rho)
COL1A2 <- rename(COL1A2, COL1A2=rho) # collagen type I alpha 2 chain
COL3A1 <- rename(COL3A1, COL3A1=rho)
COL5A2 <- rename(COL5A2, COL5A2=rho)
ILK <- rename(ILK, ILK=rho) # integrin linked kinase
CXCR4 <- rename(CXCR4, CXCR4=rho)
CXCL12 <- rename(CXCL12, CXCL12=rho) # chemokine, also known as SDF1
ITGAV <- rename(ITGAV, ITGAV=rho)
TGFBR1 <- rename(TGFBR1, TGFBR1=rho)

# combine to form one df, using inner_join to ensure mut corresponds
rho.comb <- inner_join(FN1, VIM, by= "P53.Mut") %>%
  # inner_join(., CDH2, by = "P53.Mut") %>%
  # inner_join(., ZEB1, by = "P53.Mut") %>%
  # inner_join(., ZEB2, by = "P53.Mut") %>%
  inner_join(., MMP9, by = "P53.Mut") %>%
  inner_join(., MMP2, by = "P53.Mut") %>%
  inner_join(., TWIST1, by = "P53.Mut") %>%
  inner_join(., SNAI1, by = "P53.Mut") %>%
  inner_join(., SNAI2, by = "P53.Mut") %>%
  inner_join(., ACTA2, by = "P53.Mut") %>%
  inner_join(., CLDN1, by = "P53.Mut") %>%
  inner_join(., TJP1, by = "P53.Mut") %>%
  inner_join(., CRB1, by = "P53.Mut") %>%
  inner_join(., KRT18, by = "P53.Mut") %>%
  inner_join(., DSP, by = "P53.Mut") %>%
  inner_join(., PKP1, by = "P53.Mut") %>%
  inner_join(., CDH1, by = "P53.Mut") %>% 
  ### recently added list
  
  inner_join(., COL1A1, by = "P53.Mut") %>% # collagen type 1 alpha 1 chain, valid
  inner_join(., ITGA5, by = "P53.Mut") %>% # alpha 5 subunit of integrin, valid
  inner_join(., TGFB1, by = "P53.Mut") %>%
  inner_join(., TGFB2, by = "P53.Mut") %>%
  inner_join(., TGFB3, by = "P53.Mut") %>%
  inner_join(., MMP3, by = "P53.Mut") %>%
  inner_join(., COL1A2, by = "P53.Mut") %>% # collagen type I alpha 2 chain
  inner_join(., COL3A1, by = "P53.Mut") %>%
  inner_join(., COL5A2, by = "P53.Mut") %>%
  inner_join(., ILK, by = "P53.Mut") %>% # integrin linked kinase
  inner_join(., CXCR4, by = "P53.Mut") %>%
  inner_join(., CXCL12, by = "P53.Mut") %>% # chemokine, also known as SDF1
  inner_join(., ITGAV, by = "P53.Mut") 
# ggplot geom raster wants long form  
rho.long <- gather(rho.comb, key = "Markers", value = "rho", -P53.Mut)
rho.long$Markers <- as.factor(rho.long$Markers)
# reorder EMT markers by type
# EPITHELIAL: CL1, OCC, E-CAD, PLAK, CRUMB3 CYTOKERATINS, ZO1
# MESENCH: FN1, N-CAD, MMP, TWIST1, SNAI1,2, ACTA2, VIM

rho.long$Markers <- factor(rho.long$Markers, levels = c(
  # MESENCHYMAL
  "ACTA2", 
  "COL1A1", # collagen type 1 alpha 1 chain, valid
  "COL1A2", # collagen type I alpha 2 chain
  "COL3A1",
  "COL5A2",
  "CXCR4",
  "CXCL12", # chemokine, also known as SDF1
  "FN1",
  "ILK", # integrin linked kinase
  "ITGA5", # alpha 5 subunit of integrin, valid
  "ITGAV",
  "MMP2",
  "MMP3",
  "MMP9",
  "SNAI1",
  "SNAI2",
  "TGFB1",
  "TGFB2",
  "TGFB3",
  "TWIST1",
  "VIM",
  
  # EPITHELIAL
  "CDH1",
  "CLDN1",
  "CRB1",
  "DSP",
  "KRT18", # keratin
  "PKP1",
  "TJP1"
)
)

# reorder by numeric then put WT first
# fist extract and order the mutations by gene location
p53.factors <- rho.long %>% 
  group_by(P53.Mut) %>% 
  summarize(l = length(rho))

p53.factors$position.gene <- str_sub(p53.factors$P53.Mut, start = 2, end = 4)
p53.factors <-p53.factors[order(p53.factors$position.gene),]
# extract as a vector of characters
gene.position <- as.character(p53.factors$P53.Mut)

# reorder the rho data frame with new level
rho.long$P53.Mut <- factor(rho.long$P53.Mut,
                           levels = c(gene.position))


# now we can create a heat map
NOX3.hm <- ggplot(rho.long) +
  aes(x = Markers, y = P53.Mut, fill = rho) +
  geom_raster() +
  xlab(NULL) +
  ylab(NULL) +
  scale_fill_gradient2(mid = "white", low = "royalblue4", high = "firebrick1",
                       limits = c(-1, 1)) + # limit sets the range of color to use
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 15, vjust=0.5),
        axis.text.y = element_text(size = 15)) +
  geom_text(aes(label = round(rho, 1))) 
# ggtitle("rho values of mRNA expression of various markers and NOX3")
ggsave(NOX3.hm, file = "NOX3_EMT_rho_PanCan_V1.png", dpi=900, width = 11, height = 6.5)
NOX3.hm

### CYBB-EMT correlations pan cancer with p53 ----

# this function will extract CYBB and a gene (fa) from the df = Pan.final and create a matrix 
# of correlation grouped by p53 mutation status
# and also create a column of
calc.rho <- function(fa) { # use this function to calculate the rho between two genes in a long.form data frame
  gene  <- Pan.final[Pan.final$Gene.Symbol %in% fa,] # use genes in ""
  CYBB <- Pan.final[Pan.final$Gene.Symbol %in% "CYBB",]
  co <- inner_join(x= gene, y= CYBB, by = "Patient.ID")
  co2 <- co %>% 
    select(
      Patient.ID,
      P53.Mut = P53.Mutation.x,
      gene = mRNA.Value.x,
      CYBB = mRNA.Value.y,
      -P53.Mutation.y,
      -Gene.Symbol.x,
      -Gene.Symbol.y,
      -Case.Study.x,
      -Case.Study.y
    )
  com <- co2[!duplicated(co2),]
  corr <- com %>% 
    group_by(P53.Mut) %>% 
    summarise(
      rho = cor(CYBB, gene, method = "spearman", use = "complete.obs"))
  corr <- mutate(corr, 
                 P53.Mut = as.factor(P53.Mut))
  corr <- mutate(corr, rho = round(rho, 2))
  data.frame(corr)
}
calc.rho.n <- function(fa) { # use this function to calculate the rho between two genes in a long.form data frame
  gene  <- Pan.final[Pan.final$Gene.Symbol %in% fa,] # use genes in ""
  CYBB <- Pan.final[Pan.final$Gene.Symbol %in% "CYBB",]
  co <- inner_join(x= gene, y= CYBB, by = "Patient.ID")
  co2 <- co %>% 
    select(
      Patient.ID,
      P53.Mut = P53.Mutation.x,
      gene = mRNA.Value.x,
      CYBB = mRNA.Value.y,
      -P53.Mutation.y,
      -Gene.Symbol.x,
      -Gene.Symbol.y,
      -Case.Study.x,
      -Case.Study.y
    )
  com <- co2[!duplicated(co2),]
  corr <- com %>% 
    group_by(P53.Mut) %>% 
    summarise(
      rho = cor(CYBB, gene, method = "spearman", use = "complete.obs"),
      n = length(gene))
  corr <- mutate(corr, 
                 P53.Mut = as.factor(P53.Mut))
  corr <- mutate(corr, rho = round(rho, 2))
  data.frame(corr)
}

# frequency of p53 mutation in selected patients
# P53_com.sub will then be merged with pan-cancer set to obtain p53 mutation for each patient
Pan.select.p53 <- inner_join(P53.final, Pan.RNA, by = "Patient.ID")

p53.freq <- Pan.select.p53 %>% 
  group_by(P53.Mutation) %>% 
  summarise(n = length(P53.Mutation))


# Now, Name the mutations of interest 
mut.interest <- c(
  "R175H",
  "R248Q",
  "R273H",
  # "R280K",
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
  "G248S",
  "R158L",
  "C175F",
  "H179R",
  "G245S",
  # "R175G",
  "WT"
)

Pan.final <- Pan.select.p53[Pan.select.p53$P53.Mutation %in% mut.interest,]
# this is a good place to clean up the environment
#rm(list=setdiff(ls(), c("Pan.final", "calc.rho", "EMT.Genes", "Pan.RNA")))

# invidually calculating and creating object (learn how to loop better)
FN1.CYBB.n <- calc.rho.n("FN1")
FN1 <- calc.rho("FN1")
VIM <- calc.rho("VIM")
CDH2 <- calc.rho("CDH2")
ZEB1 <- calc.rho("ZEB1")
ZEB2 <- calc.rho("ZEB2")
MMP9 <- calc.rho("MMP9")
MMP2 <- calc.rho("MMP2")
TWIST1 <- calc.rho("TWIST1")
SNAI1 <- calc.rho("SNAI1")
SNAI2 <- calc.rho("SNAI2")
ACTA2 <- calc.rho("ACTA2")
CLDN1 <- calc.rho("CLDN1")
TJP1 <- calc.rho("TJP1")
CRB1 <- calc.rho("CRB1")
KRT18 <- calc.rho("KRT18")
DSP <- calc.rho("DSP")
PKP1 <- calc.rho("PKP1")
CDH1 <- calc.rho("CDH1")
CYBB <- calc.rho("CYBB")
### recently added list

COL1A1 <- calc.rho("COL1A1") # collagen type 1 alpha 1 chain, valid
ITGA5 <- calc.rho("ITGA5") # alpha 5 subunit of integrin, valid
TGFB1 <- calc.rho("TGFB1")
TGFB2 <- calc.rho("TGFB2")
TGFB3 <- calc.rho("TGFB3")
MMP3 <- calc.rho("MMP3")
WNT5A <- calc.rho("WNT5A")
WNT5B <- calc.rho("WNT5B")
COL1A2 <- calc.rho("COL1A2") # collagen type I alpha 2 chain
COL3A1 <- calc.rho("COL3A1")
COL5A2 <- calc.rho("COL5A2")
ILK <- calc.rho("ILK") # integrin linked kinase
CXCR4 <- calc.rho("CXCR4")
CXCL12 <- calc.rho("CXCL12") # chemokine, also known as SDF1
ITGAV <- calc.rho("ITGAV")
TGFBR1 <- calc.rho("TGFBR1")
TJP1 <- calc.rho("TJP1") # zonuli1
CRB1 <- calc.rho("CRB1") # cell polarity complex crumbs1, valid
KRT18 <- calc.rho("KRT18") # keratin 18, valid 

# rename the rho in each dataframe to its respective gene_CYBB
FN1 <- rename(FN1,FN1 = rho)
VIM <- rename(VIM,VIM = rho)
CDH2 <- rename(CDH2,CDH2 = rho)
ZEB1 <- rename(ZEB1,ZEB1 = rho)
ZEB2 <- rename(ZEB2,ZEB2 = rho)
MMP9 <- rename(MMP9,MMP9 = rho)
MMP2 <- rename(MMP2,MMP2 = rho)
TWIST1 <- rename(TWIST1,TWIST1 = rho)
SNAI1 <- rename(SNAI1,SNAI1 = rho)
SNAI2 <- rename(SNAI2,SNAI2 = rho)
ACTA2 <- rename(ACTA2,ACTA2 = rho)
CLDN1 <- rename(CLDN1,CLDN1 = rho)
TJP1 <- rename(TJP1,TJP1 = rho)
CRB1 <- rename(CRB1,CRB1 = rho)
KRT18 <- rename(KRT18,KRT18 = rho)
DSP <- rename(DSP,DSP = rho)
PKP1 <- rename(PKP1,PKP1 = rho)
CDH1 <- rename(CDH1,CDH1 = rho)
CYBB <- rename(CYBB,CYBB = rho)

# newly added

COL1A1 <- rename(COL1A1, COL1A1=rho) # collagen type 1 alpha 1 chain, valid
ITGA5 <- rename(ITGA5, ITGA5=rho) # alpha 5 subunit of integrin, valid
TGFB1 <- rename(TGFB1, TGFB1=rho)
TGFB2 <- rename(TGFB2, TGFB2=rho)
TGFB3 <- rename(TGFB3, TGFB3=rho)
MMP3 <- rename(MMP3, MMP3=rho)
WNT5A <- rename(WNT5A, WNT5A=rho)
WNT5B <- rename(WNT5B, WNT5B=rho)
COL1A2 <- rename(COL1A2, COL1A2=rho) # collagen type I alpha 2 chain
COL3A1 <- rename(COL3A1, COL3A1=rho)
COL5A2 <- rename(COL5A2, COL5A2=rho)
ILK <- rename(ILK, ILK=rho) # integrin linked kinase
CXCR4 <- rename(CXCR4, CXCR4=rho)
CXCL12 <- rename(CXCL12, CXCL12=rho) # chemokine, also known as SDF1
ITGAV <- rename(ITGAV, ITGAV=rho)
TGFBR1 <- rename(TGFBR1, TGFBR1=rho)

# combine to form one df, using inner_join to ensure mut corresponds
rho.comb <- inner_join(FN1, VIM, by= "P53.Mut") %>%
  # inner_join(., CDH2, by = "P53.Mut") %>%
  # inner_join(., ZEB1, by = "P53.Mut") %>%
  # inner_join(., ZEB2, by = "P53.Mut") %>%
  inner_join(., MMP9, by = "P53.Mut") %>%
  inner_join(., MMP2, by = "P53.Mut") %>%
  inner_join(., TWIST1, by = "P53.Mut") %>%
  inner_join(., SNAI1, by = "P53.Mut") %>%
  inner_join(., SNAI2, by = "P53.Mut") %>%
  inner_join(., ACTA2, by = "P53.Mut") %>%
  inner_join(., CLDN1, by = "P53.Mut") %>%
  inner_join(., TJP1, by = "P53.Mut") %>%
  inner_join(., CRB1, by = "P53.Mut") %>%
  inner_join(., KRT18, by = "P53.Mut") %>%
  inner_join(., DSP, by = "P53.Mut") %>%
  inner_join(., PKP1, by = "P53.Mut") %>%
  inner_join(., CDH1, by = "P53.Mut") %>% 
  ### recently added list
  
  inner_join(., COL1A1, by = "P53.Mut") %>% # collagen type 1 alpha 1 chain, valid
  inner_join(., ITGA5, by = "P53.Mut") %>% # alpha 5 subunit of integrin, valid
  inner_join(., TGFB1, by = "P53.Mut") %>%
  inner_join(., TGFB2, by = "P53.Mut") %>%
  inner_join(., TGFB3, by = "P53.Mut") %>%
  inner_join(., MMP3, by = "P53.Mut") %>%
  inner_join(., COL1A2, by = "P53.Mut") %>% # collagen type I alpha 2 chain
  inner_join(., COL3A1, by = "P53.Mut") %>%
  inner_join(., COL5A2, by = "P53.Mut") %>%
  inner_join(., ILK, by = "P53.Mut") %>% # integrin linked kinase
  inner_join(., CXCR4, by = "P53.Mut") %>%
  inner_join(., CXCL12, by = "P53.Mut") %>% # chemokine, also known as SDF1
  inner_join(., ITGAV, by = "P53.Mut") 
# ggplot geom raster wants long form  
rho.long <- gather(rho.comb, key = "Markers", value = "rho", -P53.Mut)
rho.long$Markers <- as.factor(rho.long$Markers)
# reorder EMT markers by type
# EPITHELIAL: CL1, OCC, E-CAD, PLAK, CRUMB3 CYTOKERATINS, ZO1
# MESENCH: FN1, N-CAD, MMP, TWIST1, SNAI1,2, ACTA2, VIM

rho.long$Markers <- factor(rho.long$Markers, levels = c(
  # MESENCHYMAL
  "ACTA2", 
  "COL1A1", # collagen type 1 alpha 1 chain, valid
  "COL1A2", # collagen type I alpha 2 chain
  "COL3A1",
  "COL5A2",
  "CXCR4",
  "CXCL12", # chemokine, also known as SDF1
  "FN1",
  "ILK", # integrin linked kinase
  "ITGA5", # alpha 5 subunit of integrin, valid
  "ITGAV",
  "MMP2",
  "MMP3",
  "MMP9",
  "SNAI1",
  "SNAI2",
  "TGFB1",
  "TGFB2",
  "TGFB3",
  "TWIST1",
  "VIM",
  
  # EPITHELIAL
  "CDH1",
  "CLDN1",
  "CRB1",
  "DSP",
  "KRT18", # keratin
  "PKP1",
  "TJP1"
)
)

# reorder by numeric then put WT first
# fist extract and order the mutations by gene location
p53.factors <- rho.long %>% 
  group_by(P53.Mut) %>% 
  summarize(l = length(rho))

p53.factors$position.gene <- str_sub(p53.factors$P53.Mut, start = 2, end = 4)
p53.factors <-p53.factors[order(p53.factors$position.gene),]
# extract as a vector of characters
gene.position <- as.character(p53.factors$P53.Mut)

# reorder the rho data frame with new level
rho.long$P53.Mut <- factor(rho.long$P53.Mut,
                           levels = c(gene.position))


# now we can create a heat map
CYBB.hm <- ggplot(rho.long) +
  aes(x = Markers, y = P53.Mut, fill = rho) +
  geom_raster() +
  xlab(NULL) +
  ylab(NULL) +
  scale_fill_gradient2(mid = "white", low = "royalblue4", high = "firebrick1",
                       limits = c(-1, 1)) + # limit sets the range of color to use
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 15, vjust=0.5),
        axis.text.y = element_text(size = 15)) +
  geom_text(aes(label = round(rho, 1))) 
# ggtitle("rho values of mRNA expression of various markers and CYBB")
ggsave(CYBB.hm, file = "CYBB_EMT_rho_PanCan_V1.png", dpi=900, width = 11, height = 6.5)
CYBB.hm



### NOX5-EMT correlations pan cancer with p53 ----

# this function will extract NOX5 and a gene (fa) from the df = Pan.final and create a matrix 
# of correlation grouped by p53 mutation status
# and also create a column of
calc.rho <- function(fa) { # use this function to calculate the rho between two genes in a long.form data frame
  gene  <- Pan.final[Pan.final$Gene.Symbol %in% fa,] # use genes in ""
  NOX5 <- Pan.final[Pan.final$Gene.Symbol %in% "NOX5",]
  co <- inner_join(x= gene, y= NOX5, by = "Patient.ID")
  co2 <- co %>% 
    select(
      Patient.ID,
      P53.Mut = P53.Mutation.x,
      gene = mRNA.Value.x,
      NOX5 = mRNA.Value.y,
      -P53.Mutation.y,
      -Gene.Symbol.x,
      -Gene.Symbol.y,
      -Case.Study.x,
      -Case.Study.y
    )
  com <- co2[!duplicated(co2),]
  corr <- com %>% 
    group_by(P53.Mut) %>% 
    summarise(
      rho = cor(NOX5, gene, method = "spearman", use = "complete.obs"))
  corr <- mutate(corr, 
                 P53.Mut = as.factor(P53.Mut))
  corr <- mutate(corr, rho = round(rho, 2))
  data.frame(corr)
}
calc.rho.n <- function(fa) { # use this function to calculate the rho between two genes in a long.form data frame
  gene  <- Pan.final[Pan.final$Gene.Symbol %in% fa,] # use genes in ""
  NOX5 <- Pan.final[Pan.final$Gene.Symbol %in% "NOX5",]
  co <- inner_join(x= gene, y= NOX5, by = "Patient.ID")
  co2 <- co %>% 
    select(
      Patient.ID,
      P53.Mut = P53.Mutation.x,
      gene = mRNA.Value.x,
      NOX5 = mRNA.Value.y,
      -P53.Mutation.y,
      -Gene.Symbol.x,
      -Gene.Symbol.y,
      -Case.Study.x,
      -Case.Study.y
    )
  com <- co2[!duplicated(co2),]
  corr <- com %>% 
    group_by(P53.Mut) %>% 
    summarise(
      rho = cor(NOX5, gene, method = "spearman", use = "complete.obs"),
      n = length(gene))
  corr <- mutate(corr, 
                 P53.Mut = as.factor(P53.Mut))
  corr <- mutate(corr, rho = round(rho, 2))
  data.frame(corr)
}

# frequency of p53 mutation in selected patients
# P53_com.sub will then be merged with pan-cancer set to obtain p53 mutation for each patient
Pan.select.p53 <- inner_join(P53.final, Pan.RNA, by = "Patient.ID")

p53.freq <- Pan.select.p53 %>% 
  group_by(P53.Mutation) %>% 
  summarise(n = length(P53.Mutation))


# Now, Name the mutations of interest 
mut.interest <- c(
  "R175H",
  "R248Q",
  "R273H",
  # "R280K",
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
  "G248S",
  "R158L",
  "C175F",
  "H179R",
  "G245S",
  # "R175G",
  "WT"
)

Pan.final <- Pan.select.p53[Pan.select.p53$P53.Mutation %in% mut.interest,]
# this is a good place to clean up the environment
#rm(list=setdiff(ls(), c("Pan.final", "calc.rho", "EMT.Genes", "Pan.RNA")))

# invidually calculating and creating object (learn how to loop better)
FN1.NOX5.n <- calc.rho.n("FN1")
FN1 <- calc.rho("FN1")
VIM <- calc.rho("VIM")
CDH2 <- calc.rho("CDH2")
ZEB1 <- calc.rho("ZEB1")
ZEB2 <- calc.rho("ZEB2")
MMP9 <- calc.rho("MMP9")
MMP2 <- calc.rho("MMP2")
TWIST1 <- calc.rho("TWIST1")
SNAI1 <- calc.rho("SNAI1")
SNAI2 <- calc.rho("SNAI2")
ACTA2 <- calc.rho("ACTA2")
CLDN1 <- calc.rho("CLDN1")
TJP1 <- calc.rho("TJP1")
CRB1 <- calc.rho("CRB1")
KRT18 <- calc.rho("KRT18")
DSP <- calc.rho("DSP")
PKP1 <- calc.rho("PKP1")
CDH1 <- calc.rho("CDH1")
NOX5 <- calc.rho("NOX5")
### recently added list

COL1A1 <- calc.rho("COL1A1") # collagen type 1 alpha 1 chain, valid
ITGA5 <- calc.rho("ITGA5") # alpha 5 subunit of integrin, valid
TGFB1 <- calc.rho("TGFB1")
TGFB2 <- calc.rho("TGFB2")
TGFB3 <- calc.rho("TGFB3")
MMP3 <- calc.rho("MMP3")
WNT5A <- calc.rho("WNT5A")
WNT5B <- calc.rho("WNT5B")
COL1A2 <- calc.rho("COL1A2") # collagen type I alpha 2 chain
COL3A1 <- calc.rho("COL3A1")
COL5A2 <- calc.rho("COL5A2")
ILK <- calc.rho("ILK") # integrin linked kinase
CXCR4 <- calc.rho("CXCR4")
CXCL12 <- calc.rho("CXCL12") # chemokine, also known as SDF1
ITGAV <- calc.rho("ITGAV")
TGFBR1 <- calc.rho("TGFBR1")
TJP1 <- calc.rho("TJP1") # zonuli1
CRB1 <- calc.rho("CRB1") # cell polarity complex crumbs1, valid
KRT18 <- calc.rho("KRT18") # keratin 18, valid 

# rename the rho in each dataframe to its respective gene_NOX5
FN1 <- rename(FN1,FN1 = rho)
VIM <- rename(VIM,VIM = rho)
CDH2 <- rename(CDH2,CDH2 = rho)
ZEB1 <- rename(ZEB1,ZEB1 = rho)
ZEB2 <- rename(ZEB2,ZEB2 = rho)
MMP9 <- rename(MMP9,MMP9 = rho)
MMP2 <- rename(MMP2,MMP2 = rho)
TWIST1 <- rename(TWIST1,TWIST1 = rho)
SNAI1 <- rename(SNAI1,SNAI1 = rho)
SNAI2 <- rename(SNAI2,SNAI2 = rho)
ACTA2 <- rename(ACTA2,ACTA2 = rho)
CLDN1 <- rename(CLDN1,CLDN1 = rho)
TJP1 <- rename(TJP1,TJP1 = rho)
CRB1 <- rename(CRB1,CRB1 = rho)
KRT18 <- rename(KRT18,KRT18 = rho)
DSP <- rename(DSP,DSP = rho)
PKP1 <- rename(PKP1,PKP1 = rho)
CDH1 <- rename(CDH1,CDH1 = rho)
NOX5 <- rename(NOX5,NOX5 = rho)

# newly added

COL1A1 <- rename(COL1A1, COL1A1=rho) # collagen type 1 alpha 1 chain, valid
ITGA5 <- rename(ITGA5, ITGA5=rho) # alpha 5 subunit of integrin, valid
TGFB1 <- rename(TGFB1, TGFB1=rho)
TGFB2 <- rename(TGFB2, TGFB2=rho)
TGFB3 <- rename(TGFB3, TGFB3=rho)
MMP3 <- rename(MMP3, MMP3=rho)
WNT5A <- rename(WNT5A, WNT5A=rho)
WNT5B <- rename(WNT5B, WNT5B=rho)
COL1A2 <- rename(COL1A2, COL1A2=rho) # collagen type I alpha 2 chain
COL3A1 <- rename(COL3A1, COL3A1=rho)
COL5A2 <- rename(COL5A2, COL5A2=rho)
ILK <- rename(ILK, ILK=rho) # integrin linked kinase
CXCR4 <- rename(CXCR4, CXCR4=rho)
CXCL12 <- rename(CXCL12, CXCL12=rho) # chemokine, also known as SDF1
ITGAV <- rename(ITGAV, ITGAV=rho)
TGFBR1 <- rename(TGFBR1, TGFBR1=rho)

# combine to form one df, using inner_join to ensure mut corresponds
rho.comb <- inner_join(FN1, VIM, by= "P53.Mut") %>%
  # inner_join(., CDH2, by = "P53.Mut") %>%
  # inner_join(., ZEB1, by = "P53.Mut") %>%
  # inner_join(., ZEB2, by = "P53.Mut") %>%
  inner_join(., MMP9, by = "P53.Mut") %>%
  inner_join(., MMP2, by = "P53.Mut") %>%
  inner_join(., TWIST1, by = "P53.Mut") %>%
  inner_join(., SNAI1, by = "P53.Mut") %>%
  inner_join(., SNAI2, by = "P53.Mut") %>%
  inner_join(., ACTA2, by = "P53.Mut") %>%
  inner_join(., CLDN1, by = "P53.Mut") %>%
  inner_join(., TJP1, by = "P53.Mut") %>%
  inner_join(., CRB1, by = "P53.Mut") %>%
  inner_join(., KRT18, by = "P53.Mut") %>%
  inner_join(., DSP, by = "P53.Mut") %>%
  inner_join(., PKP1, by = "P53.Mut") %>%
  inner_join(., CDH1, by = "P53.Mut") %>% 
  ### recently added list
  
  inner_join(., COL1A1, by = "P53.Mut") %>% # collagen type 1 alpha 1 chain, valid
  inner_join(., ITGA5, by = "P53.Mut") %>% # alpha 5 subunit of integrin, valid
  inner_join(., TGFB1, by = "P53.Mut") %>%
  inner_join(., TGFB2, by = "P53.Mut") %>%
  inner_join(., TGFB3, by = "P53.Mut") %>%
  inner_join(., MMP3, by = "P53.Mut") %>%
  inner_join(., COL1A2, by = "P53.Mut") %>% # collagen type I alpha 2 chain
  inner_join(., COL3A1, by = "P53.Mut") %>%
  inner_join(., COL5A2, by = "P53.Mut") %>%
  inner_join(., ILK, by = "P53.Mut") %>% # integrin linked kinase
  inner_join(., CXCR4, by = "P53.Mut") %>%
  inner_join(., CXCL12, by = "P53.Mut") %>% # chemokine, also known as SDF1
  inner_join(., ITGAV, by = "P53.Mut") 
# ggplot geom raster wants long form  
rho.long <- gather(rho.comb, key = "Markers", value = "rho", -P53.Mut)
rho.long$Markers <- as.factor(rho.long$Markers)
# reorder EMT markers by type
# EPITHELIAL: CL1, OCC, E-CAD, PLAK, CRUMB3 CYTOKERATINS, ZO1
# MESENCH: FN1, N-CAD, MMP, TWIST1, SNAI1,2, ACTA2, VIM

rho.long$Markers <- factor(rho.long$Markers, levels = c(
  # MESENCHYMAL
  "ACTA2", 
  "COL1A1", # collagen type 1 alpha 1 chain, valid
  "COL1A2", # collagen type I alpha 2 chain
  "COL3A1",
  "COL5A2",
  "CXCR4",
  "CXCL12", # chemokine, also known as SDF1
  "FN1",
  "ILK", # integrin linked kinase
  "ITGA5", # alpha 5 subunit of integrin, valid
  "ITGAV",
  "MMP2",
  "MMP3",
  "MMP9",
  "SNAI1",
  "SNAI2",
  "TGFB1",
  "TGFB2",
  "TGFB3",
  "TWIST1",
  "VIM",
  
  # EPITHELIAL
  "CDH1",
  "CLDN1",
  "CRB1",
  "DSP",
  "KRT18", # keratin
  "PKP1",
  "TJP1"
)
)

# reorder by numeric then put WT first
# fist extract and order the mutations by gene location
p53.factors <- rho.long %>% 
  group_by(P53.Mut) %>% 
  summarize(l = length(rho))

p53.factors$position.gene <- str_sub(p53.factors$P53.Mut, start = 2, end = 4)
p53.factors <-p53.factors[order(p53.factors$position.gene),]
# extract as a vector of characters
gene.position <- as.character(p53.factors$P53.Mut)

# reorder the rho data frame with new level
rho.long$P53.Mut <- factor(rho.long$P53.Mut,
                           levels = c(gene.position))


# now we can create a heat map
NOX5.hm <- ggplot(rho.long) +
  aes(x = Markers, y = P53.Mut, fill = rho) +
  geom_raster() +
  xlab(NULL) +
  ylab(NULL) +
  scale_fill_gradient2(mid = "white", low = "royalblue4", high = "firebrick1",
                       limits = c(-1, 1)) + # limit sets the range of color to use
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 15, vjust=0.5),
        axis.text.y = element_text(size = 15)) +
  geom_text(aes(label = round(rho, 1))) 
# ggtitle("rho values of mRNA expression of various markers and NOX5")
ggsave(NOX5.hm, file = "NOX5_EMT_rho_PanCan_V1.png", dpi=900, width = 11, height = 6.5)
NOX5.hm



### CYBA-EMT correlations pan cancer with p53 ----

# this function will extract CYBA and a gene (fa) from the df = Pan.final and create a matrix 
# of correlation grouped by p53 mutation status
# and also create a column of
calc.rho <- function(fa) { # use this function to calculate the rho between two genes in a long.form data frame
  gene  <- Pan.final[Pan.final$Gene.Symbol %in% fa,] # use genes in ""
  CYBA <- Pan.final[Pan.final$Gene.Symbol %in% "CYBA",]
  co <- inner_join(x= gene, y= CYBA, by = "Patient.ID")
  co2 <- co %>% 
    select(
      Patient.ID,
      P53.Mut = P53.Mutation.x,
      gene = mRNA.Value.x,
      CYBA = mRNA.Value.y,
      -P53.Mutation.y,
      -Gene.Symbol.x,
      -Gene.Symbol.y,
      -Case.Study.x,
      -Case.Study.y
    )
  com <- co2[!duplicated(co2),]
  corr <- com %>% 
    group_by(P53.Mut) %>% 
    summarise(
      rho = cor(CYBA, gene, method = "spearman", use = "complete.obs"))
  corr <- mutate(corr, 
                 P53.Mut = as.factor(P53.Mut))
  corr <- mutate(corr, rho = round(rho, 2))
  data.frame(corr)
}
calc.rho.n <- function(fa) { # use this function to calculate the rho between two genes in a long.form data frame
  gene  <- Pan.final[Pan.final$Gene.Symbol %in% fa,] # use genes in ""
  CYBA <- Pan.final[Pan.final$Gene.Symbol %in% "CYBA",]
  co <- inner_join(x= gene, y= CYBA, by = "Patient.ID")
  co2 <- co %>% 
    select(
      Patient.ID,
      P53.Mut = P53.Mutation.x,
      gene = mRNA.Value.x,
      CYBA = mRNA.Value.y,
      -P53.Mutation.y,
      -Gene.Symbol.x,
      -Gene.Symbol.y,
      -Case.Study.x,
      -Case.Study.y
    )
  com <- co2[!duplicated(co2),]
  corr <- com %>% 
    group_by(P53.Mut) %>% 
    summarise(
      rho = cor(CYBA, gene, method = "spearman", use = "complete.obs"),
      n = length(gene))
  corr <- mutate(corr, 
                 P53.Mut = as.factor(P53.Mut))
  corr <- mutate(corr, rho = round(rho, 2))
  data.frame(corr)
}

# frequency of p53 mutation in selected patients
# P53_com.sub will then be merged with pan-cancer set to obtain p53 mutation for each patient
Pan.select.p53 <- inner_join(P53.final, Pan.RNA, by = "Patient.ID")

p53.freq <- Pan.select.p53 %>% 
  group_by(P53.Mutation) %>% 
  summarise(n = length(P53.Mutation))


# Now, Name the mutations of interest 
mut.interest <- c(
  "R175H",
  "R248Q",
  "R273H",
  # "R280K",
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
  "G248S",
  "R158L",
  "C175F",
  "H179R",
  "G245S",
  # "R175G",
  "WT"
)

Pan.final <- Pan.select.p53[Pan.select.p53$P53.Mutation %in% mut.interest,]
# this is a good place to clean up the environment
#rm(list=setdiff(ls(), c("Pan.final", "calc.rho", "EMT.Genes", "Pan.RNA")))

# invidually calculating and creating object (learn how to loop better)
FN1.CYBA.n <- calc.rho.n("FN1")
FN1 <- calc.rho("FN1")
VIM <- calc.rho("VIM")
CDH2 <- calc.rho("CDH2")
ZEB1 <- calc.rho("ZEB1")
ZEB2 <- calc.rho("ZEB2")
MMP9 <- calc.rho("MMP9")
MMP2 <- calc.rho("MMP2")
TWIST1 <- calc.rho("TWIST1")
SNAI1 <- calc.rho("SNAI1")
SNAI2 <- calc.rho("SNAI2")
ACTA2 <- calc.rho("ACTA2")
CLDN1 <- calc.rho("CLDN1")
TJP1 <- calc.rho("TJP1")
CRB1 <- calc.rho("CRB1")
KRT18 <- calc.rho("KRT18")
DSP <- calc.rho("DSP")
PKP1 <- calc.rho("PKP1")
CDH1 <- calc.rho("CDH1")
CYBA <- calc.rho("CYBA")
### recently added list

COL1A1 <- calc.rho("COL1A1") # collagen type 1 alpha 1 chain, valid
ITGA5 <- calc.rho("ITGA5") # alpha 5 subunit of integrin, valid
TGFB1 <- calc.rho("TGFB1")
TGFB2 <- calc.rho("TGFB2")
TGFB3 <- calc.rho("TGFB3")
MMP3 <- calc.rho("MMP3")
WNT5A <- calc.rho("WNT5A")
WNT5B <- calc.rho("WNT5B")
COL1A2 <- calc.rho("COL1A2") # collagen type I alpha 2 chain
COL3A1 <- calc.rho("COL3A1")
COL5A2 <- calc.rho("COL5A2")
ILK <- calc.rho("ILK") # integrin linked kinase
CXCR4 <- calc.rho("CXCR4")
CXCL12 <- calc.rho("CXCL12") # chemokine, also known as SDF1
ITGAV <- calc.rho("ITGAV")
TGFBR1 <- calc.rho("TGFBR1")
TJP1 <- calc.rho("TJP1") # zonuli1
CRB1 <- calc.rho("CRB1") # cell polarity complex crumbs1, valid
KRT18 <- calc.rho("KRT18") # keratin 18, valid 

# rename the rho in each dataframe to its respective gene_CYBA
FN1 <- rename(FN1,FN1 = rho)
VIM <- rename(VIM,VIM = rho)
CDH2 <- rename(CDH2,CDH2 = rho)
ZEB1 <- rename(ZEB1,ZEB1 = rho)
ZEB2 <- rename(ZEB2,ZEB2 = rho)
MMP9 <- rename(MMP9,MMP9 = rho)
MMP2 <- rename(MMP2,MMP2 = rho)
TWIST1 <- rename(TWIST1,TWIST1 = rho)
SNAI1 <- rename(SNAI1,SNAI1 = rho)
SNAI2 <- rename(SNAI2,SNAI2 = rho)
ACTA2 <- rename(ACTA2,ACTA2 = rho)
CLDN1 <- rename(CLDN1,CLDN1 = rho)
TJP1 <- rename(TJP1,TJP1 = rho)
CRB1 <- rename(CRB1,CRB1 = rho)
KRT18 <- rename(KRT18,KRT18 = rho)
DSP <- rename(DSP,DSP = rho)
PKP1 <- rename(PKP1,PKP1 = rho)
CDH1 <- rename(CDH1,CDH1 = rho)
CYBA <- rename(CYBA,CYBA = rho)

# newly added

COL1A1 <- rename(COL1A1, COL1A1=rho) # collagen type 1 alpha 1 chain, valid
ITGA5 <- rename(ITGA5, ITGA5=rho) # alpha 5 subunit of integrin, valid
TGFB1 <- rename(TGFB1, TGFB1=rho)
TGFB2 <- rename(TGFB2, TGFB2=rho)
TGFB3 <- rename(TGFB3, TGFB3=rho)
MMP3 <- rename(MMP3, MMP3=rho)
WNT5A <- rename(WNT5A, WNT5A=rho)
WNT5B <- rename(WNT5B, WNT5B=rho)
COL1A2 <- rename(COL1A2, COL1A2=rho) # collagen type I alpha 2 chain
COL3A1 <- rename(COL3A1, COL3A1=rho)
COL5A2 <- rename(COL5A2, COL5A2=rho)
ILK <- rename(ILK, ILK=rho) # integrin linked kinase
CXCR4 <- rename(CXCR4, CXCR4=rho)
CXCL12 <- rename(CXCL12, CXCL12=rho) # chemokine, also known as SDF1
ITGAV <- rename(ITGAV, ITGAV=rho)
TGFBR1 <- rename(TGFBR1, TGFBR1=rho)

# combine to form one df, using inner_join to ensure mut corresponds
rho.comb <- inner_join(FN1, VIM, by= "P53.Mut") %>%
  # inner_join(., CDH2, by = "P53.Mut") %>%
  # inner_join(., ZEB1, by = "P53.Mut") %>%
  # inner_join(., ZEB2, by = "P53.Mut") %>%
  inner_join(., MMP9, by = "P53.Mut") %>%
  inner_join(., MMP2, by = "P53.Mut") %>%
  inner_join(., TWIST1, by = "P53.Mut") %>%
  inner_join(., SNAI1, by = "P53.Mut") %>%
  inner_join(., SNAI2, by = "P53.Mut") %>%
  inner_join(., ACTA2, by = "P53.Mut") %>%
  inner_join(., CLDN1, by = "P53.Mut") %>%
  inner_join(., TJP1, by = "P53.Mut") %>%
  inner_join(., CRB1, by = "P53.Mut") %>%
  inner_join(., KRT18, by = "P53.Mut") %>%
  inner_join(., DSP, by = "P53.Mut") %>%
  inner_join(., PKP1, by = "P53.Mut") %>%
  inner_join(., CDH1, by = "P53.Mut") %>% 
  ### recently added list
  
  inner_join(., COL1A1, by = "P53.Mut") %>% # collagen type 1 alpha 1 chain, valid
  inner_join(., ITGA5, by = "P53.Mut") %>% # alpha 5 subunit of integrin, valid
  inner_join(., TGFB1, by = "P53.Mut") %>%
  inner_join(., TGFB2, by = "P53.Mut") %>%
  inner_join(., TGFB3, by = "P53.Mut") %>%
  inner_join(., MMP3, by = "P53.Mut") %>%
  inner_join(., COL1A2, by = "P53.Mut") %>% # collagen type I alpha 2 chain
  inner_join(., COL3A1, by = "P53.Mut") %>%
  inner_join(., COL5A2, by = "P53.Mut") %>%
  inner_join(., ILK, by = "P53.Mut") %>% # integrin linked kinase
  inner_join(., CXCR4, by = "P53.Mut") %>%
  inner_join(., CXCL12, by = "P53.Mut") %>% # chemokine, also known as SDF1
  inner_join(., ITGAV, by = "P53.Mut") 
# ggplot geom raster wants long form  
rho.long <- gather(rho.comb, key = "Markers", value = "rho", -P53.Mut)
rho.long$Markers <- as.factor(rho.long$Markers)
# reorder EMT markers by type
# EPITHELIAL: CL1, OCC, E-CAD, PLAK, CRUMB3 CYTOKERATINS, ZO1
# MESENCH: FN1, N-CAD, MMP, TWIST1, SNAI1,2, ACTA2, VIM

rho.long$Markers <- factor(rho.long$Markers, levels = c(
  # MESENCHYMAL
  "ACTA2", 
  "COL1A1", # collagen type 1 alpha 1 chain, valid
  "COL1A2", # collagen type I alpha 2 chain
  "COL3A1",
  "COL5A2",
  "CXCR4",
  "CXCL12", # chemokine, also known as SDF1
  "FN1",
  "ILK", # integrin linked kinase
  "ITGA5", # alpha 5 subunit of integrin, valid
  "ITGAV",
  "MMP2",
  "MMP3",
  "MMP9",
  "SNAI1",
  "SNAI2",
  "TGFB1",
  "TGFB2",
  "TGFB3",
  "TWIST1",
  "VIM",
  
  # EPITHELIAL
  "CDH1",
  "CLDN1",
  "CRB1",
  "DSP",
  "KRT18", # keratin
  "PKP1",
  "TJP1"
)
)


# reorder by numeric then put WT first
# fist extract and order the mutations by gene location
p53.factors <- rho.long %>% 
  group_by(P53.Mut) %>% 
  summarize(l = length(rho))

p53.factors$position.gene <- str_sub(p53.factors$P53.Mut, start = 2, end = 4)
p53.factors <-p53.factors[order(p53.factors$position.gene),]
# extract as a vector of characters
gene.position <- as.character(p53.factors$P53.Mut)

# reorder the rho data frame with new level
rho.long$P53.Mut <- factor(rho.long$P53.Mut,
                           levels = c(gene.position))


# now we can create a heat map
CYBA.hm <- ggplot(rho.long) +
  aes(x = Markers, y = P53.Mut, fill = rho) +
  geom_raster() +
  xlab(NULL) +
  ylab(NULL) +
  scale_fill_gradient2(mid = "white", low = "royalblue4", high = "firebrick1",
                       limits = c(-1, 1)) + # limit sets the range of color to use
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 15, vjust=0.5),
        axis.text.y = element_text(size = 15)) +
  geom_text(aes(label = round(rho, 1))) 
# ggtitle("rho values of mRNA expression of various markers and CYBA")
ggsave(CYBA.hm, file = "CYBA_EMT_rho_PanCan_V1.png", dpi=900, width = 11, height = 6.5)
CYBA.hm



### Pan-Cancer NOX4 mRNA expression by p53 mutation status ----
NOX4.mRNA <- Pan.select.p53[Pan.select.p53$Gene.Symbol %in% "NOX4", ]
NOX4.p53 <- inner_join(NOX4.mRNA, P53.final, by = "Patient.ID")
NOX4.u <- NOX4.p53[!duplicated(NOX4.p53),]

mut.interest <- c(
  "R175H",
  "R248Q",
  "R273H",
  "R280W",
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
  "G248S",
  "R158L",
  "H179R",
  "G245S",
  "R175G",
  "R280W",
  "WT"
)  

NOX4.final <- NOX4.u[NOX4.u$P53.Mutation.y %in% mut.interest,]


# remove duplicates
NOX4.final <- select(NOX4.final, Patient.ID, mut = P53.Mutation.y, mRNA.Value)
NOX4.final <- NOX4.final[!duplicated(NOX4.final),]

# reorder by median
# calculate median of each
p53cohort <- NOX4.final %>% 
  group_by(mut) %>% 
  summarize(median = median(mRNA.Value),
            n = length(mut))
p53cohort.ordered <- p53cohort[order(p53cohort$median),]
ord <- as.character(p53cohort.ordered$mut)
ord

# we also want WT to be at the end, using the moveMe Function
moveMe <- function(invec, movecommand) {
  movecommand <- lapply(strsplit(strsplit(movecommand, ";")[[1]], ",|\\s+"), 
                        function(x) x[x != ""])
  movelist <- lapply(movecommand, function(x) {
    Where <- x[which(x %in% c("before", "after", "first", "last")):length(x)]
    ToMove <- setdiff(x, Where)
    list(ToMove, Where)
  })
  myVec <- invec
  for (i in seq_along(movelist)) {
    temp <- setdiff(myVec, movelist[[i]][[1]])
    A <- movelist[[i]][[2]][1]
    if (A %in% c("before", "after")) {
      ba <- movelist[[i]][[2]][2]
      if (A == "before") {
        after <- match(ba, temp)-1
      } else if (A == "after") {
        after <- match(ba, temp)
      }    
    } else if (A == "first") {
      after <- 0
    } else if (A == "last") {
      after <- length(myVec)
    }
    myVec <- append(temp, values = movelist[[i]][[1]], after = after)
  }
  myVec
} #credit to mrdwab on GitHub

ord2 <- moveMe(ord, "WT first")
ord2
# we can manually reorder them by using levels 
NOX4.final$mut <- factor(NOX4.final$mut,
                         levels = c(
                           ord2
                         ))


NOX4.pan.mRNA.v <- ggplot(NOX4.final) +
  aes(x = mut, y = log(mRNA.Value)) +
  geom_sina(alpha = 1) +
  geom_boxplot(alpha = 0.5) +
  geom_boxplot(outlier.shape = NA, lwd = 1.5) +
  geom_sina(color = "royalblue", alpha = 0.8) +
  ylab("Relative Fold Change in NOX4 mRNA") +
  xlab(NULL) +
  ggtitle(cancer) +
  theme(
    axis.text.x = element_text(
      angle = 30,
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
    panel.background = element_rect(fill = "white", color = "black", size = 4),
    axis.ticks = element_line(size = 2)
  ) +
  geom_hline(yintercept = median(log(NOX4.final$mRNA.Value[NOX4.final$mut == "WT"])),
             color = "red",
             size = 1.5) 
NOX4.pan.mRNA.v

ggsave(NOX4.pan.mRNA.v, file = "PanCan_NOX4_mRNA_by_p53_Mutations_violin.pdf", dpi=900, width = 13, height = 6.5)

### Patient-cohort characteristic summary ----
# https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tissue-source-site-codes 
# we will extract tissue type using the TSS (TCGA-XX) code (where XX denotes tissue)
# I previously copied the code table into a cvs file

library(readr)
TSS <- read_csv("Tissue_source_site_code.csv")
# rename first column as Code
colnames(TSS)[1] <- "Code"
# now extract TSS from data set
NOX4.final$Code <- str_sub(NOX4.final$Patient.ID, start = 6, end = 7)

# R recognizes NA as a empty value but in this case NA is Uterine Carcinosarcoma(DUKE)
# so we change both code NA and empty NA to something else
NOX4.final$Code[NOX4.final$Code == 'NA'] <- "DUKE"
TSS$Code[is.na(TSS$Code)] <- "DUKE"

# now assign cancer type to our data set
TSS.final <- left_join(NOX4.final, TSS, by = "Code")

summary.of.cohort <- TSS.final %>% 
  group_by(Study) %>% 
  summarise(n = length(Study))

write.csv(summary.of.cohort, file = "pancan_dataset_composition.csv")

### Example NOX4-EMT scatter plot faceted by p53 mutations ---- 

# P53_com.sub will then be merged with pan-cancer set to obtain p53 mutation for each patient
Pan.select.p53 <- inner_join(P53.final, Pan.RNA, by = "Patient.ID")

p53.freq <- Pan.select.p53 %>% 
  group_by(P53.Mutation) %>% 
  summarise(n = length(P53.Mutation))

# Now, Name the mutations of interest 
mut.interest <- c(
  "R175H",
  "R248Q",
  "R273H",
  # "R280K",
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
  "G248S",
  "R158L",
  "C175F",
  "H179R",
  "G245S",
  # "R175G",
  "WT"
)

N4 <- Pan.final[Pan.final$Gene.Symbol == "NOX4",]
N4 <- select(N4, Patient.ID, P53.Mutation, NOX4 = mRNA.Value)
F1 <- Pan.final[Pan.final$Gene.Symbol == "FN1",]
F1 <- select(F1, Patient.ID, P53.Mutation, FN1 = mRNA.Value)

com <- inner_join(N4, F1, by = "Patient.ID")
com <- select(com, Patient.ID, mut = P53.Mutation.x, NOX4, FN1)
com <- com[!duplicated(com),]
com$mut <- as.character(com$mut)

#calculate rho by p53 mut
rho <- com %>% 
  group_by(mut) %>% 
  summarise(rho = (cor(NOX4, FN1, method = "spearman")))
rho$rho <- round(rho$rho, 2) # round by 2 decimals
# calculate n for each group
counts <- com %>% 
  group_by(mut) %>% 
  summarise(n = length(FN1))


rho.demo <- ggplot(com) + aes(x= log(NOX4), y= log(FN1)) +
  geom_point()+
  geom_smooth(method = "loess") +
  xlab("NOX4 Expression Fold Change")+
  ylab("FN1 Expression Fold Change") +
  facet_wrap(~mut, nrow = 4) +
  geom_text(data = rho, size = 3, # this adds the rho value by p53 mut
            aes(label = (paste("rho", rho, sep = "="))), 
            x = 7, y = 6 ) +
  geom_text(data = counts, size = 3, # this adds the rho value by p53 mut
            aes(label = (paste("n", n, sep = "="))), 
            x = 7, y = 5 ) +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.background = element_rect(fill = "white", color = "black", size = 1),
    strip.text.x = element_text(size=9),
    strip.background = element_rect(colour="black", fill= "white", size = 1)) 

rho.demo

ggsave(rho.demo, file = "rho_demo.png", dpi=900, width = 10, height = 6.5)

length(com$Patient.ID)


### Pan-Cancer Survival download ----
BLCA <- dl.surv("BLCA") # bladder car
BRCA <- dl.surv("BRCA") # breast car
CESC <- dl.surv("CESC")
COAD <- dl.surv("COAD") # colon adeno
COADREAD <- dl.surv("COADREAD") # colorectal adenocarcinoma
ESCA <- dl.surv("ESCA")

gc() # garbage collection

GBM  <- dl.surv("GBM")  # gliblastoma multiform
HNSC <- dl.surv("HNSC") # head and neck sq carc
KIRC <- dl.surv("KIRC") # kidney renal clear cell car
KIRP <- dl.surv("KIRP") # kidney renal pap cell car

gc()

LGG  <- dl.surv("LGG")
LIHC <- dl.surv("LIHC") # liver hep car
LUSC <- dl.surv("LUSC") # lung squ cell car
LUAD <- dl.surv("LUAD") # lung adeno
MESO <- dl.surv("MESO")
OV   <- dl.surv("OV")   # ovarian ser cyst
PRAD <- dl.surv("PRAD") # prostate adeno
PAAD <- dl.surv("PAAD")

gc() 

PCPG <- dl.surv("PCPG")
READ <- dl.surv("READ")
SARC <- dl.surv("SARC") # sarcoma
#SKCM <- dl.surv("SKCM")
STAD <- dl.surv("STAD") # stomach adeno
THCA <- dl.surv("THCA") # thyroid carcinoma
gc()
THCA <- dl.surv("THCA")
UCEC <- dl.surv("UCEC") # uterine corpus endometrial carcinoma
UCS  <- dl.surv("UCS")

gc()

#rbind for PanCan
Pan.surv <- rbind( BLCA, # bladder car
                   BRCA, # breast car
                   CESC, # cervical
                   COAD, # colon adeno
                   ESCA, # esophogeal 
                   GBM,  # gliblastoma multiform
                   HNSC, # head and neck sq carc
                   KIRC, # kidney renal clear cell car
                   KIRP, # kidney renal pap cell car
                   LGG,# brain lower grade glioma
                   LIHC, # liver hep car
                   LUSC, # lung squ cell car
                   LUAD, # lung adeno
                   MESO,# mesothelioma
                   OV,   # ovarian ser cyst
                   PRAD, # prostate adeno
                   PAAD, # pancreatic
                   # PCPG, # pheochyromocytoma
                   READ, # rectum
                   SARC, # sarcoma
                   # SKCM,
                   STAD, # stomach adeno
                   THCA, # thyroid carcinoma
                   UCEC, # uterine corpus endometrial carcinoma
                   UCS
)

### DEPRECIATED ### Patient Survival by TP53 designated with high or low NOX4 exp ----
# Pan.surv <- Pan.surv[!duplicated(Pan.surv),]
# muts <- select(NOX4.final, Patient.ID, mut)
# muts$Patient.ID <- str_sub(muts$Patient.ID, start = 1, end = 12)
# dup <- muts[(duplicated(muts)),]
# 
# muts <- muts[!duplicated(muts),]
# # we note that there are duplicates b/c there are patients with >1 tumor samples
# # we will remove the duplicated one, since these patients' samples are all tumors
# # and that (data not shown) the tumors carry same mutations
# surv <- inner_join(muts, Pan.surv, by = "Patient.ID")
# head(surv)
# # now we assign high or low status based on mutation status
# # first create an empty column called assign
# surv$assign <- NA
# # Higher group: c("R175H", "V157F", "R273H", R157L", "G245D", "H193R", "H179R")
# surv$assign[surv$mut == "R175H"] <- "Mutant" # find cells in assign col where mut matches this, and assign it as Mutant
# surv$assign[surv$mut == "V157F"] <- "Mutant"
# surv$assign[surv$mut == "R273H"] <- "Mutant"
# surv$assign[surv$mut == "R157L"] <- "Mutant"
# surv$assign[surv$mut == "G245D"] <- "Mutant"
# surv$assign[surv$mut == "H193R"] <- "Mutant"
# surv$assign[surv$mut == "H179R"] <- "Mutant"
# surv$assign[surv$mut == "R158L"] <- "Mutant"
# surv$assign[surv$mut == "R273C"] <- "Low"
# surv$assign[surv$mut == "R249S"] <- "Low"
# surv$assign[surv$mut == "R175G"] <- "Low"
# surv$assign[surv$mut == "R245S"] <- "Low"
# surv$assign[surv$mut == "G245S"] <- "Low"
# surv$assign[surv$mut == "R273L"] <- "Low"
# surv$assign[surv$mut == "R282W"] <- "Low"
# surv$assign[surv$mut == "R248Q"] <- "Low"
# surv$assign[surv$mut == "Y220C"] <- "Average"
# surv$assign[surv$mut == "R248W"] <- "Average"
# surv$assign[surv$mut == "WT"] <- "WT"
# 
# ##### there are na. rewrite this code
# summary(is.na(surv$assign))
# View(surv)
# 
# surv$vital_status <- as.numeric(surv$vital_status)
# # now we plot it using survminer
# require("survival")
# fit <- survfit(Surv(time, vital_status) ~ assign, data = surv)
# 
# ggsurvplot(fit,
#            risk.table = T,
#            conf.int = F,
#            pval=T)
# #survival_groupedby_p53_assigned_nox4_high_low
# 
# ### PanCan Patient survival by NOX4 expression levels
# ### PanCan Patient survival by housekeeping gene expression levels
### Patient Survival by high/low NOX4, mutant or WT-p53 and overall survival ----
Pan.surv <- Pan.surv[!duplicated(Pan.surv),]
muts <- select(NOX4.final, Patient.ID, mut, NOX4 = mRNA.Value) # include nox4
muts$Patient.ID <- str_sub(muts$Patient.ID, start = 1, end = 12)
dup <- muts[(duplicated(muts)),]

muts <- muts[!duplicated(muts),]
# we note that there are duplicates b/c there are patients with >1 tumor samples
# we will remove the duplicated one, since these patients' samples are all tumors
# and that (data not shown) the tumors carry same mutations

## Overall survival ## 
surv <- inner_join(muts, Pan.surv, by = "Patient.ID")
head(surv)

surv.os <- select(surv, status = vital_status, time, NOX4)
# determine optimium cut off as high or low
surv_cutpoint(surv.os, time = "time", event = "status", variable = "NOX4")
surv.os$Assign[surv.os$NOX4 <= 7.3239] <- "low"
surv.os$Assign[surv.os$NOX4 > 7.3239] <- "high"

surv.fit <- survfit(Surv(time, as.numeric(status)) ~ Assign, data = surv.os)
plot <- ggsurvplot(surv.fit, data= surv.os, 
                   risk.table = F,       # show risk table.
                   pval = TRUE,
                   size = 1.5,
                   break.time.by = 1500,
                   fontsize = 30,
                   xlab= "Time in Days",
                   ggtheme = theme_light(),
                   surv.median.line = "hv",
                   title = "Overall Survival by NOX4",
                   font.x = 0,
                   font.y = 0,
                   font.tickslab = 20)

plot
ggsave(print(plot), file = "OS_NOX4.png", dpi = 900)


## mutants only ###
surv <- inner_join(muts, Pan.surv, by = "Patient.ID")
head(surv)

surv <- surv[!surv$mut == "WT",]

surv.os <- select(surv, status = vital_status, time, NOX4)
# determine optimium cut off as high or low
surv_cutpoint(surv.os, time = "time", event = "status", variable = "NOX4")
surv.os$Assign[surv.os$NOX4 <= 5.8501] <- "low"
surv.os$Assign[surv.os$NOX4 > 5.8501] <- "high"

surv.fit <- survfit(Surv(time, as.numeric(status)) ~ Assign, data = surv.os)
plot <- ggsurvplot(surv.fit, data= surv.os, 
                   risk.table = F,       # show risk table.
                   pval = TRUE,
                   size = 1.5,
                   break.time.by = 1500,
                   fontsize = 30,
                   xlab= "Time in Days",
                   ggtheme = theme_light(),
                   surv.median.line = "hv",
                   title = "Mutants only Survival by NOX4",
                   font.x = 0,
                   font.y = 0,
                   font.tickslab = 20)

plot
ggsave(print(plot), file = "Mut_NOX4.png", dpi = 900)

## WTs only ###
surv <- inner_join(muts, Pan.surv, by = "Patient.ID")
head(surv)

surv <- surv[surv$mut == "WT",]

surv.os <- select(surv, status = vital_status, time, NOX4)
# determine optimium cut off as high or low
surv_cutpoint(surv.os, time = "time", event = "status", variable = "NOX4")
surv.os$Assign[surv.os$NOX4 <= 195.5592] <- "low"
surv.os$Assign[surv.os$NOX4 > 195.5592] <- "high"

surv.fit <- survfit(Surv(time, as.numeric(status)) ~ Assign, data = surv.os)
plot <- ggsurvplot(surv.fit, data= surv.os, 
                   risk.table = F,       # show risk table.
                   pval = TRUE,
                   size = 1.5,
                   break.time.by = 1500,
                   fontsize = 30,
                   xlab= "Time in Days",
                   ggtheme = theme_light(),
                   surv.median.line = "hv",
                   title = "WT only Survival by NOX4",
                   font.x = 0,
                   font.y = 0,
                   font.tickslab = 20)

plot
ggsave(print(plot), file = "WT_NOX4.png", dpi = 900)
