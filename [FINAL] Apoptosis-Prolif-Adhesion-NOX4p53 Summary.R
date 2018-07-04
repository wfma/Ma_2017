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
library(readxl)
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
    CNA_SNP = F, # copy number alterations in somatic cells provided by segmented sequecing
    CNA_Seq = F, # copy number alterations provided by NGS sequences
    CNA_CGH = F, # copy number alternations provided by CGH platform
    CNV_SNP = F, # copy number alterantion in germline cells
    # Methylation = T, # methylation provided by array platform
    # RPPA = T, # reverse phase protein array expression
    RNAseq2_Gene_Norm = TRUE,
    # normalized count
    fileSizeLimit = 99999,
    # getUUIDs = T,
    destdir = "FireHose Data",
    forceDownload = T
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

### Chooosing Interested gene expression----
EMT.Genes <- c(
  "FN1",
  "VIM",
  "CDH2",
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
  "MIF", # this and other genes not used are for another project
  "GAPDH",
  "TLR4",# keratin 18, valid 
  # below are genes to screen for in last section
  "BAX", # from structural biology of p53 Joerger et al.x
  "BBC3",
  "PMAIP1",
  "PCNA",
  "GADD45A",
  "PAK1",
  "RRM2B",
  "CDKN1A", #P21
  "SERPINB5",
  "THBS1",
  "ANXA5",
  "ABL1", #senensece marker by sabioscience
  "AKT1",
  "ALDH1A3",
  "E2F1",
  "HRAS",
  "IGFBP3",
  "MYC",
  "CALR",
  "SOD1",
  "SOD2",
  "PRKCD",
  "MAPK14",
  "PTEN",
  "BCL2L11",
  "BID",
  "BAD",
  "BLK",
  "BIK",
  "BCL2L1",
  "BAG1",
  "BAG2",
  "BAG3",
  "BAG4",
  "BAG5",
  "CASP1",
  "CASP2",
  "CASP3",
  "CASP4",
  "CASP5",
  "CASP6",
  "CASP7",
  "CASP8",
  "CASP9",
  "CASP10",
  "CASP11",
  "CASP12",
  "CASP13",
  "CASP14",
  "TP53",
  "ATM",
  "BMI1",
  "CCND1",
  "CCNE1",
  "CDK2",
  "CDK4",
  "CDK6",
  "CDK1",
  "CDKN1A",
  "CDKN2A",
  "CDKN2D",
  "CHEK1",
  "CHEK2",
  "E2F3",
  "ETS1",
  "ETS2",
  "MDM2",
  "RB1",
  "RBL2",
  "ATM", #DNA DAMAGE SENSENCE INITIATORS
  "NBN",
  "TERF2",
  "TERT",
  "TP53BP1",
  "PRKCD",
  "PTPRJ", # aka DEP1; Althubiti et al., 2014 Characterization of novel markers of senescence and their prognostic potential in cancer
  "LAT2", # aka NTAL
  "SLC9A3R1", # EBP50
  "STX4",
  "VAMP3",
  "B2M",# B2MG
  "LANCL1",
  "VPS26A",
  "PLD3",
  "ACHY", # senesence markers from Carnero, 2013
  "PPP1A",
  "CSN2",
  "BRF1",
  "PGM1",
  "IGFBP3",
  "IGFBP7",
  "SERPINE1",
  "MAPK14",
  "MAP2K6",
  "SMURF2",
  # angiogenesis 
  
  "ANG", # angiogensis
  "ANGPT1", # angiopoietin1
  "ANGPT2",
  "ANPEP", # alanyl aminopeptidase
  "TYMP", # thymidine phosphorylase 
  "FGF1", # fibroblast growth factor1
  "FGF2",
  "VEGFD", # vascular endothelial growth factor D
  "FLT1", # fms related tyrosine kinase1 interacts with VEGFRA/b
  "KDR", # encodes one of the two receptors of VEGF
  "NRP2", # neuropilin2
  "PGF", # placental growth factor
  "VEGFA",
  "VEGFB",
  "VEGFC",
  "ADGRB1", # bai1, bines to p53 and induced by wtp53
  "ANGPTL1",# angiopoietin like 4
  "HIF1A",
  #BaseExcisionRepair(BER)
  "APEX1", "APEX2", "CCNO", "LIG3", "MPG", "MUTYH", "NEIL1", "NEIL2", "NEIL3", "NTHL1", "OGG1", "PARP1", "PARP2", "PARP3", "POLB", "SMUG1", "TDG", "UNG", "XRCC1",
  
  #NucleotideExcisionRepair(NER)
  "ATXN3", "BRIP1", "CCNH", "CDK7", "DDB1", "DDB2", "ERCC1", "ERCC2(XPD)", "ERCC3(XPB)", "ERCC4", "ERCC5", "ERCC6", "ERCC8", "LIG1", "MMS19", "PNKP", "POLL", "RAD23A", "RAD23B", "RPA1", "RPA3", "SLK", "XAB2", "XPA", "XPC",
  
  #MismatchRepair(MMR)
  "MLH1", "MLH3", "MSH2", "MSH3", "MSH4", "MSH5", "MSH6", "PMS1", "PMS2", "POLD3", "TREX1",
  
  #Double-StrandBreak(DSB)Repair
  "BRCA1", "BRCA2", "DMC1", "FEN1", "LIG4", "MRE11A", "PRKDC", "RAD21", "RAD50", "RAD51", "RAD51C", "RAD51B(RAD51L1)", "RAD51D(RAD51L3)", "RAD52", "RAD54L", "XRCC2", "XRCC3", "XRCC4", "XRCC5", "XRCC6",
  
  #OtherDNARepairGenes
  "ATM", "ATR", "EXO1", "MGMT(AGT)", "RAD18", "RFC1", "TOP3A", "TOP3B", "XRCC6BP1", # angiogenesis 
  
  "ANG", # angiogensis
  "ANGPT1", # angiopoietin1
  "ANGPT2",
  "ANPEP", # alanyl aminopeptidase
  "TYMP", # thymidine phosphorylase 
  "FGF1", # fibroblast growth factor1
  "FGF2",
  "VEGFD", # vascular endothelial growth factor D
  "FLT1", # fms related tyrosine kinase1 interacts with VEGFRA/b
  "KDR", # encodes one of the two receptors of VEGF
  "NRP2", # neuropilin2
  "PGF", # placental growth factor
  "VEGFA",
  "VEGFB",
  "VEGFC",
  "ADGRB1", # bai1, bines to p53 and induced by wtp53
  "ANGPTL1",# angiopoietin like 4
  "HIF1A",
  #BaseExcisionRepair(BER)
  "APEX1", "APEX2", "CCNO", "LIG3", "MPG", "MUTYH", "NEIL1", "NEIL2", "NEIL3", "NTHL1", "OGG1", "PARP1", "PARP2", "PARP3", "POLB", "SMUG1", "TDG", "UNG", "XRCC1",
  
  #NucleotideExcisionRepair(NER)
  "ATXN3", "BRIP1", "CCNH", "CDK7", "DDB1", "DDB2", "ERCC1", "ERCC2", "ERCC3", "ERCC4", "ERCC5", "ERCC6", "ERCC8", "LIG1", "MMS19", "PNKP", "POLL", "RAD23A", "RAD23B", "RPA1", "RPA3", "SLK", "XAB2", "XPA", "XPC",
  
  #MismatchRepair(MMR)
  "MLH1", "MLH3", "MSH2", "MSH3", "MSH4", "MSH5", "MSH6", "PMS1", "PMS2", "POLD3", "TREX1",
  
  #Double-StrandBreak(DSB)Repair
  "BRCA1", "BRCA2", "DMC1", "FEN1", "LIG4", "MRE11A", "PRKDC", "RAD21", "RAD50", "RAD51", "RAD51C", "RAD51B", "RAD51D", "RAD52", "RAD54L", "XRCC2", "XRCC3", "XRCC4", "XRCC5", "XRCC6",
  
  #OtherDNARepairGenes
  "ATM", "ATR", "EXO1", "MGMT", "RAD18", "RFC1", "TOP3A", "TOP3B", "XRCC6BP1"
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
### NOX4-adhesions/proliferation pan cancer with p53 ----

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


# to check the n of each gene 
ct <- Pan.final %>% 
  group_by(Gene.Symbol) %>% 
  summarise(n = length (Gene.Symbol)) # counts n

# extract those that do not have zero n's as a vector
apo.genes <- c(  # below are genes to screen for in last section
  "BAX", # from structural biology of p53 Joerger et al.x
  "BBC3",
  "PMAIP1",
  "PCNA",
  "GADD45A",
  "CDK1",
  "PAK1",
  "RRM2B",
  "CDKN1A", #P21
  "SERPINB5",
  "THBS1",
  "ANXA5",
  "ABL1", #senensece marker by sabioscience
  "AKT1",
  "ALDH1A3",
  "E2F1",
  "HRAS",
  "IGFBP3",
  "MYC",
  "CALR",
  "SOD1",
  "SOD2",
  "PRKCD",
  "MAPK14",
  "PTEN",
  "BCL2L11",
  "BID",
  "BAD",
  "BLK",
  "BIK",
  "BCL2L1",
  "BAG1",
  "BAG2",
  "BAG3",
  "BAG4",
  "BAG5",
  "CASP1",
  "CASP2",
  "CASP3",
  "CASP4",
  "CASP5",
  "CASP6",
  "CASP7",
  "CASP8",
  "CASP9",
  "CASP10",
  "CASP11",
  "CASP12",
  "CASP13",
  "CASP14",
  "TP53",
  "ATM",
  "BMI1",
  "CCND1",
  "CCNE1",
  "CDK2",
  "CDK4",
  "CDK6",
  "CDKN1A",
  "CDKN2A",
  "CDKN2D",
  "CHEK1",
  "CHEK2",
  "E2F3",
  "ETS1",
  "ETS2",
  "MDM2",
  "RB1",
  "RBL2",
  "ATM", #DNA DAMAGE SENSENCE INITIATORS
  "NBN",
  "TERF2",
  "TERT",
  "TP53BP1",
  "PRKCD",
  "PTPRJ", # aka DEP1; Althubiti et al., 2014 Characterization of novel markers of senescence and their prognostic potential in cancer
  "LAT2", # aka NTAL 
  "SLC9A3R1", # EBP50 
  "STX4",
  "VAMP3", 
  "B2M",# B2MG
  "LANCL1", 
  "VPS26A", 
  "PLD3",
  "ACHY", # senesence markers from Carnero, 2013
  "PPP1A",
  "CSN2",
  "BRF1",
  "PGM1",
  "IGFBP3",
  "IGFBP7",
  "SERPINE1",
  "MAPK14",
  "MAP2K6",
  "SMURF2"
)

apopsen <- ct$Gene.Symbol[ct$n != 0] # select ones that do have gene values

n.WT <- Pan.final[Pan.final$P53.Mutation =="WT",] # select wt patients
ct.WT <- n.WT %>%  # count wt n
  group_by(Gene.Symbol) %>% 
  summarise(n = length (Gene.Symbol))

n.Mut <- Pan.final[Pan.final$P53.Mutation != "WT",] # select Mut patients
ct.Mut <- n.Mut %>%  # count wt n
  group_by(Gene.Symbol) %>% 
  summarise(n = length (Gene.Symbol))

## seeded straight from Pan.final and filtering for only 2 genes so there are a lot more N than filtering for all at once.
DrawComparisonGraph <- function(bravo) {
  Bravo.expr <- Pan.final[Pan.final$Gene.Symbol == bravo,]
  NOX4.expr <- Pan.final[Pan.final$Gene.Symbol == "NOX4",]
  
  # identify each as WT or Mutant carrier
  NOX4.expr$Status <- NA # first create an empty column
  NOX4.expr$Status[NOX4.expr$P53.Mutation == "WT"] <- "WT"
  NOX4.expr$Status[NOX4.expr$P53.Mutation != "WT"] <- "Mut" # if NOT wt then assign "mut"
  
  NOX4.expr <- select(NOX4.expr, -P53.Mutation, -Case.Study, - Gene.Symbol, NOX4.expr = mRNA.Value, Status)
  
  Bravo.expr$Status <- NA # first create an empty column
  Bravo.expr$Status[Bravo.expr$P53.Mutation == "WT"] <- "WT"
  Bravo.expr$Status[Bravo.expr$P53.Mutation != "WT"] <- "Mut" # if NOT wt then assign "mut"
  
  Bravo.expr <- select(Bravo.expr, -P53.Mutation, -Case.Study, - Gene.Symbol, Bravo.expr = mRNA.Value, -Status)
  
  
  head(NOX4.expr)
  head(Bravo.expr)
  
  combi <- inner_join(NOX4.expr, Bravo.expr, by = "Patient.ID")
  head(combi)
  
  graphix <- 
    ggplot(combi) + 
    aes (x = log(NOX4.expr), y = log(Bravo.expr)) +
    geom_point(aes(color = Status), alpha = 0.1, shape = 3) +
    geom_smooth(aes(color= Status), method = lm, size = 1) +
    xlab("NOX4 Fold Changes") +
    ylab(paste(bravo, "Fold Changes", sep = " ")) +
    ggtitle(bravo) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, vjust=0.5),
          axis.text.y = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          axis.title.x = element_text(size = 10),
          panel.grid.minor.y=element_blank(),
          panel.grid.minor.x=element_blank(),
          panel.grid.major.x=element_blank(),
          panel.grid.major.y=element_blank(),
          panel.background = element_rect(fill = "white", color = "black", size = 3),
          axis.ticks = element_line(size=2))
  
  ggsave(graphix, file = paste(bravo, ".png", sep = ""), dpi=600, width = 5, height = 3)
  graphix
  WT <- combi[combi$Status == "WT",]
  Mut <- combi[combi$Status == "Mut",]
  
  # get rho and n
  rho.wt <- cor(WT$NOX4.expr, WT$Bravo.expr, method = "spearman")
  n.wt <- length(WT$Patient.ID)
  rho.mut <- cor(Mut$NOX4.expr, Mut$Bravo.expr, method = "spearman")
  n.mut <- length(Mut$Patient.ID)
  
  # print rho and n
  print(paste("WT rho ", rho.wt, sep = "= " ))
  print(paste("WT n ", n.wt, sep = "= "))
  print(paste("Mut rho ", rho.mut, sep = "= " ))
  print(paste("Mut n ", n.mut, sep = "= "))
  
}
# now run the loop
lapply(apopsen, DrawComparisonGraph)


DrawComparisonGraph("CDK1")


#### get NOX4 and a function of wt/mut p53 level
DrawComparisonGraph.with.p53 <- function(bravo) {
  Bravo.expr <- Pan.final[Pan.final$Gene.Symbol == bravo,]
  TP53.expr <- Pan.final[Pan.final$Gene.Symbol == "TP53",]
  
  # identify each as WT or Mutant carrier
  TP53.expr$Status <- NA # first create an empty column
  TP53.expr$Status[TP53.expr$P53.Mutation == "WT"] <- "WT"
  TP53.expr$Status[TP53.expr$P53.Mutation != "WT"] <- "Mut" # if NOT wt then assign "mut"
  
  TP53.expr <- select(TP53.expr, -P53.Mutation, -Case.Study, - Gene.Symbol, TP53.expr = mRNA.Value, Status)
  
  Bravo.expr$Status <- NA # first create an empty column
  Bravo.expr$Status[Bravo.expr$P53.Mutation == "WT"] <- "WT"
  Bravo.expr$Status[Bravo.expr$P53.Mutation != "WT"] <- "Mut" # if NOT wt then assign "mut"
  
  Bravo.expr <- select(Bravo.expr, -P53.Mutation, -Case.Study, - Gene.Symbol, Bravo.expr = mRNA.Value, -Status)
  
  
  head(TP53.expr)
  head(Bravo.expr)
  
  combi <- inner_join(TP53.expr, Bravo.expr, by = "Patient.ID")
  head(combi)
  
  graphix <- 
    ggplot(combi) + 
    aes (x = log(TP53.expr), y = log(Bravo.expr)) +
    geom_point(aes(color = Status), alpha = 0.1, shape = 3) +
    geom_smooth(aes(color= Status), method = lm, size = 1) +
    xlab("TP53 Fold Changes") +
    ylab(paste(bravo, "Fold Changes", sep = " ")) +
    ggtitle(bravo) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, vjust=0.5),
          axis.text.y = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          axis.title.x = element_text(size = 10),
          panel.grid.minor.y=element_blank(),
          panel.grid.minor.x=element_blank(),
          panel.grid.major.x=element_blank(),
          panel.grid.major.y=element_blank(),
          panel.background = element_rect(fill = "white", color = "black", size = 3),
          axis.ticks = element_line(size=2))
  
  ggsave(graphix, file = paste(bravo, "_p53.png", sep = ""), dpi=600, width = 5, height = 3)
  graphix
  WT <- combi[combi$Status == "WT",]
  Mut <- combi[combi$Status == "Mut",]
  
  # get rho and n
  rho.wt <- cor(WT$TP53.expr, WT$Bravo.expr, method = "spearman")
  n.wt <- length(WT$Patient.ID)
  rho.mut <- cor(Mut$TP53.expr, Mut$Bravo.expr, method = "spearman")
  n.mut <- length(Mut$Patient.ID)
  
  # print rho and n
  print(paste("WT rho ", rho.wt, sep = "= " ))
  print(paste("WT n ", n.wt, sep = "= "))
  print(paste("Mut rho ", rho.mut, sep = "= " ))
  print(paste("Mut n ", n.mut, sep = "= "))
  
}
# now run the loop
lapply(apopsen, DrawComparisonGraph.with.p53)
DrawComparisonGraph.with.p53("NOX4")

#### To make summary table of Pan.final here ----
TSS <- read_excel("~/Desktop/R_Cache/Tissue_source_site_code.xlsx")
head(TSS)

Pan.final$Code <- str_sub(Pan.final$Patient.ID, start = 6, end = 7)
TSS.pan.final <- inner_join(Pan.final, TSS, by = "Code") 
  
TSS.pan.final.2 <- select(TSS.pan.final, Patient.ID, Gene.Symbol, Study)
View(TSS.pan.final.2)

TSS.pan.final.raw <- TSS.pan.final.2 %>% 
  group_by(Gene.Symbol) %>% 
  summarise(n.study = length (Patient.ID) )

TSS.pan.final.count <- TSS.pan.final.2 %>% 
  group_by(Gene.Symbol, Study) %>% 
  summarise(n.study = length (Patient.ID) ) # listed by genes
  
View(TSS.pan.final.spread)
