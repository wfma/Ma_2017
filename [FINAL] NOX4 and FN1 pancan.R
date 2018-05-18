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
Pan.final

NOX4 <- Pan.final[Pan.final$Gene.Symbol %in% "NOX4",]
NOX4 <- select(NOX4, 
               NOX4.mRNA = mRNA.Value,
               Patient.ID,
               P53.Mutation)

FN1  <- Pan.final[Pan.final$Gene.Symbol %in% "FN1",]
FN1 <- select(FN1, 
               FN1.mRNA = mRNA.Value,
               Patient.ID,
               P53.Mutation)

combiFN1NOX4 <- inner_join(FN1, NOX4, by = "Patient.ID")
length(combiFN1NOX4$FN1.mRNA)

FN1NOX4plot <- ggplot(combiFN1NOX4) + aes( x = log(NOX4.mRNA), y = log(FN1.mRNA)) +
  geom_point(alpha = 0.3) + geom_smooth(method = "lm", color = "orange", size = 2) +
  xlab("NOX4 Fold Change") +
  ylab("FN1 Fold Change")+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 20, vjust=0.5),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 23),
        axis.title.x = element_text(size = 23),
        panel.grid.minor.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.background = element_rect(fill = "white", color = "black", size = 4),
        axis.ticks = element_line(size=2))

cor(combiFN1NOX4$FN1.mRNA, combiFN1NOX4$NOX4.mRNA, method = "spearman")
modellm <- lm(log(combiFN1NOX4$FN1.mRNA) ~ log(combiFN1NOX4$NOX4.mRNA))
summary(modellm)
n.wt <- length(combiFN1NOX4$Patient.ID)

ggsave(FN1NOX4plot, file = "FN1NOX4plot.png", width = 7, height = 6, dpi = 600)

