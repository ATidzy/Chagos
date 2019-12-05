# Required packages
# require("readtext")
# (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager", dependencies=TRUE)
#  install.packages("beepr", dependencies=TRUE)
#  BiocManager::install("dada2", dependencies=TRUE)
#  BiocManager::install("DECIPHER")
#  BiocManager::install("phyloseq")
#  BiocManager::install("decontam")
#  BiocManager::install("vegan")
#  install.packages("tidyverse")
#  install.packages("ggplot2")
#  installpackages("stringr")
#  install.packages("beepr")
# install.packages("readtext", dependencies = TRUE) 
library(DECIPHER); packageVersion("DECIPHER")
library(dada2); packageVersion("dada2")
library(beepr)
library(phyloseq)
library(decontam)
library(Biostrings); packageVersion("Biostrings")
library(ggplot2)
library(vegan)
library(tidyverse)
theme_set(theme_bw())
# 
# sessionInfo()

path <- "./raw_data/"
list.files(path)

# Forward   and reverse fastq filenames have format: 
#SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))


# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 3)
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])


# Place filtered files in filtered subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))


names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) 
head(out)
# On Windows always set multithread=FALSE
# multithread provides multiple threads of execution concurrently, 
# which is not supported by the operating system.

#error model
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

?file
?open.connection

#Reload point
filtFs <- list.files(path = "./raw_data/filtered", pattern = "_F_filt.fastq.gz",full.names = TRUE)
filtRs <- list.files(path = "./raw_data/filtered", pattern = "_R_filt.fastq.gz",full.names = TRUE)
# ?get
# filtFs <- file.path(path= "./raw_data/filtered", pattern= "_F_filt.fastq.gz")
# filtFs <- get(list.files(path = "./raw_data/filtered", pattern = "_F_filt.fastq.gz"))

# Dereplicate
derepFs <- derepFastq(filtFs, verbose=TRUE,n = 1e5)
derepRs <- derepFastq(filtRs, verbose=TRUE)


# Name the derep-class objects by sample name
names(derepFs) <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
names(derepRs) <- sapply(strsplit(basename(filtRs), "_"), `[`, 1)


# apply the core sample inference algorithm 
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)


#Inspecting the returned dada-class object:
dadaFs[[1]]
dadaRs[[1]]

#constructing the merged “contig” sequences
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers[[1]])


#Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)


# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))


#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", 
                                    multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)


# rename seqtab object samples
seqtab.df <- as.data.frame(seqtab.nochim)
row.names(seqtab.df) <- map(strsplit(row.names(seqtab.df), "_"),1)


# re-load point for non-chimeric sequences
write.csv(x=seqtab.df, sep=",", file ="./seqtab.csv")
saveRDS(taxa,"./seqtab.csv")


#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), 
               rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
track <- as.data.frame(track)

plot(track$nonchim)


#Assign taxonomy
#seqtab.nochim <- read.csv("./seqtab.csv")
# Create a DNAStringSet from the ASVs
dna <- DNAStringSet(getSequences(seqtab.nochim)) 
load("../Chagos/Chagos_Project/Silva/SILVA_SSU_r132_March2018.RData") #TRAINING SET
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) 
#processors=1 because multithread is not supported on OS. Otherwise use NULL
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") 


# Convert the output object of class "Taxa" to a matrix 
#analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)
#add 16s species
# Inspect taxonomic assignments 
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
saveRDS(taxa,"taxa.RDS")

# re-load point for taxonomy
taxa <- readRDS("taxa.RDS")

#Import into phyloseq:


# read in metadata
library(readxl)
meta = read_excel("metadata/Chagos_PlateMpas_MetaData_1_Geoff.xlsx")

# subset meta and seqtab to match
in.meta = which(names(seqtab.nochim[,1]) %in% meta$`Pooled library ID (from GIS)`== TRUE)
seqtab.nochim = seqtab.nochim[in.meta,]

in.seqtab = which(meta$`Pooled library ID (from GIS)`%in% names(seqtab.nochim[,1]))
meta = meta[in.seqtab,]

# re-order
meta = meta[order(meta$`Pooled library ID (from GIS)`),]
row.names(meta) <- meta$`Pooled library ID (from GIS)`
# Check 
#identical(row.names(seqtab.nochim), meta$`Pooled library ID (from GIS)`)
row.names(seqtab.nochim)
meta$`Pooled library ID (from GIS)` <-  as.character(meta$`Pooled library ID (from GIS)`)

currentorder = order(as.numeric(names(seqtab.nochim[,1])))
seqtab.nochim=seqtab.nochim[currentorder,]
names(seqtab.nochim[,1])

identical(row.names(seqtab.nochim), meta$`Pooled library ID (from GIS)`)

meta$Negatives = (meta$salinity == "Blank")
# Remove contaminants ####
# find contaminants
contams = isContaminant(seqtab.nochim, neg = meta$Negatives, normalize = TRUE)
table(contams$contaminant)
write.csv(contams, file = "likely_contaminants.csv", row.names = TRUE)

# remove them
seqtab.nochim = seqtab.nochim[,(which(contams$contaminant != TRUE))]
seqtab.nochim = seqtab.nochim[meta$Negatives == FALSE,]
meta = meta[meta$Negatives == FALSE,]

# subset meta and seqtab to match
in.meta = which(names(seqtab.nochim[,1]) %in% meta$`Pooled library ID (from GIS)` == TRUE)
seqtab.nochim = seqtab.nochim[in.meta,]

in.seqtab = which(meta$`Pooled library ID (from GIS)` %in% names(seqtab.nochim[,1]))
meta = meta[in.seqtab,]

#re-order
meta = meta[order(meta$`Pooled library ID (from GIS)`),]
row.names(meta) <- meta$`Pooled library ID (from GIS)`
# Check
identical(row.names(seqtab.nochim), meta$`Pooled library ID (from GIS)`)


# make phyloseq object ####
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(meta), 
               tax_table(taxa))
# save it
saveRDS(ps, "phyloseq_object_16S.RDS")