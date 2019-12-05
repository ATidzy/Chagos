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
#library(DECIPHER); packageVersion("DECIPHER")
library(dada2); packageVersion("dada2")
#library(beepr)
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
#plotQualityProfile(fnFs[1:2])
#plotQualityProfile(fnRs[1:2])


# Place filtered files in filtered subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))


names(filtFs) <- sample.names
names(filtRs) <- sample.names

#out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250,160),
maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
compress=TRUE, multithread=TRUE) 
head(out)
# On Windows always set multithread=FALSE
# multithread provides multiple threads of execution concurrently, 
# which is not supported by the operating system.

#error model
#errF <- learnErrors(filtFs, multithread=TRUE)
#errR <- learnErrors(filtRs, multithread=TRUE)
#plotErrors(errF, nominalQ=TRUE)
#plotErrors(errR, nominalQ=TRUE)

#Reload point
filtFs <- list.files(path = "./raw_data/filtered", pattern = "_F_filt.fastq.gz",full.names = TRUE)
filtRs <- list.files(path = "./raw_data/filtered", pattern = "_R_filt.fastq.gz",full.names = TRUE)
# ?get
# filtFs <- file.path(path= "./raw_data/filtered", pattern= "_F_filt.fastq.gz")
# filtFs <- get(list.files(path = "./raw_data/filtered", pattern = "_F_filt.fastq.gz"))

# Dereplicate
derepFs <- derepFastq(filtFs, verbose=TRUE,n = 1e5)
saveRDS(derepFs,"./derepFs.RDS")

derepRs <- derepFastq(filtRs, verbose=TRUE)
saveRDS(derepRs,"./derepRs.RDS")
