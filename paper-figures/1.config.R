# Also install data.table, tidyr
library(R.utils)
# adding it at the bottom may brak igraph plotting
library(ggiraph)  
library(stringr)
library(dplyr)
library(ggplot2)
library(data.table)
library(ggnetwork)
library("xlsx")
set.seed(1)
options(stringsAsFactors = F)


###################################################### user defined paths ###############################
DROPBOX.PATH =  "/Users/bognasmug/MGG Dropbox/"
PROJECT.PATH = sprintf("%s/Projects/divRBP/",DROPBOX.PATH)
# DROPBOX.PATH =  "/Users/rmostowy/MGG Dropbox/"

###################################################### user defined parameters ###############################
# parameters related to functional annotation
# exploring how annotation depends on HMM-HMM comparison parameters
# for Fig S1a
MINIMUM.COVERAGES.FOR.ANNOTATION.EXPLORATION = seq(0.05, 1, 0.05)
MAIN.PROBS.FOR.ANNOTATION = c(50, 95)
# for fig S1b and  S2
MINIMUM.PROBABILITIES.FOR.ANNOTATION.EXPLORATION = seq(5, 100, 5)
MAIN.COVS.FOR.ANNOTATION = c(0.1, 0.3, 0.8)
# default parameters for functonal annotation
MINIMUM.ALIGNMENT.LENGTH.FOR.ANNOTATION = 30
DEFAULT.NINIMUM.COV.FOR.ANNOTATION = 0.8
DEFAULT.MINIMUM.PROB.FOR.ANNOTATION = 95
DEFAULT.NAXIMUM.EVAL.FOR.ANNOTATION = 10^5
# number of proteins/sequences needed to be in a functional category to be included in our analyses
MIN.NUM.PROT = 20
MIN.NUM.PHROG.SEQ.PER.ANNOT = 500
# fir fig S2
MIN.NUM.SHARED.SEQ.FOR.SHARED.ANNOT.GRAPH = 20

# parameters related to domain detectiom
NINIMUM.QCOV.FOR.DOMAIN = 0
NINIMUM.SCOV.FOR.DOMAIN = 0.7
MINIMUM.PROB.FOR.DOMAIN = 95
MAXIMUM.EVAL.FOR.DOMAIN = 10^(-5)
MIN.NUM.REPRSEQS.FOR.ECOD.PRESENCE = 5

# parameters related to pairwise similarity and seq. mosaicism
MINIMUM.PROB.FOR.PAIRWISE.HIT = 95
MINIMUM.COV.FOR.PROTEIN.SIMILARITY = 0.5
MIN.HIT.LENGTH = 50
MINIMUM.PIDENT.FOR.PAIRWISE.HIT = 30
MIN.NUM.REPRSEQ.PAIRS.FOR.BETWEEN.FUNCTION.MOSAICISM = 20

# min.prop. proteins in class with any ECOD hit (FigS8)
MIN.PROP.PROTS.WITH.ECOD.HIT = 0.2

###################################################### sort all vectors in case they are input in bad order ############################
MAIN.COVS.FOR.ANNOTATION = sort(MAIN.COVS.FOR.ANNOTATION)
MAIN.PROBS.FOR.ANNOTATION = sort(MAIN.PROBS.FOR.ANNOTATION)
MINIMUM.PROBABILITIES.FOR.ANNOTATION.EXPLORATION = sort(MINIMUM.PROBABILITIES.FOR.ANNOTATION.EXPLORATION)
MAIN.COVS.FOR.ANNOTATION = sort(MAIN.COVS.FOR.ANNOTATION)
######################################################  create all necessary paths ###################################################### 
source(sprintf("%scode/protein-mosaicism/paper-figures/helpers.R", PROJECT.PATH))
DATA.PATH = sprintf("%s/phage-pp-workdir-refseq-hhblits/output/", PROJECT.PATH)

# TO DO: CHANGE IT BACK TO MAIN CLUSTER PROJECT ONCE IT IS CORRECT
#FAMILIES.FILEPATH = "/Users/bognasmug/MGG Dropbox/Projects/divRBP/rafals-figures/data/families/family-table.txt"  #"/Users/bognasmug/MGG Dropbox/Bogna Smug/Projects/divRBP/phage-pp-workdir-refseq-hhblits/output/prot-families/families/families/dataset-full/family-table.txt"
FAMILIES.RAW.FILEPATH = sprintf("%sprot-families/families/dataset-full/mcl/repr-hits-pairwise-prob95-cov80-mcl.out", DATA.PATH)
#PROFILE.SIMILARITY.TABLE = sprintf("%sprot-families/families/dataset-full/repr-hits-pairwise-prob50.csv", DATA.PATH)
PROFILE.SIMILARITY.TABLE = "/Users/bognasmug/MGG Dropbox/Projects/divRBP/rafals-figures/data/families/table_qcov-scov_all.txt"


# all vs all / recent HGT pairs
HHALIGN_RECENT_MOSAICISM_PATH = sprintf("%sprot-families/all-by-all/hhblits/table-hhr-hgt-hhalign.txt", DATA.PATH)

# data from mmseq clustes from which the representative sequences were selected and the lengths of representative sequences
CLUSTERING_RESULTS_PATH = sprintf("%sprot-families/representative/clustering.tsv", DATA.PATH)
PROTEIN_NAMES_MAPPING_PATH = sprintf("%sprot-families/representative/name-table.txt", DATA.PATH)
REPR.SEQ.LENGTH.FILENAME = sprintf("%sprot-families/representative/repr-seqs-lengths.txt", DATA.PATH)
# data from HMM-HMM comparison results and annotation metadata
PHROG.HHPRED.HITS.PATH = sprintf("%sprot-families/functional/hhblits-phrogs.txt", DATA.PATH)
ECOD.DOMAIN.HITS.PATH = sprintf("%sprot-families/functional/hhblits-ecod.txt", DATA.PATH)
# downloaded from: http://prodata.swmed.edu/ecod/complete/distribution:
ECOD.DOMAIN.DESCRIPTION.FILEPATH = sprintf("%sDatabases/ECOD/ECOD/ecod.develop283.domains.txt", DROPBOX.PATH)
PHROG.TABLE.PATH = sprintf("%sDatabases/Phrogs/phrog_annot_v4.tsv", DROPBOX.PATH)
MANUAL.PHROG.CLASS.MAPPING = sprintf("%sDatabases/Phrogs/custom/v3_phrogs-table-rafal-3_12.xlsx", DROPBOX.PATH) 
# where to output tables and figures
OUTPUT.DATA.PATH = sprintf("%spaper-figures/tables/", PROJECT.PATH)
OUTPUT.FIGURES.PATH = sprintf("%spaper-figures/", PROJECT.PATH)
dir.create(OUTPUT.FIGURES.PATH, recursive = TRUE)
dir.create(OUTPUT.DATA.PATH, recursive = TRUE)



######################################################## visual parameters ################################################################
# theme for shiny-like plots
theme.no.verical = theme(
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.major.y = element_blank(),
  panel.grid.minor.y = element_blank()
)
# color and label mapping
PHROG.COLOR.MAP = c("integration and excision" = "#FF7F00", 
                  "head and packaging" = "#C7E9C0", 
                  "transcription regulation" = "#984EA3", 
                  "connector" = "#41AB5D",
                  "tail" = "#006D2C",
                  "lysis" = "#E41A1C", 
                  "DNA, RNA and nucleotide metabolism" = "#377EB8", 
                  "other" = "#999999",
                  "moron, auxiliary metabolic gene and host takeover" = "#F781BF")

phrogs.class.order = c("head and packaging", "connector", "tail", 
                        "DNA, RNA and nucleotide metabolism", "integration and excision",
                        "moron, auxiliary metabolic gene and host takeover", "transcription regulation",
                        "lysis", "other")

phrog.class.labels = c(
  'integration and excision' = "integration\nand excision",
  'head and packaging' = "head\nand\npackaging",
  "connector" = "connector",
  "tail" = "tail",
  "lysis" = "lysis",
  'transcription regulation' = "transcription\nregulation",
  'DNA, RNA and nucleotide metabolism' = "DNA, RNA and\nnucleotide\nmetabolism",
  "other" = "other",
  'moron, auxiliary metabolic gene and host takeover' = "metabolic"#"moron, auxiliary\nmetabolic gene\nand host takeover"
)
