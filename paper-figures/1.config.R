# Also install data.table, tidyr
library(R.utils)
# adding it at the bottom may brak igraph plotting
library(ggiraph)  
library(stringr)
library(dplyr)
library(ggplot2)
library(data.table)
library(ggnetwork)
set.seed(1)
options(stringsAsFactors = F)


DROPBOX.PATH =  "/Users/bognasmug/MGG Dropbox/"
# TO DO: CHANGE THAT WHEN WE HAVE NEW DATA FROM THE CLUSTER
#FAMILIES.FILEPATH = "/Users/bognasmug/MGG Dropbox/Bogna Smug/Projects/divRBP/phage-pp-workdir-refseq-hhblits/output/prot-families/families/families/dataset-full/family-table.txt"
FAMILIES.RAW.FILEPATH = "/Users/bognasmug/MGG Dropbox/Projects/divRBP/phage-pp-workdir-refseq-hhblits/output/prot-families/families/dataset-full/mcl/repr-hits-pairwise-prob95-mcl.out"
# DROPBOX.PATH =  "/Users/rmostowy/MGG Dropbox/"
# FAMILIES.FILEPATH <- "/Users/rmostowy/MGG Dropbox/Rafal Mostowy/Projects/Data/divRBP-large/families/family-table.txt"


PROJECT.PATH = sprintf("%s/Projects/divRBP/",DROPBOX.PATH)
source(sprintf("%spaper-figures/scripts/helpers.R", PROJECT.PATH))
DATA.PATH = sprintf("%s/phage-pp-workdir-refseq-hhblits/output/", PROJECT.PATH)
PROFILE.SIMILARITY.TABLE = sprintf("%sprot-families/families/dataset-full/repr-hits-pairwise-prob50.csv", DATA.PATH)
PHROG.TABLE.PATH = sprintf("%sDatabases/Phrogs/phrog_annot_v4.tsv", DROPBOX.PATH)
PHROG.HHPRED.HITS.PATH = sprintf("%sprot-families/functional/hhblits-phrogs.txt", DATA.PATH)
# downloaded from: http://prodata.swmed.edu/ecod/complete/distribution:
ECOD.DOMAIN.DESCRIPTION.FILEPATH = sprintf("%sDatabases/ECOD/ECOD/ecod.develop283.domains.txt", DROPBOX.PATH)
MANUAL.PHROG.CLASS.MAPPING = sprintf("%sDatabases/Phrogs/custom/v3beta_phrogs-table-rafal-27_09.xlsx", DROPBOX.PATH) 
GENE.ONTOLOGY.FILE.PATH =  sprintf("%sDatabases/GeneOntology-org/releases/2021-07-02/go-basic.obo", DROPBOX.PATH)
GENOME.METADATA.PATH = sprintf("%sDatabases/NCBI/Virus/download_2021-12-03/sequences.csv", DROPBOX.PATH)

# where to output tables and figures
OUTPUT.DATA.PATH = sprintf("%spaper-figures/tables/", PROJECT.PATH)
OUTPUT.FIGURES.PATH = sprintf("%spaper-figures/", PROJECT.PATH)
dir.create(OUTPUT.FIGURES.PATH, recursive = TRUE)
dir.create(OUTPUT.DATA.PATH, recursive = TRUE)



MINIMUM.COVERAGES.FOR.ANNOTATION = seq(0.05, 1, 0.05)
MINIMUM.PROBABILITIES.FOR.ANNOTATION = seq(5, 100, 5)
NAXIMUM.EVALS.FOR.ANNOTATION = 10^seq(-10,5,1)
MINIMUM.ALIGNMENT.LENGTH.FOR.ANNOTATION = 30

DEFAULT.NINIMUM.COV.FOR.ANNOTATION = 0.8
DEFAULT.MINIMUM.PROB.FOR.ANNOTATION = 95
DEFAULT.NAXIMUM.EVAL.FOR.ANNOTATION = 10^5
MAIN.PROBS.FOR.ANNOTATION = c(50, 95)
MAIN.COVS.FOR.ANNOTATION = c(0.1, 0.8)



NINIMUM.QCOV.FOR.DOMAIN = 0
NINIMUM.SCOV.FOR.DOMAIN = 0.7
MINIMUM.PROB.FOR.DOMAIN = 95
MINIMUM.PROB.FOR.DOMAIN.RELAXED = 70
MAXIMUM.EVAL.FOR.DOMAIN = 10^(-5)
MIN.NUM.PROT = 20
MIN.NUM.PHROG.SEQ.PER.ANNOT = 500


# f
MINIMUM.PROB.FOR.PAIRWISE.HIT = 95
# two same
MINIMUM.COV.FOR.PROTEIN.SIMILARITY = 0.3
MIN.HIT.LENGTH = 30
MINIMUM.PIDENT.FOR.PAIRWISE.HIT = 30

#MIN.PERC.PROTEINS.WITH.H.DOMAIN = 20

MIN.NUM.PROT.WITH.2.X.DOMAINS = 10
MIN.PROB.COND.FOR.AMBIGUOUS.ANNOTATION = 0.3
MIN.NUM.ECOD.MOSAIC.PAIRS.BETWEEN.FUNCTIONS = 3
MIN.NUM.ECOD.MOSAIC.PAIR.FAMILIES.BETWEEN.FUNCTIONS = 3
MIN.NUM.ECOD.MOSAIC.PAIR.FAMILIES.WITHIN.FUNCTIONS = 3
SHARED.DOMAIN.LEVEL = "t_id"

FIG.WIDTH = 20
FIG.HEIGHT = 10

theme.no.verical = theme(
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.major.y = element_blank(),
  panel.grid.minor.y = element_blank()
)
text.size= 10


PHROG.COLOR.MAP = c("integration and excision" = "#FF7F00", 
                  "head and packaging" = "#C7E9C0", 
                  "transcription regulation" = "#984EA3", 
                  "connector" = "#41AB5D",
                  "tail" = "#006D2C",
                  "lysis" = "#E41A1C", 
                  "DNA, RNA and nucleotide metabolism" = "#377EB8", 
                  "other" = "#999999",
                  "moron, auxiliary metabolic gene and host takeover" = "#F781BF")

