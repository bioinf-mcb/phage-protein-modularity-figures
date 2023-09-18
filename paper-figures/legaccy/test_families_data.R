# rename Rafal's table so that it looks like the table from Krzysiek's pipeline
protein.similarity.data.raw = read.csv("/Users/bognasmug/MGG Dropbox/Projects/divRBP/rafals-figures/data/families/table_qcov-scov_all.txt", header = TRUE)
protein.similarity.data.raw.new = protein.similarity.data.raw %>% 
  select(qname = query,sname = subject,prob,scov_min = scov.min,scov_max = scov.max,
         qcov_min = qcov.min, qcov_max = qcov.max,pident,max_cov = max.cov,min_cov = min.cov)
write.csv(protein.similarity.data.raw.new, 
      "/Users/bognasmug/MGG Dropbox/Projects/divRBP/minimal-input-data/python_pipeline/output/prot-families/families/dataset-full/repr-hits-pairwise-prob50.csv", 
      sep = ",", row.names = FALSE)

# test
# families
#"/Users/bognasmug/MGG Dropbox//Projects/divRBP/minimal-input-data/other/family-table.txt"
FAMILIES.FILEPATH = "/Users/bognasmug/MGG Dropbox/Projects/divRBP/rafals-figures/data/families/family-table.txt"
FAMILIES.RAW.FILEPATH = "/Users/bognasmug/MGG Dropbox/Projects/divRBP/rafals-figures/data/families/mcl-clusters-I2.txt" #"/Users/bognasmug/Downloads/out.seq.mci.I20"
families.raw = readLines(FAMILIES.RAW.FILEPATH) %>%
  stringi::stri_split_lines() %>%
  lapply(FUN = function(x) {unlist(strsplit(x, split = "\t"))})
names(families.raw) = paste0("fam", 1:length(families.raw))
families = families.raw %>% stack()
names(families) = c("qname", "family")
repr.seq.lengths = data.table::fread(file = REPR.SEQ.LENGTH.PATH)
qnames.with.no.family = setdiff(repr.seq.lengths$name, families$qname)
families.singletons = data.frame(qname = qnames.with.no.family) %>% mutate(family = paste0("fam", 1+length(families.raw):length(families.raw)+length(qnames.with.no.family)))
families = rbind(families, families.singletons) %>%
  filter(!(qname %in% unaccepted.genes$qname))

families2 = data.table::fread(FAMILIES.FILEPATH, header = TRUE) %>%
  select(family, qname = members) %>%
  filter(!(qname %in% unaccepted.genes$qname))

aa = families %>% as.data.frame() %>% group_by(family) %>% summarise(n = n())
aa2 = families2 %>% group_by(family) %>% summarise(n = n())
f1=aa %>% filter(n < 6 & n > 2) %>% head(1) %>% pull(family)
ff = families %>% filter(family == f1)
q = families %>% filter(qname == ff$qname[1])
f2 = families2 %>% filter(qname == q$qname)
ff2 = families2 %>% filter(family == f2$family)
ff
ff2

ggplot(aa) + geom_histogram(aes(x=n))
ggplot(aa2) + geom_histogram(aes(x=n))




FAMILIES.RAW.FILEPATH = sprintf("%spython_pipeline/output/prot-families/families/dataset-full/mcl/repr-hits-pairwise-prob95-cov80-mcl.out", DATA.PATH)
families.raw = readLines(FAMILIES.RAW.FILEPATH) %>%
  stringi::stri_split_lines() %>%
  lapply(FUN = function(x) {unlist(strsplit(x, split = "\t"))})
names(families.raw) = paste0("fam", 1:length(families.raw))
families = families.raw %>% stack()
names(families) = c("qname", "family")
repr.seq.lengths = data.table::fread(file = REPR.SEQ.LENGTH.PATH)
qnames.with.no.family = setdiff(repr.seq.lengths$name, families$qname)
families.singletons = data.frame(qname = qnames.with.no.family) %>% mutate(family = paste0("fam", 1+length(families.raw):length(families.raw)+length(qnames.with.no.family)))
families = rbind(families, families.singletons) %>%
  filter(!(qname %in% unaccepted.genes$qname))

FAMILIES.RAW.FILEPATH2 = sprintf("%spython_pipeline/output/prot-families/families/dataset-full/mcl/repr-hits-pairwise-prob95-cov80-mcl-OLD.out", DATA.PATH)
families.raw2= readLines(FAMILIES.RAW.FILEPATH2) %>%
  stringi::stri_split_lines() %>%
  lapply(FUN = function(x) {unlist(strsplit(x, split = "\t"))})
names(families.raw2) = paste0("fam", 1:length(families.raw2))
families2 = families.raw2 %>% stack()
names(families2) = c("qname", "family")
qnames.with.no.family2 = setdiff(repr.seq.lengths$name, families2$qname)
families.singletons2 = data.frame(qname = qnames.with.no.family2) %>% mutate(family2 = paste0("fam", 1+length(families.raw2):length(families.raw2)+length(qnames.with.no.family2)))
families2 = rbind(families2, families.singletons2) %>%
  filter(!(qname %in% unaccepted.genes$qname))




PROFILE.SIMILARITY.TABLE = sprintf("%spython_pipeline/output/prot-families/families/dataset-full/repr-hits-pairwise-prob50.csv", DATA.PATH)
PROFILE.SIMILARITY.TABLE2 = sprintf("%spython_pipeline/output/prot-families/families/dataset-full/repr-hits-pairwise-prob50_OLD.csv", DATA.PATH)
protein.similarity.data.raw = read.csv(PROFILE.SIMILARITY.TABLE,header = TRUE) 
protein.similarity.data.raw2 = read.csv(PROFILE.SIMILARITY.TABLE2,header = TRUE) 