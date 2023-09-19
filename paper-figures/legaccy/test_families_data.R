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



# those seem the same

QCOVDSCVNEW = "/Users/bognasmug/MGG Dropbox/Projects/divRBP/phage-pp-workdir-refseq-hhblits/output/prot-families/families/dataset-full/repr-hits-pairwise-prob50.csv"
QCOVSCOVRAFAL = "/Users/bognasmug/MGG Dropbox/Projects/divRBP/minimal-input-data/python_pipeline/output/prot-families/families/dataset-full/repr-hits-pairwise-prob50.csv"
protein.similarity.data.raw.new = read.csv(QCOVDSCVNEW,header = TRUE) 
protein.similarity.data.raw.Rafal = read.csv(QCOVSCOVRAFAL,header = TRUE) 

p1 = protein.similarity.data.raw.new %>% arrange(qname, sname)
p2 = protein.similarity.data.raw.Rafal %>% arrange(qname, sname)

pfulljoin = full_join(p1, p2, by = c('qname', 'sname'))
# no numeric diference
summary(pfulljoin$prob.x - pfulljoin$prob.y)
summary(pfulljoin$pident.x - pfulljoin$pident.y)
summary(pfulljoin$scov_min.x - pfulljoin$scov_min.y)
summary(pfulljoin$scov_max.x - pfulljoin$scov_max.y)
summary(pfulljoin$qcov_min.x - pfulljoin$qcov_min.y)
summary(pfulljoin$qcov_max.x - pfulljoin$qcov_max.y)
summary(pfulljoin$min_cov.x - pfulljoin$min_cov.y)
summary(pfulljoin$max_cov.x - pfulljoin$max_cov.y)


p1notp2 = pfulljoin %>% filter(is.na(prob.y))
p2notp1 = pfulljoin %>% filter(is.na(prob.x))
p2a =  p2 %>% inner_join(p2notp1 %>% select(qname, sname)) %>% arrange(qname, sname)
p1a =  p1 %>% inner_join(p1notp2 %>% select(qname, sname)) %>% arrange(qname, sname)

p2notp1inv = p2a %>% rename(qname.new = sname, sname.new = qname, qcov.new_min = scov_min, qcov.new_max = scov_max,  scov.new_max = qcov_max,  scov.new_min = qcov_min) %>% 
  rename(qname = qname.new, sname = sname.new, qcov_min = qcov.new_min, qcov_max = qcov.new_max,  scov_max = scov.new_max ,  scov_min = scov.new_min) %>%
  select(qname, sname,prob, scov_min, scov_max, qcov_min, qcov_max, pident, max_cov, min_cov) %>%
  arrange(qname, sname)

ss = inner_join(p1a, p2notp1inv, by = c("qname", "sname"))
nrow(p1a) 
nrow(ss) 
summary(ss$prob.x - ss$prob.y)
summary(ss$pident.x - ss$pident.y)
summary(ss$scov_min.x - ss$scov_min.y)
summary(ss$scov_max.x - ss$scov_max.y)
summary(ss$qcov_min.x - ss$qcov_min.y)
summary(ss$qcov_max.x - ss$qcov_max.y)
summary(ss$min_cov.x - ss$min_cov.y)
summary(ss$max_cov.x - ss$max_cov.y)

FAM.NEW = "/Users/bognasmug/MGG Dropbox/Projects/divRBP/minimal-input-data/python_pipeline/output/prot-families/families/dataset-full/mcl/repr-hits-pairwise-prob95-cov80-mcl-I20.out"
FAM.R = "/Users/bognasmug/MGG Dropbox/Projects/divRBP/minimal-input-data/python_pipeline/output/prot-families/families/dataset-full/mcl/repr-hits-pairwise-prob95-cov80-mcl-RAFAL.out"

families.raw = readLines(FAM.NEW) %>%
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


families.raw2= readLines(FAM.R) %>%
  stringi::stri_split_lines() %>%
  lapply(FUN = function(x) {unlist(strsplit(x, split = "\t"))})
names(families.raw2) = paste0("fam", 1:length(families.raw2))
families2 = families.raw2 %>% stack()
names(families2) = c("qname", "family")
qnames.with.no.family2 = setdiff(repr.seq.lengths$name, families2$qname)
families.singletons2 = data.frame(qname = qnames.with.no.family2) %>% mutate(family2 = paste0("fam", 1+length(families.raw2):length(families.raw2)+length(qnames.with.no.family2)))
families2 = rbind(families2, families.singletons2) %>%
  filter(!(qname %in% unaccepted.genes$qname))


ggplot(families %>% group_by(family) %>% summarise(n = n_distinct(qname))) +
  geom_histogram(aes(x=n)) +
  xlim(c(0,100))
ggplot(families2 %>% group_by(family) %>% summarise(n = n_distinct(qname))) +
  geom_histogram(aes(x=n))+
  xlim(c(0,100))



n_distinct(families$family)
n_distinct(families2$family)

MCLINPUTKRZYSIEK = "/Users/bognasmug/MGG Dropbox/Projects/divRBP/minimal-input-data/python_pipeline/output/prot-families/families/dataset-full/mcl/repr-hits-pairwise-prob50-cov80-mcl-in.abc"
aaa = read.table(MCLINPUTKRZYSIEK)
MCLIINPUTRAFAL = "/Users/bognasmug/MGG Dropbox/Projects/divRBP/rafals-figures/data/families/prot-network-fam.txt"
aaa2 = read.table(MCLIINPUTRAFAL)
aaa2 = aaa2[2:nrow(aaa2),]
aaa = aaa %>% arrange(V1, V2)
aaa2 = aaa2 %>% arrange(V1, V2)

aR = aaa2  %>%
  mutate(reprseq.pair = CreatePair(V1, V2, collapse = "&"))
aK = aaa  %>%
  mutate(reprseq.pair = CreatePair(V1, V2, collapse = "&"))

bb = anti_join(aR %>% distinct(V3, reprseq.pair) %>% mutate(V3 = as.numeric(V3)), aK %>% distinct(V3, reprseq.pair) %>% mutate(V3 = as.numeric(V3)), by = "reprseq.pair")
bb %>% filter(V3 < 0.99) %>% head()
aaa %>% filter(V2 == "reprseq000014" & V1 == "reprseq012215")
aaa2 %>% filter(V1 == "reprseq000014" & V2 == "reprseq012215")

protein.similarity.data.raw.new %>% filter(qname == "reprseq000014" & sname == "reprseq012215" )
protein.similarity.data.raw.Rafal %>% filter(qname == "reprseq000014" & sname == "reprseq012215" )
