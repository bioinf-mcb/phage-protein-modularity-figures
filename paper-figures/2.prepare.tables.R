
######################################### ANNOTATIONS ########################################################
# read in our hhr hits to PHROGS, PHROG id tpo annotation mapping, and our custom mapping between old and new categories and annotations
phrogs.annotation.table = data.table::fread(PHROG.TABLE.PATH) %>% 
  select(phrog, annot, category)
hhr.phrog.hits.raw = data.table::fread(file = PHROG.HHPRED.HITS.PATH)
phrog.classes.manual.mapping = readxl::read_excel(MANUAL.PHROG.CLASS.MAPPING) %>%
  select(annot = funct.orig, category = class.orig, funct.new, class.new, no.seqs, no.phrogs) %>%
  group_by(funct.new, class.new) %>%
  mutate(no.seqs = sum(no.seqs)) %>%
  ungroup() %>%
  mutate(include = if_else(no.seqs >= MIN.NUM.PHROG.SEQ.PER.ANNOT, "TRUE", "FALSE")) %>%
  select(annot, category, class.new = class.new, funct.new = funct.new, include) 

# parse the hr table to get phrid ID, hit length, scov, qcov
hhr.phrog.hits.parsed = Parse.Hhr.Raw.Result(hhr.phrog.hits.raw) %>%
  select(qname, sname, qcov, scov, prob, eval, hit.length) %>%
  filter(hit.length >= MINIMUM.ALIGNMENT.LENGTH.FOR.ANNOTATION) 
hhr.phrog.hits.parsed$phrog = as.numeric(gsub("phrog_", "", hhr.phrog.hits.parsed$sname))

hhr.phrogs.annotated = hhr.phrog.hits.parsed %>% 
  left_join(phrogs.annotation.table, by = "phrog") %>%
  left_join(phrog.classes.manual.mapping, by = c("annot", "category")) %>%
  distinct(qname, qcov, scov, prob, eval, annotation = funct.new, category = class.new, include) %>%
  filter(!is.na(annotation) & category != "unknown function")

data.table::fwrite(hhr.phrogs.annotated, file = sprintf("%shhr.phrogs.annotated.txt",OUTPUT.DATA.PATH))
rm(hhr.phrog.hits.parsed)
rm(hhr.phrog.hits.raw)

######################################## ANNOTATIONS ###################################################
annotation.indexes = hhr.phrogs.annotated %>% 
  distinct(annotation, category) %>% 
  arrange(category, annotation) %>%
  mutate(annotation.index = row_number())

relaxely.annotated.proteins.including.multi.annot = GetPhrogAnnotationsData(hhr.phrogs.annotated, annotation.indexes, min.cov = MAIN.COVS.FOR.ANNOTATION[1], min.prob = DEFAULT.MINIMUM.PROB.FOR.ANNOTATION, min.eval = DEFAULT.NAXIMUM.EVAL.FOR.ANNOTATION)
mid.annotated.proteins.including.multi.annot = GetPhrogAnnotationsData(hhr.phrogs.annotated, annotation.indexes, min.cov = MAIN.COVS.FOR.ANNOTATION[2], min.prob = DEFAULT.MINIMUM.PROB.FOR.ANNOTATION, min.eval = DEFAULT.NAXIMUM.EVAL.FOR.ANNOTATION)
surely.annotated.proteins.including.multi.annot = GetPhrogAnnotationsData(hhr.phrogs.annotated, annotation.indexes, min.cov = MAIN.COVS.FOR.ANNOTATION[3], min.prob = DEFAULT.MINIMUM.PROB.FOR.ANNOTATION, min.eval = DEFAULT.NAXIMUM.EVAL.FOR.ANNOTATION)

data.table::fwrite(surely.annotated.proteins.including.multi.annot, file = sprintf("%ssurely.annotated.proteins.including.multi.annot",OUTPUT.DATA.PATH))
data.table::fwrite(relaxely.annotated.proteins.including.multi.annot, file = sprintf("%srelaxely.annotated.proteins.including.multi.annot",OUTPUT.DATA.PATH))
data.table::fwrite(mid.annotated.proteins.including.multi.annot, file = sprintf("%smid.annotated.proteins.including.multi.annot",OUTPUT.DATA.PATH))

# We only include categories that have MIN.NUM.PROT (including the proteins with multiple annotation) as defined at high coverage threshold
included.category.sizes = surely.annotated.proteins.including.multi.annot %>%
  filter(include == "TRUE") %>%
  group_by(annotation.index, annotation, category) %>%
  summarise(num.seq.per.annotation = n_distinct(qname[num.annots == 1]),
            num.seq.per.annotation.including.multi.annot = n_distinct(qname)) %>%
  ungroup()  %>%
  filter(num.seq.per.annotation.including.multi.annot >= MIN.NUM.PROT) %>%
  select(-num.seq.per.annotation.including.multi.annot)
data.table::fwrite(included.category.sizes, file = sprintf("%sincluded.category.sizes",OUTPUT.DATA.PATH))


########################################### ECOD DOMAIN HITS ####################################
# read hits to ECOD and the mapping between ecod id and t_names, h_names etc
classes = rep("character", 16)
ecod.domain.descriptions = read.delim(ECOD.DOMAIN.DESCRIPTION.FILEPATH, 
                                      skip = 4, header = TRUE, 
                                      colClasses = classes, sep = "\t") %>%
  select(X.uid, ecod_domain_id, arch_name, x_name, h_name, t_name, f_name, f_id) 

ecod.domain.hits.all = data.table::fread(file = ECOD.DOMAIN.HITS.PATH) %>%
  Parse.Hhr.Raw.Result() %>% 
  filter(  eval <= MAXIMUM.EVAL.FOR.DOMAIN &
           qcov >= NINIMUM.QCOV.FOR.DOMAIN &
           scov >= NINIMUM.SCOV.FOR.DOMAIN &
           prob >= MINIMUM.PROB.FOR.DOMAIN) %>%
  rename(domain = sname) %>%
  data.table::as.data.table()

# Note that some f_name or t_name can be found within different x groups so we need to work on full paths to a specific F,T,H,X
ecod.annotation.map = Get.Ecod.Annotation.Mapping(distinct_domains = ecod.domain.hits.all %>% distinct(domain) %>% pull(), 
                                                   ecod.domain.descriptions = ecod.domain.descriptions) 

ecod_id_to_names_map = Get.Ecod.Id.To.Name.Map(ecod.annotation.map) 
data.table::fwrite(ecod_id_to_names_map, file = sprintf("%secod_id_to_names_map.txt",OUTPUT.DATA.PATH))

ecod.domains.hits = ecod.domain.hits.all  %>%
  left_join(ecod.annotation.map, by = c("domain"))  
# THIS MAY BE COMMENTED OUT
#ecod.domains.hits = ecod.domains.hits[!grepl("DEAD", ecod.domains.hits$f_name)]
data.table::fwrite(ecod.domains.hits, file = sprintf("%secod.domains.hits.txt",OUTPUT.DATA.PATH))
rm(ecod.domain.hits.raw)


###################################### PAIRWISE HITS AND SEQUENCE MOSAICISM ##################################
seq.lengths = data.table::fread(REPR.SEQ.LENGTH.FILENAME)
protein.similarity.data.raw = read.csv(PROFILE.SIMILARITY.TABLE,header = TRUE) 
protein.similarity.data = protein.similarity.data.raw %>%
  # TODO: delete when we use the original file
  select(qname = query, sname = subject, prob, scov, qcov, pident) %>%
  left_join(seq.lengths %>% select(qname = name, qlength = length)) %>%
  left_join(seq.lengths %>% select(sname = name, slength = length))  %>%
  mutate(q.hit.length = qlength*qcov,
         s.hit.length = slength*scov) %>%
  as.data.table()
protein.similarity.data[, hit.length := pmin(s.hit.length, q.hit.length)]
protein.similarity.data[, max.cov := pmax(scov, qcov, na.rm = TRUE)]
protein.similarity.data[, min.hit.len := pmin(q.hit.length, s.hit.length, na.rm = TRUE)]

protein.similarity.data = protein.similarity.data %>%
  mutate(not.similar = max.cov <= MINIMUM.COV.FOR.PROTEIN.SIMILARITY,
         share.any.fragment = (prob >= MINIMUM.PROB.FOR.PAIRWISE.HIT/100 & 
                                 min.hit.len >= MIN.HIT.LENGTH)) %>%
  mutate(share.a.fragment = share.any.fragment & pident >= MINIMUM.PIDENT.FOR.PAIRWISE.HIT/100,
         share.a.fragment.pident10 = share.any.fragment & pident >= 0.1,
         share.a.fragment.pident30 = share.any.fragment & pident >= 0.3,
         share.a.fragment.pident50 = share.any.fragment & pident >= 0.5,
         share.a.fragment.pident70 = share.any.fragment & pident >= 0.7,
         share.a.fragment.pident90 = share.any.fragment & pident >= 0.9) %>%
  mutate(mosaic = share.a.fragment & not.similar,
         mosaic.pident10 = share.a.fragment.pident10 & not.similar,
         mosaic.pident30 = share.a.fragment.pident30 & not.similar,
         mosaic.pident50 = share.a.fragment.pident50 & not.similar,
         mosaic.pident70 = share.a.fragment.pident70 & not.similar,
         mosaic.pident90 = share.a.fragment.pident90 & not.similar)
data.table::fwrite(protein.similarity.data, file = sprintf("%sprotein.similarity.data.txt",OUTPUT.DATA.PATH))



###################################### OTHER ##################################
# raw pairwise hits
hhr.table.filename = sprintf("%s/prot-families/all-by-all/hhblits/table-hhr.txt", DATA.PATH)
table.hhr = data.table::fread(hhr.table.filename) %>%
  data.table::setkey()


# families
#families.raw = readLines(FAMILIES.RAW.FILEPATH) %>%
#  stringi::stri_split_lines() %>%
#  lapply(FUN = function(x) {unlist(strsplit(x, split = "\t"))})
#names(families.raw) = paste0("fam", 1:length(families.raw))
#families = families.raw %>% stack()
#names(families) = c("qname", "family")
#repr.seq.lengths = data.table::fread(file = sprintf("%sprot-families/representative/repr-seqs-lengths.txt", DATA.PATH))
#qnames.with.no.family = setdiff(repr.seq.lengths$name, families$qname)
#families.singletons = data.frame(qname = qnames.with.no.family) %>% mutate(family = paste0("fam", 1+length(families.raw):length(families.raw)+length(qnames.with.no.family)))
#families = rbind(families, families.singletons)
families = data.table::fread(FAMILIES.FILEPATH) %>% select(qname = members, family)
data.table::fwrite(families, file = sprintf("%sfamilies.txt", OUTPUT.DATA.PATH))


# cluste sizes
clustering = read.table(CLUSTERING_RESULTS_PATH)
names(clustering) = c('cluster', 'seq')
prot.names = read.table(PROTEIN_NAMES_MAPPING_PATH, sep = ",", header = TRUE) 
names(prot.names) = c("qname", "cluster")
cluster.sizes = clustering %>%
  left_join(prot.names, by = 'cluster') %>%
  group_by(qname) %>%
  summarise(n.prot.in.reprseq = n_distinct(seq))
data.table::fwrite(cluster.sizes, file = sprintf("%snum.prot.in.reprseq.txt", OUTPUT.DATA.PATH))



# recent mosaic pairs: read hhalign results and check if they are indeed mosaic
# there is always just one hit per (qname, sname)
recent.mosaicism.hhr = data.table::fread(file = HHALIGN_RECENT_MOSAICISM_PATH) %>%
  rowwise() %>%
  mutate(pair = paste(sort(c(qname, sname)), collapse = " AND ", sep = "")) %>%
  as.data.table()
multiple.hits.data = recent.mosaicism.hhr %>% group_by(pair) %>% summarise(n = n()) %>% filter(n > 1)
if (nrow(multiple.hits.data) > 0 ) {print("There are multiple hits per pairin hhalign data.")}

recent.mosaicism.hhr = recent.mosaicism.hhr[,q.hit.length := (as.numeric(qend) - as.numeric(qstart) + 1)]
recent.mosaicism.hhr = recent.mosaicism.hhr[,s.hit.length := (as.numeric(send) - as.numeric(sstart) + 1)]
recent.mosaicism.hhr = recent.mosaicism.hhr[,qcov := q.hit.length/as.numeric(qlength)]
recent.mosaicism.hhr = recent.mosaicism.hhr[,scov := s.hit.length/as.numeric(slength)]
recent.mosaicism.hhr = recent.mosaicism.hhr[, max.cov := pmax(qcov,scov)]
recent.mosaicism.hhr = recent.mosaicism.hhr[, prob := prob]
recent.mosaicism.hhr = recent.mosaicism.hhr[, pident := pident/100]
recent.mosaicism.hhr[, min.hit.len := pmin(q.hit.length, s.hit.length, na.rm = TRUE)]
recent.mosaicism.hhr = recent.mosaicism.hhr %>%
  mutate(not.similar = max.cov <= MINIMUM.COV.FOR.PROTEIN.SIMILARITY,
         share.any.fragment = (prob >= MINIMUM.PROB.FOR.PAIRWISE.HIT & 
                              min.hit.len >= MIN.HIT.LENGTH)) %>%
mutate(share.a.fragment.pident10 = share.any.fragment & pident >= 0.1,
        share.a.fragment.pident30 = share.any.fragment & pident >= 0.3,
        share.a.fragment.pident50 = share.any.fragment & pident >= 0.5,
        share.a.fragment.pident70 = share.any.fragment & pident >= 0.7,
        share.a.fragment.pident90 = share.any.fragment & pident >= 0.9) %>%
  mutate(
         mosaic.pident10 = share.a.fragment.pident10 & not.similar,
         mosaic.pident30 = share.a.fragment.pident30 & not.similar,
         mosaic.pident50 = share.a.fragment.pident50 & not.similar,
         mosaic.pident70 = share.a.fragment.pident70 & not.similar,
         mosaic.pident90 = share.a.fragment.pident90 & not.similar)
data.table::fwrite(recent.mosaicism.hhr, file = sprintf("%srecent.mosaicism.hhr.txt",OUTPUT.DATA.PATH))

