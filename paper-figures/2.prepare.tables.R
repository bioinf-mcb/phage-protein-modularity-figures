
######################################### ANNOTATIONS ######################################################## Uniclust and GOs
phrogs.annotation.table = data.table::fread(PHROG.TABLE.PATH) %>% select(phrog, annot, category)
hhr.phrog.hits.raw = data.table::fread(file = PHROG.HHPRED.HITS.PATH)
phrog.classes.manual.mapping = readxl::read_excel(MANUAL.PHROG.CLASS.MAPPING) %>%
  select(annot = funct.orig, category = class.orig, structural.part = funct.new, category.corrected = class.new, no.seqs, no.phrogs) %>%
  mutate(structural = if_else(category.corrected %in% c("head and packaging", "tail", "connector"), "Y", "N")) %>%
  group_by(structural.part, category.corrected) %>%
  mutate(no.seqs = sum(no.seqs)) %>%
  ungroup() %>%
  mutate(include = if_else(no.seqs >= MIN.NUM.PHROG.SEQ.PER.ANNOT, "TRUE", "FALSE")) %>%
  select(annot, category, category.corrected = category.corrected, structural.part = structural.part, structural, include) 

hhr.phrog.hits.parsed = Parse.Hhr.Raw.Result(hhr.phrog.hits.raw) %>%
  select(qname, sname, qcov,scov,prob,eval,hit.length, annotation) %>%
  filter(hit.length >= MINIMUM.ALIGNMENT.LENGTH.FOR.ANNOTATION) 

hhr.phrog.hits.parsed$phrog = as.numeric(gsub("phrog_", "", hhr.phrog.hits.parsed$sname))

hhr.phrogs.annotated = hhr.phrog.hits.parsed %>% 
  left_join(phrogs.annotation.table, by = "phrog") %>%
  left_join(phrog.classes.manual.mapping, by = c("annot", "category")) %>%
  filter(!is.na(annot) & category != "unknown function") %>%
  #distinct(qname, annotation = structural.part, category = category.corrected, structural, include) %>% 
  distinct(qname, qcov, scov, prob, eval, annotation = structural.part, category = category.corrected, structural, include)

data.table::fwrite(hhr.phrogs.annotated, file = sprintf("%shhr.phrogs.annotated.txt",OUTPUT.DATA.PATH))
rm(hhr.phrog.hits.parsed)
rm(hhr.phrog.hits.raw)

######################################## ANNOTATIONS ###################################################

i=0
hhr.phrogs.list = list()
for (this.cov in MAIN.COVS.FOR.ANNOTATION) {
  i = i+1
  hhr.phrogs.list[[i]] = Get.significant.phrog.hits(hhr.phrogs.annotated, 
                                                    min.prob.threshold = DEFAULT.MINIMUM.PROB.FOR.ANNOTATION, 
                                                    min.qcov.threshold = this.cov, 
                                                    min.scov.threshold = this.cov, 
                                                    max.eval.threshold = DEFAULT.NAXIMUM.EVAL.FOR.ANNOTATION) %>%
    distinct(qname, annotation, category, structural, include) %>%
    mutate(annotation.coverage = this.cov)
}

hhr.phrogs = do.call('rbind', hhr.phrogs.list) %>% 
  mutate(annotation.coverage = factor(annotation.coverage, levels = MAIN.COVS.FOR.ANNOTATION))
data.table::fwrite(hhr.phrogs, file = sprintf("%shhr.phrogs",OUTPUT.DATA.PATH))

# now only look at annotations that we want to include
annotation.indexes = hhr.phrogs %>% 
  distinct(annotation, category) %>% 
  arrange(category, annotation) %>%
  mutate(annotation.index = row_number())

annotated.proteins.raw = hhr.phrogs %>%
  filter(!is.na(annotation)) %>%
  left_join(annotation.indexes, by = c("annotation", "category")) %>% 
  select(qname, annotation.index, annotation, category, structural, annotation.coverage, include) 

relaxely.annotated.proteins.including.multi.annot = annotated.proteins.raw %>%
  filter(annotation.coverage == min(MAIN.COVS.FOR.ANNOTATION)) %>%
  select(qname, annotation.index, annotation, category, structural, annotation.coverage, include) %>%
  group_by(qname) %>%
  # number of annotations including the excluded ones
  mutate(num.annots = n_distinct(annotation)) %>%
  ungroup() 

mid.annotated.proteins.including.multi.annot = annotated.proteins.raw %>%
  filter(annotation.coverage == 0.3) %>%
  select(qname, annotation.index, annotation, category, structural, annotation.coverage, include) %>%
  group_by(qname) %>%
  # number of annotations including the excluded ones
  mutate(num.annots = n_distinct(annotation)) %>%
  ungroup() 

surely.annotated.proteins.including.multi.annot = annotated.proteins.raw %>%
  filter(annotation.coverage == DEFAULT.NINIMUM.COV.FOR.ANNOTATION) %>%
  select(qname, annotation.index, annotation, category, structural, annotation.coverage, include) %>%
  group_by(qname) %>%
  mutate(num.annots = n_distinct(annotation)) %>%
  ungroup() 
data.table::fwrite(surely.annotated.proteins.including.multi.annot, file = sprintf("%ssurely.annotated.proteins.including.multi.annot",OUTPUT.DATA.PATH))
data.table::fwrite(relaxely.annotated.proteins.including.multi.annot, file = sprintf("%srelaxely.annotated.proteins.including.multi.annot",OUTPUT.DATA.PATH))
data.table::fwrite(mid.annotated.proteins.including.multi.annot, file = sprintf("%smid.annotated.proteins.including.multi.annot",OUTPUT.DATA.PATH))

# We only include categories that have MIN.NUM.PROT including the proteins with multiple annotation
included.category.sizes = rbind(surely.annotated.proteins.including.multi.annot, relaxely.annotated.proteins.including.multi.annot) %>%
  group_by(annotation.index, annotation, category, structural, annotation.coverage, include) %>%
  summarise(num.seq.per.annotation = n_distinct(qname[num.annots == 1]),
            num.seq.per.annotation.including.multi.annot = n_distinct(qname)) %>%
  ungroup()  %>%
  filter(num.seq.per.annotation.including.multi.annot >= MIN.NUM.PROT & include == "TRUE") %>%
  select(-num.seq.per.annotation.including.multi.annot)
data.table::fwrite(included.category.sizes, file = sprintf("%sincluded.category.sizes",OUTPUT.DATA.PATH))


######################################### DOMAINS #######################################################
########################################### ECOD DOMAIN HITS ####################################
# TODO: THIS MAY BE COMMENTED OUT
#ecod.domains.hits = ecod.domains.hits[!grepl("DEAD", ecod.domains.hits$f_name)]

classes = rep("character", 16)
ecod.domain.descriptions = read.delim(ECOD.DOMAIN.DESCRIPTION.FILEPATH, 
                                      skip = 4, header = TRUE, 
                                      colClasses = classes, sep = "\t") %>%
  select(X.uid, ecod_domain_id, arch_name, x_name, h_name, t_name, f_name, f_id) 





ecod.domain.hits.raw = data.table::fread(file = ECOD.DOMAIN.HITS.PATH)

ecod.domain.hits.parsed = Parse.Hhr.Raw.Result(ecod.domain.hits.raw) %>% 
  filter(  eval <= MAXIMUM.EVAL.FOR.DOMAIN &
           qcov >= NINIMUM.QCOV.FOR.DOMAIN &
           scov >= NINIMUM.SCOV.FOR.DOMAIN)

ecod.domain.hits.all = rbind(ecod.domain.hits.parsed %>% 
  filter(prob >= MINIMUM.PROB.FOR.DOMAIN) %>%
    mutate(minimum.prob = MINIMUM.PROB.FOR.DOMAIN),
  ecod.domain.hits.parsed %>% 
    filter(prob >= MINIMUM.PROB.FOR.DOMAIN.RELAXED) %>%
    mutate(minimum.prob = MINIMUM.PROB.FOR.DOMAIN.RELAXED)) %>%
  rename(domain = sname) %>%
  data.table::as.data.table()

# Note that some f_name or t_name can be found within different x groups so we need to work on full paths to a specific F,T,H,X
ecod.annotation.map = Get.Ecod.Annotation.Mapping.Numeric(distinct_domains = ecod.domain.hits.all %>% distinct(domain) %>% pull(), 
                                              ecod.domain.descriptions) 

ecod_id_to_names_map = 
  ecod.annotation.map %>% 
  group_by(x_id) %>% 
  summarise(name = paste(sort(unique(x_name)), collapse = " | ")) %>% 
  select(id = x_id, name) %>%
  mutate(level = "X")  %>%
  rbind(ecod.annotation.map %>% 
          group_by(h_id) %>% 
          summarise(name = paste(sort(unique(h_name)), collapse = " | ")) %>% 
          select(id = h_id, name) %>%
          mutate(level = "H")) %>%
  rbind(ecod.annotation.map %>% 
          group_by(t_id) %>% 
          summarise(name = paste(sort(unique(t_name)), collapse = " | ")) %>% 
          select(id = t_id, name) %>%
          mutate(level = "T") ) %>%
  rbind(ecod.annotation.map %>% 
          group_by(f_id) %>% 
          summarise(name = paste(sort(unique(f_name)), collapse = " | ")) %>% 
          select(id = f_id, name) %>%
          mutate(level = "F") )


no_h_name_map = ecod.annotation.map %>% 
  filter(h_name == "NO_H_NAME") %>% 
  group_by(h_id) %>% 
  summarise(h_name = paste(unique(t_name), collapse = " | ", sep = "")) %>%
  select(id = h_id, name = h_name) %>%
  mutate(level = "H")

ecod_id_to_names_map = rbind(ecod_id_to_names_map %>% filter(name != "NO_H_NAME"),
                             no_h_name_map)
data.table::fwrite(ecod_id_to_names_map, file = sprintf("%secod_id_to_names_map.txt",OUTPUT.DATA.PATH))


ecod.domains.hits = ecod.domain.hits.all  %>%
  left_join(ecod.annotation.map, by = c("domain"))  #%>%
  # delete the ones where family is unknown
  #filter(!is.na(f_ind))
  
  
data.table::fwrite(ecod.domains.hits, file = sprintf("%secod.domains.hits.txt",OUTPUT.DATA.PATH))
rm(ecod.domain.hits.raw)


###################################### PAIRWISE HITS ##################################
seq.lengths = data.table::fread(REPR.SEQ.LENGTH.FILENAME)
protein.similarity.data.raw = read.csv(PROFILE.SIMILARITY.TABLE,
                                   header = TRUE) 
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

hhr.table.filename = sprintf("%s/prot-families/all-by-all/hhblits/table-hhr.txt", DATA.PATH)
table.hhr = data.table::fread(hhr.table.filename) %>%
  data.table::setkey()


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
