
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





######################################### DOMAINS #######################################################
classes = rep("character", 16)
ecod.domain.descriptions = read.delim(ECOD.DOMAIN.DESCRIPTION.FILEPATH, 
                                      skip = 4, header = TRUE, 
                                      colClasses = classes, sep = "\t") %>%
  select(X.uid, ecod_domain_id, arch_name, x_name, h_name, t_name, f_name, f_id) 


ecod.domain.hits.raw = data.table::fread(
  file = sprintf("%sprot-families/functional/hhblits-ecod.txt", DATA.PATH))

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
  left_join(ecod.annotation.map, by = c("domain"))  
data.table::fwrite(ecod.domains.hits, file = sprintf("%secod.domains.hits.txt",OUTPUT.DATA.PATH))
rm(ecod.domain.hits.raw)


###################################### PAIRWISE HITS ##################################
hhr.table.filename = sprintf("%s/prot-families/all-by-all/hhblits/table-hhr.txt", DATA.PATH)
table.hhr = data.table::fread(hhr.table.filename) %>%
  data.table::setkey()


hiting.protein.lengths = table.hhr %>%
  distinct(qname, sname, qlength, slength)

protein.similarity.data = read.csv(PROFILE.SIMILARITY.TABLE,
                                   header = TRUE) %>%
  mutate(similar.proteins = 
           qcov >= MINIMUM.COV.FOR.PROTEIN.SIMILARITY | 
           scov >= MINIMUM.COV.FOR.PROTEIN.SIMILARITY) %>%
  left_join(hiting.protein.lengths)



pairwise.hits.data = table.hhr %>% 
  filter(prob >= MINIMUM.PROB.FOR.PAIRWISE.HIT) %>%
  mutate(hit.length = qend - qstart + 1) %>%
  filter(hit.length > MIN.HIT.LENGTH) %>%
  left_join(protein.similarity.data %>% select(qname, sname, similar.proteins))


#seq.mosaicism.data = table.hhr %>%
#  filter(prob >= MINIMUM.PROB.FOR.PAIRWISE.HIT &
#         pident >= MINIMUM.PIDENT.FOR.PAIRWISE.HIT) %>%
#  mutate(hit.length = qend - qstart + 1) %>%
#  filter(hit.length > MINIMUM.HIT.LENGTH) %>%
#  group_by(qname, sname)
#  mutate(share.a.frgment = TRUE) %>%
#  left_join(protein.similarity.data %>% select(qname, sname, similar.proteins))

  
data.table::fwrite(pairwise.hits.data, file = sprintf("%spairwise.hits.data.txt",OUTPUT.DATA.PATH))
rm(table.hhr)

#families = data.table::fread(FAMILIES.FILEPATH) %>% select(qname = members, family)
families.raw = readLines(FAMILIES.RAW.FILEPATH) %>%
  stringi::stri_split_lines() %>%
  lapply(FUN = function(x) {unlist(strsplit(x, split = "\t"))})
names(families.raw) = paste0("fam", 1:length(families.raw))
families = families.raw %>% stack()
names(families) = c("qname", "family")
data.table::fwrite(families, file = sprintf("%sfamilies.txt", OUTPUT.DATA.PATH))
