# check if cds are ok
aaletters = str_split("ABCDEFGHIKLMNPQRSTUVWYZ", "") %>% unlist()
aalettersstr = sprintf("[%s]", paste0(c(aaletters, tolower(aaletters)), sep = "", collapse = ""))
all.cds = seqinr::read.fasta(ALL.CDS.PATH, as.string= TRUE)
name.table = read.table(NAME_TABLE_PATH, sep = ",", header = TRUE)

num.unaccepted.chars.in.sequences = lapply(all.cds, CheckAAcode) 
num.unaccepted.chars.in.sequences.df = num.unaccepted.chars.in.sequences %>% stack()
names(num.unaccepted.chars.in.sequences.df) = c("num.unaccepted.chars", "ncbi.id")
which(num.unaccepted.chars.in.sequences.df > 0) %>% length()
annotated.reprseq.with.unallowed.chars = num.unaccepted.chars.in.sequences.df %>%
  left_join(name.table %>% select(ncbi.id, qname = repr.name)) 
# qnames where the any member of the cluster  has too many "X" in the sequnce
unaccepted.genes = annotated.reprseq.with.unallowed.chars  %>% 
  filter(num.unaccepted.chars.in.sequences > ACCEPTED.NUM.UNKNOWN.AA) 
data.table::fwrite(unaccepted.genes, file = sprintf("%sunaccepted.genes.txt",OUTPUT.DATA.PATH))



######################################### ANNOTATIONS ########################################################
# families
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

#families2 = data.table::fread(FAMILIES.FILEPATH, header = TRUE) %>%
#  select(family, qname = members) %>%
#  filter(!(qname %in% unaccepted.genes$qname))
data.table::fwrite(families, file = sprintf("%sfamilies.txt", OUTPUT.DATA.PATH))




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
  filter(hit.length >= MINIMUM.ALIGNMENT.LENGTH.FOR.ANNOTATION) %>%
  filter(!(qname %in% unaccepted.genes$qname))
hhr.phrog.hits.parsed$phrog = as.numeric(gsub("phrog_", "", hhr.phrog.hits.parsed$sname))

hhr.phrogs.annotated = hhr.phrog.hits.parsed %>% 
  left_join(phrogs.annotation.table, by = "phrog") %>%
  left_join(phrog.classes.manual.mapping, by = c("annot", "category")) %>%
  distinct(qname, qcov, scov, prob, eval, annotation = funct.new, category = class.new, include) %>%
  filter(!is.na(annotation) & category != "unknown function") 

data.table::fwrite(hhr.phrogs.annotated, file = sprintf("%shhr.phrogs.annotated.txt",OUTPUT.DATA.PATH))
rm(hhr.phrog.hits.parsed)
rm(hhr.phrog.hits.raw)

# Add anti-defense descriptions
anti_defense_descriptions =read.xlsx(ANTI_DEFENSE_DESCRIPTION_PATH, sheetName = "Anti-defense & related", endRow = 246) %>%
  rowwise() %>%
  mutate(accession.short = GetAccessionShort(Accession))

refseq.hits.antidefense = data.table::fread(file = ANTI_DEFENSE_DATA_TABLE_PATH) %>%
  Parse.Hhr.Raw.Result() %>%
  select(qname.new = sname, sname.new = qname , qcov.new = scov, scov.new = qcov, prob, eval, hit.length) %>%
  select(qname = qname.new, sname = sname.new, qcov = qcov.new, scov = scov.new, prob, eval, hit.length) %>%
  filter(hit.length >= MINIMUM.ALIGNMENT.LENGTH.FOR.ANNOTATION) %>%
  left_join(anti_defense_descriptions %>% select(sname = accession.short, annotation = annotation)) %>%  
  mutate(category = "antidefense", include = TRUE) %>%
  distinct(qname, qcov, scov, prob, eval, annotation, category, include) %>%
  filter(!(qname %in% unaccepted.genes$qname))

######################################## ANNOTATIONS ###################################################
hhr.annotated = hhr.phrogs.annotated %>% 
  rbind(refseq.hits.antidefense)

manual.indices = data.frame(
  annotation = c("unknown", "antidefense"),
  category = c("unknown", "antidefense"),
  annotation.index = c(UNKNOWN.ANNOTATION.INDEX, GENERAL.ANTIDEFENSE.ANNOTATION.INDEX))

annotation.indices = 
  rbind(
    hhr.annotated %>%
      distinct(annotation, category) %>% 
      arrange(category, annotation) %>%
      mutate(annotation.index = row_number()),
    manual.indices)

relaxely.annotated.proteins.including.multi.annot = GetPhrogAnnotationsData(hhr.annotated, annotation.indices, min.cov = MAIN.COVS.FOR.ANNOTATION[1], min.prob = DEFAULT.MINIMUM.PROB.FOR.ANNOTATION, min.eval = DEFAULT.NAXIMUM.EVAL.FOR.ANNOTATION)
mid.annotated.proteins.including.multi.annot = GetPhrogAnnotationsData(hhr.annotated, annotation.indices, min.cov = MAIN.COVS.FOR.ANNOTATION[2], min.prob = DEFAULT.MINIMUM.PROB.FOR.ANNOTATION, min.eval = DEFAULT.NAXIMUM.EVAL.FOR.ANNOTATION)
surely.annotated.proteins.including.multi.annot = GetPhrogAnnotationsData(hhr.annotated, annotation.indices, min.cov = MAIN.COVS.FOR.ANNOTATION[3], min.prob = DEFAULT.MINIMUM.PROB.FOR.ANNOTATION, min.eval = DEFAULT.NAXIMUM.EVAL.FOR.ANNOTATION)

data.table::fwrite(annotation.indices, file = sprintf("%sannotation.indices",OUTPUT.DATA.PATH))
data.table::fwrite(hhr.annotated, file = sprintf("%shhr.annotated",OUTPUT.DATA.PATH))
data.table::fwrite(surely.annotated.proteins.including.multi.annot, file = sprintf("%ssurely.annotated.proteins.including.multi.annot",OUTPUT.DATA.PATH))
data.table::fwrite(relaxely.annotated.proteins.including.multi.annot, file = sprintf("%srelaxely.annotated.proteins.including.multi.annot",OUTPUT.DATA.PATH))
data.table::fwrite(mid.annotated.proteins.including.multi.annot, file = sprintf("%smid.annotated.proteins.including.multi.annot",OUTPUT.DATA.PATH))


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
  left_join(ecod.annotation.map, by = c("domain"))  %>%
  filter(!(qname %in% unaccepted.genes$qname))
# THIS MAY BE COMMENTED OUT
#ecod.domains.hits = ecod.domains.hits[!grepl("DEAD", ecod.domains.hits$f_name)]
data.table::fwrite(ecod.domains.hits, file = sprintf("%secod.domains.hits.txt",OUTPUT.DATA.PATH))
rm(ecod.domain.hits.raw)


###################################### PAIRWISE HITS AND SEQUENCE MOSAICISM ##################################
seq.lengths = data.table::fread(REPR.SEQ.LENGTH.FILENAME)
protein.similarity.data.raw = read.csv(PROFILE.SIMILARITY.TABLE,header = TRUE) %>%
  #rename(qname = query, sname = subject)
  rename(scov.min = scov_min, qcov.min = qcov_min, scov.max = scov_max, qcov.max = qcov_max, max.cov = max_cov, min.cov = min_cov) 
protein.similarity.data = protein.similarity.data.raw %>%
  select(qname, sname, prob, pident, scov.min, qcov.min, scov.max, qcov.max, max.cov, min.cov) %>%
  left_join(seq.lengths %>% select(qname = name, qlength = length)) %>%
  left_join(seq.lengths %>% select(sname = name, slength = length))  %>%
  mutate(q.hit.length = qlength*qcov.min,
         s.hit.length = slength*scov.min) %>%
  filter(!(qname %in% unaccepted.genes$qname) & !(sname %in% unaccepted.genes$qname)) %>%
  as.data.table()
#protein.similarity.data[, hit.length := pmin(s.hit.length, q.hit.length)]
#protein.similarity.data[, max.cov := pmax(scov, qcov, na.rm = TRUE)]
protein.similarity.data[, min.hit.len := pmin(q.hit.length, s.hit.length, na.rm = TRUE)]
protein.similarity.data = protein.similarity.data %>%
  # max.cov is max(max.scov, max.qcov)
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

protein.similarity.data.pident.above.30.pc = protein.similarity.data %>% 
  filter(pident >= 0.3) %>%
  mutate(reprseq.pair = CreatePair(qname, sname, collapse = "&"))
data.table::fwrite(protein.similarity.data.pident.above.30.pc, file = sprintf("%sprotein.similarity.data.pident.above.30.pc.txt",OUTPUT.DATA.PATH))


# cluste sizes: 
cluster.sizes = name.table %>% select(cluster = cluster.name, qname = repr.name, ncbi.id) %>%
  filter(!(qname %in% unaccepted.genes$qname)) %>%
  group_by(qname) %>%
  summarise(n.prot.in.reprseq = n_distinct(ncbi.id)) 
data.table::fwrite(cluster.sizes, file = sprintf("%snum.prot.in.reprseq.txt", OUTPUT.DATA.PATH))



############################# process metadata ##########################
# take NCBI genus and predicted lifestyle
metadata = read.csv(METADATA_PATH) %>% select(Accession, Molecule, Genus, Family,Order, Class, Phylum, Kingdom,Baltimore.Group, bacphlip_virulent_score, Jumbophage, Genus_ICTV_38, Family_ICTV38, Host) %>%
  mutate(lifestyle = if_else(bacphlip_virulent_score >= MIN.VIRULENT.BACPHLIP.SCORE, "virulent",
                             if_else(bacphlip_virulent_score <= MAX.TEMPERATE.BACPHLIP.SCORE, "temperate", "unknown"))) %>%
  mutate(GeneticMaterial = if_else(grepl("DNA", Molecule), "DNA", if_else(grepl("RNA", Molecule),"RNA", NA)))
# manually found RNA phages missing from inphrared
# from missing phages: 11 RNA phages: NC_012091, NC_042071 NC_042072 NC_042069 NC_055057 NC_042068 NC_012092 NC_042073. NC_004301 NC_055058 NC_001417 (found manually)
# about 15 obsolete DNA genomes
# metadata$GeneticMaterial[metadata$Accession %in% c('NC_012091', 'NC_042071', 'NC_042072', 'NC_042069', 'NC_055057', 'NC_042068', 'NC_012092', 'NC_042073', 'NC_004301', 'NC_055058', 'NC_001417')] = "RNA"


# note: some Genuses are unclassified and some lifestyles are "unknown"
reprseq.metadata.table = name.table %>%
  mutate(Accession =sub("\\..*", "", ncbi.id)) %>%
  mutate(unaccepted.gene = repr.name %in% unaccepted.genes) %>%
  filter(!unaccepted.gene) %>%
  left_join(metadata, by = "Accession") %>%
  group_by(repr.name) %>%
  summarise(num.genomes.with.this.rHMM = n_distinct(Accession),
            genus = unique.unclassified.rm(Genus),
            family = unique.unclassified.rm(Family),
            rder = unique.unclassified.rm(Order),
            class = unique.unclassified.rm(Class),
            kingdom = unique.unclassified.rm(Kingdom),
            phylum = unique.unclassified.rm(Phylum),
            baltimore.group = unique.unclassified.rm(Baltimore.Group),
            lifestyle = unique.unclassified.rm(lifestyle),
            jumbophage = unique.unclassified.rm(Jumbophage),
            genus_ICTV38 = unique.unclassified.rm(Genus_ICTV_38),
            family_ICTV38 = unique.unclassified.rm(Family_ICTV38),
            host = unique.unclassified.rm(Host),
            molecule = unique.unclassified.rm(GeneticMaterial)) %>%
  ungroup()


data.table::fwrite(reprseq.metadata.table, file = sprintf("%sreprseq.metadata.table_p09.txt", OUTPUT.DATA.PATH))
#data.table::fwrite(missing.genomes, file = sprintf("%smissing.genomes.txt", OUTPUT.DATA.PATH))
data.table::fwrite(metadata, file = sprintf("%smetadata.txt", OUTPUT.DATA.PATH))



