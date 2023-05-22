no.self.hit.qnames = c('reprseq010712',
                       'reprseq011861',
                       'reprseq020849',
                       'reprseq038307',
                       'reprseq046493', 
                       'reprseq053180',
                       'reprseq061119',
                       'reprseq061358',
                       'reprseq065673',
                       'reprseq068528',
                       'reprseq070520',
                       'reprseq073233',
                       'reprseq077124', 
                       'reprseq085833',
                       'reprseq090243',
                       'reprseq094782',
                       'reprseq095733',
                       'reprseq111039',
                       'reprseq112537',
                       'reprseq113881',
                       'reprseq118405',
                       'reprseq124736',
                       'reprseq126575',
                       'reprseq128529')


no.self.hits.and.with.repeats.qnames = c('reprseq073233', 'reprseq077124', 'reprseq126575', 'reprseq053180',
'reprseq085833', 'reprseq020849', 'reprseq128529', 'reprseq065673',
'reprseq070520', 'reprseq010712', 'reprseq090243', 'reprseq111039',
'reprseq112537')


hhr.phrogs.annotated = data.table::fread(sprintf("%shhr.phrogs.annotated.txt",OUTPUT.DATA.PATH))
relaxely.annotated.proteins.including.multi.annot =data.table::fread(sprintf("%srelaxely.annotated.proteins.including.multi.annot",OUTPUT.DATA.PATH)) %>%
  left_join(families)
anots = hhr.phrogs.annotated %>% filter(qname %in% no.self.hit.qnames) %>% group_by(qname) %>% arrange(qcov, scov) %>% mutate(row = row_number()) %>% filter(row == 1 )
df = data.frame(qname = no.self.hit.qnames) %>%
  left_join(anots) %>%
  mutate(hhrepid = qname %in% no.self.hits.and.with.repeats.qnames)


hhsearch.filename = "/Users/bognasmug/MGG Dropbox/Projects/divRBP/phage-pp-workdir-refseq-hhblits-DEBUG/table-hhr-hhsearch-DEBUG-3160.csv"
hhblits.filename = "/Users/bognasmug/MGG Dropbox/Projects/divRBP/phage-pp-workdir-refseq-hhblits-DEBUG/table-hhr-hhblits-DEBUG-3160.csv"
hhsearch.table = read.csv(hhsearch.filename, 
                     sep = ",",
                     header = TRUE) %>%
  mutate(hhrepid = hhrepid == "True")

hhblits.table = read.csv(hhblits.filename, 
                     sep = ",",
                     header = TRUE)  %>%
  mutate(hhrepid = hhrepid == "True")

library(dplyr)
library(ggplot2)

self.hits.hhsearch = hhsearch.table %>%
  filter(qname == sname) %>%
  distinct(qname, sname) %>%
  nrow()

self.hits.hhblits = hhblits.table %>%
  filter(qname == sname) %>%
  distinct(qname, sname) %>%
  nrow()

# note that sometimes we have more than one hit for the same pair!!
table.both.summary = hhsearch.table %>% 
  distinct(qname, sname, eval, hhrepid) %>%
  mutate(method = "hhsearch") %>%
  rbind(hhblits.table %>% distinct(qname, sname, eval, hhrepid) %>% mutate(method = "hhblits")) %>%
  group_by(qname, sname, hhrepid) %>%
  summarise(method = paste(sort(unique(method)), collapse = "&", sep = ""),
            e = if_else(mean(eval) < 10^(-1), "<10^(-1)", ">=10^(-1)"),
            hhrepid = any(hhrepid)) %>%
  ungroup() 

ggplot(data = table.both.summary %>% group_by(method,e, hhrepid) %>% summarise(n = n()) %>% as.data.frame()) +
  geom_col(aes(x = method, y = n, fill = e, alpha = hhrepid), position = position_dodge2()) + ylab("Num pairs") 


pairs.absent.in.hhblits.raw = table.both.summary %>% 
  filter(method == "hhsearch") %>%
  inner_join(hhsearch.table) %>%
  rowwise() %>%
  mutate(hit.length = qend - qstart + 1,
         log10eval = log10(eval)) %>%
  distinct(qname, sname, hit.length, pident, prob, log10eval, bitscore) 

pairs.absent.in.hhblits = pairs.absent.in.hhblits.raw%>%
  tidyr::gather(key = "stat", value = "value", hit.length, pident, prob, log10eval, bitscore)

  


ggplot(pairs.absent.in.hhblits) +
  geom_histogram(aes(x=value, fill = stat)) + 
  facet_grid(. ~ stat, scales = "free") +
  ggtitle("Hits missing from hhblits")

pairs.absent.in.hhblits %>% filter(stat == "log10eval.times100") %>% pull(value) %>% summary()


hhsearch.table.with.suffixes = hhsearch.table
colnames(hhsearch.table.with.suffixes)<-paste(colnames(hhsearch.table.with.suffixes),"hhsearch",sep="_")
hhblits.table.with.suffixes = hhblits.table
colnames(hhblits.table.with.suffixes)<-paste(colnames(hhblits.table.with.suffixes),"hhblits",sep="_")

# note that sometimes we have more than one hit for the same pair!!
# so when we join by qname, sname we may get multiple pairs, some of which  are different
# so we overestmate the differences between hits
discreapacies.in.common.hits =  table.both.summary %>% 
  filter(method == "hhblits&hhsearch") %>%
  inner_join(hhsearch.table.with.suffixes %>% rename(qname = qname_hhsearch,
                                                     sname = sname_hhsearch), by = c("qname", "sname")) %>%
  inner_join(hhblits.table.with.suffixes %>% rename(qname = qname_hhblits,
                                                    sname = sname_hhblits), by = c("qname", "sname")) %>%
  mutate(diff.prob = prob_hhsearch - prob_hhblits,
         diff.log10.eval = log10(eval_hhsearch) - log10(eval_hhblits),
         diff.hit.length = qend_hhsearch - qstart_hhsearch - qend_hhblits + qstart_hhblits,
         diff.bitscore = bitscore_hhsearch - bitscore_hhblits,
         diff.pident = pident_hhsearch - pident_hhblits,
         e_hhblits = if_else(eval_hhblits < 10^(-3), "<10^(-3)", ">=10^(-3)"),
         e_hhsearch = if_else(eval_hhsearch < 10^(-3), "<10^(-3)", ">=10^(-3)")) %>%
  mutate(roughly.same.hit = if_else(condition = 
           abs(diff.prob) < 5 & 
           (abs(diff.log10.eval) < 0.2 | log10(eval_hhsearch) < -10) & 
           abs(diff.hit.length) < 5 &
           abs(diff.bitscore) < 5 &
           abs(diff.pident) < 3,
           true = "roughly.same.hit",
           false = if_else(condition = abs(diff.prob) > 30 | 
                            (abs(diff.log10.eval) > 2 & log10(eval_hhsearch) > -5) | 
                            abs(diff.hit.length) > 50 |
                            abs(diff.bitscore) > 100 |
                            abs(diff.pident) > 20,
                          true = "very.different.hit",
                          false = "somewhat.different.hit")))



ggplot(discreapacies.in.common.hits #%>%
         #mutate(log10e_hhsearch = log10(eval_hhsearch)
       ,
       aes(x = prob_hhsearch, y = pident_hhsearch, col = roughly.same.hit)) +
  geom_point() +
  facet_grid(. ~ e_hhsearch, scales = "free") +
  ggtitle(" Hits found by both methods")

ggplot(discreapacies.in.common.hits %>%
         group_by(roughly.same.hit, e_hhsearch, hhrepid) %>%
         summarise(n = n())) +
  geom_col(aes(x = roughly.same.hit, y = n, fill = e_hhsearch, alpha = hhrepid), position = position_dodge2()) + ylab("Num pairs") + xlab("") +
  ggtitle("Are hits between hhsearch and hhblits similar?")

discreapacies.in.common.hits.long.rsh = discreapacies.in.common.hits %>%
  distinct(qname, sname, diff.prob, diff.log10.eval, diff.hit.length, diff.bitscore, diff.pident) %>%
  tidyr::gather(key = "stat", value = "value", -qname, -sname)


ggplot(discreapacies.in.common.hits.long.rsh) +
  geom_histogram(aes(x = value, fill = stat)) +
  facet_grid(. ~ stat, scales = "free") +
  ggtitle("Discrepacies between somewhat different hits found by both methods")#+ +


discreapacies.in.common.hits.long.high.prob= discreapacies.in.common.hits %>%
  filter(prob_hhsearch > 95) %>%
  distinct(qname, sname, diff.prob, diff.log10.eval, diff.hit.length, diff.bitscore, diff.pident) %>%
  tidyr::gather(key = "stat", value = "value", -qname, -sname)


ggplot(discreapacies.in.common.hits.long.high.prob) +
  geom_histogram(aes(x = value, fill = stat)) +
  facet_grid(. ~ stat, scales = "free") +
  xlim(c(-100,100)) +
  ggtitle("Discrepacies between somewhat different hits found by both methods [prob > 95%]")#+ +

discreapacies.in.common.hits.long.low.prob= discreapacies.in.common.hits %>%
  filter(prob_hhsearch < 95) %>%
  distinct(qname, sname, diff.prob, diff.log10.eval, diff.hit.length, diff.bitscore, diff.pident) %>%
  tidyr::gather(key = "stat", value = "value", -qname, -sname)


ggplot(discreapacies.in.common.hits.long.low.prob) +
  geom_histogram(aes(x = value, fill = stat)) +
  facet_grid(. ~ stat, scales = "free") +
  xlim(c(-100,100)) +
  ggtitle("Discrepacies between somewhat different hits found by both methods [prob < 95%]")#+ +


very.different.hits = discreapacies.in.common.hits %>%
  filter(roughly.same.hit == "very.different.hit")


# how often hhsearch gives a high probability when it is low by hhblits?
aa = very.different.hits %>%
  filter(prob_hhsearch > 90 &
           prob_hhblits < 50)
bb = pairs.absent.in.hhblits.raw %>%
  filter(prob > 90) 
# num problematic cases
nrow(aa) + nrow(bb)

# what are the worrying cases?
aa = discreapacies.in.common.hits %>%
  filter(bitscore_hhsearch > 100 &
           prob_hhblits < 50)
bb = pairs.absent.in.hhblits.raw %>%
  filter(bitscore > 100) 
cc = hhsearch.table %>%
  filter(bitscore > 100)
# num problematic cases
(nrow(aa) + nrow(bb))/nrow(cc)


aa = discreapacies.in.common.hits %>%
  filter(bitscore_hhsearch > 200 &
           prob_hhblits < 50)
bb = pairs.absent.in.hhblits.raw %>%
  filter(bitscore > 200) 
(nrow(aa) + nrow(bb))/nrow(cc)
  
  #geom_jitter(aes(x = stat, y = value, fill = stat))