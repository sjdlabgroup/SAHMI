library(optparse)
library(stringr)
library(ShortRead)
library(dplyr)
library(Matrix)
library(data.table)
library(tidyverse)

option_list = list(
  make_option(c("--fa1"), action="store", help = "path to first fa1"),
  make_option(c("--fa2"), action="store", help = "path to first fa2"),
  make_option(c("--taxa"), action="store", help = "tsv containing taxa to extract"),
  make_option(c("--out_path"), action="store", help = "output path to save files"),
  make_option(c("--sample_name"), action="store", help = "sample name"),
  make_option(c("--kraken_report"), action="store", help = "path to kraken report"),
  make_option(c("--mpa_report"), action="store", help = "path to mpa report"),
  make_option(c("--nbarcodes"), action="store", default = 4, help = "minimum number of barcodes on which to run sckmer correlation"),
  make_option(c("--p"), action="store", default = 0.05, help = "sckmer correlations p-value threshold to keep taxa"),
  make_option(c("--cb_len"), action="store", default = 16, help = "cell barcode length"),
  make_option(c("--umi_len"), action="store", default = 10, help = "umi length"),
  make_option(c("--nFilter"), action="store", default = 130, help = "filter reads with >n of one nucleotide"),
  make_option(c("--movingAverage"), action="store", default = 25, help = "window for moving average to find barcode cutoff on elbow plot")
)
opt = parse_args(OptionParser(option_list = option_list))

# get barcodes, umis, and tax-ids
print(paste('Started extracting barcode data from fastq files for', opt$sample_name))

reads = readFasta(opt$fa1) 
filter <- polynFilter(threshold=opt$nFilter, nuc=c("A","T","G","C")) %>% compose()
reads = reads[filter(reads)]
sequences = sread(reads)
headers = ShortRead::id(reads)
barcode = substr(sequences, 1, opt$cb_len)
umi = substr(sequences, opt$cb_len+1, opt$cb_len + opt$umi_len) 
taxid = gsub('.*taxid\\|', '', headers) 

taxa = read.table(opt$taxa)

kr = read.delim(opt$kraken_report, header = F)
kr = kr[-c(1:2), ] %>% mutate(V8 = trimws(V8)) %>% mutate(V8 = str_replace_all(V8, '[^[:alnum:]]', '_'))
kr.names = kr$V8 %>% as.character() %>% trimws() %>% str_replace_all('_', ' ') %>% str_squish()
mpa = read.delim(opt$mpa_report, header = F)
mpa$taxid = NA

for(i in 2:nrow(mpa)){
  t_names = mpa[i,1] %>% as.character() %>% 
    strsplit('\\|') %>% 
    unlist() %>% 
    str_remove('.*__') %>% 
    str_replace_all('[^[:alnum:]]', '_') 
  mpa$taxid[i] = paste0('*', paste(kr$V7[match(t_names, kr$V8)], collapse = '*'), '*')
}

x = paste0('\\*', taxa$V1, '\\*', collapse = '|')
ii = str_which(mpa$taxid, x)
x = paste0('\\*', taxa$V1, '\\*.*')
ii = sapply(x, function(y) str_extract(mpa$taxid[ii], y) %>% strsplit('\\*') %>% unlist() %>% unique()) %>% unlist() %>% unique()

i = which(taxid %in% ii)

barcode = barcode[i]
umi = umi[i]
taxid = taxid[i]

bc = cbind(barcode, umi, taxid)
s.bc = bc %>% unique() %>% data.frame()
s.bc$umi = 1
s.bc = s.bc %>% group_by(barcode, taxid) %>% summarize(umi = sum(umi), .groups = 'keep') %>% arrange(desc(umi))
rm(bc)

write.table(s.bc, file = paste0(opt$out_path, opt$sample_name, '.all.barcodes.txt'),
            quote = F, sep='\t', row.names = F, col.names = T)

s.mat = sparseMatrix(as.integer(s.bc$barcode), as.integer(s.bc$taxid), x=s.bc$umi)
colnames(s.mat) = levels(s.bc$taxid)
rownames(s.mat) = levels(s.bc$barcode)
s.mat = t(s.mat)

print(paste('Started classifying reads for', opt$sample_name))

# count parent classifications for each read
df = list()
counter = 0
for(i in 1:nrow(s.mat)){
  print(i)
  if(length(which(kr$V7 == rownames(s.mat)[i])) == 0){next}
  
  tax = mpa[which(kr$V7 == rownames(s.mat)[i]), 1] %>%
    as.character() %>% 
    strsplit('\\|') %>% 
    unlist() %>% 
    data.frame() %>% 
    separate(col = 1, c('rank', 'name'), sep = '__') %>% 
    subset(rank %in% c('k','p','c','o','f','g','s')) %>%
    mutate(name = str_replace_all(name, '[^[:alnum:]]', ' ') %>% str_squish())
  
  if(nrow(tax) == 0){next}
  
  id = c()
  for(j in 1:nrow(tax)){id[j] = kr$V7[kr.names == tax$name[j]]}
  tax$id = id
  
  row = s.mat[i,]
  row = row[row>0]
  for(j in 1:nrow(tax)){
    counter = counter + 1
    tax2 = tax[1:j, ]
    df[[counter]] = tibble(barcode = names(row), 
                           counts = row,
                           taxid = tax2$id[j], 
                           rank = tax2$rank[j], 
                           kingdom = tax2$name[tax2$rank == 'k'][1],
                           phylum = tax2$name[tax2$rank == 'p'][1],
                           class = tax2$name[tax2$rank == 'c'][1],
                           order = tax2$name[tax2$rank == 'o'][1],
                           family = tax2$name[tax2$rank == 'f'][1],
                           genus = tax2$name[tax2$rank == 'g'][1],
                           species = tax2$name[tax2$rank == 's'][1])
  }
}

df = rbindlist(df, use.names = T)

df = df %>% 
  group_by(barcode, taxid, rank, kingdom, phylum, class, order, family, genus, species) %>% 
  summarize(counts = sum(counts), .groups = 'keep') %>% 
  arrange(desc(counts))

# # remove empty barcodes
# moving.average <- function(x, n = opt$movingAverage){stats::filter(x, rep(1 / n, n), sides = 2)}
# bc.depth = colSums(s.mat) %>% sort(decreasing = T)
# slope = bc.depth %>% moving.average(n = 50) %>% diff(na.rm = T)
# n_bc = which(abs(slope) < 10^-3)[1]
# 
# s.mat = s.mat[, names(bc.depth)[1:n_bc]]
# ind = which(rowSums(s.mat) == 0)
# if(length(ind)>0){s.mat = s.mat[-ind, ]}

# save
write.table(df, file = paste0(opt$out_path, opt$sample_name, '.counts.txt'),
            quote = F, sep='\t', row.names = F, col.names = T)

print('Finished')