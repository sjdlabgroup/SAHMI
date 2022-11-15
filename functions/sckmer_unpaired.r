library(optparse)
library(tidyverse)
library(data.table)
library(dplyr)
library(stringr)
library(ShortRead)

option_list = list(
  make_option(c("--sample_name"), action="store", help = "sample name"),
  make_option(c("--fa1"), action="store", help = "path to fasta 1 file"),
  make_option(c("--microbiome_output_file"), action="store", help = "path to microbiome output file"),
  make_option(c("--cb_len"), action="store", default=16, help = "nucleutide length of cell barcodes"),
  make_option(c("--umi_len"), action="store", default=10, help = "nucleutide length of umis"),
  make_option(c("--kraken_report"), action="store", help = "path to standard kraken report"),
  make_option(c("--mpa_report"), action="store", help = "path to standard kraken mpa report"),
  make_option(c("--ranks"), action="store", default = c('G', 'S'), help = "taxa ranks to analyze"),
  make_option(c("--host"), action="store", default = 9606, help = "host taxid to exclude"),
  make_option(c("--min_frac"), action="store", default = 0.5, help = "minimum fraction of kmers directly assigned to taxid to use read"),
  make_option(c("--kmer_len"), action="store", default = 35, help = "Kraken kmer length"),
  make_option(c("--nsample"), action="store", default = 1000, help = "Max number of barcodes to sample per taxa"),
  make_option(c("--out_path"), action="store", help = "output path")
)
opt = parse_args(OptionParser(option_list = option_list))

# opt = list(fa1 = '/projects/sd948/bassel/datasets/liao-covid/kraken2uniq/SRR11537951_1.fa',
#            fa2 = '/projects/sd948/bassel/datasets/liao-covid/kraken2uniq/SRR11537951_2.fa',
#            microbiome_output_file = '/projects/sd948/bassel/datasets/liao-covid/kraken2uniq/SRR11537951.microbiome.output.txt',
#            cb_len = 16,
#            umi_len = 10,
#            kraken_report = '/projects/sd948/bassel/datasets/liao-covid/kraken2uniq/SRR11537951.kraken.report.txt',
#            mpa_report = '/projects/sd948/bassel/datasets/liao-covid/kraken2uniq/SRR11537951.kraken.report.mpa.txt',
#            ranks = c('G', 'S'),
#            host = 9606, 
#             min_frac = 0.5,
#            kmer_len = 35,
#            nsample = 1000)

reads1 = readFasta(opt$fa1)
sequences1 = sread(reads1) 
# reads2 = readFasta(opt$fa2)
# sequences2 = sread(reads2) 

headers = ShortRead::id(reads1)
barcode = substr(sequences1, 1, opt$cb_len)
umi = substr(sequences1, opt$cb_len+1, opt$cb_len + opt$umi_len)
taxid = gsub('.*taxid\\|', '', headers)
id = gsub('\\s.*', '', headers)

# headers2 = ShortRead::id(reads2)
# id2 = gsub('\\s.*', '', headers2)

# i = intersect(id, id2)
# ii = which(id %in% i)
# barcode = barcode[ii]
# umi = umi[ii]
# taxid = taxid[ii]
# id = id[ii]
# sequences1 = sequences1[ii]

# ii = which(id2 %in% i)
# sequences2 = sequences2[ii]

kr = read.delim(opt$kraken_report, header = F)
kr = kr[-c(1:2), ] %>% mutate(V8 = trimws(V8)) %>% mutate(V8 = str_replace_all(V8, '[^[:alnum:]]', '_'))
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

microbiome_output_file = read.delim(opt$microbiome_output_file, header = F)

microbiome_output_file = microbiome_output_file %>% 
  select(-V1) %>% 
  separate(V3, into = c('name', 'taxid'), sep = '\\(taxid') %>% 
  mutate(taxid = str_remove(taxid, '\\)') %>% trimws(),
         name = trimws(name))

tx = kr$V7[kr$V6 %in% opt$ranks] %>% setdiff(opt$host)
tx = microbiome_output_file$taxid[microbiome_output_file$taxid %in% tx] %>% unique()


barcode_kmer = list()
counter = 0 
for(taxa in tx){
  counter = counter + 1
  cat(paste('\r', 'taxa processed:', round(counter/length(tx)*100, 3), '%'))
  
  lin = str_subset(mpa$taxid, paste0('\\*', taxa, '\\*')) %>% 
    str_extract(paste0('\\*', taxa, '\\*.*')) %>%
    str_remove('^\\*') %>% 
    str_remove('\\*$') %>% 
    str_split('\\*') %>% 
    unlist() %>% 
    as.numeric() %>% 
    unique()
  
  full.lin = str_subset(mpa$taxid, paste0('\\*', taxa, '\\*')) %>% 
    str_remove('^\\*') %>% 
    str_remove('\\*$') %>% 
    str_split('\\*') %>% 
    unlist() %>% 
    as.numeric() %>% 
    unique()
  
  out = subset(microbiome_output_file, taxid %in% lin) %>% separate(V5, into = c('r1', 'r2'), sep = '\\|\\:\\|') 
  out$r1[str_which(out$r1, paste0(' ', opt$host, ':'))] = NA
  out$r2[str_which(out$r2, paste0(' ', opt$host, ':'))] = NA
  out$r1[out$r1 == ''] = NA; out$r2[out$r2 == ''] = NA
  out = subset(out, !is.na(r1) | !is.na(r2))
  out$r1 = trimws(out$r1)
  out$r2 = trimws(out$r2)
  
  
  if(nrow(out) == 0){next}
  
  i = which(id %in% out$V2)
  seq = data.frame(r1 = sequences1[i] %>% as.character())
  barcode.x = barcode[i]
  umi.x = umi[i]
  
  tax.df = list()
  if(nrow(out) > opt$nsample){n = sample(nrow(out), opt$nsample)} else {n=1:nrow(out)}
  counter2 = 0
  
  for(i in n){
    # cat(paste('\r', 'barcodes processed:', round(counter2/length(n)*100, 3), '%   '))
    for(mate in c('r1', 'r2')){
      
      r = data.frame(pos = out[[mate]][i] %>% strsplit('\\s') %>% unlist()) %>% 
        separate(pos, into = c('taxid', 'nkmer'), sep = ':', convert = T) 
      # if(any(r$taxid %in% c(0, full.lin) == F)){next}
      
      r = r %>% 
        mutate(fkmer = nkmer/sum(nkmer),
               nt_start = cumsum(nkmer) - nkmer + 1,
               nt_end = cumsum(nkmer) + opt$kmer_len - 1) %>% 
        mutate(nt_len = nt_end - nt_start + 1) %>% 
        ungroup() %>% 
        # subset(taxid == taxa)
        subset(taxid %in% c(0,full.lin))
      
      if(sum(r$fkmer) < opt$min_frac){next}
      
      counter2 = counter2 + 1
      
      if(nrow(r) > 0){
        kmer = c()
        for(k in 1:nrow(r)){
          for(m in 1:r$nkmer[k]){
            kmer = c(kmer, substr(seq[[mate]][i], r$nt_start[k] + m - 1, r$nt_start[k] + m + opt$kmer_len - 2))
          }
        }
        tax.df[[counter2]] = data.frame(barcode = barcode.x[i], taxid = taxa, k = kmer, n = sum(r$nt_len[r$taxid %in% lin]))
      } else {tax.df[[counter2]] = data.frame(barcode = barcode.x[i], taxid = taxa, k = NA, n = NA)}
    }
  }
  
  if(length(tax.df) == 0){next}
  
  tax.df = bind_rows(tax.df) %>%
    tibble() %>% 
    subset(!is.na(k)) %>%
    group_by(barcode, taxid) %>%
    summarize(kmer = length(k),
              uniq = length(unique(k)),
              .groups = 'keep')
  
  barcode_kmer[[counter]] = tax.df
  # cat('\n')
}

barcode_kmer = rbindlist(barcode_kmer)

# c = barcode_kmer %>%
#   subset(kmer > 1) %>%
#   group_by(taxid) %>%
#   mutate(nn = n()) %>%
#   subset(nn > 2) %>%
#   group_by(taxid) %>%
#   summarize(r = cor.test(kmer, uniq, method = 'spearman', use = 'pairwise.complete')$estimate,
#             p = cor.test(kmer, uniq, method = 'spearman', use = 'pairwise.complete')$p.value) %>%
#   mutate(p = p.adjust(p)) 

write.table(barcode_kmer, file = paste0(opt$out_path, opt$sample_name, '.sckmer.txt'), quote = F, row.names = F)

paste('Done')



