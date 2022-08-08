library(optparse)
library(tidyverse)
library(dplyr)
library(stringr)

option_list = list(
  make_option(c("--sample_name"), action="store", help = "sample name"),
  make_option(c("--fq"), action="store", help = "path to fastq file"),
  make_option(c("--kraken_report"), action="store", help = "path to standard kraken report"),
  make_option(c("--mpa_report"), action="store", help = "path to standard kraken mpa report"),
  make_option(c("--out_path"), action="store", help = "output path"),
  make_option(c("--keep_original"), action="store", default=T, help ="delete original fastq file? T/F"),
  make_option(c("--ntaxid"), action="store", default=7000, help = "number of taxids to extract at a time")
)
opt = parse_args(OptionParser(option_list = option_list))

kr = read.delim(opt$kraken_report, header = F)
kr = kr[-c(1:2), ]
mpa = read.delim(opt$mpa_report, header = F)
n = str_which(mpa$V1, 'k__Bacteria|k__Fungi|k__Viruses')
taxid = kr$V7[n]
taxid.list = split(taxid, ceiling(seq_along(taxid)/opt$ntaxid))

if(file.exists(paste0(opt$out_path, opt$sample_name, '_line_numbers.txt'))){
  system(paste0('rm ', opt$out_path, opt$sample_name, '_line_numbers.txt'))
}

for(i in 1:length(taxid.list)){
  print(paste('Finding reads', i, '/', length(taxid.list)))
  
  taxid = paste0("taxid|", taxid.list[[i]], collapse = "\\|") 
  taxid = paste0("'", taxid, "'")
  str = paste0("grep -wn ", taxid, " ", opt$fq, " | grep -Eo '^[^:]+' >> ", opt$out_path, opt$sample_name, "_line_numbers.txt")
  system(str)
}

h = read.delim(paste0(opt$out_path, opt$sample_name, '_line_numbers.txt'), header = F)
r = h$V1+1
d = data.frame(h=h$V1,r=r) %>% rownames_to_column('n') %>% pivot_longer(-n)
write.table(d$value, file = paste0(opt$out_path, opt$sample_name, '_line_numbers.txt'), row.names = F, col.names = F)

print('Extracting reads')
str = paste0("awk 'NR==FNR{ a[$1]; next }FNR in a' ", opt$out_path, opt$sample_name, '_line_numbers.txt ', opt$fq, " > ", opt$out_path, opt$sample_name, ".fa")
system(str)

str = paste0("sed -i 's/@/>/g' ", opt$out_path, opt$sample_name, '.fa')
system(str)

str = paste0(opt$out_path, opt$sample_name, '_line_numbers.txt')
system(paste('rm', str))

if(opt$keep_original == F){
  system(paste('rm', fq))
}

print('Done')

