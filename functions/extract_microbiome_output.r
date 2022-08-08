library(optparse)
library(tidyverse)
library(dplyr)
library(stringr)

option_list = list(
  make_option(c("--sample_name"), action="store", help = "sample name"),
  make_option(c("--output_file"), action="store", help = "path to kraken output file"),
  make_option(c("--kraken_report"), action="store", help = "path to kraken report"),
  make_option(c("--mpa_report"), action="store", help = "path to standard kraken mpa report"),
  make_option(c("--out_path"), action="store", help = "output path"),
  make_option(c("--keep_original"), action="store", default=T, help ="delete original fastq file? T/F"),
  make_option(c("--ntaxid"), action="store", default = 8000, help = "number of taxids to extract at a time")
)
opt = parse_args(OptionParser(option_list = option_list))

kr = read.delim(opt$kraken_report, header = F)
kr = kr[-c(1:2), ]
mpa = read.delim(opt$mpa_report, header = F)
n = str_which(mpa$V1, 'k__Bacteria|k__Fungi|k__Viruses')
taxid = kr$V7[n]
taxid.list = split(taxid, ceiling(seq_along(taxid)/opt$ntaxid))

if(file.exists(paste0(opt$out_path, opt$sample_name, '.microbiome.output.txt'))){
  system(paste0('rm ', opt$out_path, opt$sample_name, '.microbiome.output.txt'))
}

for(i in 1:length(taxid.list)){
  print(paste('Extracting output data', i, '/', length(taxid.list)))
  
  taxid = paste0("(taxid ", taxid.list[[i]], ")", collapse = "\\|")
  taxid = paste0("'", taxid, "'")
  str = paste0("grep -w ", taxid, " ", opt$output_file, " >> ", opt$out_path, opt$sample_name, ".microbiome.output.txt")
  system(str)
}

if(opt$keep_original == F){
  system(paste('rm', opt$output_file))
}

print('Done')
