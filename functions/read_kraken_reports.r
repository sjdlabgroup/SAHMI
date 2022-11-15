read_kraken_reports = function(files, sample_names = NULL, study_name = NULL, min_reads = 2, min_uniq = 2){
  require(data.table)
  require(tibble)
  
  if(is.null(sample_names)){sample_names = files}
  if(is.null(study_name)){study_name = NA}
  if(length(study_name) == 1){study_name = rep(study_name, length(files))}
  
  df = list()
  n = 0
  for(i in 1:length(files)){
    if(round(i/length(files)*100, 2) > n){n = round(i/length(files)*100, 2); cat(paste0('\r',n,'% done   '))}
    x = read.delim(files[i], header = F)
    x$V8 = trimws(x$V8)
    total_reads = x$V2[1] + x$V2[2]
    n_microbiome_reads = sum(x$V2[x$V8 %in% c('Bacteria', 'Fungi', 'Viruses')])
    df[[i]] = data.frame(study = study_name[i], sample = sample_names[i],
                         rank = x$V6, taxid = x$V7, name = x$V8, 
                         reads = x$V2, min = x$V4, uniq = x$V5, 
                         rpm = x$V2/total_reads*10^6,
                         rpmm = x$V2/n_microbiome_reads*10^6)
  }
  df = rbindlist(df) %>% tibble()
  cat('\n')
  
  df = subset(df, reads >= min_reads & uniq >= min_uniq) 
  df
}