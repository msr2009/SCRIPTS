read.fused.data <- function(input_table){
  colnames = c('readID', 'sequence', 'length', 'paired_mismatch', 'paired_match', 'gap', 'unresolvable', 'mutrun', 
               'read_avgquality', 'read1_chastity', 'read1_Ncount', 'read2_avgquality', 'read2_chastity', 'read2_Ncount',
               'read1_mut_quality', 'read2_mut_quality', 'fused_quality')
  dat <- read.table(input_table, sep='\t', col.names = colnames)
  return(dat)
}