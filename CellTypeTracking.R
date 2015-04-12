setwd("E:/Dropbox/workingspace/CellTypeTracker")
#assume input file are pileup files

#function that generate all the unique combination of alleles(genotype) of given ploidy
generate.genotype.permutation=function(ploidy)
{
  #alelles=c('A','C','G','T')
  return(unique(t(combn(c(rep('A',ploidy),rep('C',ploidy),rep('G',ploidy),rep('T',ploidy)),ploidy))))
  
}

site.likelihood=function(seq,qual,ploidy)
{
  genotype.list=generate.genotype.permutation(ploidy)
}
overall.likelihood=function()
{
  
}



ploidy.estimate=function()
{
  
}
