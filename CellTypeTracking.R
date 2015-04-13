#setwd("E:/Dropbox/workingspace/CellTypeTracker")
setwd("C:/Users/PC/Dropbox/workingspace/CellTypeTracker")
#assume input file are pileup files

#allele frequency prior
allele.A.freq=0.25
allele.C.freq=0.25
allele.G.freq=0.25
allele.T.freq=0.25

countCharOccurrences <- function(char, s) {
  s2 <- gsub(char,"",s)
  return (nchar(s) - nchar(s2))
}
#prior genotype likelihood
genotype.prior=function(genotype)
{
  numA=sum(countCharOccurrences("A",genotype))
  numC=sum(countCharOccurrences("C",genotype))
  numG=sum(countCharOccurrences("G",genotype))
  numT=sum(countCharOccurrences("T",genotype))
  total=length(genotype)
  prior=(factorial(total)/(factorial(numA)*factorial(numC)*factorial(numG)*factorial(numT)))*(allele.A.freq^numA)*(allele.C.freq^numC)*(allele.G.freq^numG)*(allele.T.freq^numT)
  return(prior)
}

#conditional prob
conditional.prob=function(base,genotype,error.status)
{
  numA=sum(countCharOccurrences("A",genotype))
  numC=sum(countCharOccurrences("C",genotype))
  numG=sum(countCharOccurrences("G",genotype))
  numT=sum(countCharOccurrences("T",genotype))
  num=c(numA,numC,numG,numT)
  names(num)=c("A","C","T","G")
  total=length(genotype)
  if(error.status==1)
  {
    if(num[base]!=0)
    {
      
    }
  }
  else
  {
      
  }
}
#function that generate all the unique combination of alleles(genotype) of given ploidy
generate.genotype.permutation=function(ploidy)
{
  #alelles=c('A','C','G','T')
  return(unique(t(combn(c(rep('A',ploidy),rep('C',ploidy),rep('G',ploidy),rep('T',ploidy)),ploidy))))
  
}

#multiply all the observed reads for a single site
single.observation.likelihood(base,quality,alpha,geno1,geno2)
{
  #if error happened
  likelihood=10^(-quality/10)*((1-alpha)*conditional.prob(base,geno1,1)+alpha*conditional.prob(base,geno2,1))
    
  #if error not happened
  likelihood=likelihood+(1-10^(-quality/10))*((1-alpha)*conditional.prob(base,geno1,0)+alpha*conditional.prob(base,geno2,0))
  return(likelihood)
}

#calculate genotype likelihood for single site, ploidy1 for celltype of interest, ploidy2 for the rest smoothies
site.likelihood=function(seq,qual,ploidy1,ploidy2,alpha)
{
  genotype.list1=generate.genotype.permutation(ploidy1)
  genotype.list2=generate.genotype.permutation(ploidy2)
  likelihood=0;
  for(geno1 in genotype.list1)
  {
    tmplike2=0;
    for(geno2 in genotype.list2)
    {
      loglik=0;
      for( i in 1:length(seq))#assume seq is of format c("A","A","A","A")
      {
        loglik=loglik+single.observation.likelihood(seq[i],qual[i],alpha,geno1,geno2)
      }
      tmplike2=tmplike2+e^loglik*genotype.prior(geno2)
    }
    likelihood=likelihood+tmplike2*genotype.prior(geno1)
  }
  return(likelihood)
}
overall.likelihood=function()
{
  
}



ploidy.estimate=function()
{
  
}
