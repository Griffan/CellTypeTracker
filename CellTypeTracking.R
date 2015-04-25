#setwd("E:/Dropbox/workingspace/CellTypeTracker")
setwd("C:/Users/PC/Dropbox/workingspace/CellTypeTracker")
#assume input file are pileup files

#allele frequency prior
allele.A.freq=0.25
allele.C.freq=0.25
allele.G.freq=0.25
allele.T.freq=0.25

countCharOccurrences <- function(char, s) {
  #s2 <- gsub(char,"",s)
  #return (nchar(s) - nchar(s2))
  return(s %in% char)
}
#prior genotype likelihood
genotype.prior=function(genotype,af)
{
  numA=af["A"]
  numC=af["C"]
  numG=af["G"]
  numT=af["T"]
  total=length(genotype)
  prior=(factorial(total)/(factorial(numA)*factorial(numC)*factorial(numG)*factorial(numT)))*(allele.A.freq^numA)*(allele.C.freq^numC)*(allele.G.freq^numG)*(allele.T.freq^numT)
  return(prior)
}

#conditional prob
conditional.prob=function(base,genotype,error.status,af)
{
#   base="T"
#   genotype=c("A","A","C","G")
#   error.status=1
#   numA=sum(countCharOccurrences("A",genotype))
#   numC=sum(countCharOccurrences("C",genotype))
#   numG=sum(countCharOccurrences("G",genotype))
#   numT=sum(countCharOccurrences("T",genotype))
#   num=c(numA,numC,numG,numT)
  num=af
 
  total=length(genotype)
  #k=sum(c(1,1,1,1)[c(num["A"]>0,num["C"]>0,num["G"]>0,num["T"]>0)])# how many kinds of alleles
  if(error.status==1)
  {
      #return(num["A"]/total*1/3+num["C"]/total*1/3+num["G"]/total*1/3+num["T"]/total*1/3-num[base]/total*1/3)
      return(sum(num/total*1/3)-num[base]/total*1/3)
  }
  else# no error event
  {
      
    return(num[base]/total)
  }
}
#function that generate all the unique combination of alleles(genotype) of given ploidy
generate.genotype.permutation=function(ploidy)
{
  #alelles=c('A','C','G','T')
  return(unique(t(combn(c(rep('A',ploidy),rep('C',ploidy),rep('G',ploidy),rep('T',ploidy)),ploidy))))
  
}

#marginal likelihood over all error status
single.observation.likelihood=function(base,quality,alpha,geno1,geno2,af1,af2)
{
  #print(quality)
  quality=strtoi(charToRaw(quality),16L)
  #print(quality)
  #if error happened
  likelihood=10^(0-quality/10)*((1-alpha)*conditional.prob(base,geno1,1,af1)+alpha*conditional.prob(base,geno2,1,af2))
    
  #if error not happened
  likelihood=likelihood+(1-10^(0-quality/10))*((1-alpha)*conditional.prob(base,geno1,0,af1)+alpha*conditional.prob(base,geno2,0,af2))
  #print(likelihood)
  return(log(likelihood))
}

#calculate genotype likelihood for single site, ploidy1 for celltype of interest, ploidy2 for the rest smoothies
site.likelihood=function(seq,qual,ploidy1,ploidy2,alpha)
{
  genotype.list1=generate.genotype.permutation(ploidy1)
  genotype.list2=generate.genotype.permutation(ploidy2)
  likelihood=0;
  for(p in 1:nrow(genotype.list1))
  {
    geno1=genotype.list1[p,]
    numA=sum(countCharOccurrences("A",geno1))
    numC=sum(countCharOccurrences("C",geno1))
    numG=sum(countCharOccurrences("G",geno1))
    numT=sum(countCharOccurrences("T",geno1))
    af1=c(numA,numC,numG,numT)
    names(af1)=c("A","C","G","T")
    tmplike2=0;
    for(q in 1:nrow(genotype.list2))
    {
      geno2=genotype.list2[q,]
      numA=sum(countCharOccurrences("A",geno2))
      numC=sum(countCharOccurrences("C",geno2))
      numG=sum(countCharOccurrences("G",geno2))
      numT=sum(countCharOccurrences("T",geno2))
      af2=c(numA,numC,numG,numT)
      names(af2)=c("A","C","G","T")
      loglik=0;
      for( i in 1:length(seq))#assume seq is of format c("A","A","A","A")
      {
        loglik=loglik+single.observation.likelihood(seq[i],qual[i],alpha,geno1,geno2,af1,af2)
      }
      #print(loglik)
      tmplike2=tmplike2+exp(loglik)*genotype.prior(geno2,af2)
    }
    #print(tmplike2)
    likelihood=likelihood+tmplike2*genotype.prior(geno1,af1)
  }
  #print(likelihood)
  return(log(likelihood))
}

#data has the pileup format: Chrom\tPos\tSNPID\tDepth\tBases\tQualities
overall.likelihood=function(data,ploidy1,ploidy2,alpha)
{
  loglikelihood=0;
  for(i in 1:length(data[,1]))
  {
    seq=as.vector(data[i,5])
    seq=strsplit(toupper(seq),"")[[1]]
    qual=as.vector(data[i,6])
    qual=strsplit(qual,"")[[1]]
    #print(seq[1])
    loglikelihood=loglikelihood+site.likelihood(seq,qual,ploidy1,ploidy2,alpha)
  }
  return(0-loglikelihood)
}


pileup=read.table(gzfile("C:/Users/PC/Desktop/HG00553.large_insert.Pileup.gz"),header=F,comment.char="",quote="\"")
#pileup=read.table(gzfile("C:/Users/Administrator/Desktop/HG00553.large_insert.Pileup.gz"),header=F,comment.char="",quote="\"")
pileup=read.table(gzfile("C:/Users/PC/Desktop/contamin_out.Pileup.gz"),header=F,comment.char="",quote="\"")

overall.likelihood(pileup[1:1000,],2,2,0.1)

ptm <- proc.time()
alpha=runif(1,0,0.5)
res <- optim( alpha, overall.likelihood, data=pileup, ploidy1=2, ploidy2=2,method="Brent",lower=0,upper=0.5)
pt2=proc.time() - ptm



ploidy.estimate=function()
{
  
}
