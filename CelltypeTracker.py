# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 12:30:59 2015

@author: PC
"""
"""
os.chdir('C:/Users/PC/Dropbox/workingspace/CellTypeTracker')
"""
import numpy as np
import math
from scipy.optimize import brent,minimize_scalar


allele_A_freq=0.25
allele_C_freq=0.25
allele_G_freq=0.25
allele_T_freq=0.25
table={'A':0,'C':1,'G':2,'T':3}
def genotype_prior(genotype):
#  numA=af['A']
#  numC=af['C']
#  numG=af['G']
#  numT=af['T']

  return (math.factorial(sum(genotype))/(math.factorial(genotype[0])*math.factorial(genotype[1])*math.factorial(genotype[2])*math.factorial(genotype[3])))*(allele_A_freq**genotype[0])*(allele_C_freq**genotype[1])*(allele_G_freq**genotype[2])*(allele_T_freq**genotype[3])

  
def conditional_prob(base,error_status,af):
    total=sum(af)+0.0
#k=sum(c(1,1,1,1)[c(num["A"]>0,num["C"]>0,num["G"]>0,num["T"]>0)])# how many kinds of alleles
    if error_status==1:
        return (1-af[base]/total)*1./3
    else:# no error event  
        return af[base]/total

import itertools
def generate_genotype_permutation(ploidy):
    perm={}
    for a in list(itertools.product(*[['A','C','G','T']]*ploidy)):
        perm[''.join(sorted(a))]=1
    return perm
    
geno_list_table={}
geno_list_allele_freq={}
geno_prior_table={}
for i in range(1,6):
    geno_list_table[i]=generate_genotype_permutation(i)
    for j in geno_list_table[i].keys():
        geno_list_allele_freq[j]=[j.count('A'),j.count('C'),j.count('G'),j.count('T')]
        geno_prior_table[j]=genotype_prior(geno_list_allele_freq[j])

def single_observation_likelihood(base,quality,alpha,af1,af2):


  #print(quality)
  quality=ord(quality)
  #print(quality)
  #if error happened
  likelihood=10**(0-quality/10)*((1-alpha)*conditional_prob(table[base],1,af1)+alpha*conditional_prob(table[base],1,af2))
    
  #if error not happened
  likelihood=likelihood+(1-10**(0-quality/10))*((1-alpha)*conditional_prob(table[base],0,af1)+alpha*conditional_prob(table[base],0,af2))
  #print(likelihood)
  return math.log(likelihood)

def site_likelihood(seq,qual,ploidy1,ploidy2,alpha):
  genotype_list1=generate_genotype_permutation(ploidy1)
  genotype_list2=generate_genotype_permutation(ploidy2)
  likelihood=0;
  for p in range(len(genotype_list1)):
    geno1=genotype_list1.keys()[p]
    af1=[geno1.count('A'),geno1.count('C'),geno1.count('G'),geno1.count('T')]
    tmplike2=0;
    for q in range(len(genotype_list2)):
      geno2=genotype_list2.keys()[q]
      af2=[geno2.count('A'),geno2.count('C'),geno2.count('G'),geno2.count('T')]
      loglik=0;
      for i in range(len(seq)):#assume seq is of format 'AAAAAAAAA'
        if(seq[i]=='N'):
            continue
        loglik=loglik+single_observation_likelihood(seq[i],qual[i],alpha,af1,af2)
      #print(loglik)
      tmplike2=tmplike2+math.exp(loglik)*genotype_prior(af2)
    
    #print(tmplike2)
    likelihood=likelihood+tmplike2*genotype_prior(af1)
  
  #print(likelihood)
  return math.log(likelihood)
  
  
def site_likelihood_fast(seq,qual,ploidy1,ploidy2,alpha):
  genotype_list1=geno_list_table[ploidy1]
  genotype_list2=geno_list_table[ploidy2]
  likelihood=0;
  for p in range(len(genotype_list1)):
    geno1=genotype_list1.keys()[p]
    af1=geno_list_allele_freq[geno1]
    tmplike2=0;
    for q in range(len(genotype_list2)):
      geno2=genotype_list2.keys()[q]
      af2=geno_list_allele_freq[geno2]
      loglik=0;
      for i in range(len(seq)):#assume seq is of format 'AAAAAAAAA'
        if(seq[i]=='N'):
            continue
        loglik=loglik+single_observation_likelihood(seq[i],qual[i],alpha,af1,af2)
      #print(loglik)
      tmplike2=tmplike2+math.exp(loglik)*geno_prior_table[geno2]
    
    #print(tmplike2)
    likelihood=likelihood+tmplike2*geno_prior_table[geno1]
  
  #print(likelihood)
  return math.log(likelihood)
  
def site_likelihood_wrapper(iterator):
  if(iterator[2]<6 and iterator[3]<6):
      return (site_likelihood_fast(*iterator))
  else:
      return (site_likelihood(*iterator))
  

import time
def overall_likelihood(alpha,data,ploidy1,ploidy2,pool):
  if(alpha<0 or alpha>0.5):
      return sys.maxint
  start = time.time()
  loglikelihood=0;

#  for i in range(0,num,4):
#    seq=[data[i][4].upper(),
#         data[i+1][4].upper(),
#         data[i+2][4].upper(),
#         data[i+3][4].upper()
#         ]
#    qual=[data[i][5],
#          data[i+1][5],
#          data[i+2][5],
#          data[i+3][5]
#          ]#print(seq)
#    #loglikelihood=loglikelihood+site_likelihood(seq,qual,ploidy1,ploidy2,alpha)
#    results=pool.map(site_likelihood_wrapper,zip(seq,qual,itertools.repeat(ploidy1),itertools.repeat(ploidy2),itertools.repeat(alpha)))
#    loglikelihood=loglikelihood+sum(results)
  
  seq=map(lambda meth:meth[4].upper(),data)
  qual=map(lambda meth:meth[5],data)
  results=pool.map(site_likelihood_wrapper,zip(seq,qual,itertools.repeat(ploidy1),itertools.repeat(ploidy2),itertools.repeat(alpha)))
  loglikelihood=sum(results)  

  print "time 1 eplsed:",time.time() - start,"s\n","alpha:",alpha,"\tlik:",-loglikelihood,"\n"
  return -loglikelihood
  
  
  
def overall_likelihood2(data,ploidy1,ploidy2,alpha):
  start = time.time()
  loglikelihood=0;
  num=len(data)

  for i in range(0,num):
    seq=data[i][4].upper()
    qual=data[i][5]
    loglikelihood=loglikelihood+site_likelihood(seq,qual,ploidy1,ploidy2,alpha)
  print "time 2 eplsed:",time.time() - start,"s"
  return -loglikelihood
  
if __name__ == '__main__':  
    from multiprocessing import Pool
    import sys
    
#    filename="C:/Users/Administrator/Desktop/HG00553.large_insert.Pileup.gz"
    filename=sys.argv[1]
    data=np.genfromtxt(filename,comments='~',dtype={'names':('chr','pos','ID','depth','seq','qual','strand','coordinate'),'formats':('i1','i32','S1','i1','S128','S128','S128','S128')})
    p=Pool(7)
#    plotdata=[]
#    alpha= range(0,100,1)
#    for a in alpha:
#        plotdata.append(overall_likelihood(data,2,2,a/100.,p))
    a=0.3
    #res = minimize_scalar(overall_likelihood,a,args=(data,2,2,p), method='brent',bounds=(0,0.5),bracket=(0,0.5))
    res = brent(overall_likelihood,args=(data,2,2,p))
    p.close()
    p.join()
    print res
   # import matplotlib.pyplot as plt
    
#    plt.plot(alpha,plotdata)
#    plt.savefig('C:/Users/Administrator/Documents/Python Scripts/foo.pdf')
#    plt.close()

    #print overall_likelihood2(data[1:1000],2,2,0.1)
    