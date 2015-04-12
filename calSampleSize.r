f=function(x){factorial(x[4])/(factorial(x[1])*factorial(x[2])*factorial(x[3])*factorial(x[4]-x[1]-x[2]-x[3]))*0.1^x[1]*0.1^x[2]*0.1^x[3]*(1-0.3)^(x[4]-x[1]-x[2]-x[3])}
integ=function(x)
{
  p=0
  total=x
  for (i in 1:(total-2))
  {
    for( j in 1:((total-1)-i))
    {
      for( k in 1:(total-i-j))
      {
        #mat=rbind(mat,c(i,j,k))
        p=p+f(c(i,j,k,total))
      }
    }
  }
return(p)
}
adaptIntegrate(f,lowerLimit = c(1,1,1,20),upperLimit = c(2,2,2,20))