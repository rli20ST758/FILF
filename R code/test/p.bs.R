#############################
# Written by: Stephanie Chen (stchen3@ncsu.edu)
# Purpose: Calculate p-value for the observed Tn from the null approximation
# Updated: Aug 4, 2018

p.bs<-function(stat,bs.stat){
  p<-mean(stat<=bs.stat)
  list(p=p,mean=mean(bs.stat),var=var(bs.stat))
}

