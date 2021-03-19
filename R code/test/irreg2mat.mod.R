#############################
# Written by: Stephanie Chen (stchen3@ncsu.edu)
# Purpose: Transform functional data in "long" format to "wide" grid, with NA for missing obs

irreg2mat.mod<-function(data,times){
  grid<-matrix(NA,nrow=length(unique(data$.id)),ncol=length(times))
  inds<-cbind(data$.id,match(data$.index,times))
  grid[inds]<-data$.value
  return(grid)
}
