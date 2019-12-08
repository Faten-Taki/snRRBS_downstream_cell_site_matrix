
#This script uses the cell_site_matrix.csv outputted from the processing steps of single nuclear RRBS

#the input file merged all matrix files for all 2 experiment Runs


#merge all files from 2nd experiment



#Open R

R

#set working directory

setwd("/path/to/files/on/local/machine/cell_site_matrix/")



#Files have each row as one CpG site and each column as 1 cell



#define merging function

multmerge = function(mypath){
  filenames=list.files(path=mypath, full.names=TRUE)
  datalist = lapply(filenames, function(x){read.csv(file=x,header=T)})
  Reduce(function(x,y) {merge(x,y,all = TRUE)}, datalist)
}


#apply function and merge all matrices in the path provided


full_data = multmerge("/path/to/files/on/local/machine/cell_site_matrix/")


write.table(full_data, file="/path/to/files/on/local/machine/cell_site_matrix/merged_matrix_2_batches.txt", sep="\t", row.names=FALSE, quote=FALSE)


#To exit R type:
quit(save="no")



