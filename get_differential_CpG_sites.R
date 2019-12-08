
#The below commands are in R and Unix (interchangeably)

#this document will enable the identification of differentially methylated CpG sites between 2 groups called: Pos and Neg

#The input file is the file we got from the document (merge_cell_site_matrices.R) in this repository (snRRBS_downstream_cell_site_matrix). This file has all matrices for 2 experimental runs merged into 1 table as text file. The file has all sites and cells (including those that are outliers). There were 288 Neg cells and 288 Pos cells. Total number of columns with the sites column is 577.




#Open R

R

#set working directory 

setwd("/path/to/files/on/local/machine/cell_site_matrix/")



#read text file in R

data = read.table(file = "merged_matrix_2_batches.txt", header = TRUE)

#transform values to numerical type

data[,2:577] = apply(data[,2:577], 2, function(x) as.numeric(as.character(x)));

#Remove outliers. Outliers were cells that had less than 50,000 CpG(s) detected, =< 98% bisulfite conversion, and had reads more than the 99th percentile of all cells. This information is in a separate excel document that include the summary information per cell. So, we can track which cell is an outlier in excel, keep track of it, and use that number to subset the dataframe (data) to remove outliers.


data_no_outlier = data.frame(cbind(data[,1:18], data[,20:23], data[,25:114], data[,116:154], data[,156:165], data[,167:172], data[,174:191], data[,193:296], data[,298:313], data[,315:320], data[,322:323], data[,325:326], data[,328:337], data[,387:400], data[,402:447], data[,449:480], data[,483:488], data[,490], data[,493], data[,495:496], data[,498], data[,500:501], data[,505:532], data[,534:535], data[,537:539], data[,541:546], data[,548:567], data[,569:577]))

#organize the columns by cell type, such that all Neg cells are first, followed by all Pos cells. 

data_no_outlier_organized = data.frame(cbind(data_no_outlier[,1:95], data_no_outlier[,187:325], data_no_outlier[,96:186], data_no_outlier[,326:498]))


#get the site names and download for your record to your folder (on local machine or server, whereever you are working). This will be used later.

sites = data$site


write.table(sites, file="/path/to/files/on/local/machine/cell_site_matrix/sites_2_batches.txt", sep="\t", row.names=FALSE, quote=FALSE)



~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~


##in Unix terminal

#organize the "sites_2_batches.txt" file to make it compatible for intersectBed down the line

awk 'BEGIN{OFS="\t"} {print $1, $2, $2}' sites_2_batches.txt > sites_2_batches_cleaned.txt
 
 
 
~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~



##back in R


sites_2_batches_cleaned = read.table("sites_2_batches_cleaned.txt", header = FALSE)   



#prepare data_2 file for intersectBed

data_2 = data.frame(cbind(sites_2_batches_cleaned,data_no_outlier_organized))


#label the columns

names(data_2) =c("chr","start","end","site","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos")


#get the sum of "1"s (methylation) per cell

data_2$Neg_meth = rowSums(data_2[,5:237], na.rm = TRUE) #NA values are ignored

data_2$Pos_meth = rowSums(data_2[,238:501], na.rm = TRUE) #NA values are ignored




#get the length of "detected" sites per cell

data_2$Neg_length = apply(data_2[,5:237],1,function(x) length(which(!is.na(x))))

data_2$Pos_length = apply(data_2[,238:501],1,function(x) length(which(!is.na(x))))

total_length = apply(data_2[,5:501],1,function(x) length(which(!is.na(x)))) #average number of cells within which a CpG is detected = 28.92865

mean(total_length)



#get the zeros of the detected sites per cell

data_2$Neg_zero = data_2$Neg_length - data_2$Neg_meth

data_2$Pos_zero = data_2$Pos_length - data_2$Pos_meth



#save what you have so far in case something crashes (optional)

write.table(data_2, file = "2_batches_before_barnard.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')



#let's subset to make the data smaller before we get differential sites

df = data.frame(cbind(data_2$Neg_meth,data_2$Pos_meth,data_2$Neg_zero,data_2$Pos_zero))

#attach the sites with the subsetted dataframe

df = cbind(data_2$site,df)

#save for your records (optional)

write.table(df, file = "df_for_barnard.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')





#I use Barnard test to get the differentially methylated sites b/t Neg and Pos cells

#load the library Barnard

library(Barnard)

#read the table that will be used for this test

df = read.table("df_for_barnard.txt")


#run Barnard and save the p.values into an object. This will take a long time (took 72 hours)

Barnard.p.val <- sapply(seq_len(nrow(df[,2:5])), function(x) {
  b = barnard.test(df[x,2],df[x,3],df[x,4],df[x,5])$p.value[2]})
  

#add the computed p.values to the dataframe
  
df$Barnard.p.val <- Barnard.p.val


#name the columns

names(df) = c("sites_2","Neg_1", "Pos_1", "Neg_0","Pos_0","p.val")


#get FDR (optional)

df$Barnard.fdr = p.adjust(df$p.val, method = "fdr", n = length(df$p.val))

#save the results

write.table(df, file = "df_after_barnard.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')

