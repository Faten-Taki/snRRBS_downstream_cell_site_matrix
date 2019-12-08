
#This file is for plotting the CpGs in the 2 cell groups (Neg and Pos) using t-SNE

#the files outputted from get_differential_CpG_sites.R will be used




#in R 

R

#read the tables containing the results from Barnard test as well as the table containing all individual cells as columns

df = read.table("df_after_barnard.txt") 

data_2 = read.table("2_batches_before_barnard.txt")


#label the columns

names(data_2) =c("chr","start","end","site","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Neg","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Pos","Neg_meth", "Pos_meth", "Neg_length","Pos_length","Neg_zero","Pos_zero")


 #label the columns in the table containing results after Barnard test
 
 names(df) = c("sites_2","Neg_1", "Pos_1", "Neg_0","Pos_0","p.val","FDR")


#concatenate df to data_2 with all the data points

data_3 = data.frame(cbind(data_2,df))

#get the number of cells in which each CpG is detected. (meaning per row)

data_3$sum.cells.per.CpG <- rowSums(data_3[, c(509, 510, 511, 512)])

#get the average number of cells in which a CpG is detected

avg_cells_per_CpG = mean(data_3$sum.cells.per.CpG) #optional (for your record)

#get the number of Neg cells in which each CpG is detected

data_3$sum.neg.cells.per.CpG = rowSums(data_3[, c(509, 511)])

#get the average number of Neg cells in which a CpG is detected #optional (for your record)

avg_neg.cells_per_CpG = mean(data_3$sum.neg.cells.per.CpG) 

#get the number of Pos cells in which a CpG is detected

data_3$sum.pos.cells.per.CpG = rowSums(data_3[, c(510, 512)])

#get the average number of Pos cells in which a CpG is detected

avg_pos.cells_per_CpG = mean(data_3$sum.pos.cells.per.CpG) #optional (for your record)





#subset data_3 with sites detected in at least 3 Pos cells and 3 Neg cells


data_3_ordered_subset_A = subset(data_3, Neg_length >= 3)

#save (optional) for your record

write.table(data_3_ordered_subset_A, file = "sites_detected_in_3_cells_per_Neg_2_batches.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')


data_3_ordered_subset_B = subset(data_3_ordered_subset_A, Pos_length >= 3) #XX sites detected in at least 3 Neg cells and 3 Pos cells


#optional get FDR using only the filtered sites

data_3_ordered_subset_B$Barnard.fdr = p.adjust(data_3_ordered_subset_B$p.val, method = "fdr", n = length(data_3_ordered_subset_B$p.val))

#optional count how many sites are differential

count_p <- length(which(data_3_ordered_subset_B$p.val < 0.05)) #X sites with p<0.05 detected in >=3 Neg and >=3 Pos cells


#save this table in a file locally

write.table(data_3_ordered_subset_B, file = "sites_detected_in_3_cells_per_group_2_batches.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')
 






#prepare to overlap with coordinates of interest (it is always good to do that, to prepare it for any intersectBed down the line)

#order data based on chr name and then site position

data_3_subset_B_chr_order = data_3_ordered_subset_B[order(data_3_ordered_subset_B$chr, -data_3_ordered_subset_B$start),]

#subset the data to get the sites with p.val < 0.05

data_3_subset_B_chr_order_p0.05 = subset(data_3_subset_B_chr_order, p.val < 0.05)


#take a look at the column names for your record

names(data_3_subset_B_chr_order_p0.05)
# [1] "chr"                   "start"                 "end"
# [4] "site"                  "Neg"                   "Neg.1"
# [7] "Neg.2"                 "Neg.3"                 "Neg.4"
# [10] "Neg.5"                 "Neg.6"                 "Neg.7"
# [13] "Neg.8"                 "Neg.9"                 "Neg.10"
# [16] "Neg.11"                "Neg.12"                "Neg.13"
# [19] "Neg.14"                "Neg.15"                "Neg.16"
# [22] "Neg.17"                "Neg.18"                "Neg.19"
# [25] "Neg.20"                "Neg.21"                "Neg.22"
# [28] "Neg.23"                "Neg.24"                "Neg.25"
# [31] "Neg.26"                "Neg.27"                "Neg.28"
# [34] "Neg.29"                "Neg.30"                "Neg.31"
# [37] "Neg.32"                "Neg.33"                "Neg.34"
# [40] "Neg.35"                "Neg.36"                "Neg.37"
# [43] "Neg.38"                "Neg.39"                "Neg.40"
# [46] "Neg.41"                "Neg.42"                "Neg.43"
# [49] "Neg.44"                "Neg.45"                "Neg.46"
# [52] "Neg.47"                "Neg.48"                "Neg.49"
# [55] "Neg.50"                "Neg.51"                "Neg.52"
# [58] "Neg.53"                "Neg.54"                "Neg.55"
# [61] "Neg.56"                "Neg.57"                "Neg.58"
# [64] "Neg.59"                "Neg.60"                "Neg.61"
# [67] "Neg.62"                "Neg.63"                "Neg.64"
# [70] "Neg.65"                "Neg.66"                "Neg.67"
# [73] "Neg.68"                "Neg.69"                "Neg.70"
# [76] "Neg.71"                "Neg.72"                "Neg.73"
# [79] "Neg.74"                "Neg.75"                "Neg.76"
# [82] "Neg.77"                "Neg.78"                "Neg.79"
# [85] "Neg.80"                "Neg.81"                "Neg.82"
# [88] "Neg.83"                "Neg.84"                "Neg.85"
# [91] "Neg.86"                "Neg.87"                "Neg.88"
# [94] "Neg.89"                "Neg.90"                "Neg.91"
# [97] "Neg.92"                "Neg.93"                "Neg.94"
#[100] "Neg.95"                "Neg.96"                "Neg.97"
#[103] "Neg.98"                "Neg.99"                "Neg.100"
#[106] "Neg.101"               "Neg.102"               "Neg.103"
#[109] "Neg.104"               "Neg.105"               "Neg.106"
#[112] "Neg.107"               "Neg.108"               "Neg.109"
#[115] "Neg.110"               "Neg.111"               "Neg.112"
#[118] "Neg.113"               "Neg.114"               "Neg.115"
#[121] "Neg.116"               "Neg.117"               "Neg.118"
#[124] "Neg.119"               "Neg.120"               "Neg.121"
#[127] "Neg.122"               "Neg.123"               "Neg.124"
#[130] "Neg.125"               "Neg.126"               "Neg.127"
#[133] "Neg.128"               "Neg.129"               "Neg.130"
#[136] "Neg.131"               "Neg.132"               "Neg.133"
#[139] "Neg.134"               "Neg.135"               "Neg.136"
#[142] "Neg.137"               "Neg.138"               "Neg.139"
#[145] "Neg.140"               "Neg.141"               "Neg.142"
#[148] "Neg.143"               "Neg.144"               "Neg.145"
#[151] "Neg.146"               "Neg.147"               "Neg.148"
#[154] "Neg.149"               "Neg.150"               "Neg.151"
#[157] "Neg.152"               "Neg.153"               "Neg.154"
#[160] "Neg.155"               "Neg.156"               "Neg.157"
#[163] "Neg.158"               "Neg.159"               "Neg.160"
#[166] "Neg.161"               "Neg.162"               "Neg.163"
#[169] "Neg.164"               "Neg.165"               "Neg.166"
#[172] "Neg.167"               "Neg.168"               "Neg.169"
#[175] "Neg.170"               "Neg.171"               "Neg.172"
#[178] "Neg.173"               "Neg.174"               "Neg.175"
#[181] "Neg.176"               "Neg.177"               "Neg.178"
#[184] "Neg.179"               "Neg.180"               "Neg.181"
#[187] "Neg.182"               "Neg.183"               "Neg.184"
#[190] "Neg.185"               "Neg.186"               "Neg.187"
#[193] "Neg.188"               "Neg.189"               "Neg.190"
#[196] "Neg.191"               "Neg.192"               "Neg.193"
#[199] "Neg.194"               "Neg.195"               "Neg.196"
#[202] "Neg.197"               "Neg.198"               "Neg.199"
#[205] "Neg.200"               "Neg.201"               "Neg.202"
#[208] "Neg.203"               "Neg.204"               "Neg.205"
#[211] "Neg.206"               "Neg.207"               "Neg.208"
#[214] "Neg.209"               "Neg.210"               "Neg.211"
#[217] "Neg.212"               "Neg.213"               "Neg.214"
#[220] "Neg.215"               "Neg.216"               "Neg.217"
#[223] "Neg.218"               "Neg.219"               "Neg.220"
#[226] "Neg.221"               "Neg.222"               "Neg.223"
#[229] "Neg.224"               "Neg.225"               "Neg.226"
#[232] "Neg.227"               "Neg.228"               "Neg.229"
#[235] "Neg.230"               "Neg.231"               "Neg.232"
#[238] "Pos"                   "Pos.1"                 "Pos.2"
#[241] "Pos.3"                 "Pos.4"                 "Pos.5"
#[244] "Pos.6"                 "Pos.7"                 "Pos.8"
#[247] "Pos.9"                 "Pos.10"                "Pos.11"
#[250] "Pos.12"                "Pos.13"                "Pos.14"
#[253] "Pos.15"                "Pos.16"                "Pos.17"
#[256] "Pos.18"                "Pos.19"                "Pos.20"
#[259] "Pos.21"                "Pos.22"                "Pos.23"
#[262] "Pos.24"                "Pos.25"                "Pos.26"
#[265] "Pos.27"                "Pos.28"                "Pos.29"
#[268] "Pos.30"                "Pos.31"                "Pos.32"
#[271] "Pos.33"                "Pos.34"                "Pos.35"
#[274] "Pos.36"                "Pos.37"                "Pos.38"
#[277] "Pos.39"                "Pos.40"                "Pos.41"
#[280] "Pos.42"                "Pos.43"                "Pos.44"
#[283] "Pos.45"                "Pos.46"                "Pos.47"
#[286] "Pos.48"                "Pos.49"                "Pos.50"
#[289] "Pos.51"                "Pos.52"                "Pos.53"
#[292] "Pos.54"                "Pos.55"                "Pos.56"
#[295] "Pos.57"                "Pos.58"                "Pos.59"
#[298] "Pos.60"                "Pos.61"                "Pos.62"
#[301] "Pos.63"                "Pos.64"                "Pos.65"
#[304] "Pos.66"                "Pos.67"                "Pos.68"
#[307] "Pos.69"                "Pos.70"                "Pos.71"
#[310] "Pos.72"                "Pos.73"                "Pos.74"
#[313] "Pos.75"                "Pos.76"                "Pos.77"
#[316] "Pos.78"                "Pos.79"                "Pos.80"
#[319] "Pos.81"                "Pos.82"                "Pos.83"
#[322] "Pos.84"                "Pos.85"                "Pos.86"
#[325] "Pos.87"                "Pos.88"                "Pos.89"
#[328] "Pos.90"                "Pos.91"                "Pos.92"
#[331] "Pos.93"                "Pos.94"                "Pos.95"
#[334] "Pos.96"                "Pos.97"                "Pos.98"
#[337] "Pos.99"                "Pos.100"               "Pos.101"
#[340] "Pos.102"               "Pos.103"               "Pos.104"
#[343] "Pos.105"               "Pos.106"               "Pos.107"
#[346] "Pos.108"               "Pos.109"               "Pos.110"
#[349] "Pos.111"               "Pos.112"               "Pos.113"
#[352] "Pos.114"               "Pos.115"               "Pos.116"
#[355] "Pos.117"               "Pos.118"               "Pos.119"
#[358] "Pos.120"               "Pos.121"               "Pos.122"
#[361] "Pos.123"               "Pos.124"               "Pos.125"
#[364] "Pos.126"               "Pos.127"               "Pos.128"
#[367] "Pos.129"               "Pos.130"               "Pos.131"
#[370] "Pos.132"               "Pos.133"               "Pos.134"
#[373] "Pos.135"               "Pos.136"               "Pos.137"
#[376] "Pos.138"               "Pos.139"               "Pos.140"
#[379] "Pos.141"               "Pos.142"               "Pos.143"
#[382] "Pos.144"               "Pos.145"               "Pos.146"
#[385] "Pos.147"               "Pos.148"               "Pos.149"
#[388] "Pos.150"               "Pos.151"               "Pos.152"
#[391] "Pos.153"               "Pos.154"               "Pos.155"
#[394] "Pos.156"               "Pos.157"               "Pos.158"
#[397] "Pos.159"               "Pos.160"               "Pos.161"
#[400] "Pos.162"               "Pos.163"               "Pos.164"
#[403] "Pos.165"               "Pos.166"               "Pos.167"
#[406] "Pos.168"               "Pos.169"               "Pos.170"
#[409] "Pos.171"               "Pos.172"               "Pos.173"
#[412] "Pos.174"               "Pos.175"               "Pos.176"
#[415] "Pos.177"               "Pos.178"               "Pos.179"
#[418] "Pos.180"               "Pos.181"               "Pos.182"
#[421] "Pos.183"               "Pos.184"               "Pos.185"
#[424] "Pos.186"               "Pos.187"               "Pos.188"
#[427] "Pos.189"               "Pos.190"               "Pos.191"
#[430] "Pos.192"               "Pos.193"               "Pos.194"
#[433] "Pos.195"               "Pos.196"               "Pos.197"
#[436] "Pos.198"               "Pos.199"               "Pos.200"
#[439] "Pos.201"               "Pos.202"               "Pos.203"
#[442] "Pos.204"               "Pos.205"               "Pos.206"
#[445] "Pos.207"               "Pos.208"               "Pos.209"
#[448] "Pos.210"               "Pos.211"               "Pos.212"
#[451] "Pos.213"               "Pos.214"               "Pos.215"
#[454] "Pos.216"               "Pos.217"               "Pos.218"
#[457] "Pos.219"               "Pos.220"               "Pos.221"
#[460] "Pos.222"               "Pos.223"               "Pos.224"
#[463] "Pos.225"               "Pos.226"               "Pos.227"
#[466] "Pos.228"               "Pos.229"               "Pos.230"
#[469] "Pos.231"               "Pos.232"               "Pos.233"
#[472] "Pos.234"               "Pos.235"               "Pos.236"
#[475] "Pos.237"               "Pos.238"               "Pos.239"
#[478] "Pos.240"               "Pos.241"               "Pos.242"
#[481] "Pos.243"               "Pos.244"               "Pos.245"
#[484] "Pos.246"               "Pos.247"               "Pos.248"
#[487] "Pos.249"               "Pos.250"               "Pos.251"
#[490] "Pos.252"               "Pos.253"               "Pos.254"
#[493] "Pos.255"               "Pos.256"               "Pos.257"
#[496] "Pos.258"               "Pos.259"               "Pos.260"
#[499] "Pos.261"               "Pos.262"               "Pos.263"
#[502] "Neg_meth"              "Pos_meth"              "Neg_length"
#[505] "Pos_length"            "Neg_zero"              "Pos_zero"
#[508] "sites_2"               "Neg_1"                 "Pos_1"
#[511] "Neg_0"                 "Pos_0"                 "p.val"
#[514] "FDR"                   "sum.cells.per.CpG"     "sum.neg.cells.per.CpG"
#[517] "sum.pos.cells.per.CpG" "Barnard.fdr"



#save the results locally

write.table(data_3_subset_B_chr_order_p0.05, file = "p0.05_sites_subset_in_3_cells_per_group_2_batchs.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')








#visualize the cells (as a function of their methylation status across the sites) using t-SNE


#tsNE for sites_detected_in_3_cells_per_group.txt (all sites detected in >= 3 cells per group (cFos(-) and cFos(+))


#load library

library(Rtsne)

#read the file that includes all sites detected in at least 3 cells per group

sites = read.table("sites_detected_in_3_cells_per_group_2_batches.txt")

#function to replace missing values with the column mean (for each cell)

NA2mean <- function(x) replace(x, is.na(x), mean(x, na.rm = TRUE))

#apply function on dataset

sites_NA2mean = replace(sites[,5:501], TRUE, lapply(sites[,5:501], NA2mean))

#attach the site label to the transformed dataset

sites_NA2mean = data.frame(cbind(sites[,4],sites_NA2mean))

#get the unique values

sites_NA2mean_unique = unique(sites_NA2mean)

#transform dataset into a matrix

DE_sites_m = as.matrix(DE_sites_NA2mean[,2:498])


#transpose the matrix

z = t(DE_sites_m)

#normalize the dataset (from the Rtsne library)

y = normalize_input(z)


#get the t-SNE 

Rtsne = Rtsne(y, check_duplicates = FALSE)


#label the cell groups with different colors

cells =c("black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red")



#plot the t-SNE in a color-coded way

plot(Rtsne$Y,col=cells, asp=1)












#visulize t-SNE for  differentially methylated sites_detected_in_3_cells_per_group.txt (DE detected in >= 3 cells per group (cFos(-) and cFos(+))


#read data file

DE_sites = read.table("p0.05_sites_subset_in_3_cells_per_group_2_batchs.txt")

#define function to replace missing values with the mean per column (for each cell)

NA2mean <- function(x) replace(x, is.na(x), mean(x, na.rm = TRUE))

#apply function on dataset

DE_sites_NA2mean = replace(DE_sites[,5:501], TRUE, lapply(DE_sites[,5:501], NA2mean))

#attach the sites label to the transformed dataset

DE_sites_NA2mean = data.frame(cbind(DE_sites[,4],DE_sites_NA2mean))

#get unique values

DE_sites_NA2mean_unique = unique(DE_sites_NA2mean)

#transform dataset into matrix

DE_sites_m = as.matrix(DE_sites_NA2mean[,2:498])


#transform matrix

z = t(DE_sites_m)

#normalize input

y = normalize_input(z)


#get t-SNE

Rtsne = Rtsne(y, check_duplicates = FALSE)


#label the cell groups with different colors

cells =c("black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red")


#plot t-SNE

plot(Rtsne$Y,col=cells, asp=1)



