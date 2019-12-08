
#after you plot and visualize the cells accroding to their methylation status, you may decide to cluster those cells

#this is downstream of the Rtsne done in "visualize_using_tSNE.R"


##the below script is from the link before

#https://jmonlong.github.io/Hippocamplus/2018/02/13/tsne-and-clustering/

#load libraries

library(ggplot2)
library(dplyr)
library(magrittr)
library(ggrepel)
library(Rtsne)


#use files outputted from the "get_differential_CpG_sites.R"

#the below is the same as the "visualize_using_tSNE.R"

#read the file

sites = read.table("sites_detected_in_3_cells_per_group_2_batches.txt")

#function to replace missing values with the mean of the column (for each cell)

NA2mean <- function(x) replace(x, is.na(x), mean(x, na.rm = TRUE))

#apply the function

sites_NA2mean = replace(sites[,5:501], TRUE, lapply(sites[,5:501], NA2mean))

#attach the sites to the dataset

sites_NA2mean = data.frame(cbind(sites[,4],sites_NA2mean))

#get unique values

sites_NA2mean_unique = unique(sites_NA2mean)

#transform dataset into matrix

sites_m = as.matrix(sites_NA2mean[,2:498])



#transpose matrix

z = t(DE_sites_m)

#normalize input

y = normalize_input(z)


#get t-SNE

Rtsne = Rtsne(y, check_duplicates = FALSE)


#label cell groups in different colors

cells =c("black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red")


#label cells

cell_name = c("N1","N2","N3","N4","N5","N6","N7","N8","N9","N10","N11","N12","N13","N14","N15","N16","N17","N18","N19","N20","N21","N22","N23","N24","N25","N26","N27","N28","N29","N30","N31","N32","N33","N34","N35","N36","N37","N38","N39","N40","N41","N42","N43","N44","N45","N46","N47","N48","N49","N50","N51","N52","N53","N54","N55","N56","N57","N58","N59","N60","N61","N62","N63","N64","N65","N66","N67","N68","N69","N70","N71","N72","N73","N74","N75","N76","N77","N78","N79","N80","N81","N82","N83","N84","N85","N86","N87","N88","N89","N90","N91","N92","N93","N94","N95","N96","N97","N98","N99","N100","N101","N102","N103","N104","N105","N106","N107","N108","N109","N110","N111","N112","N113","N114","N115","N116","N117","N118","N119","N120","N121","N122","N123","N124","N125","N126","N127","N128","N129","N130","N131","N132","N133","N134","N135","N136","N137","N138","N139","N140","N141","N142","N143","N144","N145","N146","N147","N148","N149","N150","N151","N152","N153","N154","N155","N156","N157","N158","N159","N160","N161","N162","N163","N164","N165","N166","N167","N168","N169","N170","N171","N172","N173","N174","N175","N176","N177","N178","N179","N180","N181","N182","N183","N184","N185","N186","N187","N188","N189","N190","N191","N192","N193","N194","N195","N196","N197","N198","N199","N200","N201","N202","N203","N204","N205","N206","N207","N208","N209","N210","N211","N212","N213","N214","N215","N216","N217","N218","N219","N220","N221","N222","N223","N224","N225","N226","N227","N228","N229","N230","N231","N232","N233","P1","P2","P3","P4","P5","P6","P7","P8","P9","P10","P11","P12","P13","P14","P15","P16","P17","P18","P19","P20","P21","P22","P23","P24","P25","P26","P27","P28","P29","P30","P31","P32","P33","P34","P35","P36","P37","P38","P39","P40","P41","P42","P43","P44","P45","P46","P47","P48","P49","P50","P51","P52","P53","P54","P55","P56","P57","P58","P59","P60","P61","P62","P63","P64","P65","P66","P67","P68","P69","P70","P71","P72","P73","P74","P75","P76","P77","P78","P79","P80","P81","P82","P83","P84","P85","P86","P87","P88","P89","P90","P91","P92","P93","P94","P95","P96","P97","P98","P99","P100","P101","P102","P103","P104","P105","P106","P107","P108","P109","P110","P111","P112","P113","P114","P115","P116","P117","P118","P119","P120","P121","P122","P123","P124","P125","P126","P127","P128","P129","P130","P131","P132","P133","P134","P135","P136","P137","P138","P139","P140","P141","P142","P143","P144","P145","P146","P147","P148","P149","P150","P151","P152","P153","P154","P155","P156","P157","P158","P159","P160","P161","P162","P163","P164","P165","P166","P167","P168","P169","P170","P171","P172","P173","P174","P175","P176","P177","P178","P179","P180","P181","P182","P183","P184","P185","P186","P187","P188","P189","P190","P191","P192","P193","P194","P195","P196","P197","P198","P199","P200","P201","P202","P203","P204","P205","P206","P207","P208","P209","P210","P211","P212","P213","P214","P215","P216","P217","P218","P219","P220","P221","P222","P223","P224","P225","P226","P227","P228","P229","P230","P231","P232","P233","P234","P235","P236","P237","P238","P239","P240","P241","P242","P243","P244","P245","P246","P247","P248","P249","P250","P251","P252","P253","P254","P255","P256","P257","P258","P259","P260","P261","P262","P263","P264")


#plot t-SNE

plot(Rtsne$Y,col=cells, asp=1)






###cluster the cells

#get the hierarchical clustering using the t-SNE as input

hc.norm = hclust(dist(Rtsne$Y))

info.norm = tibble(truth = factor(cell_name))

info.norm %<>% mutate(tsne1 = Rtsne$Y[, 1], tsne2 = Rtsne$Y[, 2])

#here, based on the visualization from t-SNE, I decided to cluster my cells into 4 clusters

info.norm$hclust = factor(cutree(hc.norm, 4))

hc.norm.cent = info.norm %>% group_by(hclust) %>% select(tsne1, 
    tsne2) %>% summarize_all(mean)
	

#plot clustering

ggplot(info.norm, aes(x = tsne1, y = tsne2, colour = hclust)) + 
geom_point(alpha = 0.3) + theme_bw() + geom_label_repel(aes(label = hclust), 
data = hc.norm.cent) + guides(colour = FALSE) +  ggtitle("Linkage criterion: Complete")
	

#get the cells per cluster and save it locally
	
write.table(info.norm, "all_sites_cell_clusters.txt", quote = FALSE, sep = '\t')


  



