install.packages("EnvStats")
install.packages("matrixStats")
#install.packages("factoextra")
library(EnvStats)
library(matrixStats)
library(dplyr)


setwd("/Users/manasarapuru/Desktop/")

#Data is being loaded in 
data = read.csv("corrected_exprs_data.csv", header = T, sep = ",")
#Skipping header
rownames(data) = data[,1]
data = data[,-1]

#Function to identify if sum is greater than log 2 15 or not
greater_than <- function(x,s1) { #function to count specific values
  x1 <- x[s1]
  x1 <- as.numeric(x1)
  required = sum(x1 > log2(15))
  return(required)
}

#4.1: filter for rows whose 20% of expression data is above log(2)15
#Totals up expression level for each gene across samples
count = apply(data, 1, greater_than, s1 = names(data)) # apply function to all rows
data1 = cbind(data, count)  # column bind the value for each row
data1 <- data1[data1$count > (20/100)*ncol(data),]  # get rows that meet the first filter
data2  = subset(data1, select = -count)   #delete the count column
#data2 contains genes that pass the first filter

#4.2 
my_varTest <- function(x,s1) { #function to count specific values
  x1 <- x[s1]
  x1 <- as.numeric(x1)
  required = varTest(x1, alternative = "two.sided", conf.level = 0.95, sigma.squared = 1)
  return(required$p.value)
}

#gets a list of pvalues using results from second filter
p_values = apply(data2, 1, my_varTest, s1 = names(data2))   #get a list of p-values.
data3 = cbind(data2, p_values)   # bind all p_values
#Filters for rows that have pvalues greater than 0.01
data4 = data3[data3$p_values < 0.01, ]
#Gets rid of pvalue column
data5  = subset(data4, select = -p_values)

#data 5 is data that passed through first 2 filters
write.csv(data5,file='noise_filter_passed_two_filters.csv')


#4.3 coefficient of variation
data6 = data.frame(data5)
#calculates row-wise Standard Deviation
data6$Std_Deviation = rowSds(as.matrix(data6))  
#calculates row-wise Means
data6$Mean = rowMeans(as.matrix(data6[,-ncol(data6)]))
# calculates Coefficient of Variance for each row
data6$coe_of_variace = data6$Std_Deviation / data6$Mean
#select genes with Coef. of Variace > 0.186
data6 = data6[data6$coe_of_variace > 0.186, ]
#Removes the additional Std. Deviation, mean and coefficient of variance rows
data7 = subset(data6, select = c(-Std_Deviation, -Mean, -coe_of_variace))

#data7 is the final data frame to use after Filter1, Filter2, Filter3.
write.csv(data7,file='noise_filter_passed_all_filters.csv')

#5.1: Hierarchy cluster of samples - Use hierarchy clustering to cluster the samples
transf <- t(data7)

#Used standard distance function to generate matrix and passed it into hclus function
clusters <- hclust(dist(transf))
#plot(clusters)


#5.2: Cut into 2 groups 
cluster_num <- cutree(clusters, 2)
table(cluster_num)
#turn into data frame, add the cluster labels to transformed table
cluster_cut_table <-as.data.frame(cluster_num)
cluster_labelled = cbind(transf, cluster_cut_table)

#5.3 create heat map of data
heat_map_table <- data7

#The samples in the heat map are being named as the geo_accession name found in the metadata file
names(heat_map_table) <- sub("_.*", "", names(heat_map_table))

#Retransform data before heatmap
heat_map_table<-t(heat_map_table)

#Load the metadata file, from 
meta = read.csv("proj_metadata.csv", header = T, sep = ",")

#get the geo accession number and subtype into another datafra
keeps <- c("geo_accession","SixSubtypesClassification")
meta_name_type = meta[keeps]

#Map meta data to actual table with this function
maptodata<-function(annote_file){
  colorsVector = ifelse(annote_file["SixSubtypesClassification"]=="C3", 
                        "red", ifelse(annote_file["SixSubtypesClassification"]!="C3", 
                        "blue","green"))
  return(colorsVector)
}
labels = maptodata(meta_name_type)
heatmap(heat_map_table, RowSideColors=labels)
rcluster_labelled = data.frame(t(cluster_labelled))

my_tTest = function(x,s1,s2) { #function get p value and t statistic from the t.test() function
  x1 = x[s1]
  x2 = x[s2]
  x1 = as.numeric(x1)
  x2 = as.numeric(x2)
  required = t.test(x1, x2, alternative="two.sided")
  p_value = required$p.value
  t_stats = required$statistic
  combine = c(p_value, t_stats)
}

#a data set for each cluster is being subsetted here
a = rcluster_labelled[,(rcluster_labelled[nrow(rcluster_labelled), ]) == 1]
b = rcluster_labelled[,(rcluster_labelled[nrow(rcluster_labelled), ]) == 2]
data8 = rcluster_labelled[!(row.names(rcluster_labelled) %in% "cluster_num"),]


results = apply(data8, 1 ,my_tTest, s1= names(a), s2 = names(b))
data9 = cbind(data8, data.frame(t(results)))
names(data9)[names(data9) == "V1"] <- "p_value"
names(data9)[names(data9) == "t"] <- "t_statistics"

#filter for p adjusted value greater than 0.05
data9$p_adjusted = p.adjust(data9$p_value)
data10 = data.frame(data9[data9$p_adjusted < 0.05,])

#Data 10 is data 6 after t_test and after filtering with p adjusted value 
data11 = subset(data10, select =c(-p_value, -t_statistics, -p_adjusted))
write.csv(data11,file='clustering_passed_all_filters.csv')
write.csv(data10[, 135:137], file='clustering_all_filters_statistics.csv' ,
          row.names = T)

#Data 5 res
results_data5 = apply(data5, 1 ,my_tTest, s1= names(a), s2 = names(b))
data12 = cbind(data5, data.frame(t(results_data5)))
names(data12)[names(data12) == "V1"] <- "p_value"
names(data12)[names(data12) == "t"] <- "t_statistics"

data12$p_adjusted = p.adjust(data12$p_value)
data13 = data.frame(data12[data12$p_adjusted < 0.05,])

#Data 14 is data second filter data analysis after t_test and after filtering with p adjusted value 
data14 = subset(data13, select =c(-p_value, -t_statistics, -p_adjusted))
write.csv(data14,file='clustering_passed_two_filters.csv')
write.csv(data13[, 135:137], file='clustering_passed_two_filters_statistics.csv' ,
          row.names = T)