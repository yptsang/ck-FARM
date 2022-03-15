# ck-FARM
An R Package to Discover Big Data Association for Business Intelligence  
  
The proposed algorithm, namely ck-FARM, is deployed in the R programming environment with applying five existing   libraries, including Matrix, arules, ca, ggplots, and factoextra. Thus, the statistical significance between process parameters can be examined, while fuzzy association rules (FARs) can be automatically built. Based on the fuzzy relationship of process parameters, corresponding business management/re-engineering strategies can be established to strengthen business intelligence capabilities.
  
Example of using the package:  
library(Matrix)  
library(arules)  
library(ca)  
library(ggplot2)  
library(factoextra)  

data  <- read.table(file="filePath", header=TRUE)  

n <- dim(data)[1]  
p <- dim(data)[2]  
max_clu <- 4  
lambda <- 0.7  
alpha <- numeric(p)+n/50  

corrana_data(data)  
mem_fun_result(data)  
result <- main_func(data, alpha, lambda)  
