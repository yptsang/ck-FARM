library(Matrix)
library(arules) 
library(ca) 
library(ggplot2)
library(factoextra) 


corr_ana <- function(p1, p2){
  ## To test the corresponding analysis between two vectors
  n <- length(p1)
  result1 <- cluster(p1, max_clu)
  k1 <- result1[1]; c1 <- result1[-1]
  result2 <- cluster(p2, max_clu)
  k2 <- result2[1]; c2 <- result2[-1]
  count <- matrix(0, k1, k2)
  rej <- 0
  for (i in 1:(k1)){
    for (j in 1:(k2)){
        c <- 1
        for (l in 1:n){
          if( (p1[l] >= c1[i])&(p1[l] <= c1[i+1])&(p2[l] >= c2[j])&(p2[l] <= c2[j+1]))
            {c <- c+1}
          
        }
        count[i,j] <- c
    }
  }
  
  res <- ca(count)
  eig <- get_eigenvalue(res)
  tr <- sum(eig$eig)     
  ch2 <- tr*sum(as.matrix(count))
  df <- (nrow(count)-1)*(ncol(count)-1)
  pval <- pchisq(ch2, df=df, lower.tail=FALSE)

  if(pval < 0.05  ) {rej <- 1}
  return(rej)
}

corrana_data <- function(data){
  n <- dim(data)[1]; m <- dim(data)[2]
  corr <- matrix(0, m, m)
 
  for (i in 1:(m-1)){
     p1 <- data[, i]
     for (j in (i+1):m){
        p2 <- data[, j]
        corr[i, j] <- corr_ana(p1, p2)
     }
  }
  return(corr)
}



cluster <- function(x, max_clu){
  rep_time <- max_clu-1
  ch <- numeric(rep_time)
  n <- length(x)
  centers <- matrix(0, nrow <- rep_time, ncol <- max_clu)
  for (i in 1:rep_time){
    k <- i + 1
    center_ini <- numeric(k)
    center_ini[1] <- min(x)
    center_ini[k] <- max(x)
    if (k > 2){
       for (j in 2:(k-1)){
          center_ini[j] <- quantile(x, j/k); u=0
          while ((center_ini[j] == center_ini[j-1]) &(u < (k-1))) {
             u <- u+1
             center_ini[j] <- (center_ini[j]+quantile(x, (j+u)/k)) / 2
          }
          if (center_ini[j] == center_ini[k]) center_ini[j] <- (center_ini[j]+center_ini[j -1]) / 2
       }
    }
    
    cc <- try(kmeans(x, center_ini, iter.max=10, nstart=1))
    if ('try-error' %in% class(cc)) next else{
      kmeans_cluster <- kmeans(x, center_ini, iter.max=10, nstart=1)
      centers[i, 1:k] <- kmeans_cluster$centers
      ch[i] <- kmeans_cluster$betweenss / (k - 1) 
      ch[i] <- ch[i]/( kmeans_cluster$tot.withinss/ (n - k))
    }
  }
  max_k <- which(ch == max(ch)) + 1
  center <- centers[max_k-1, 1:(max_k)]
  center <- sort(center)
  result <- c(max_k, min(x), center,  max(x))
  return(result)
}


mem_fun <- function(cluster_result){
  ## cluster_result is a vector
  k <- cluster_result[1]
  min_x <- cluster_result[2]
  max_x <- cluster_result[k + 3]
  centers <- cluster_result[3:(k+2)]   
  mem_zigzag <- matrix(0, k, 3)

  if (k ==2){
     mem_zigzag[1, ] <- c(min_x, centers)
     mem_zigzag[2, ] <- c(centers, max_x)
  }
  if (k > 2){
     centers_new <- c(min_x, centers, max_x)
     for (i in 1:k){
         mem_zigzag[i, ] <- centers_new[i:(i+2)] 
     }   
  }
  return(mem_zigzag)
}

memval1 <- function(x, end_point){
  ## x is a number; end_point are three endpoints in the lines
  ## The left  side of the membership functions
  y <- 0
{
  if ((x >= end_point[1]) &(x <= end_point[2]) )  {
    y <- 1
  }
  if ((x >= end_point[2]) &(x <= end_point[3]) )    {
    y <- (end_point[3] - x)/(end_point[3] - end_point[2]) 
  }
}
  return(y)
}

memval2 <- function(x, end_point){
 ## x is a one-dimensional observed data; end_point are three endpoints in the lines
 ## The middle part of the functions, the triangular shape
  y <- 0
  if ((x >= end_point[1]) &(x <= end_point[2]) )   {y <- (x - end_point[1])/(end_point[2] - end_point[1]) }
  if ((x >= end_point[2]) &(x <= end_point[3]) )   {y <- (end_point[3] - x)/(end_point[3] - end_point[2]) }
  return(y)
}

memval3 <- function(x, end_point){
 ## x is a one-dimensional observed data; end_point are three endpoints in the lines
 ## The right side of the membership functions
  y <- 0
  if ((x >= end_point[1]) &(x <= end_point[2]) )    {y <- (x - end_point[1])/(end_point[2] - end_point[1]) }
  if ((x >= end_point[2]) &(x <= end_point[3]) )    {y <- 1 }
  return(y)
}



mem_fun_val <- function(x){
  ## x here is a vector; e.g, x is the observed data in one process
  n <- length(x)
  cluster_result <- cluster(x, max_clu)
  k <- cluster_result[1]  
  mem_zigzag <- mem_fun(cluster_result)
  mem_value <- matrix(0, n, k)

  for (i in 1:n){
     mem_value[i, 1] <- memval1(x[i], mem_zigzag[1, ])
     mem_value[i, k] <- memval3(x[i], mem_zigzag[k, ])

     if (k > 2){
       for (j in 2:(k-1)){
         mem_value[i, j] <- memval2(x[i], mem_zigzag[j, ])
       }
     }
   }
  return(mem_value)
}

terms_count <-  function(data){
  m <- dim(data)[2]
  terms_num <- numeric(m) 
  
  for (i in 1:m){
    terms_num[i] <- cluster(data[,i], max_clu)[1]
  }
  return(terms_num)
}

mem_fun_result <-  function(data){
  m <- dim(data)[2]
  terms <- terms_count(data)
  k <- max(terms)
  mem_fun <- matrix(0, m, k+3)
 
  for (i in 1:m){
    u <- cluster(data[,i], max_clu)
    kk <- length(u)
    mem_fun[i, 1:kk] <- u
  }
  return( mem_fun)
}


mem_data <- function(data){
  ## To calculate the membership value of the data for each term in each process

  n <- dim(data)[1]; m <- dim(data)[2] 
  terms_num <- terms_count(data)
  mem_val <- array(0, dim <- c(m, max(terms_num), n))
  for (j in 1:m){
    mem_val[j, 1:terms_num[j], ] <- t(mem_fun_val(data[, j]))
  }

   return(mem_val)
}



item_set1 <- function(mem_val, alpha){
  ## L1; to find the 1-itemset
   n <- dim(mem_val)[3];  m <- length(alpha)
   terms_num <- terms_count(data)
   count <- matrix(0, m, max(terms_num))
   for (i in 1:m){ 
      for (j in 1:terms_num[i]){
          a <- mem_val[i, j, ]
          count[i, j] <- sum(a)
      }
   }   
   k_sum <- 0
   itemset1 <- matrix(0, m, max(terms_num))
   for (i in 1:m){
     ind1 <- which(count[i, ] > alpha[i])      
     itemset1[i, ind1] <- 1
   }
   return(itemset1)
}



itemset_value <- function(mem_val, index1){
  ## index here is a vector
  ## mem_val is a array; not a matrix
  ## The output y is also a vector

  n <- dim(mem_val)[3]; m <- dim(mem_val)[2]
  L <- length(index1)/2
  y <- numeric(n)
  
  for(i in 1:n){
    u <- 1
    for (j in 1:L){
      u <- u * mem_val[index1[2*j-1], index1[2*j], i] 
    }
    y[i] <- u
  }
  return(y)
}




item_set2 <- function(mem_val, itemset1, alpha){
  ## L2; to find the 2-itemset
  n <- dim(mem_val)[3];  m <- length(alpha)
  terms_num <- terms_count(data)
  itemset_new <- matrix(0, m, max(terms_num))

  k <- numeric(m)
  for (i in 1:m){
    k[i] <- length(which(itemset1[i,] > 0))
  }
  process_item <- which(k > 0)
  k1 <- k[process_item]
  m1 <- length(process_item)
  item_num <- 0
  for (l in 1:m1){
     item_num <- item_num + sum(k1[-(1:l)])*k1[l]
  }   
  index <- matrix(0, item_num, 4)
  index_new <- matrix(0, item_num, 4)
  yy <- numeric(item_num)
  u <- 0; cc <- 0
  if (m1 > 2){
    for (i in 1:(m1-1)){
      u1 <- process_item[i]
      for (i2 in 1:k[u1]){
        u2 <- which(itemset1[u1, ]>0)[i2]
        for (j in (i+1):m1){
          u3 <- process_item[j]
          for (j2 in 1:k[u3]){
             u4 <- which(itemset1[u3, ]>0)[j2]
             u <- u+1
             index[u, ] <- c(u1, u2, u3, u4)

             yy[u] <- sum(itemset_value(mem_val, index[u, ]))

             if (yy[u] > alpha[u2]*alpha[u4]/n){
               cc <- cc+1                    
               index_new[cc, ] <- c(u1, u2, u3, u4); 
              # set1 <- which(itemset_new[u1, ] > 0)
              # set2 <- which(itemset_new[u3, ] > 0)
              # itemset_new[u1, ] <- c(itemset_new[u1, set1], u2)
              # itemset_new[u3, ] <- c(itemset_new[u1, set2], u4)
                          
            }
          }
        }
      }      
    }
  }  
  index2 <-  index_new[1:cc, ]
  return(index2)
}



unique_index <- function(index){
  n <- dim(index)[1]; m <- dim(index)[2]/2
  index.set <- matrix(0, n*m, 2)
  for (j in 1:m){
    index.set[((j - 1)*n + 1) :(j*n), ] <- index[, (2*j-1):(2*j)]
  }
  unique_set <- index.set[!duplicated(index.set), ]
  return(unique_set)
}

rowbind <- function(index){
  n <- dim(index)[1]; m <- dim(index)[2]
  a <- matrix(0, 1, n*m)
  for (i in 1:n){
    a[(1+i*m-m):(i*m)] <- index[i, ]
  }
  return(a)
} 

index <- function(a1, a2){
## This function is to find the row number of a1 in matrix a2
## the column of a1 and a2 are the same
  n1 <- dim(a1)[1]; m1 <- dim(a1)[2]
  n2 <- dim(a2)[1];
  ind <- numeric(n1) 
  for (i in 1:n1){
     for (j in 1:n2){
        if (all(a1[i, ] == a2[j, ])){
          ind[i] <- j
        }
     }
  }
  return(ind)
}

vect_odd <- function(a){
## To select the odd index in vector in a 
  n <- length(a)
  b <- numeric(ceiling(n/2))
  for (i in 1:(ceiling(n/2))){
     b[i] <- a[2*i-1]
  }
  return(b)
}

vect_same <- function(a1, a2){
## To test whether the elements in a1 and a2 are the same, not considering the order
  y <- 0 
  m <- length(a1)/2
  z <- matrix(0, m, m)

  for (i in 1:m){
     for (j in 1:m){
        if (all(a1[(2*i-1):(2*i)] == a2[(2*j-1):(2*j)])){
           z[i, j] <- 1
        }
     }    
  }
  if (min(apply(z, 1, max)) == 1) y <- 1
  return(y)
}




unique_matrix <- function(item_set){
   n <- dim(item_set)[1]; m <- dim(item_set)[2]
   
   c <- 0
   u <- matrix(0, n, n)
   d <- numeric(n)
   for (i in 1:(n-1)){
      u1 <- item_set[i, ]
      ind1 <-  vect_odd(u1)
      if (length(unique(ind1)) < m/2) {
          c <- c + 1;
          d[c] <- i
      }
     else{
        for (j in (i+1):n){
           v1 <- item_set[j, ] 
           u[i, j] <- vect_same(u1, v1)
           if (u[i, j] == 1) {
              c <- c + 1;
              d[c] <- j
           }             
        } 
     }      
   } 
   u1 <- item_set[n, ]
   ind1 <-  vect_odd(u1)
   if (length(unique(ind1)) < m/2) {
        c <- c + 1;
        d[c] <- n
   }

   ind <- unique(d[1:c])
   item_new <- item_set[-ind, ]
   return(item_new)
}



item_set3 <- function(mem_val, itemset,  alpha){
  index_new <- itemset
  samplesize <- dim(mem_val)[3]
  n <- dim(index_new)[1]; m <- dim(index_new)[2]  
  a <- matrix(0, m/2, 2)
  
  process <- unique_index(index_new)
  k <- dim(process)[1]  
  ind <- c(1:k)
  item_set <- matrix(0, n*(k - m/2), m+2)

  t <- 0
  for (i in 1:n){    
     for (j in 1:(m/2)){
        a[j, ] <- index_new[i, ((j-1)*2 + 1 ):(2*j)]        
     }
     ind1 <- index(a, process)
     if (process[max(ind1), 1] < max(process[, 1]) ){
        ind2 <- which(process[, 1] > process[max(ind1), 1])
        k1 <- length(ind2)
        for (l in 1: k1){
           t <- t+1
           item_set[t, ] <- c(index_new[i, ], process[ind2[l],])
        }
     }
  }
  #L_set <- unique_matrix(item_set[1:t, ])
  L_set <- item_set[1:t, ]
  
  item_new <- L_set
  n1 <- dim(L_set)[1]
  yy <- numeric(n1)
  cc <- 0
  for (k in 1:n1){
     yy[k] <- sum(itemset_value(mem_val, L_set[k, ]))
     ind3 <- vect_odd(L_set[k, ])
     if (yy[k] > prod(alpha[ind3])/samplesize){
          cc <- cc+1                    
          item_new[cc, ] <- L_set[k, ] 
     } 
  }
  if (cc > 0){
    index2 <-  L_set[1:cc, ]
  }
  if (cc == 0) {index2 <- 0}
  return(index2)
}






asso_rules <- function(mem_val, index, lambda){
## Here index is a vector
## index is the rule
  n <- dim(mem_val)[1]; 
  n2 <- length(index)/2
  rules <- matrix(0, n2, n2*2)
  c <- 0
  for (i in 1:n2){
    u <- (1+(i-1)*2):(2*i)
    ind1 <- index[-u]
    if( (sum(itemset_value(mem_val, ind1))>0) & (sum(itemset_value(mem_val, index))/sum(itemset_value(mem_val, ind1))> lambda) ){
       c <- c+1
       rules[c, ] <- c(index[-u], index[u])      
    }
  }
  if (c == 0) {result <- 0}
  if (c > 0) {result <- rules[1:c, ]}
  return(result)
}


asso_rules_set <- function(itemset, mem_val, lambda){
   n1 <- dim(itemset)[1]; m1 <-  dim(itemset)[2]
   rules <- matrix(0, n1*m1, m1)
   c <- 0
   for (i in 1:n1){
      index <- itemset[i, ]
      u <- asso_rules(mem_val, index, lambda)
      
      if (max(u) > 0){   
          if (is.null(dim(u)[1])) { 
             rules[c+1, ] <- u
             c <- c + 1
          }
          else{      
             rules[(c+1):(c+dim(u)[1]), ] <- u
             c <- c + dim(u)[1]
          }
      }
   }
   if (c == 0) {result <- 0}
   if (c > 0) {result <- rules[1:c, ]}
   return(result)
}
  

rule_length <- function(result){
## To record the index of the obtained rules
 
   n1 <- dim(result)[1]; m <- dim(result)[2]
   ind1 <- which(apply(result, 1, min) == 0)

   rule_len<- matrix(0, nrow <- n1, ncol <- 2)
   rule_len[, 1] <- rep(1:n1)
   for (i in 1: n1){
      rule_len[i, 2] <- max(which(result[i, ] >0))  
   }
   return(rule_len)
}





rule_interest <- function(result1, proc_interest){
## This fucntion is to select the interested rules
  n1 <- dim(result1)[1]
  result_new <- result1
  k <- 0
  for (i in 1:n1){
    l1 <- length(which(result1[i, ] > 0))
    ind1 <- which(result1[i, ] > 0)[l1-1]
    if (result1[i, ind1] == proc_interest){
      k <- k+1
      result_new[k, ] <- result1[i, ]      
    }  
  }
  if (k == 0) {result_new <- 0} else{result_new <- result_new[1:k, ]}
  return(result_new)
}
  


rule_index <- function(result){
## To record the index of the obtained rules
 
   n1 <- dim(result)[1]; m <- dim(result)[2]
   ind1 <- which(apply(result, 1, min) == 0)

   rule_ind <- matrix(0, nrow <- n1, ncol <- 2)
   rule_ind[, 1] <- rep(1:n1)
   for (i in 1: n1){
      rule_ind[i, 2] <- max(which(result[i, ] >0))  
   }
   return(rule_ind)
}


condifence <- function(data, alpha, lambda, result){
## This function is to calculate the confidnece levele of the proposed algorithm
## result <- main_func(data, alpha, lambda)
   n <- dim(data)[1]
   mem_result <- mem_fun_result(data)
   
   
   n1 <- dim(result)[1]; m <- dim(result)[2]
 

   rule_ind <- rule_index(result)  

   rule_cluster <- rule_ind[, 2][!duplicated(rule_ind[, 2])]
   rule_cluster_num <- length(rule_ind[, 2][!duplicated(rule_ind[, 2])])
 
   ind_rule_cluster <- numeric(rule_cluster_num + 1)

   for (i in 1: rule_cluster_num){
     ind_rule_cluster[i] <- min(which(rule_ind[, 2]== rule_cluster[i]))
   }
   ind_rule_cluster[rule_cluster_num + 1] <- n1+1

  
  p_support <- matrix(0, n1, 3)

  t <- 0
  for (i in 2:(rule_cluster_num+1)){
     index1 <- ind_rule_cluster[i-1]:(ind_rule_cluster[i] - 1)
     a <- result[index1, 1:rule_cluster[i-1]]
   
     a_if <- result[index1, 1:(rule_cluster[i-1]-2)]
     a_then <- result[index1, (rule_cluster[i-1] - 1):(rule_cluster[i-1])]
  
     for (j in 1:length(index1)){
        t <- t + 1
        rule <- a[j, ]      
        if_condition <- matrix(0, length(rule)/2 - 1, 2);      
        for (s in 1:(length(rule)/2 - 1)){
           if_condition[s, 1] <- rule[2*s-1]
           if_condition[s, 2] <- rule[2*s]
        }
        then_condition <- a_then[j, ]
        items <- dim(if_condition)[1]
        endpoint <- matrix(0, items+1, 3)         
        count1 <- numeric(n)
        count2 <- numeric(n)

        for (k in 1:n){      
           record <- numeric(items + 1)

           for (kk in 1:items){
              process_if <- if_condition[kk, 1]
              term_if <- if_condition[kk, 2]; 
              endpoint[kk, ] <- mem_result[if_condition[kk, 1], (term_if + 1):(term_if + 3)]

              if ((data[k, process_if] >= endpoint[kk, 1]) & (data[k, process_if] <= endpoint[kk, 3])){
                 record[kk] <- 1
              }        
            }
        
            process_then <- then_condition[1]
            term_then <- then_condition[2]
            endpoint[items+1, ] <- mem_result[then_condition[1], (term_then + 1):(term_then + 3)]
            if ((data[k, process_then] >= endpoint[items+1, 1]) & (data[k, process_then] <= endpoint[items+1, 3])){
                 record[items+1] <- 1
            }  

           count1[k] <- prod(record[-(items+1)])
           count2[k] <- prod(record)
        }
       p_support[t, 1] <- sum(count1); p_support[t, 2] <- sum(count2)
       p_support[t, 3] <- sum(count2)/sum(count1)
     }
  }
  return(p_support)
}


main_func <- function(data, alpha, lambda){

  ## x is the observed data
  ## alpha is a vector, and lambda is a number

  n <- dim(data)[1]; m <- dim(data)[2] 

  mem_val <- mem_data(data)
  L1_set <- item_set1(mem_val, alpha) 
  L2_set <- item_set2(mem_val, L1_set, alpha)
  rules <- numeric(2)

  u <- 0
  if (max(L2_set) > 0){
     result2 <- asso_rules_set(L2_set, mem_val, lambda)
     if (max(result2) > 0) {
       u <- u+ dim(result2)[1]; u1 <- dim(result2)[2]  
       rules <- result2    
     }
  
     itemset <- item_set3(mem_val, L2_set, alpha)
     if (max(itemset) > 0) {
        result <- asso_rules_set(itemset, mem_val, lambda)
        rules2 <- matrix(0, nrow <- u, ncol <- dim(result)[2])
        rules2[, 1:u1] <- rules
        if (max(result) > 0) {
            u <- u+ dim(result)[1]; u1 <- dim(result)[2]
            rules <- rbind(rules2, result)         
        }       
     }

    kk <- 0
    while(max(itemset) > 0){
       kk <- kk+1
       itemset_new <- itemset
       itemset <- item_set3(mem_val, itemset_new, alpha)
       if (max(itemset) != 0) {
         result_new <- asso_rules_set(itemset, mem_val, lambda)     

         if (max(result_new) > 0) {
            u <- u + dim(result_new)[1]; u1 <- dim(rules)[2]
            rules2 <- matrix(0, nrow <- u, ncol <- dim(result_new)[2])
            rules2[1:dim(rules)[1], 1:u1] <- rules
            rules2[(dim(rules)[1]+1):u, ] <- result_new

            rules <- rules2
            result <- result_new          
         }
       }
    }
  } 
 return(rules)
}

#Read the data file in txt
data  <- read.table(file="/data/data_v3.txt", header=FALSE)

#Settings of the parameters
data <- as.matrix(data)
n <- dim(data)[1]
p <- dim(data)[2]

proc_interest <- 1

max_clu <- 4
lambda <- 0.7
alpha <- numeric(p)+n/50

#Execution of the proposed functions

corrana_data(data)
mem_result <- mem_fun_result(data)

time0 = proc.time()
result1 <- main_func(data, alpha, lambda)
result <- rule_interest(result1, proc_interest)
c <- condifence(data, alpha, lambda, result)
proc.time()-time0

result
c