
### Trajectory analysis on MDA clones ###

# Functions at: https://github.com/chenweng1991/redeemR/blob/master/R/BuidTree.R



#' Tomerge_v2
#'
#' This function is to quickly merge two dataframe by rownames, but can choose to leave A or B all information
#' @param A  dataframe A
#' @param B  dataframe B
#' @return  return a data frame with merged information
#' @export
#' @examples
#' Gettopgenes(XX.ob,number)
Tomerge_v2<-function(A,B,leavex=T,leavey=F){
  mergeAB<-merge(A,B,by="row.names",all.x=leavex,all.y=leavey)
  row.names(mergeAB)<-mergeAB[,1]
  mergeAB<-mergeAB[,-1]
  return(mergeAB)
}



#' Filter genes
#' 
#' This function is used to filter out genes with zero-expression in too many samples
#' 
#' @param 
#' 
#' @return 
filter_genes <- function(input_expr_df, clones = "GBC", perc_thr = 10){
  n_clones <- length(input_expr_df[, clones])
  n_clones_thr <- ceiling((n_clones / 100) * perc_thr)
  
  n_zeros <- colSums(input_expr_df[,14:dim(input_expr_df)[2]] == 0) 
  features <- names(n_zeros)[!n_zeros>n_clones_thr] # Retain only genes with non-zero expr in at least 90% of clones
  
  out_expr_df <- input_expr_df[, c(colnames(input_expr_df)[1:13], features)] 
  return(out_expr_df)
}


#' Convert to unique data.frame 
#' 
#' @param  
#' @return 
convert_results_to_df <- function(LinOut.result, x_cond = F){
  if(!x_cond){
    LinOut.result.df <- merge(LinOut.result$slopes, LinOut.result$ps, by = "row.names") %>%
      merge(., LinOut.result$qs, by.x = "Row.names", by.y = "row.names", all = T)
    colnames(LinOut.result.df) <- c("genes", "slopes", "ps", "qs")
    
    return(LinOut.result.df)    
  } else if(x_cond){
    for(exp_cond in names(LinOut.result)){
      LinOut.result.cond <- LinOut.result[[exp_cond]]
      LinOut.result.cond.df <- merge(LinOut.result.cond$slopes, LinOut.result.cond$ps, by = "row.names") %>%
        merge(., LinOut.result.cond$qs, by.x = "Row.names", by.y = "row.names", all = T)
      colnames(LinOut.result.cond.df) <- c("genes", "slopes", "ps", "qs")
      
      LinOut.result[[exp_cond]] <- LinOut.result.cond.df
    }
    return(LinOut.result)
  }
}


#' Run_Lin_regression
#' 
#' Run Poisson regression on clone-level gene expression data
#' 
#' @param LinOut input data.frame with clone-level gene expression data plus metastatic potential
#' @param regress_factor  
#' @param n.cores  =8
#' @param qval compute q-values additionally to p-values. Sometimes q-values computation is not successful. Keep p-values then.
#' @param tot_UMIs include total n. of UMIs x clone instead of median. This allows to have a functioning poi model.
#' @export
#' @import foreach doParallel 
Run_Lin_regression<-function(LinOut, 
                             regress_factor = c("met_potential"),
                             n.cores = 8, 
                             qval = F, 
                             tot_UMIs = F){ 
  library(dplyr)
  library(foreach)
  library(doParallel)
  library(qvalue)
  
  genes=colnames(LinOut)[14:ncol(LinOut)]
  
  # To allow the model to correct for clone-size, compute total n. of counts x clone
  if(tot_UMIs == T){
    clonesize.umi <- LinOut[,14:dim(LinOut)[2]] %>% rowSums 
    LinOut <- LinOut %>% mutate(tot_UMIs = clonesize.umi) %>% 
      relocate(., tot_UMIs, .after = nUMIs)
  }
  
  # used to set up a parallel computing environment
  my.cluster <- parallel::makeCluster(n.cores) # This command creates a cluster of worker processes for parallel computation. n.cores specifies the number of CPU cores you want to use for parallel computation.
  print(my.cluster)
  doParallel::registerDoParallel(cl = my.cluster) # This command registers the parallel backend for the foreach package, allowing you to execute foreach loops in parallel.
  
  res<-foreach(gene=genes[1:length(genes)]) %dopar%{
    ps<-c()
    slopes<-c()
    for(i in 1:length(regress_factor)){ # For each factor for which we want to perform poi regression (i.e. for us, only clonal fitness)
      if(tot_UMIs == T){
        f<-as.formula(paste(gene,"~",regress_factor,"+ log(tot_UMIs)")[i])
      }else{
        f<-as.formula(paste(gene,"~",regress_factor,"+ log(nUMIs)")[i])
      }
      md<-lm(f,data=LinOut)
      md.summary<-summary(md)
      slope<-md.summary$coefficients[2,1] # Save slope from fitted model
      p<-md.summary$coefficients[2,4] # Save p-value
      slopes<-c(slopes,slope) # Store each slope, of model fitted on gene ~ each regress_factor
      ps<-c(ps,p) # Store each p-value
    }
    return(list(Gene=gene,ps=ps,slopes=slopes)) # res will be a list of lists (?)
  }
  parallel::stopCluster(cl = my.cluster)
  
  if(length(regress_factor) > 1){
    genes<-sapply(res,function(x){x[[1]]})
    ps<-sapply(res,function(x){x[[2]]}) %>% t %>% as.data.frame 
    slopes<-sapply(res,function(x){x[[3]]}) %>% t %>% as.data.frame    
  } else if(length(regress_factor) == 1){
    genes<-sapply(res,function(x){x[[1]]})
    ps<-sapply(res,function(x){x[[2]]}) %>% as.data.frame 
    slopes<-sapply(res,function(x){x[[3]]}) %>% as.data.frame
  }
  
  row.names(ps)<-genes
  row.names(slopes)<-genes
  colnames(ps)<-regress_factor
  colnames(slopes)<-regress_factor
  
  # Save either q-values of p-values based on choice 
  if(qval == F){
    return(list(ps=ps,slopes=slopes))  
  }else{
    qs <- apply(ps, 2, function(x){p.adjust(x, method = "BH")})  %>% as.data.frame
    return(list(ps=ps,qs=qs,slopes=slopes))
  }
}



#' Run_Lin_regression_poi
#' 
#' Run Poisson regression on clone-level gene expression data
#' 
#' @param LinOut input data.frame with clone-level gene expression data plus metastatic potential
#' @param regress_factor  
#' @param n.cores  =8
#' @param qval compute q-values additionally to p-values. Sometimes q-values computation is not successful. Keep p-values then.
#' @param tot_UMIs include total n. of UMIs x clone instead of median. This allows to have a functioning poi model.
#' @export
#' @import foreach doParallel 
Run_Lin_regression_poi<-function(LinOut, 
                                 regress_factor = c("met_potential"),
                                 n.cores = 8, 
                                 qval = F, 
                                 tot_UMIs = F){ 
  library(dplyr)
  library(foreach)
  library(doParallel)
  library(qvalue)
  
  genes=colnames(LinOut)[14:ncol(LinOut)]
  
  # To allow the model to correct for clone-size, compute total n. of counts x clone
  if(tot_UMIs == T){
    clonesize.umi <- LinOut[,14:dim(LinOut)[2]] %>% rowSums 
    LinOut <- LinOut %>% mutate(tot_UMIs = clonesize.umi) %>% 
      relocate(., tot_UMIs, .after = nUMIs)
  }
  
  # used to set up a parallel computing environment
  my.cluster <- parallel::makeCluster(n.cores) # This command creates a cluster of worker processes for parallel computation. n.cores specifies the number of CPU cores you want to use for parallel computation.
  print(my.cluster)
  doParallel::registerDoParallel(cl = my.cluster) # This command registers the parallel backend for the foreach package, allowing you to execute foreach loops in parallel.
  
  res<-foreach(gene=genes[1:length(genes)]) %dopar%{
    ps<-c()
    slopes<-c()
    for(i in 1:length(regress_factor)){ # For each factor for which we want to perform poi regression (i.e. for us, only clonal fitness)
      if(tot_UMIs == T){
        f<-as.formula(paste(gene,"~",regress_factor,"+ log(tot_UMIs)")[i])
      }else{
        f<-as.formula(paste(gene,"~",regress_factor,"+ log(nUMIs)")[i])
      }
      md<-glm(f,data=LinOut,family=poisson(link="log"))
      md.summary<-summary(md)
      slope<-md.summary$coefficients[2,1] # Save slope from fitted model
      p<-md.summary$coefficients[2,4] # Save p-value
      slopes<-c(slopes,slope) # Store each slope, of model fitted on gene ~ each regress_factor
      ps<-c(ps,p) # Store each p-value
    }
    return(list(Gene=gene,ps=ps,slopes=slopes)) # res will be a list of lists (?)
  }
  parallel::stopCluster(cl = my.cluster)
  
  if(length(regress_factor) > 1){
    genes<-sapply(res,function(x){x[[1]]})
    ps<-sapply(res,function(x){x[[2]]}) %>% t %>% as.data.frame 
    slopes<-sapply(res,function(x){x[[3]]}) %>% t %>% as.data.frame    
  } else if(length(regress_factor) == 1){
    genes<-sapply(res,function(x){x[[1]]})
    ps<-sapply(res,function(x){x[[2]]}) %>% as.data.frame 
    slopes<-sapply(res,function(x){x[[3]]}) %>% as.data.frame
  }
  
  row.names(ps)<-genes
  row.names(slopes)<-genes
  colnames(ps)<-regress_factor
  colnames(slopes)<-regress_factor
  
  # Save either q-values of p-values based on choice 
  if(qval == F){
    return(list(ps=ps,slopes=slopes))  
  }else{
    qs <- apply(ps, 2, function(x){p.adjust(x, method = "BH")})  %>% as.data.frame
    return(list(ps=ps,qs=qs,slopes=slopes))
  }
}



#' Run_Lin_regression_nb
#' 
#' Run Negative Binomial regression on clone-level gene expression data
#' Negative Binomial reg. is more suitable then Poisson reg. when overdispersion is present
#' 
#' @param LinOut input data.frame with clone-level gene expression data plus metastatic potential
#' @param regress_factor  
#' @param n.cores  =8
#' @param qval compute q-values additionally to p-values. Sometimes q-values computation is not successful. Keep p-values then.
#' @param tot_UMIs include total n. of UMIs x clone instead of median. This allows to have a functioning poi model.
#' @export
#' @import foreach doParallel 
Run_Lin_regression_nb<-function(LinOut, 
                                regress_factor = c("met_potential"),
                                n.cores = 8,
                                qval = F, 
                                tot_UMIs = T){
  library(dplyr)
  library(foreach)
  library(doParallel)
  library(qvalue)
  library(MASS)
  
  genes=colnames(LinOut)[14:ncol(LinOut)]
  
  # To allow the model to correct for clone-size, compute total n. of counts x clone
  if(tot_UMIs == T){
    clonesize.umi <- LinOut[,14:dim(LinOut)[2]] %>% rowSums 
    LinOut <- LinOut %>% mutate(tot_UMIs = clonesize.umi) %>% 
      relocate(., tot_UMIs, .after = nUMIs)
  }
  
  # used to set up a parallel computing environment
  my.cluster <- parallel::makeCluster(n.cores) # This command creates a cluster of worker processes for parallel computation. n.cores specifies the number of CPU cores you want to use for parallel computation.
  print(my.cluster)
  doParallel::registerDoParallel(cl = my.cluster) # This command registers the parallel backend for the foreach package, allowing you to execute foreach loops in parallel.
  
  res<-foreach(gene=genes[1:length(genes)]) %dopar%{
    ps<-c()
    slopes<-c()
    for(i in 1:length(regress_factor)){ # For each factor for which we want to perform poi regression (i.e. for us, only clonal fitness)
      if(tot_UMIs == T){
        f<-as.formula(paste(gene,"~",regress_factor,"+ log(tot_UMIs)")[i])
      }else{
        f<-as.formula(paste(gene,"~",regress_factor,"+ log(nUMIs)")[i])
      }
      md<-MASS::glm.nb(f,data=LinOut)
      md.summary<-summary(md)
      slope<-md.summary$coefficients[2,1] # Save slope from fitted model
      p<-md.summary$coefficients[2,4] # Save p-value
      slopes<-c(slopes,slope) # Store each slope, of model fitted on gene ~ each regress_factor
      ps<-c(ps,p) # Store each p-value
    }
    return(list(Gene=gene,ps=ps,slopes=slopes)) # res will be a list of lists (?)
  }
  parallel::stopCluster(cl = my.cluster)
  
  if(length(regress_factor) > 1){
    genes<-sapply(res,function(x){x[[1]]})
    ps<-sapply(res,function(x){x[[2]]}) %>% t %>% as.data.frame 
    slopes<-sapply(res,function(x){x[[3]]}) %>% t %>% as.data.frame    
  } else if(length(regress_factor) == 1){
    genes<-sapply(res,function(x){x[[1]]})
    ps<-sapply(res,function(x){x[[2]]}) %>% as.data.frame 
    slopes<-sapply(res,function(x){x[[3]]}) %>% as.data.frame
  }
  
  row.names(ps)<-genes
  row.names(slopes)<-genes
  colnames(ps)<-regress_factor
  colnames(slopes)<-regress_factor
  
  # Save either q-values of p-values based on choice 
  if(qval == F){
    return(list(ps=ps,slopes=slopes))  
  }else{
    qs <- apply(ps, 2, function(x){p.adjust(x, method = "BH")})  %>% as.data.frame
    return(list(ps=ps,qs=qs,slopes=slopes))
  }
}



#' Volcano plot of results
#' 
#' Plot a Volcano with q-value ~ slopes of genes resulting from poi regression
#' 
#' @param 
#' @return 
PlotLinRegress_Vocano <- function(LinOut.result.df, slot = "qs", pline = 0.001, qline = 0.05){
  library(ggpubr)
  library(purrr)
  
  datatoplot <- LinOut.result.df
  ifelse(slot == "qs", datatoplot <- datatoplot %>% mutate(score = -log(qs)*abs(slopes)), 
         datatoplot <- datatoplot %>% mutate(score = -log(ps)*abs(slopes)))
  
  #Label <- subset(datatoplot, qs < 0.2) %>% .[order(.$score,decreasing=T),] %>% .[1:80,] %>% .[, "genes"]
  #slope.L <- subset(datatoplot, qs < 0.2)$slopes %>% min
  #slope.R <- subset(datatoplot,qs < 0.2)$slopes %>% max
  genes_of_int <- c("PAEP", "CST1")
  Label <- datatoplot[order(datatoplot$score,decreasing=T),] %>% .[1:20,] %>% .[, "genes"]
  Label <- c(Label, genes_of_int)
  slope.L <- datatoplot$slopes %>% min
  slope.R <- datatoplot$slopes %>% max
  datatoplot$Label <- ifelse(datatoplot$genes %in% Label, datatoplot$genes, "")
  Name = "volcano"
  
  if(slot=="qs"){
    # Cap q-values to min q-value to avoid having zeros (leading to -log10(0) = Inf)
    datatoplot$qs[datatoplot$qs == 0] <- min(datatoplot$qs[datatoplot$qs != 0])
    qval.top <- -log10(datatoplot$qs) %>% max()
    p <- ggplot(datatoplot)+aes(slopes, -log10(qs), label = Label)+
      geom_point()+
      geom_hline(yintercept=-log10(qline),linetype=2)+
      geom_vline(xintercept=c(-0.02,0.02),linetype=2)+
      xlim(slope.L-0.1,slope.R+0.1)+
      ylim(0, qval.top+10)+
      geom_text_repel(force=20,size=6,max.overlaps=100)+theme_pubr()+ggtitle(Name)
  }else{
    # Cap p-values to min p-value to avoid having zeros (leading to -log10(0) = Inf)
    datatoplot$ps[datatoplot$ps == 0] <- min(datatoplot$ps[datatoplot$ps != 0])
    pval.top <- -log10(datatoplot$ps) %>% max()
    p <- ggplot(datatoplot)+aes(slopes,-log10(ps),label=Label)+
      geom_point()+
      geom_hline(yintercept=-log10(pline),linetype=2)+
      geom_vline(xintercept=c(-0.02,0.02),linetype=2)+
      xlim(slope.L-0.1,slope.R+0.1)+
      ylim(0, pval.top+10)+
      geom_text_repel(force=20,size=6,max.overlaps=100)+theme_pubr()+ggtitle(Name)
  }
}


