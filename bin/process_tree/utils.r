# Code
library(ape)
library(phangorn)
library(tidyverse)
library(foreach)
library(doParallel)
library(doMC)
library(pryr)


##


split_profile <- function(profile){
  as.integer(unlist(strsplit(profile,split="")))
}


##


add_derived_profile_info <- function(profile_df, samples=sprintf("s%s",0:(nchar(profile_df$profile[1])-1))){

  base=rep(0,nchar(profile_df$profile[1]))
  samples_private=sapply(1:length(base),function(i){this=base;this[i]=1;paste(this,collapse="")})
  missing_samples=setdiff(samples_private,profile_df$profile)
  if(length(missing_samples)>0){
    profile_df=rbind(profile_df,data.frame(profile=missing_samples,edge_length=0))
  }
  profile_df$mut_count=nchar(profile_df$profile)-nchar(gsub("1","",profile_df$profile))
  profile_df$profile_int=lapply(profile_df$profile,split_profile)
  profile_df$id=1:dim(profile_df)[1]
  profile_df$label=profile_df$profile
  idx=which(profile_df$mut_count==1)
  profile_df$label[idx]=sapply(idx,function(i) samples[which(profile_df$profile_int[[i]]==1)])
  
  return(profile_df)

}


##


get_ancestral_nodes <- function(node,edge,exclude_root=TRUE){

  idx=which(edge[,2]==node)
  parents=node ##Include the node
  while(length(idx)>0){
    if(length(idx)>1){
      stop("multiple parents!")
    }
    parent=edge[idx,1]
    parents=c(parents,parent)
    idx=which(edge[,2]==parent)
  }
  if(exclude_root){
    parents[-length(parents)] ##The last node is the root.
  }else{
    parents
  }

}



##


reconstruct_genotype_summary <- function(phylo){
  
  dat <- phylo$edge
  samples <- phylo$tip.label
  N <- length(samples)
  zeros <- rep(0,N)
  branches <- lapply(1:dim(dat)[1],function(i) c())                     # List. One element per phylo$edge row (edge). 1:n_edges numbers.
  map_node <- match(1:max(dat[,2]),dat[,2])                         # Map each node to its branch (i.e., the branch where the node is a child)
  for(i in 1:N){                                                    # Each tip (index) is mapped to each of its parents' branch
    parents <- get_ancestral_nodes(i, edge=dat, exclude_root=TRUE)
    for(j in parents){
      branches[[map_node[j]]] <- append(branches[[map_node[j]]], i)
    }
  }
  profiles <- sapply(branches, function(x){ tmp=zeros; tmp[x]=1; paste(tmp, collapse="") })
  df <- data.frame(profile=profiles, edge_length=phylo$edge.length, stringsAsFactors=FALSE)
  df <- add_derived_profile_info(df, phylo$tip.label)
  
  return(df)

}


##


assign_variants <- function(tree, mtr, mtr.bi, depth, n.cores=8, assign.to='Child'){
  
  # Get cell profiles
  df <- reconstruct_genotype_summary(tree)
  df_profile_mtx.t <- df$profile_int %>% do.call(rbind,.) %>% t

  ## Set expected global heteroplasmy level 
  heteroplasmy <- matrix(mtr/(depth+0.000001))
  Global.heteroplasmy <- heteroplasmy[!is.na(heteroplasmy)] %>% .[.!=0] %>% median  # this is a number ~ 0.03-0.05, an expected p for modeling dropout
  
  ## Generate a binomial-based probability distribution for each cell to model the dropout
  Zero.p <- dbinom(0,depth,Global.heteroplasmy) # PROBABILITY OF OBSERVING AF=0, with a global AF ~0.05
  Zero.p[Zero.p==1] <- 0.99
  Zero.logV <- log(Zero.p)  # A matrix of log(probability of zero given w/ mutation in the cell)  Variant vs cell
  NonZero.logV <- log(1-Zero.p) # A matrix of log(probability of 1 or >=1 given w/ mutation in the cell)  Variant vs cell
  Zero.logP <- log(0.95) # A matrix of log(probability of zero given no mutation in the cell)  Variant vs cell
  NonZero.logP <- log(1-0.95) # A matrix of log(probability of 1 or >=1 given no mutation in the cell)  Variant vs cell
  
  # Calculate LogLihelihoods for each variant
  registerDoMC(n.cores)
  Loglik <- foreach (i=1:nrow(mtr.bi), .combine='rbind') %dopar%  # To loop through all variants
  {
    ## Make indicator matrices for 1-1 (Inside the clade, and is mutation) and 0-0 (outside of clade, and is not mutation) 
    x <- df_profile_mtx.t + mtr.bi[i,]
      
    x_11 <- matrix(nrow=nrow(x),ncol=ncol(x)) # cell x branch
    x_00 <- matrix(nrow=nrow(x),ncol=ncol(x)) # cell x branch
    x_11[which(x!=2)] <- 0                   
    x_11[which(x==2)] <- 1                   # Cell inside the branch with mutation 
    x_00[which(x!=0)] <- -1  
    x_00[which(x==0)] <- 1                   # Cell outside the branch and no mutation 
    x_00[which(x_00==-1)]<-0
      
    ## Make indicator matrices for 1-0 and 0-1 scenarios
    x_10 <- df_profile_mtx.t - mtr.bi[i,] 
    x_10[which(x_10==-1)] <- 0                # Cell inside the branch and no mutation 
    x_01 <- mtr.bi[i,]-df_profile_mtx.t
    x_01[which(x_01==-1)] <- 0                # Cell outside the branch with mutation 
      
    ## Compute the Loglikihood of mutation i across branches
    Loglik_11 <- x_11*NonZero.logV[i,]
    Loglik_10 <- x_10*Zero.logV[i,]
    Loglik_01 <- x_01*NonZero.logP
    Loglik_00 <- x_00*Zero.logP
    Loglik.i <- colSums(Loglik_11+Loglik_10+Loglik_01+Loglik_00)
      
    return(Loglik.i)

  }
  
  # Get final probabilities of assignment
  row.names(Loglik) <- row.names(mtr.bi)
  n <- dim(mtr.bi)[1]
  edge_ml <- apply(Loglik,1,which.max) %>% as.numeric
  p <- exp(Loglik-Loglik[(as.numeric(edge_ml)-1)*n+1:n])
  p <- p/rowSums(p)
  
  # Create the final report
  if (assign.to=="Parent"){edge.node=1} else {edge.node=2}
  var_assignment <- data.frame(
    temp=tree$edge[,edge.node][apply(p,1,which.max)], # NBB: final assignment of each variant to the CHILD NODE of each branch
    p=apply(p,1,max)                                   
  )
  var_assignment <- setNames(var_assignment, c(assign.to, 'p')) %>% arrange(desc(p))
  
  # Add variants AF and prevalence within and outside the assigned branches tips (terminal nodes, cells)
  af.m <- mtr/(depth+0.000001)
  var_assignment <- cbind(
    var_assignment,
    apply(
      var_assignment %>% rownames_to_column(var='variant'), 1, function(x) { 
        node <- as.numeric(x[assign.to])
        mut <- x['variant']
        cells <- tree$tip.label[Descendants(tree, node, 'tips')[[1]]]
        other <- setdiff(colnames(mtr.bi), cells)
        af.1 <- af.m[mut, cells] %>% mean
        p.1 <- sum(mtr.bi[mut, cells]) / length(cells)
        af.0 <- af.m[mut, other] %>% mean
        p.0 <- sum(mtr.bi[mut, other]) / length(other)
        return(c(af1=af.1, af0=af.0, p1=p.1, p0=p.0))
      }
    ) %>% 
    t %>% as.data.frame %>% mutate_all(~ifelse(is.na(.), 0, .))
  )
  
  return(var_assignment)
  
}


##


MakeAllNodes <- function(phy, results, prob.cut=0.3, assign.to='Child'){
  

  NodeVariant.summary <- results %>% filter(p>=prob.cut) %>% count(.data[[assign.to]]) %>% arrange(desc(n)) 
  EdgeTable <- phy$edge %>% as.data.frame
  colnames(EdgeTable) <- c("Parent","Child")
  Allnodes <- merge(EdgeTable, NodeVariant.summary, by=assign.to, all.x=T)    
  Allnodes$CladeSize <- Descendants(phy, Allnodes[[assign.to]]) %>% sapply(.,length)  
  Allnodes$n[is.na(Allnodes$n)]<-0                                         
  
  return(Allnodes)
  
}
 

##


cut_tree <- function(phy, results, N=1, prob.cut=0.3, MinCell=5, Dumpcut=100) {
  
  Allnodes<-MakeAllNodes(phy, results, prob.cut=prob.cut, assign.to='Child') # Get all nodes (i.e., children of a branch), their n of confidently assigned variants and clade size (i.e., their descendants)
  
  # Filter nodes whose list of ancestors assigned-variants >= N, retaining only the up-most
  SumV.df_CloneNode <- data.frame(
    Node=phy$edge[,2],
    n_muts_ancestors=Ancestors(phy, phy$edge[,2], type="all") %>% sapply(., function(x){Allnodes %>% filter(Parent %in% x) %>% .$n %>% sum})
  ) %>% subset(n_muts_ancestors>=N) 
  SumV.df_CloneNode.filtered <- SumV.df_CloneNode[Ancestors(phy, SumV.df_CloneNode$Node) %>% sapply(.,function(x){!any(x %in% SumV.df_CloneNode$Node)}),]
  SumV.df_CloneNode.filtered <- data.frame(SumV.df_CloneNode.filtered, CladeSize=sapply(Descendants(phy,SumV.df_CloneNode.filtered$Node,type="tips"),length)) %>% .[order(.$CladeSize),]
  FinalCloneNodes <-SumV.df_CloneNode.filtered[,c(1,3)] %>% .[order(.$CladeSize),]  
  
  # Small clones
  Smallclones<-c()
  while(FinalCloneNodes$CladeSize[1]<MinCell){              # Sorted in ascending order. [1] is the smallest clade
    SibNode<-Siblings(phy, FinalCloneNodes$Node[1])         # Get the sibling node of the smallest clade
    if(length(Descendants(phy,SibNode)[[1]])<100){          # If the smallest clade sibling has also a very small clade
      Parent<-data.frame(
        Node=Ancestors(phy,FinalCloneNodes$Node[1],type="parent"),
        CladeSize=sapply(Descendants(phy,Ancestors(phy,FinalCloneNodes$Node[1],type="parent")),length)
      )
      FinalCloneNodes<-rbind(FinalCloneNodes[-1,],Parent)   # Replace the smallest clade node with its parent node
      FinalCloneNodes<-FinalCloneNodes[order(FinalCloneNodes$CladeSize),]  # Re-order
    }else{                                                  # If not, remove the clade from FinalCloneNodes and add to Smallclones
      Smallclones<-rbind(Smallclones,FinalCloneNodes[1,])  
      FinalCloneNodes<-FinalCloneNodes[-1,]            
    }
  }
  FinalCloneNodes<-rbind(FinalCloneNodes,Smallclones)
  FinalCloneNodes<-FinalCloneNodes[!duplicated(FinalCloneNodes),]
  # Again, retain only the most ancestral clonal nodes
  FinalCloneNodes<-FinalCloneNodes[Ancestors(phy,FinalCloneNodes$Node) %>% sapply(.,function(x){!any(x %in% FinalCloneNodes$Node)}),]
  FinalCloneNodes<-FinalCloneNodes[order(FinalCloneNodes$CladeSize,decreasing=T),] # Now sort in decreasing order
  
  # Big clones
  Bigclones<-c()
  while(FinalCloneNodes$CladeSize[1]>100){                 # Sorted in ascending order. [1] is the largest clade. 
    node=FinalCloneNodes$Node[1]                           # Largest clade
    # X<-as.list(rep(0,nrow(FinalCloneNodes)))
    # names(X)<-FinalCloneNodes$Node
    CurrentSize<-sapply(Descendants(phy,node),length)
    if(length(Descendants(phy,node,type="children"))>2){  # If the current node is a polytomy, take the two subclades with higher cell number
      x<-Descendants(phy,node,type="children")
      x.size<-sapply(Descendants(phy,x),length)
      a<-Descendants(phy,node,type="children")[order(x.size,decreasing=T)[1]]  
      b<-Descendants(phy,node,type="children")[order(x.size,decreasing=T)[2]]
    }else{
      a<-Descendants(phy,node,type="children")[1]         # Take the two children
      b<-Descendants(phy,node,type="children")[2]
    }
    
    # test.a <- Allnodes %>% filter(Parent==a) %>% .$n %>% sum()>0        # Test if a and b are mutation supported
    # test.b <- Allnodes %>% filter(Parent==b) %>% .$n %>% sum()>0
    
    a.size<-sapply(Descendants(phy,a),length)             # Take their clade size
    b.size<-sapply(Descendants(phy,b),length)
    if(min(a.size,b.size)>MinCell){                       # If both clade sizes is > MinCell
      subreturn_small<-data.frame(Node=c(a,b)[which.min(c(a.size,b.size))],CladeSize=sapply(Descendants(phy,c(a,b)[which.min(c(a.size,b.size))]),length))
      subreturn_big<-data.frame(Node=c(a,b)[which.max(c(a.size,b.size))],CladeSize=sapply(Descendants(phy,c(a,b)[which.max(c(a.size,b.size))]),length))
      FinalCloneNodes<-FinalCloneNodes[-1,]               # Substitute the initial clade with the two subclades
      FinalCloneNodes<-rbind(FinalCloneNodes,subreturn_small)
      FinalCloneNodes<-rbind(FinalCloneNodes,subreturn_big)
    }else{
      node_sub=c(a,b)[which.max(c(a.size,b.size))]        # If one of the two is too small... Take the top one and continue going south...
      a_sub<-Descendants(phy,node_sub,type="children")[1] 
      b_sub<-Descendants(phy,node_sub,type="children")[2]
      a_sub.size<-sapply(Descendants(phy,a_sub),length)   
      b_sub.size<-sapply(Descendants(phy,b_sub),length)   
      if(min(a_sub.size,b_sub.size)>MinCell){            # If now all the two subclades are sized enough, add them to final clones
        subreturn_small.2<-data.frame(Node=c(a_sub,b_sub)[which.min(c(a_sub.size,b_sub.size))],CladeSize=sapply(Descendants(phy,c(a_sub,b_sub)[which.min(c(a_sub.size,b_sub.size))]),length))
        subreturn_big.2<-data.frame(Node=c(a_sub,b_sub)[which.max(c(a_sub.size,b_sub.size))],CladeSize=sapply(Descendants(phy,c(a_sub,b_sub)[which.max(c(a_sub.size,b_sub.size))]),length))
        FinalCloneNodes<-FinalCloneNodes[-1,]
        FinalCloneNodes<-rbind(FinalCloneNodes,subreturn_small.2)
        FinalCloneNodes<-rbind(FinalCloneNodes,subreturn_big.2)
      } else { # If there is still one below MinCell...
        Dump<-min(a_sub.size,b_sub.size)                
        node_sub=c(a_sub,b_sub)[which.max(c(a_sub.size,b_sub.size))]  # Take the biggest subclade
        while(Dump<Dumpcut){
          a_sub<-Descendants(phy,node_sub,type="children")[1]         # Take children
          b_sub<-Descendants(phy,node_sub,type="children")[2]
          a_sub.size<-sapply(Descendants(phy,a_sub),length)
          b_sub.size<-sapply(Descendants(phy,b_sub),length)
          if(min(a_sub.size,b_sub.size)<MinCell & min(a_sub.size,b_sub.size)>0){    # If there is still one very small...
            Dump<-Dump+min(a_sub.size,b_sub.size)                                   # Update Dump
            node_sub=c(a_sub,b_sub)[which.max(c(a_sub.size,b_sub.size))]            # Take the biggest to descend again
          } else {  
            subreturn_small.2<-data.frame(Node=c(a_sub,b_sub)[which.min(c(a_sub.size,b_sub.size))],CladeSize=sapply(Descendants(phy,c(a_sub,b_sub)[which.min(c(a_sub.size,b_sub.size))]),length))
            subreturn_big.2<-data.frame(Node=c(a_sub,b_sub)[which.max(c(a_sub.size,b_sub.size))],CladeSize=sapply(Descendants(phy,c(a_sub,b_sub)[which.max(c(a_sub.size,b_sub.size))]),length))
            FinalCloneNodes<-FinalCloneNodes[-1,]
            FinalCloneNodes<-rbind(FinalCloneNodes,subreturn_small.2)
            FinalCloneNodes<-rbind(FinalCloneNodes,subreturn_big.2)
            break
          }
        }
        if(Dump>=Dumpcut){
          FinalCloneNodes<-FinalCloneNodes[-1,]
          Bigclones<-rbind(Bigclones,data.frame(Node=node,CladeSize=sapply(Descendants(phy,node),length)))
        }
      }
    }  
    FinalCloneNodes<-FinalCloneNodes[order(FinalCloneNodes$CladeSize,decreasing=T),]          # Sort by most abundant
  }
  
  # Final assignment
  FinalCloneNodes<-rbind(FinalCloneNodes,Bigclones)
  FinalCloneNodes<-FinalCloneNodes[order(FinalCloneNodes$CladeSize,decreasing=T),] 
  FinalCloneNodes<-FinalCloneNodes[complete.cases(FinalCloneNodes),]
  FinalCloneNodes$Clone<-1:nrow(FinalCloneNodes)
  SumV.df_CloneNode.filtered<-SumV.df_CloneNode.filtered[order(SumV.df_CloneNode.filtered$CladeSize,decreasing=T),]
  SumV.df_CloneNode.filtered$Clone<-1:nrow(SumV.df_CloneNode.filtered)
  cell_assignment<-FinalCloneNodes %>% 
    apply(., 1, function(x){
      data.frame(
        Cell=phy$tip.label[unlist(Descendants(phy,as.numeric(x[1])))],
        Clade_merge=x[1],
        Clone_merge=x[3],
        row.names=NULL)
  }) %>% do.call(rbind,.)
  row.names(cell_assignment) <- NULL

  # Prep results
  res <- list(nodes_df=Allnodes,supp_nodes=SumV.df_CloneNode.filtered,cell_assignment=cell_assignment)
  
  return(res)
  
}


##