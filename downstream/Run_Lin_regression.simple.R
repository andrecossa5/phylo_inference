
### Run Lin, Poi or Nb regression ###

# Notebook with trajectory analysis at: https://github.com/sankaranlab/redeem_reproducibility/blob/master/Note-6%20HSC%20clonal%20output%20analysis.ipynb
# Script to run Regression at: https://github.com/sankaranlab/redeem_reproducibility/blob/master/API/Run_Lin_regression.R

library(dplyr)
library(tibble)
library(openxlsx)
library(ggrepel)
library(rlang)

SEED <- 4321
set.seed(SEED)

path_utils <- fs::path("/Users/ieo6983/Desktop/breast_alberto/BC_chemo_reproducibility/downstream/MDA/trajectory/utils_trajectory.simple.R")

path_input_table <- fs::path("/Users/ieo6983/Desktop/breast_alberto/data/MDA/agg_for_poisson.RDS")
path_results_data <- fs::path("/Users/ieo6983/Desktop/breast_alberto/results/MDA/data/")
path_results_plots <- fs::path("/Users/ieo6983/Desktop/breast_alberto/results/MDA/plots/")

for(dir in c(path_results_data, path_results_plots)){
  if(!dir.exists(dir)){
    dir.create(dir, recursive = T)
  }
}

source(path_utils)


##


# Define parameters 
mode <- "poi"   # poi (poisson), nb (negative binomial), lm (linear)
core <- as.numeric(8) 
name <- "reg_out" 
clone_column <- "GBC" # define which table's column contains clone info 
regression_factor <- "met_potential"  


##


## Load and filter data

print("Read in")
print(path_input_table)
LinOut.df<-readRDS(path_input_table)

# Filter out genes with 0-expression in more than 10% of samples
LinOut.df.filt <- filter_genes(LinOut.df, clones = clone_column)


##


## Run regression

if(mode=="lm"){
  print("run linear model")
  LinOut.result <- Run_Lin_regression(LinOut.df.filt, n.cores = core)
}else if(mode=="poi"){
  print("run poisson model")
  LinOut.result <- Run_Lin_regression_poi(LinOut.df.filt, regress_factor = regression_factor, 
                                          n.cores = core, qval = T, tot_UMIs = T)
}else if(mode=="nb"){
  # W: nb does not support vectors of zeros as response variable (genes with 0 expression across all clones)
  print("run negative-binomial model")
  LinOut.result <- Run_Lin_regression_nb(LinOut.df.filt, n.cores = core, qval = T, tot_UMIs = T)
}

# Convert results to data.frame
if("qs" %in% names(LinOut.result)){
  LinOut.result.df <- convert_results_to_df(LinOut.result)    
}else{
  LinOut.result.df <- merge(LinOut.result$slopes, LinOut.result$ps, by = "row.names")
  colnames(LinOut.result.df) <- c("genes", "slopes", "ps")  
}


##


# Plot results (volcano) & save output

#par(mfrow = c(1, 1))
#qline = 0.05; pline = 0.001
#p1 <- PlotLinRegress_Vocano(LinOut.result.df, slot = "ps", pline = 0.01)

# Save output 
#saveRDS(LinOut.result, fs::path(path_results, paste0(name, ".", mode, ".RDS")))
#LinOut.result.df %>% arrange(., qs) %>% write.xlsx(., fs::path(path_results, paste0(name, ".", mode, ".xlsx")), rowNames=T)


##


#' @FIXME:
#'Run_Lin_regression() contains a step to compute tot_UMIs that depends on how the input table is formatted.
#'Change to make generalizable, or assert that input table as the required format. 




