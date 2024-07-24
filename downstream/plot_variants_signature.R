
library(tidyverse)

# Using data: /hpcnfs/scratch/PGP/acossa/AML_MITO/phylo/phylo_old/lymphoid/sAML1/set5/
path_variants <- fs::path("/Users/ieo6983/Desktop/phylo_inference/data/phylo/phylo_old/lymphoid/sAML1/set5/input_folder/variants.csv")


##


#' Function to plot bulk level mutation signatures
#'
#' This function allows you to plot the mito mutation signatures
#' @param cell_variants a vector of variants formated as c('93_A_G''103_G_A''146_T_C'
#' @return p from ggplot2
#' @examples
#' MutationProfile.bulk(DN1CD34_1.Variants.feature.lst[[name]]$Variants
MutationProfile.bulk<-function(cell_variants){
  library(redeemR)
  library(ggplot2)
  
  # Annotate with called variants
  data(ref_all_long)
  called_variants <- strsplit(cell_variants,"_") %>% sapply(.,function(x){paste(x[1],x[2],">",x[3],sep="")})
  ref_all_long$called <- ref_all_long$variant %in% called_variants
  # Compute changes in expected/observed
  total <- dim(ref_all_long)[1]
  total_called <- sum(ref_all_long$called)
  prop_df <- ref_all_long %>% group_by(three_plot, group_change, strand) %>%
    dplyr::summarize(observed_prop_called = sum(called)/total_called, expected_prop = n()/total, n = n()) %>%
    mutate(fc_called = observed_prop_called/expected_prop)
  prop_df$change_plot <- paste0(prop_df$group_change, "_", prop_df$three_plot)
  # Visualize
  p1 <- ggplot(prop_df, aes(x = change_plot, fill = strand, y = fc_called)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme(axis.title.x=element_blank(),
          axis.text.x =element_blank())+
    scale_fill_manual(values= c("firebrick", "dodgerblue3")) +theme_bw()+
    theme(legend.position = "bottom",axis.text.x  = element_blank(),axis.ticks.x = element_blank()) +
    scale_y_continuous(expand = c(0,0)) +
    geom_hline(yintercept = 1, linetype =2, color = "black") +
    labs(x = "Change in nucleotide", y = "Substitution Rate")+
    facet_grid(.~group_change,scales = "free",space="free")
  return(p1)
}


##


# Load and format variants
variants <- read_csv(path_variants, col_names = "variant")
variants <- str_replace(variants$variant, ">", "_")

# Plot the mutation signatures
options(repr.plot.width=16, repr.plot.height=6,repr.plot.res=80)
p<-MutationProfile.bulk(variants)+ggtitle("Title")+theme(title =element_text(size=20))
print(p)


##

