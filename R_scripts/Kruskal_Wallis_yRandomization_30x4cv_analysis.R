library(rstatix)
library(tidyverse)
library(xtable)

# **************************************************  
# Do Kruskal-Wallis and Dunn test 
# per VS_method between chi values
# **************************************************  

# KRUSKAL-WALLIS

# The following function allows to obtain the `y_scrambling`
# csv file path from each protein
get_filname <- function(prot_name){
  path <- paste0('./', prot_name, 
                 '/5_Machine_Learning/y_scrambling_', 
                 prot_name,'.csv')
  return(path)
}

# Get the path files for all proteins
cdk2_yRand_file  <- get_filname('cdk2')
fxa_yRand_file   <- get_filname('fxa')
egfr_yRand_file  <- get_filname('egfr')
hsp90_yRand_file <- get_filname('hsp90')

# Read the files and combine them into a 
# vector of dataframes
df_cdk2  = read.csv(cdk2_yRand_file,  header = TRUE)
df_fxa   = read.csv(fxa_yRand_file,   header = TRUE)
df_egfr  = read.csv(egfr_yRand_file,  header = TRUE)
df_hsp90 = read.csv(hsp90_yRand_file, header = TRUE)

# Define the dataframes to iterate over them
# and apply the statistical analysis
dataFrames = list('cdk2' = df_cdk2,
                  'fxa'  = df_fxa,
                  'egfr' = df_egfr,
                  'hsp90'= df_hsp90)

# Define the metrics and the SBVS methods to evaluate
metrics    = unique(df_cdk2$metric)
vs_methods = unique(df_cdk2$vs_method)

# Compute the Kruskal-Walis analysis 
# per chi value, per sbvs method, per metric
tab_res <- c()
for (prot in names(dataFrames)) {
  for (m in metrics) {
    for (vs_m in vs_methods) {
      df <- dataFrames[[prot]]
      kw <- df %>%
        filter(m == metric) %>%
        filter(vs_m == vs_method) %>%
        kruskal_test(score ~ chi) %>%
        select(- c(.y., method)) %>%
        mutate(signif = case_when(
                  p < 0.001 ~ '***',
                  p < 0.01 ~ '**',
                  p < 0.05 ~ '*',
                  TRUE ~ 'NS'
                )) %>%
        cbind(Protein  = toupper(prot),
              Metric    = toupper(m), 
              VS_method = toupper(vs_m), .)
      
      tab_res <- rbind(tab_res, kw)
    }
  }
}

xtable(tab_res)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
# Dunn's TEST
vs_methods = unique(df_cdk2$vs_method)
tab_res_dunn <- c()
for (prot in names(dataFrames)) {
  for (m in metrics) {
    for (vs_m in vs_methods) {
      df <- dataFrames[[prot]]
      dunn <- df %>%
        filter(m == metric) %>%
        filter(vs_m == vs_method) %>%
        dunn_test(score ~ chi) %>%
        select(- c(.y., n2)) %>%
        filter(p.adj.signif == 'ns') %>%
        mutate(VS_method = rep(vs_m, nrow(.))) %>%
        mutate(Metric = rep(toupper(m), nrow(.)))  %>%
        mutate(Protein = rep(toupper(prot), nrow(.)))
      
      tab_res_dunn <- rbind(tab_res_dunn, dunn)
    }
  }
}
# format the table
tab_res_dunn <- tab_res_dunn %>% relocate(Protein, Metric, VS_method)
# Filter out the Dummy classifier results
tab_res_dunn <- tab_res_dunn %>%
  filter(VS_method != 'DClf')
xtable(tab_res_dunn)