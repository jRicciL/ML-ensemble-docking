library(rstatix)
library(tidyverse)
library(xtable)

# Analysis
cdk2_yRand_file = './yRand_30x4cv_cdk2.csv'
fxa_yRand_file = './yRand_30x4cv_fxa.csv'

df_cdk2 = read.csv(cdk2_yRand_file, header = TRUE)
df_fxa  = read.csv(fxa_yRand_file, header = TRUE)
dataFrames = list('cdk2' = df_cdk2, 'fxa' = df_fxa)

# Do Kruskal-Wallis and Dunn test per VS_method between chi values
metrics = unique(df_cdk2$X)
vs_methods = unique(df_cdk2$X.1)

tab_res <- c()
for (prot in names(dataFrames)) {
  for (metric in metrics) {
    for (vs_method in vs_methods) {
      df <- dataFrames[[prot]]
      kw <- df %>%
        filter(X == metric) %>%
        filter(X.1 == vs_method) %>%
        kruskal_test(score ~ Method) %>%
        select(- c(.y., method)) %>%
        mutate(signif = case_when(
                  p < 0.001 ~ '***',
                  p < 0.01 ~ '**',
                  p < 0.05 ~ '*',
                  TRUE ~ 'NS'
                )) %>%
        cbind(Protein = toupper(prot),
              Metric = toupper(metric), 
              VS_method = toupper(vs_method), .)
      
      tab_res <- rbind(tab_res, kw)
    }
  }
}

xtable(tab_res)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
# Dunn's TEST
vs_methods = c("ml_lr",   "ml_xgb",  "cs_MEAN", "cs_GEO",  "cs_MIN" )
tab_res_dunn <- c()
for (prot in names(dataFrames)) {
  for (metric in metrics) {
    for (vs_method in vs_methods) {
      df <- dataFrames[[prot]]
      dunn <- df %>%
        filter(X == metric) %>%
        filter(X.1 == vs_method) %>%
        dunn_test(score ~ Method) %>%
        select(- c(.y., n2)) %>%
        filter(p.adj.signif == 'ns') %>%
        mutate(VS_method = rep(vs_method, nrow(.))) %>%
        mutate(Metric = rep(toupper(metric), nrow(.)))  %>%
        mutate(Protein = rep(toupper(prot), nrow(.)))
      # %>%
      #   cbind(Protein = toupper(prot),
      #         Metric = toupper(metric), 
      #         VS_method = toupper(vs_method), .)
      
      tab_res_dunn <- rbind(tab_res_dunn, dunn)
    }
  }
}
tab_res_dunn <- tab_res_dunn %>% relocate(Protein, Metric, VS_method)

xtable(tab_res_dunn)
