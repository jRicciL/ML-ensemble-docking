library(rstatix)
library(ggpubr)
library(tidyr)
library(tidyverse)
library(dplyr)
library(ez)
library(xtable)

# REPEAT THIS ANALYSIS PER METRIC AND PER PROTEIN
# MEtrics: roc_auc and nef_12_Ra
# Proteins: cv30x4_cdk2.csv and cv30x4_fxa.csv

do_RM_ANOVA <- function(csv_file, metric_colname='roc_auc') {
  df = read.csv(csv_file, header = TRUE)
  
  df = df %>%
    filter(X == metric_colname) %>%
    select(-X.1, -X)
  
  # Melting the data
  df_melt <- df  %>%
    mutate(rep = factor(1:nrow(.))) %>%
    pivot_longer(cols=c(everything(), -rep), 
                 names_to='VS_method', 
                 values_to='score')
  
  # Check for outliers
  df_melt %>%
    group_by(VS_method) %>%
    identify_outliers(score)
  
  df_melt %>%
    group_by(VS_method) %>%
    shapiro_test(score)
  
  shapiro.test(df$csGEO)
  
  anova <- aov(score ~ VS_method, data=df_melt)
  print(summary(anova))
  # Perfrom the repeated measures ANOVA
  res.aov <- anova_test(data = df_melt, 
                        dv = score, 
                        wid = rep, 
                        within = VS_method)
  get_anova_table(res.aov)
  return(res.aov)
}

# Perform the analysis

# CDK2
cdk2_file = '../cdk2/5_Machine_Learning/cv30x4_cdk2.csv'
## ROC
res.aov.cdk2.roc <- do_RM_ANOVA(csv_file = cdk2_file, metric_colname = 'roc_auc')
## NEF
res.aov.cdk2.nef <- do_RM_ANOVA(csv_file = cdk2_file, metric_colname = 'nef_Ra')

# FXa
fxa_file = '../fxa/5_Machine_Learning/cv30x4_fxa.csv'
## ROC
res.aov.fxa.roc <- do_RM_ANOVA(csv_file = fxa_file, metric_colname = 'roc_auc')
## NEF
res.aov.fxa.nef <- do_RM_ANOVA(csv_file = fxa_file, metric_colname = 'nef_Ra')

# EGFR
egfr_file = '../egfr/5_Machine_Learning/cv30x4_egfr.csv'
## ROC
res.aov.egfr.roc <- do_RM_ANOVA(csv_file = egfr_file, metric_colname = 'roc_auc')
## NEF
res.aov.egfr.nef <- do_RM_ANOVA(csv_file = egfr_file, metric_colname = 'nef_Ra')

# HSP90
hsp90_file = '../hsp90/5_Machine_Learning/cv30x4_hsp90.csv'
## ROC
res.aov.hsp90.roc <- do_RM_ANOVA(csv_file = hsp90_file, metric_colname = 'roc_auc')
## NEF
res.aov.hsp90.nef <- do_RM_ANOVA(csv_file = hsp90_file, metric_colname = 'nef_Ra')


results <- list(
  res.aov.cdk2.roc, res.aov.cdk2.nef,
  res.aov.fxa.roc, res.aov.fxa.nef,
  res.aov.egfr.roc, res.aov.egfr.nef,
  res.aov.hsp90.roc, res.aov.hsp90.nef
  )

tab_res <- c()
for (res in results) {
    tab_res <- rbind(tab_res, cbind(res$ANOVA,
          res$`Mauchly's Test for Sphericity`[2:4]))
}

# Print the latex table
xtable(
  tab_res
)

# Citation
citation('rstatix')
citation('tidyverse')




