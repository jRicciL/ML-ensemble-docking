
library(reshape2)
library(PMCMR)
library(ggplot2)

p_labels = c("p < 0.001", "p < 0.01", 'p < 0.05', 'NS')
names(p_labels) <- p_labels
p_colors = c('#FA3A6A',   "#FD9EB2", "#FBDEE3", '#53CDF6')
names(p_colors) <- p_labels

plot_p_vals_heatmap <- function(df_R) {
    df_melt <- df_R %>%
        mutate(rep = factor(1:nrow(df_R))) %>%
        pivot_longer(cols=c(everything(), -rep), names_to='Method', 
                     values_to='score')

    p_values <- posthoc.friedman.nemenyi.test(formula=score ~ Method | rep, data=df_melt)$p.value
    p_values <- p_values[rev(rownames(p_values)), ]

    m_p_values <- melt(p_values, na.rm = TRUE)
    
    m_p_values$cutoffs <- cut(m_p_values$value, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), right=TRUE, labels=p_labels)

    ggplot(data = m_p_values, aes(x=Var1, y=Var2, fill=cutoffs)) + 
      geom_tile(color='white', size=1) +
       scale_fill_manual(labels = p_labels, 
                         values = p_colors,
                        name='Significance') + 
       geom_text(aes(Var1, Var2, label = round(value, 4)), color = "black", size = 3.5) +
       theme(
        # text=element_text(family="Trebuchet MS"),
             plot.title = element_text(hjust = 0.5, size=13),
             plot.subtitle = element_text(hjust = 0.5, size=11),
             panel.background = element_rect(fill = "white"),
             legend.position = c(0.85, 0.8),
             axis.title.x = element_blank(),
             axis.title.y = element_blank(),
             axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=10, color='black'),
             axis.text.y = element_text(size=10, color='black')
             )
    
    }

get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }

plot_p_vals_heatmap_kruskal <- function(df_R){
    df_R <- get_lower_tri(df_R)
    df_R <- df_R[rev(rownames(df_R)), ]
    m_p_values <- melt(as.matrix(df_R), na.rm = TRUE)
    m_p_values$cutoffs <- cut(m_p_values$value, 
                              breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), 
                              right=TRUE, labels=p_labels)
    m_p_values
    ggplot(data = m_p_values, aes(x=Var1, y=Var2, fill=cutoffs)) + 
          geom_tile(color='white', size=1) +
           scale_fill_manual(labels = p_labels, 
                             values = p_colors,
                            name='Significance') + 
           geom_text(aes(Var1, Var2, label = round(value, 4)), 
                     color = "black", size = 3.5) +
           theme(
            # text=element_text(family="Trebuchet MS"),
                 plot.title = element_text(hjust = 0.5, size=13),
                 plot.subtitle = element_text(hjust = 0.5, size=11),
                 panel.background = element_rect(fill = "white"),
                 legend.position = c(0.85, 0.8),
                 axis.title.x = element_blank(),
                 axis.title.y = element_blank(),
                 axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=10, color='black'),
                 axis.text.y = element_text(size=10, color='black')
                 )
    }