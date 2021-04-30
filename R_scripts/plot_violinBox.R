library(extrafont)
library(ggthemes)
library(ggsci)
library(plyr)
library(tidyr)
library(dplyr)
library(gghalves)
library(scales)
library(matrixStats)
library(tidyverse)

# initialize some variables
theme_set(theme_gray(base_size = 10))
full_names <- c('LR', 'GBT', 'DClf', 'csAVG', 'csGEO', 'csMIN')
cbbPalette <- c( '#785EF0', '#3F93D2', '#44AA99', '#FE6100', '#DC267F', '#FFB000')
names(cbbPalette) <- full_names
lb <- function(x) { median(x) - 1.57*IQR(x)/sqrt(length(x)) }
ub <- function(x) { median(x) + 1.57*IQR(x)/sqrt(length(x)) }

# Plot
plot_swarm_box <- function(df, cbbPalette = cbbPalette, decreasing_order = TRUE, y_label='AUC-ROC', 
                           y_min=0.4, y_max=1, dot_size=8, bin_width=0.001, base_h_line=0.5) {
    
    names_order <- names(sort(apply(df, 2, FUN=median), decreasing = decreasing_order))
    df <- df[, names_order]
    
    df_melted <- df %>%
        mutate(rep = factor(1:nrow(.))) %>%
        pivot_longer(cols=c(everything(), -rep), names_to='method', values_to='score')
    
    df_melted$method <- factor(df_melted$method, levels = names_order)
   

    sumld<- ddply(df_melted, ~method, summarise, 
                  mean = mean(score), 
                  median = median(score), 
                  lower = lb(score), 
                  upper = ub(score))

    ggplot(data = df_melted, 
           mapping = aes(x = method, 
                         y = score, 
                         fill = method)) + 
      geom_hline(yintercept= base_h_line, linetype="dashed", color="#444444") +
      geom_violin(width=01, lwd=0.2, color='black') +
      theme(text=element_text(family="Trebuchet MS")) + 
      stat_summary(fun.data = mean_sdl, 
                   fun.args = list(mult = 1), 
                   geom = "pointrange", 
                   position = position_nudge(0.5)) +
      geom_errorbar(data = sumld, aes(ymin = lower, ymax = upper, y = median), 
                    position = position_nudge(x = -0.17), width = 0) + 
      geom_point(data = sumld, aes(x = method, y = median), colour='#000000',
                    position = position_nudge(x = -0.17), size = 0.8, stroke=0.5) +
      geom_point(data = sumld, aes(x = method, y = median), colour='white',
                    position = position_nudge(x = -0.17), size = 0.7, stroke=0.1) +
        theme(legend.position = "none", panel.border = element_rect(colour = "black", fill=NA, size=1),
              panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
              panel.grid.major.y = element_line(size = 0.1, linetype = 'solid', colour = "grey"), 
              panel.grid.minor.y = element_line(size = 0.05, linetype = 'solid', colour = "darkgrey"),
              panel.grid.major.x = element_blank()
             ) + 
          labs(x = "Methods (ML/CS)", 
               y = y_label) +
      scale_y_continuous(breaks = seq(y_min, y_max, 0.1), limits = c(y_min, y_max)) +  
      scale_fill_manual(values=cbbPalette)
}

plot_violin <- function(df, cbbPalette = cbbPalette, decreasing_order = TRUE, y_label='AUC-ROC', scale='area',
                           y_min=0.4, y_max=1, dot_size=8, bin_width=0.001, base_h_line=0.5,
                           violin_width=1.1) {
    
    names_order <- names(sort(apply(df, 2, FUN=median), decreasing = decreasing_order))
    df <- df[, names_order]
    
    df_melted <- df %>%
        mutate(rep = factor(1:nrow(.))) %>%
        pivot_longer(cols=c(everything(), -rep), names_to='method', values_to='score')
    
    df_melted$method <- factor(df_melted$method, levels = names_order)
   

    sumld<- ddply(df_melted, ~method, summarise, 
                  mean = mean(score), 
                  median = median(score), 
                  lower = lb(score), 
                  upper = ub(score))

    ggplot(data = df_melted, 
           mapping = aes(x = method, 
                         y = score, 
                         fill = method)) + 
      geom_hline(yintercept= base_h_line, linetype="dotted", color="#444444", lwd=0.8) +
      geom_violin(width=violin_width, lwd=0.3, color='black', alpha=0.6, trim=TRUE, scale=scale) +
      geom_boxplot(aes(outlier.color=method), notch=TRUE, width=0.3, lwd=0.2, color='black',
                  outlier.size=0.5, outier.shape=21) +
      # theme(text=element_text(family="Trebuchet MS")) + 
      stat_summary(fun.data = mean_sdl, 
                   fun.args = list(mult = 1), 
                   geom = "pointrange", 
                   position = position_nudge(0.5)) +
      geom_errorbar(data = sumld, aes(ymin = lower, ymax = upper, y = median), 
                    width = 0) + 
      geom_point(data = sumld, aes(x = method, y = median), colour='#000000',
                    size = 0.6, stroke=0.5) +
      geom_point(data = sumld, aes(x = method, y = median), colour='white',
                    size = 0.6, stroke=0.1) +
        theme(legend.position = "none", panel.border = element_rect(colour = "black", fill=NA, size=0.6),
              panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.6, linetype = "solid"),
              panel.grid.major.y = element_line(size = 0.2, linetype = 'solid', colour = "grey"), 
              panel.grid.minor.y = element_line(size = 0.2, linetype = 'solid', colour = "lightgrey"),
              axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
              panel.grid.major.x = element_blank(),
              plot.title = element_text(hjust = 0.5, size=11)
             ) + 
          labs(x = "Methods (ML/CS)", 
               y = y_label) +
      scale_y_continuous(breaks = seq(y_min, y_max, 0.1), limits = c(y_min, y_max)) + 
      # scale_x_discrete(guide = guide_axis(n.dodge = 2)) + 
    #   scale_fill_brewer(palette = "") +
      scale_fill_manual(values=cbbPalette)
}


add_ref_values <- function(text, value, color='#888888', y_add=0.02, x=0.5, size=2.3) {
        list(
            geom_hline(yintercept= value, linetype="dashed", color=color, size=0.3),
            annotate('text', x=x, y= value + y_add, label=text, color=color,  size=size, hjust=0)
        )
    }


plot_lines <- function(df, cbbPalette=cbbPalette, y_label='AUC-ROC', y_min=0.4, y_max=1, switch_x=TRUE, 
                       line_size=1, point_size=2.2, error_dodge=0.05, error_width=1.5, error_size=1, 
                       legend.position='none', shape=21, add_ribbon=FALSE, ribbon_alpha=0.2, title.size=11,
                       ticks.text.size=8, ticks.text.angle=0, 
                       base_h_line=0.5, x_label="Percentage of shuffled labels (%)", title='',
                       include_color_scale = TRUE) {

    ggplot(data = df, 
           mapping = aes(x = index, 
                         y = mean, 
                         color = selection)) + 
        geom_hline(yintercept= base_h_line, 
                   linetype="dashed", color="#333333") +
        geom_line(size=line_size, alpha=0.7) + 
        # theme(text=element_text(family="Trebuchet MS")) + 
        {if(switch_x)scale_x_reverse()} +
        geom_errorbar(aes(ymin=mean-std, ymax=mean+std), width=error_width, size=error_size, position=position_dodge(error_dodge)) +
        {if(add_ribbon) geom_ribbon(aes(ymin=mean-std, ymax=mean+std, fill=selection), alpha=ribbon_alpha, colour=NA)} +
        # geom_point(color='black', size = point_size + 0.6, stroke = 0.5, shape=shape)+
        geom_point(aes(fill=selection), colour='black', size = point_size, stroke = 0.4, shape=shape)+
        ggtitle(title) +
        scale_y_continuous(breaks = seq(y_min, y_max, 0.1), limits = c(y_min, y_max)) + 
        theme(legend.position = legend.position, 
              legend.title=element_text(size=8),
              legend.text=element_text(size=7),
              panel.border = element_rect(colour = "black", fill=NA, size=0.6),
              panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.6, linetype = "solid"),
              panel.grid.major.y = element_line(size = 0.2, linetype = 'solid', colour = "grey"), 
              panel.grid.minor.y = element_line(size = 0.2, linetype = 'solid', colour = "lightgrey"),
              axis.title.x = element_text(size=8),
              axis.text.x=element_text(size=ticks.text.size, angle=ticks.text.angle),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              plot.title = element_text(hjust = 0.5, size=title.size)
             ) + 
        labs(x = x_label, y = y_label) +
        {if(include_color_scale) 
        scale_color_manual(values=rev(cbbPalette), name='selection')} +
        {if(include_color_scale) 
        scale_fill_manual(values=rev(cbbPalette), name='selection')}
    }
