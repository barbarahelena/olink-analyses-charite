# Volcano plot Olink
# Barbara Verhaar, b.j.verhaar@amsterdamumca.nl

# Libraries
library(tidyverse)
library(ggrepel)
library(ggpubr)
library(ggsci)

abs_log <- function(x){
    x[x==0] <- 1
    si <- sign(x)
    si * log2(si*x)
}

# Data with ranks
df2 <- readRDS("data/proteins_expression_ranks.RDS")
df_sum <- df2 %>% group_by(protein) %>% 
    summarise(rankdiff = mean[which(Treatment == "Pneumonia + CI")] - mean[which(Treatment == "Control + CI")],
              qval = q.value,
              pval = p.value,
              sig = sig,
              # mean = mean,
              # treatment = Treatment
    ) %>% filter(duplicated(protein)) 
head(df_sum)
df_sum %>% filter(protein == "Mia")
df_sum %>% filter(protein == "Gcg")

## Volcano plot
df_sum <- df_sum %>% 
    mutate(group = case_when(
    rankdiff > 0 ~ paste0("lower in CI+S.pn."),
    rankdiff < 0 ~ paste0("higher in CI+S.pn.")
), 
sigdir = case_when(
    sig != "" & group == "higher in CI+S.pn." ~ paste0("up"),
    sig != "" & group == "lower in CI+S.pn." ~ paste0("down"),
    sig == "" ~ paste0("no")
), 
    sigdir = as.factor(sigdir),
    delabel = case_when(
       sig != "" ~ paste0(protein),
       sig == "" ~ paste0("")
    ),
    delabel2 = case_when(
        protein == "Gcg" ~ paste0("Glucagon"),
        protein != "Gcg" ~ paste0("")
    )
    ) %>% 
    ungroup(.) 
df_sum
#df_sum$rank <- c(c(-53:-1), c(39:1))

set.seed(1234)


ggplot(df_sum, aes(x = rankdiff, y = -log10(pval), color = group, label = delabel)) +
    theme_Publication() +
    theme(axis.title = element_text(size = rel(0.8))) +
    geom_hline(aes(yintercept = -log10(0.05)), color = "darkgrey", linetype = "dashed") +
    geom_segment( aes(x=rankdiff, xend=rankdiff, y=0, yend=-log10(pval)), color = "grey") +
    #geom_vline(aes(xintercept = -2), color = "darkgrey", linetype = "dashed") +
    #geom_vline(aes(xintercept = 2), color = "darkgrey", linetype = "dashed") +
    #geom_rect(aes(ymin = 0, ymax = 5, xmin = -2, xmax = 2), fill = "lightgrey", color = "lightgrey") +
    geom_point(alpha = 0.6) +
    geom_text_repel(size = 3, seed = 24, box.padding = 0.5, min.segment.length = 0,
                    point.padding = 0, color = "black", fontface = "bold", force_pull = 0,
                    nudge_x = 0.05, nudge_y = 0.5, segment.color = "grey70", 
                    force = 1.5, max.overlaps = 10) +
    scale_color_manual(values = c(ggsci::pal_lancet()(2)), guide = "none") +
    #scale_x_continuous(limits = c(55, 40)) +
    labs(x = "Mean rank difference (pneumonia - control)",
         y = "-log10(p-value)") 

ggsave("results/pdf/230525_volcanoplot_rank.pdf", width = 5, height = 5, device = "pdf")    
ggsave("results/svg/230525_volcanoplot_rank.svg", width = 5, height = 5, device = "svg")
