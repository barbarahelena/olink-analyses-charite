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

# Data
df <- rio::import("data/proteins_ttest_diff.RDS")
df$log2mean <- abs_log(df$mean)
df_sum <- df %>% group_by(protein) %>%
    summarise(log2fold = log2mean[which(Treatment == "Pneumonia + CI")] - log2mean[which(Treatment == "Control + CI")],
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
    log2fold > 0 ~ paste0("lower in CI+S.pn."),
    log2fold < 0 ~ paste0("higher in CI+S.pn.")
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
ggplot(df_sum, aes(x = log2fold, y = -log10(pval), color = group, label = delabel)) +
    theme_Publication() +
    theme(axis.title = element_text(size = rel(0.8))) +
    geom_hline(aes(yintercept = -log10(0.05)), color = "darkgrey", linetype = "dashed") +
    geom_vline(aes(xintercept = -2), color = "darkgrey", linetype = "dashed") +
    geom_vline(aes(xintercept = 2), color = "darkgrey", linetype = "dashed") +
    #geom_rect(aes(ymin = 0, ymax = 5, xmin = -2, xmax = 2), fill = "lightgrey", color = "lightgrey") +
    geom_point(alpha = 0.6) +
    geom_text_repel(size = 3, seed = 24, box.padding = 0.5, min.segment.length = 0,
                    point.padding = 0, color = "black", fontface = "bold", force_pull = 0,
                    nudge_x = 0.05, nudge_y = 0.5, segment.color = "grey70",
                    force = 1.5, max.overlaps = 10) +
    scale_color_manual(values = c(ggsci::pal_lancet()(2)), guide = "none") +
    scale_x_continuous(limits = c(-3.5, 3.5), breaks = c(-5:5)) +
    labs(x = "Log2 fold change (Pneumonia+CI - Control+CI)",
         y = "-log10(p-value)")

ggsave("results/pdf/230428_volcanoplot_pval.pdf", width = 5, height = 5, device = "pdf")
ggsave("results/svg/230428_volcanoplot_pval.svg", width = 5, height = 5, device = "svg")
ggsave("results/png/230428_volcanoplot_pval.png", width = 5, height = 5, device = "png")

set.seed(1234)
ggplot(df_sum, aes(x = log2fold, y = -log10(pval), color = group, label = delabel2)) +
    theme_Publication() +
    theme(axis.title = element_text(size = rel(0.8))) +
    geom_hline(aes(yintercept = -log10(0.05)), color = "darkgrey", linetype = "dashed") +
    geom_vline(aes(xintercept = -2), color = "darkgrey", linetype = "dashed") +
    geom_vline(aes(xintercept = 2), color = "darkgrey", linetype = "dashed") +
    #geom_rect(aes(ymin = 0, ymax = 5, xmin = -2, xmax = 2), fill = "lightgrey", color = "lightgrey") +
    geom_point(alpha = 0.6) +
    geom_text_repel(size = 3, seed = 24, box.padding = 0.5, min.segment.length = 0,
                    point.padding = 0, color = "black", fontface = "bold", force_pull = 0,
                    nudge_x = 0.05, nudge_y = 0.5, segment.color = "grey70",
                    force = 1.5, max.overlaps = 10) +
    scale_color_manual(values = c(ggsci::pal_lancet()(2)), guide = "none") +
    scale_x_continuous(limits = c(-3.5, 3.5), breaks = c(-5:5)) +
    labs(x = "Log2 fold change (Pneumonia+CI - Control+CI)",
         y = "-log10(p-value)")

ggsave("results/pdf/230428_volcanoplot_gluc.pdf", width = 5, height = 5, device = "pdf")
ggsave("results/svg/230428_volcanoplot_gluc.svg", width = 5, height = 5, device = "svg")
ggsave("results/png/230428_volcanoplot_gluc.png", width = 5, height = 5, device = "png")

set.seed(1234)
ggplot(df_sum, aes(x = log2fold, y = -log10(qval), color = group, label = delabel2)) +
    theme_Publication() +
    theme(axis.title = element_text(size = rel(0.8))) +
    geom_hline(aes(yintercept = -log10(0.05)), color = "darkgrey", linetype = "dashed") +
    geom_vline(aes(xintercept = -2), color = "darkgrey", linetype = "dashed") +
    geom_vline(aes(xintercept = 2), color = "darkgrey", linetype = "dashed") +
    #geom_rect(aes(ymin = 0, ymax = 5, xmin = -2, xmax = 2), fill = "lightgrey", color = "lightgrey") +
    geom_point(alpha = 0.6) +
    geom_text_repel(size = 3, seed = 24, box.padding = 0.5, min.segment.length = 0,
                    point.padding = 0, color = "black", fontface = "bold", force_pull = 0,
                    nudge_x = 0.05, nudge_y = 0.5, segment.color = "grey70",
                    force = 1.5, max.overlaps = 10) +
    scale_color_manual(values = c(ggsci::pal_lancet()(2)), guide = "none") +
    scale_x_continuous(limits = c(-3.5, 3.5), breaks = c(-5:5)) +
    labs(x = "Log2 fold change (Pneumonia+CI - Control+CI)",
         y = "-log10(q-value)")

ggsave("results/pdf/230428_volcanoplot_qval.pdf", width = 5, height = 5, device = "pdf")
ggsave("results/svg/230428_volcanoplot_qval.svg", width = 5, height = 5, device = "svg")
ggsave("results/png/230428_volcanoplot_qval.png", width = 5, height = 5, device = "png")

