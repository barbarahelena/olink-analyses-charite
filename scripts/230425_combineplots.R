## Arrange plots on one grid
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Libraries
library(ComplexHeatmap)
library(tidyverse)
library(ggsci)
library(circlize)
library(cowplot)

theme_Publication <- function(base_size=12, base_family="sans") {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(1.0), hjust = 0.5),
                text = element_text(family = 'Helvetica'),
                panel.background = element_rect(colour = NA),
                plot.background = element_rect(colour = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(1)),
                axis.title.y = element_text(angle=90,vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text = element_text(), 
                axis.line = element_line(colour="black"),
                axis.ticks = element_line(),
                panel.grid.major = element_line(colour="#f0f0f0"),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                legend.position = "bottom",
                # legend.direction = "horizontal",
                legend.key.size= unit(0.2, "cm"),
                legend.spacing  = unit(0, "cm"),
                # legend.title = element_text(face="italic"),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold")
        ))
    
} 

## Data
df3 <- readRDS("data/proteins_expression_data.RDS") # all protein expression data, unscaled
df3_scale <- readRDS("data/proteins_expression_scaled_data.RDS") # all protein expression data (scaled to mean 0 and sd 1)
res <- readRDS("data/proteins_ttest_diff.RDS") %>% select(protein, p.value, q.value, sig) %>% filter(!duplicated(protein))
qval_sig <- "Gcg"
pval_sig <- unique(df5_sel$protein)

## Heatmap
col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("dodgerblue4", "white", "firebrick")) # color mapping for blocks heatmap
pvalue_col_fun = colorRamp2(c(0, 1, 3), c("lightgrey", "white", "firebrick1")) # color mapping for -log10(pvalue)
font_treatment <- 12 # bigger fonts because heatmap is smaller
font_proteins <- 10
width_pvalue <- 0.5
width_blocks <- 12
height_blocks <- 0.5

# matrix of all proteins (3rd column until end) for heatmap (only works with matrices)
df3_matrix <- t(as.matrix(df3_scale[,3:ncol(df3_scale)])) 
df3_matrix3 <- df3_matrix[which(rownames(df3_matrix) %in% pval_sig),] # select only sig proteins
res_p <- res %>% filter(protein %in% pval_sig)
res_p$sig <- case_when(res_p$q.value >= 0.05 ~ paste0(""),
                       res_p$q.value < 0.05 ~ paste0("*"))

(pl4 <- Heatmap(df3_matrix3, name = "Z-score", 
                col = col_fun, 
                rect_gp = gpar(col = "white", lwd = 2),
                column_title_gp = grid::gpar(fontsize = font_treatment, fontface = "bold"),
                column_split = c(rep("CI + PBS", 5), rep("CI + S.pneu",5)),  
                column_order = colnames(df3_matrix),
                show_column_dend = FALSE,
                top_annotation = HeatmapAnnotation(Treatment = anno_block(gp = gpar(fill = c("darkgrey", "dodgerblue"))), 
                                                   height = unit(0.2, "cm")),
                left_annotation = rowAnnotation(qvalue = anno_simple(-log10(res_p$q.value), 
                                                                     which = 'row',
                                                                     col = pvalue_col_fun, 
                                                                     pch = res_p$sig,
                                                                     pt_size = unit(1, "snpc")*1.0),
                                                annotation_name_side = "bottom",
                                                width = unit(width_pvalue, "cm"),
                                                annotation_name_gp = grid::gpar(fontsize = 8)),
                width = unit(width_blocks, "cm"),
                height = unit(height_blocks*nrow(df3_matrix3), "cm"),
                row_names_side = "left", 
                column_gap = unit(2, "mm"),
                row_dend_side = "right",
                row_names_gp = grid::gpar(fontsize = font_proteins)))

pl4_anno <- draw(pl4, annotation_legend_list = list(lgd_qvalue, lgd_sig))

heatmap_grob <- grid.grabExpr(draw(pl4_anno))


## Boxplots
res_box <- list()
for(a in 1:3){
    protein <- pval_sig[a]
    proteinname <- paste0(pval_sig[a])
    print(proteinname)
    df_protein <- df3_scale %>% select(Treatment, all_of(protein))
    df_protein$protein_y <- df_protein[,2]
    comp <- list(c("CI+PBS", "CI+S.pn."))
    (pl <- ggplot(data = df_protein, aes(x = Treatment, y = protein_y)) + 
            ggpubr::stat_compare_means(method = "t.test", label = "p.signif", comparisons = comp,
                                       hide.ns = TRUE, bracket.size = 0.5, size = 5) +
            geom_boxplot(aes(fill = Treatment), color = "black", outlier.shape = NA, 
                         width = 0.5, alpha = 0.9) +
            geom_jitter(color = "black", height = 0, width = 0.1) +
            scale_fill_manual(guide = "none", values = c("darkgrey","dodgerblue")) +
            ylim(NA, max(df_protein$protein_y)*1.3) +
            labs(title=proteinname, y="Protein expression (z-score)") +
            theme_Publication())
    res_box[[a]] <- pl
}

pl_boxsig <- gridExtra::grid.arrange(grobs=list(res_box[[1]], res_box[[2]], res_box[[3]]), ncol=3)

## PCA
res.pca <- prcomp(df3[,4:ncol(df3)], scale = TRUE)
evar <- get_eigenvalue(res.pca)
evar_pc12 <- c(evar$variance.percent[1], evar$variance.percent[2])
evar_pc12 <- format(round(evar_pc12, 2), 2)
x_axis <- str_c("PC1 (", evar_pc12[1], "%)")
y_axis <- str_c("PC2 (", evar_pc12[2], "%)")
df.pca <- as.data.frame(res.pca$x[,c("PC1", "PC2")])
df.pca$sampleID <- df3$sampleID
df.tot <- left_join(df.pca, df3, by = "sampleID")

(plot_pca <- df.tot %>% 
        ggplot(aes(PC1, PC2)) +
        geom_point(aes(color = Treatment), size = 2) +
        xlab(x_axis) +
        ylab(y_axis) +
        theme_Publication() +
        scale_color_manual(values = c("grey48", "dodgerblue")) +
        guides(fill = guide_legend(override.aes = list(shape = 21, size = 2))) +
        stat_ellipse(aes(color = Treatment), type = "t", level = 0.9) +
        ggtitle("Principal component analysis"))


plot_grid(heatmap_grob, plot_pca, pl_boxsig, ncol = 2, nrow = 2,
          rel_widths = c(2, 1), labels = "AUTO")
ggsave("results/pdf/plot_combi.pdf", width = 12, height = 8)
ggsave("results/svg/plot_combi.svg", width = 12, height = 8)
