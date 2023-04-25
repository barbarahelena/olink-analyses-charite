## Heatmaps for Olink assay
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Libraries
library(ComplexHeatmap)
library(tidyverse)
library(ggsci)
library(circlize)

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
df3_scale <- readRDS("data/proteins_expression_scaled_data.RDS") # all protein expression data (scaled to mean 0 and sd 1)
res <- readRDS("data/proteins_ttest_diff.RDS") %>% select(protein, p.value, q.value, sig) %>% filter(!duplicated(protein))
qval_sig <- "Gcg"
pval_sig <- unique(df5_sel$protein)


### Heatmaps ###
col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("dodgerblue4", "white", "firebrick")) # color mapping for blocks heatmap

pvalue_col_fun = colorRamp2(c(0, 1, 3), c("lightgrey", "white", "firebrick1")) # color mapping for -log10(pvalue)

df3_matrix <- t(as.matrix(df3_scale[,3:ncol(df3_scale)])) # matrix of all proteins (3rd column until end) for heatmap (only works with matrices)

(pl1 <- Heatmap(df3_matrix, name = "Z-score", 
                col = col_fun, 
                rect_gp = gpar(col = "white", lwd = 2),
                column_title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
                column_split = c(rep("CI + PBS", 5), rep("CI + S.pneu",5)),  
                column_order = colnames(df3_matrix),
                show_column_dend = FALSE,
                clustering_method_rows = "ward.D2",
                top_annotation = HeatmapAnnotation(Treatment = anno_block(gp = gpar(fill = c("darkgrey", "dodgerblue"))), 
                                                    height = unit(0.2, "cm")),
                left_annotation = rowAnnotation(pvalue = anno_simple(-log10(res$p.value), 
                                                                     which = 'row',
                                                                     col = pvalue_col_fun, 
                                                                     pch = res$sig,
                                                                     pt_size = unit(1, "snpc")*0.5
                                                                     ),
                                                annotation_name_side = "bottom",
                                                width = unit(0.5, "cm"),
                                                annotation_name_gp = grid::gpar(fontsize = 8)),
                width = unit(6, "cm"), 
                height = unit(0.25*93, "cm"), 
                row_names_side = "left", 
                column_gap = unit(2, "mm"),
                row_dend_side = "right",
                row_names_gp = grid::gpar(fontsize = 7)))

lgd_pvalue = Legend(title = "p-value", col_fun = pvalue_col_fun, at = c(0, 1, 2, 3), 
                    labels = c("1", "0.1", "0.01", "0.001"))
lgd_sig = Legend(pch = "*", type = "points", labels = "< 0.05")
pl1_anno <- draw(pl1, annotation_legend_list = list(lgd_pvalue, lgd_sig))

pdf("results/pdf/230424_heatmap_all.pdf", width = 8, height = 10)
pl1_anno
dev.off()
svg("results/svg/230424_heatmap_all.svg", width = 8, height = 10)
pl1_anno
dev.off()

# From here on, smaller heatmaps: for only Gcg, and for only sig proteins
font_treatment <- 12 # bigger fonts because heatmap is smaller
font_proteins <- 10
width_pvalue <- 0.5
width_blocks <- 12
height_blocks <- 0.5

# First heatmap for proteins sig after multiple testing (only Gcg)
df3_matrix2 <- t(df3_matrix[which(rownames(df3_matrix) %in% qval_sig),])
rownames(df3_matrix2) <- qval_sig
res_q <- res %>% 
    mutate(sig = case_when(
        q.value < 0.05 ~ paste0("*"),
        q.value >= 0.05 ~ paste0("")
    )) %>% 
    filter(q.value < 0.05)
(pl2 <- Heatmap(df3_matrix2, name = "Z-score", 
                col = col_fun, 
                rect_gp = gpar(col = "white", lwd = 2),
                column_title_gp = grid::gpar(fontsize = font_treatment, fontface = "bold"),
                column_split = c(rep("CI + PBS", 5), rep("CI + S.pneu",5)),  
                column_order = colnames(df3_matrix),
                show_column_dend = FALSE,
                top_annotation = HeatmapAnnotation(Treatment = anno_block(gp = gpar(fill = c("darkgrey", "dodgerblue"))), 
                                                   height = unit(0.2, "cm")),
                left_annotation = rowAnnotation(qvalue = anno_simple(-log10(res_q$q.value), 
                                                                     which = 'row',
                                                                     col = pvalue_col_fun, 
                                                                     pch = res_q$sig),
                                                annotation_name_side = "bottom",
                                                width = unit(width_pvalue, "cm"),
                                                annotation_name_gp = grid::gpar(fontsize = 8)),
                width = unit(width_blocks, "cm"),
                height = unit(height_blocks*nrow(df3_matrix2), "cm"),
                row_names_side = "left", 
                column_gap = unit(2, "mm"),
                row_dend_side = "right",
                row_names_gp = grid::gpar(fontsize = font_proteins)))

lgd_qvalue = Legend(title = "q-value", col_fun = pvalue_col_fun, at = c(0, 1, 2, 3), 
                    labels = c("1", "0.1", "0.01", "0.001"))
pl2_anno <- draw(pl2, annotation_legend_list = list(lgd_qvalue, lgd_sig))

pdf("results/pdf/230424_heatmap_Gcg.pdf", width = 7, height = 3)
pl2_anno
dev.off()
svg("results/svg/230424_heatmap_Gcg.svg", width = 7, height = 3)
pl2_anno
dev.off()

## Now heatmap with all significant proteins BEFORE multiple testing corr + p value
df3_matrix3 <- df3_matrix[which(rownames(df3_matrix) %in% pval_sig),]
rownames(df3_matrix3) <- pval_sig
res_p <- res %>% filter(protein %in% pval_sig)
(pl3 <- Heatmap(df3_matrix3, name = "Z-score", 
                col = col_fun, 
                rect_gp = gpar(col = "white", lwd = 2),
                column_title_gp = grid::gpar(fontsize = font_treatment, fontface = "bold"),
                column_split = c(rep("CI + PBS", 5), rep("CI + S.pneu",5)),  
                column_order = colnames(df3_matrix),
                show_column_dend = FALSE,
                top_annotation = HeatmapAnnotation(Treatment = anno_block(gp = gpar(fill = c("darkgrey", "dodgerblue"))), 
                                                   height = unit(0.2, "cm")),
                left_annotation = rowAnnotation(pvalue = anno_simple(-log10(res_p$p.value), 
                                                                     which = 'row',
                                                                     col = pvalue_col_fun, 
                                                                     pch = res_p$sig,
                                                                     pt_size = unit(1, "snpc")*0.65
                                                                     ),
                                                annotation_name_side = "bottom",
                                                width = unit(width_pvalue, "cm"),
                                                annotation_name_gp = grid::gpar(fontsize = 8)),
                width = unit(width_blocks, "cm"),
                height = unit(height_blocks*nrow(df3_matrix3), "cm"),
                row_names_side = "left", 
                column_gap = unit(2, "mm"),
                row_dend_side = "right",
                row_names_gp = grid::gpar(fontsize = font_proteins)))

pl3_anno <- draw(pl3, annotation_legend_list = list(lgd_pvalue, lgd_sig))

pdf("results/pdf/230424_heatmap_sig.pdf", width = 8, height = 4)
pl3_anno
dev.off()
svg("results/svg/230424_heatmap_sig.svg", width = 8, height = 4)
pl3_anno
dev.off()

## Now heatmap with all significant proteins BEFORE multiple testing corr + q value (=FDR corrected)
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
                                                                     pt_size = unit(1, "snpc")*1),
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

pdf("results/pdf/230424_heatmap_sig_q.pdf", width = 8, height = 4)
pl4_anno
dev.off()
svg("results/svg/230424_heatmap_sig_q.svg", width = 8, height = 4)
pl4_anno
dev.off()


## Now heatmap with all significant proteins BEFORE multiple testing corr WITHOUT p or q value
(pl5 <- Heatmap(df3_matrix3, name = "Z-score", 
                col = col_fun, 
                rect_gp = gpar(col = "white", lwd = 2),
                column_title_gp = grid::gpar(fontsize = font_treatment, fontface = "bold"),
                column_split = c(rep("CI + PBS", 5), rep("CI + S.pneu",5)),  
                column_order = colnames(df3_matrix),
                show_column_dend = FALSE,
                top_annotation = HeatmapAnnotation(Treatment = anno_block(gp = gpar(fill = c("darkgrey", "dodgerblue"))), 
                                                   height = unit(0.2, "cm")),
                width = unit(width_blocks, "cm"),
                height = unit(height_blocks*nrow(df3_matrix3), "cm"),
                row_names_side = "left", 
                column_gap = unit(2, "mm"),
                row_dend_side = "right",
                row_names_gp = grid::gpar(fontsize = font_proteins)))

pdf("results/pdf/230424_heatmap_sig_noval.pdf", width = 8, height = 4)
pl5
dev.off()
svg("results/svg/230424_heatmap_sig_noval.svg", width = 8, height = 4)
pl5
dev.off()
