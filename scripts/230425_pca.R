## PCA for Olink assay
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Libraries
library(tidyverse)
library(factoextra)
library(ggrepel)

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

### Data ###
df3 <- rio::import("data/proteins_expression_data.csv", format = "csv2")

### Principal component analysis ###
res.pca <- prcomp(df3[,3:ncol(df3)], center = TRUE, scale = TRUE)
fviz_eig(res.pca)
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
        scale_color_manual(values = c("black", "red")) +
        guides(fill = guide_legend(override.aes = list(shape = 21, size = 2))) +
        stat_ellipse(aes(color = Treatment), type = "t", level = 0.9) +
        ggtitle("Principal component analysis"))
ggsave("results/pdf/pca_plot_unlabeled.pdf", width = 5, height = 5)
ggsave("results/svg/pca_plot_unlabeled.svg", width = 5, height = 5)

(plot_pca <- df.tot %>% 
        ggplot(aes(PC1, PC2)) +
        geom_point(aes(color = Treatment), size = 2) +
        xlab(x_axis) +
        ylab(y_axis) +
        theme_Publication() +
        scale_color_manual(values = c("black", "red")) +
        guides(fill = guide_legend(override.aes = list(shape = 21, size = 2))) +
        stat_ellipse(aes(color = Treatment), type = "t", level = 0.9))
ggsave("results/pdf/pca_plot_notitle.pdf", width = 5, height = 5)
ggsave("results/svg/pca_plot_notitle.svg", width = 5, height = 5)

(plot_pca_labs <- df.tot %>% 
        ggplot(aes(PC1, PC2)) +
        geom_point(aes(color = Treatment), size = 2) +
        xlab(x_axis) +
        ylab(y_axis) +
        geom_text_repel(aes(label = sampleID)) +
        theme_Publication() +
        scale_color_manual(values = c("black", "red")) +
        guides(fill = guide_legend(override.aes = list(shape = 21, size = 2))) +
        stat_ellipse(aes(color = Treatment), type = "t", level = 0.9) +
        ggtitle("Principal component analysis"))
ggsave("results/pdf/pca_plot_labeled.pdf", width = 5, height = 5)
ggsave("results/svg/pca_plot_labeled.svg", width = 5, height = 5)

