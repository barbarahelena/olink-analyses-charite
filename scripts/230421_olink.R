## Data processing for Olink assay
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Libraries
library(tidyverse)

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
samples <- rio::import("data/Q-01893_Olink-Sample Manifest1.xlsx")
samples$sampleID <- samples$SampleID
samples <- samples %>% select(sampleID, Treatment)
df <- rio::import("data/Q-01893_Ramezani Rad_NPX.xlsx", skip = 3)
names(df)
df <- df[,1:93]

proteins <- df[c(1:2, 17:19),]
proteins <- as.data.frame(t(as.matrix(proteins)))
dim(proteins)
head(proteins)
names(proteins)
colnames(proteins) <- c(proteins$`1`[1], proteins$`2`[1], 
                        proteins$`17`[1], proteins$`18`[1],
                        proteins$`19`[1])
proteins <- proteins[2:nrow(proteins),]
proteins$protein <- rownames(proteins)
head(proteins)
proteins <- proteins %>% mutate(across(3:4, as.numeric))

df2 <- df[4:15,]
rownames(df2) <- df2$Assay
df2 <- df2[,-1]
#df2 <- as.data.frame(t(as.matrix(df2)))
df2$sampleID <- rownames(df2)
df2 <- df2 %>% 
    #filter(!str_detect(sampleID, "CONTROL")) %>% 
    mutate(across(1:(ncol(.)-1), as.numeric)) %>% 
    full_join(samples, ., by = "sampleID")

## Check if proteins fall under limit of detection (LOD)
res_lod <- c()
for(a in 3:94){
    proteinname <- colnames(df2)[a]
    values <- df2[,proteinname]
    values_controls <- df2 %>% filter(str_detect(sampleID, "CONTROL")) %>% select(all_of(proteinname))
    values_trA <- df2 %>% filter(str_detect(Treatment, "PBS")) %>% select(all_of(proteinname))
    values_trB <- df2 %>% filter(str_detect(Treatment, "S.pn.")) %>% select(all_of(proteinname))
    lod <- proteins$LOD[which(proteins$protein == proteinname)]
    no_lod <- sum(values < (lod))
    no_lod_controls <- sum(values_controls < lod)
    no_lod_trA <- sum(values_trA < lod)
    no_lod_trB <- sum(values_trB < lod)
    row_lod <- c(proteinname, no_lod, no_lod_controls, no_lod_trA, no_lod_trB)
    names(row_lod) <- c("protein", "under_lod", "under_lod_controls", "under_lod_PBS", "under_lod_Strep")
    res_lod <- rbind(res_lod, row_lod)
    proteinname <- NULL
}

res_lod <- as.data.frame(res_lod)
rownames(res_lod) <- res_lod$protein
res_lod <- res_lod %>% mutate(across(2:5, as.numeric))
res_lod_long <- res_lod %>% pivot_longer(., cols = 3:5, names_to = "group", values_to = "number") %>% 
    mutate(group = str_remove(group, "under_lod_")) %>% 
    arrange(under_lod) %>% 
    mutate(protein = as.factor(protein),
           protein = forcats::fct_inorder(protein)) %>% 
    filter(under_lod != 0)

ggplot(data = res_lod_long, aes(x = protein, y = number, fill = group)) +
    geom_bar(position = "stack", stat = "identity") +
    scale_fill_manual(values = c("firebrick", "darkgrey", "dodgerblue")) +
    scale_y_continuous(breaks = 1:12) +
    geom_hline(aes(yintercept = 12), linetype = "dashed") +
    coord_flip() +
    theme_Publication() +
    labs(y = "number below LOD", title = "Number of proteins below LOD per group")
ggsave("results/proteins_belowLOD.pdf", width = 5, height = 8)

# How many proteins have less than 5 missings for the groups of interest (PBS + Strep)
res_lod %>% filter((under_lod_PBS + under_lod_Strep) < 5) %>% nrow(.) # 71
protein_lod <- res_lod %>% filter((under_lod_PBS + under_lod_Strep) < 5) %>% select(protein)

# 
df3 <- df2 %>% select(sampleID, Treatment, all_of(protein_lod$protein)) %>% 
    filter(!str_detect(sampleID, "CONTROL")) %>% 
    mutate(Treatment = as.factor(Treatment)) %>% 
    arrange(Treatment)

df3_scale <- df3 %>% mutate(across(3:ncol(.), ~scale(.x)))

df4 <- df3
res <- c()
for(a in 3:73){
    df4$prot <- df4[,a]
    t.res <- t.test(df4$prot ~ df4$Treatment)
    res.row <- c(colnames(df4)[a], t.res$p.value)
    names(res.row) <- c("protein", "p.value")
    res <- rbind(res, res.row)
    df4$prot <- NULL
}

res <- as.data.frame(res)
rownames(res) <- res$protein
res <- res %>% mutate(protein = as.factor(protein), p.value = as.numeric(p.value))
res$q.value <- p.adjust(res$p.value, method = "fdr")
nrow(res %>% filter(p.value < 0.05))
nrow(res %>% filter(q.value < 0.05))

res <- res %>% mutate(sig = case_when(
    p.value >= 0.05 ~ paste0(""),
    p.value < 0.0001 ~ paste0("****"),
    p.value < 0.001 ~ paste0("***"),
    p.value < 0.01 ~ paste0("**"),
    p.value < 0.05 ~ paste0("*")
))

df.sum <- df4 %>% group_by(Treatment) %>% summarise(across(2:72, 
                                                    list(mean = ~mean(.x, na.rm = TRUE),
                                                         sd = ~sd(.x, na.rm = TRUE),
                                                         n = ~length(.x)), 
                                                    .names = "{.col}_{.fn}"))

df.sumlong <- df.sum %>% pivot_longer(2:ncol(.), names_to = c("protein", "stat"),
                                      names_sep = "_",
                                      values_to = "number") %>% 
    pivot_wider(id_cols = 1:2, values_from = "number", names_from = "stat")

df5 <- left_join(res, df.sumlong, by = "protein") 
df5_sel <- df5 %>% 
    filter(p.value < 0.05) %>% 
    arrange(q.value)
qval_sig <- "Gcg"
pval_sig <- unique(df5_sel$protein)

write_csv2(df5, "data/proteins_ttest_diff.csv")
write_csv2(df3, "data/proteins_expression_data.csv")
#write_csv2(df3_scale, "data/proteins_expression_scaled_data.csv")

saveRDS(df5, "data/proteins_ttest_diff.RDS")
saveRDS(df3, "data/proteins_expression_data.RDS")
saveRDS(df3_scale, "data/proteins_expression_scaled_data.RDS")
