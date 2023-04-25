## Data processing for Olink assay
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Libraries
library(tidyverse)

## Data
samples <- rio::import("data/Q-01893_Olink-Sample Manifest1.xlsx")
samples$sampleID <- as.numeric(samples$SampleID)
samples <- samples %>% select(sampleID, Treatment)
df <- rio::import("data/Q-01893_Ramezani Rad_NPX.xlsx", skip = 3)
names(df)
df <- df[,1:93]

proteins <- df[1:2,]
proteins <- as.data.frame(t(as.matrix(proteins)))
dim(proteins)
head(proteins)
names(proteins)
colnames(proteins) <- c(proteins$`1`[1], proteins$`2`[1])
proteins <- proteins[2:nrow(proteins),]
proteins$protein <- rownames(proteins)
head(proteins)

df2 <- df[4:15,]
rownames(df2) <- df2$Assay
df2 <- df2[,-1]
#df2 <- as.data.frame(t(as.matrix(df2)))
df2$sampleID <- rownames(df2)
df2 <- df2 %>% filter(!str_detect(sampleID, "CONTROL")) %>% 
    mutate(across(everything(.), as.numeric))

df3 <- left_join(samples, df2, by = "sampleID") %>% 
    mutate(Treatment = as.factor(Treatment)) %>% 
    arrange(Treatment)

df3_scale <- df3 %>% mutate(across(3:ncol(.), ~scale(.x)))

df4 <- df3
res <- c()
for(a in 3:94){
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

df.sum <- df4 %>% group_by(Treatment) %>% summarise(across(2:93, 
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
write_csv2(df3_scale, "data/proteins_expression_scaled_data.csv")

saveRDS(df5, "data/proteins_ttest_diff.RDS")
saveRDS(df3, "data/proteins_expression_data.RDS")
saveRDS(df3_scale, "data/proteins_expression_scaled_data.RDS")
