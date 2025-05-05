rm(list = ls())
library(dplyr)
library(tidyr)
library(ggplot2)
source("utility_functions_1.R")

m <- read.table("IABP_three_conventional_TMT_plex.txt",  header = T, stringsAsFactors = F, sep = "\t")
colnames(m)

m <- m[complete.cases(m),]
m <- 2^m

# to use the ABP data as 100 for relative protein abidance
m[,1:10] <- m[,1:10] / m[,2]
m[,11:21] <- m[,11:21] /  m[,12]
m[,21:30] <- m[,21:30] / m[,22]

m <- m*100

colnames(m)

groups <- c("NP",
            "ABP",
            "PAR_C1",
            "PAR_C2",
            "PAR_C3",
            "PAR_C4",
            "PAR_C5",
            "PAR_C6",
            "PAR_C7",
            "PAR_C8")

name1 <- paste(rep(groups, 3),
               rep(paste0("R", 1:3), each = 10),
               sep = "_")


colnames(m) <- name1

m1 <- m

###### peform stats ################
# Loop through the comparison groups and perform row-wise t-tests

rowwise_ttest <- function(df, group1_cols, group2_cols) {
  apply(df[, c(group1_cols, group2_cols)], 1, function(x) {
    t.test(x[1:length(group1_cols)], x[(length(group1_cols) + 1):length(x)])$p.value
  })
}

# note should be divide 
rowwise_log2FC <- function(df, group1_cols, group2_cols) {
  apply(df[, c(group1_cols, group2_cols)], 1, function(x) {
    mean(x[(length(group1_cols) + 1):length(x)]) /  mean(x[1:length(group1_cols)])
  })
}

colnames(m1)
# do the highest concentration 

par_vs_abp_p <- rowwise_ttest(m1, c("ABP_R1", "ABP_R2", "ABP_R3"), c("PAR_C8_R1", "PAR_C8_R2", "PAR_C8_R3"))

par_vs_abp_log2FC <- rowwise_log2FC(m1, c("ABP_R1", "ABP_R2", "ABP_R3"), c("PAR_C8_R1", "PAR_C8_R2", "PAR_C8_R3"))
par_vs_abp_log2FC <- log(par_vs_abp_log2FC, 2)

# define target
abp_vs_np_p  <- rowwise_ttest(m1,   c("NP_R1", "NP_R2", "NP_R3"), c("ABP_R1", "ABP_R2", "ABP_R3"))
abp_vs_np_log2FC  <- rowwise_log2FC(m1, c("NP_R1", "NP_R2", "NP_R3"),  c("ABP_R1", "ABP_R2", "ABP_R3"))
abp_vs_np_log2FC <- log(abp_vs_np_log2FC, 2)


# change the cutoff

vplot(m1, "NP", "ABP", 0.05, log2(1.2))
vplot(m1, "ABP", "PAR", 0.05, log2(1.2))




# combine the data
df_stat <- data.frame(abp_vs_np_p = abp_vs_np_p,
                      abp_vs_np_log2FC = abp_vs_np_log2FC,
                      par_vs_abp_p = par_vs_abp_p,
                      par_vs_abp_log2FC = par_vs_abp_log2FC)

m2 <- cbind(m1, df_stat)

write.table(m2, "IABP_three_conventional_TMT_with_stats.txt", sep = "\t", row.names = TRUE, col.names = NA)

# select sig for using two criteira
m2a <- filter(m2,abp_vs_np_p < 0.05, abp_vs_np_log2FC > log(1.2, 2) )
sig_par <- filter(m2a, par_vs_abp_log2FC < -log(1.2, 2), par_vs_abp_p < 0.05)

write.csv(sig_par, "Target_list_from_conventional_analysis.csv")


