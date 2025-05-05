
# within group normalization #####################

rm(list = ls())

library(dplyr)
source("utility_functions.R")
################ within group normalization before stats ###################
m <- read.csv("IABP_IABP_plex4.txt", sep = "\t",
              header = T, stringsAsFactors = F, row.names = 1)
dim(m)

colnames(m)
row.names(m)
m <- m[complete.cases(m),]

# assign sample names

groups <- c("NP",
            "ABP",
            "PAR")

name1 <- paste(rep(groups, c(3,3,4)),
               c(rep(paste0("R", 1:3), 2),
                 rep(paste0("R", 1:4), 1)),
               sep = "_")


colnames(m) <- name1

my_groups <- unique(gsub("_.*", "", colnames(m)))

my_list <- vector("list", length = length(my_groups))

for(i in 1:length(my_groups)) {
  dt <- dplyr::select(m, starts_with(my_groups[i]))
  dt1 <- global_median(dt)
  my_list[[i]] <- dt1
}

m1 <- do.call(cbind, my_list)
apply(m1, 2, median, na.rm = T)

par(mfrow = c(1,2), mar = c(4.5, 5.5 ,3,3))
box_plot(m, "Before within-group norm")
box_plot(m1, "After within-group norm")
dim(m1)

table(is.na(m1))

###### peform stats ################
# Loop through the comparison groups and perform row-wise t-tests

rowwise_ttest <- function(df, group1_cols, group2_cols) {
  apply(df[, c(group1_cols, group2_cols)], 1, function(x) {
    t.test(x[1:length(group1_cols)], x[(length(group1_cols) + 1):length(x)])$p.value
  })
}


rowwise_log2FC <- function(df, group1_cols, group2_cols) {
  apply(df[, c(group1_cols, group2_cols)], 1, function(x) {
    mean(x[(length(group1_cols) + 1):length(x)]) - mean(x[1:length(group1_cols)])
  })
}

colnames(m1)

par_vs_abp_p <- rowwise_ttest(m1, c("ABP_R1", "ABP_R2", "ABP_R3"), c("PAR_R1", "PAR_R2", "PAR_R3", "PAR_R4"))

par_vs_abp_log2FC <- rowwise_log2FC(m1, c("ABP_R1", "ABP_R2", "ABP_R3"), c("PAR_R1", "PAR_R2", "PAR_R3", "PAR_R4"))


# define target
abp_vs_np_p  <- rowwise_ttest(m1,   c("NP_R1", "NP_R2", "NP_R3"), c("ABP_R1", "ABP_R2", "ABP_R3"))
abp_vs_np_log2FC  <- rowwise_log2FC(m1, c("NP_R1", "NP_R2", "NP_R3"),  c("ABP_R1", "ABP_R2", "ABP_R3"))

dev.off

vplot(m1, "NP", "ABP")
vplot(m1, "ABP", "PAR")


# change the cutoff

vplot(m1, "NP", "ABP", 0.05, log2(2))
vplot(m1, "ABP", "PAR", 0.05, log2(1.2))



# output


# combine the data
df_stat <- data.frame(abp_vs_np_p = abp_vs_np_p,
                      abp_vs_np_log2FC = abp_vs_np_log2FC,
                      par_vs_abp_p = par_vs_abp_p,
                      par_vs_abp_log2FC = par_vs_abp_log2FC)

m2 <- cbind(m1, df_stat)

write.table(m2, "IABP_IABP_plex4_with_stats.txt", sep = "\t", row.names = TRUE, col.names = NA)

# select sig for using two criteira
m2a <- filter(m2,abp_vs_np_p < 0.05, abp_vs_np_log2FC > log(1.2, 2) )
sig_par <- filter(m2a, par_vs_abp_log2FC < -log(1.2, 2), par_vs_abp_p < 0.05)

write.csv(sig_par, "Target_list_from_IABP_analysis.csv")
