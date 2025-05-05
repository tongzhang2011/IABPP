library(dplyr)
library(tidyr)
library(ggplot2)

source("plotting functions.R")

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

### calculate the mean for each conditions among the three replciates

# Create a vector of condition numbers
conditions <- sprintf("%02d", 1:10)

for (condition in conditions) {
  # Identify columns related to the current condition
  cols <- grep(paste0("_", condition, "$"), names(m), value = TRUE)
  
  # Add mean column to the dataframe
  m <- m %>%
    mutate(!!paste0("mean_", condition) := rowMeans(select(., all_of(cols)), na.rm = TRUE))
}


m <- select(m, starts_with("mean"))
colnames(m)
m <- select(m, -mean_01)

# m <- filter(m, mean_10 < 100) # last channel < 100


abu_plot_area_v2(1)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
strong_target <- c("ACOT1",
                   "MGLL",
                   "LCAT",
                   "EST1C",
                   "CHLE",
                   "EST1D")
length(strong_target)

row.names(m)


table(gsub("_MOUSE", "", 
           gsub(".*\\|", "", row.names(m))) %in% strong_target)



strong_index <- which((gsub("_MOUSE", "", 
                            gsub(".*\\|", "", row.names(m))) %in% strong_target) == TRUE)



# for figure 4
library(gridExtra)
grid.arrange(abu_plot_area_v2(strong_index[6]),
             abu_plot_area_v2(strong_index[1]),
             abu_plot_area_v2(strong_index[2]), ncol = 3)


# for supplemental plot
m <- m[my_plot_list,]

# reorder the proteins based on area
m <- m[match(my_plot_list, row.names(m)), ]



plot_list <- list()
for (i in 1:15) {
  plot_list[[i]]  <- abu_plot_area_v2(i)
}

grid.arrange(
  grobs = plot_list,    # List of plots
  ncol = 5              # Number of columns in the grid
)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 





# Define the number of plots you want to create
num_plots <- 41

# Specify the file name for the PDF
pdf("41_plots.pdf", width = 8, height = 6)

# Loop to create and save plots
for (i in 1:num_plots) {
  # Generate a plot
  p <- abu_plot_area_v2(i)
  
  # Print the plot to the PDF
  print(p)
}

dev.off()



abu_plot_area_v2(strong_index[5])



my_area <-  numeric(dim(m)[1])
for (i in 1: dim(m)[1]) {
  my_area[i] <-  area_calc(i)
}

m1 <- m
m1$area <- my_area
m1 <- m1[order(m1$area),]
write.csv(m1, "all_41_proteins_with_area_information.csv", row.names = T)

m1 <- filter(m1, area < 0.8) # 15 proteins

my_plot_list <- row.names(m1)
