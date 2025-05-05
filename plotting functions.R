library(dplyr)
library(tidyr)
library(ggplot2)

m <- read.table("IABP_three_conventional_TMT_plex.txt",  header = T, stringsAsFactors = F, sep = "\t")
colnames(m)

abu_plot <- function(i) {
  title <- gsub("_MOUSE", "", 
                gsub(".*\\|", "", row.names(m)[i]))
  
  y_min <- min(as.numeric(m[i,1:30]), na.rm = T)
  y_max <- max(as.numeric(m[i,1:30]), na.rm = T)
  
  plot(as.numeric(m[i,1:10]), type = "l", col = "black",
       ylim = c(y_min, y_max),
       xaxt = "n", xlab='',
       ylab = "Relative Abundance",
       main = title)
  lines(1:10, as.numeric(m[i,11:20]),  col = "blue")
  lines(1:10, as.numeric(m[i,21:30]), col = "red")
  
  points(1:10, as.numeric(m[i, 1:10]),  col = "black")
  points(1:10, as.numeric(m[i,11:20]),  col = "blue")
  points(1:10, as.numeric(m[i,21:30]), col = "red")
  
  text(x = 1:10,
       y = par("usr")[3] - 0.1*(y_max - y_min),
       labels = c("NP",
                  "ABP",
                  "Par, C1",
                  "Par, C2",
                  "Par, C3",
                  "Par, C4",
                  "Par, C5",
                  "Par, C6",
                  "Par, C7",
                  "Par, C8"),
       xpd = NA,
       ## Rotate the labels by 35 degrees.
       srt = 35,
       cex = 0.75)
}



abu_plot_area <- function(i) {
  x <- 0:8 # x axis 
  y <- as.numeric(m[i,1:9])
  # to close the polygon, need additional three points
  p1.x <- 8
  p1.y <- y[length(y)]
  
  p2.x <- 0
  p2.y <- p1.y <- y[length(y)]
  
  p3.x <- 0
  p3.y <- y[1]
  
  x = c(x, p1.x, p2.x, p3.x)
  y = c(y, p1.y, p2.y, p3.y)
  
  n <- length(x)
  
  # Shoelace formula to calculate the area
  area <- abs(sum(x * c(y[-1], y[1]) - y * c(x[-1], x[1])) / 2)
  relative_AUC = sprintf("%.2f", area / 900)
  
  # finish area calculation
  
  
  # Create a data frame for plotting
  df <- data.frame(x = x, y = y)
  
  
  ### some parameters for plotting ##############
  y_min <- min(as.numeric(m[i,1:9]), na.rm = T)
  y_max <- max(as.numeric(m[i,1:9]), na.rm = T)
  
  title <- gsub("_MOUSE", "", 
                gsub(".*\\|", "", row.names(m)[i]))
  
  data_point_color <- c(rep("black",9), "black", "white", "black")
  
  mylables <- c("ABP",
                paste0("Con ", 1:8))
  ### some parameters for plotting ##############
  
  
  # Create the plot using ggplot2
  p1 <- ggplot(df, aes(x = x, y = y)) +
    ylim(c(0,y_max))+
    geom_polygon(fill = "lightblue", color = "white") +
    geom_point(color = data_point_color)+
    labs(x = "", 
         y = "Relative Abundance",
         title = title) + 
    theme_minimal() +
    theme(axis.title.x = element_text(size = 12, color = "black"),
          axis.title.y = element_text(size = 12, color = "black"),
          axis.text.x = element_text(size = 10, color = "black"),
          axis.text.y = element_text(size = 10, color = "black"))+
    geom_abline(intercept = c(y_min, 100), slope = 0, color = "blue", linetype = "dashed")+
    scale_x_continuous(breaks = c(0:8),
                       labels = mylables)+
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
          axis.text.x = element_text(angle = 90, 
                                     color = "black"),
          plot.title = element_text(hjust = 0.5,
                                    size = 13, face = 2,
                                    color = "black"))+
    annotate("text", x=7, y=95, 
             label= paste0("AUC = ", relative_AUC))
  
  return(p1)
}


# data m will be used for plotting
# i <- 1

abu_plot_area_v2 <- function(i) {
  x <- 0:8 # x axis 
  y <- as.numeric(m[i,1:9])
  # to close the polygon, need additional three points
  p1.x <- 8
  p1.y <- 0
  
  p2.x <- 0
  p2.y <- 0
  
  p3.x <- 0
  p3.y <- y[1]
  
  x = c(x, p1.x, p2.x, p3.x)
  y = c(y, p1.y, p2.y, p3.y)
  
  n <- length(x)
  
  # Shoelace formula to calculate the area
  area <- abs(sum(x * c(y[-1], y[1]) - y * c(x[-1], x[1])) / 2)
  relative_AUC = sprintf("%.2f", area / 900)
  
  # finish area calculation
  
  
  # Create a data frame for plotting
  df <- data.frame(x = x, y = y)
  
  
  ### some parameters for plotting ##############
  y_min <- min(as.numeric(m[i,1:9]), na.rm = T)
  y_max <- max(as.numeric(m[i,1:9]), na.rm = T)
  
  title <- gsub("_MOUSE", "", 
                gsub(".*\\|", "", row.names(m)[i]))
  
  data_point_color <- c(rep("black",9), "white", "white", "black")
  
  # mylables <- c("ABP",
  #               paste0("Con ", 1:8))
  
  
  mylables <- c("0 (ABP)",
                "0.01",
                "0.05",
                "0.1",
                "0.2",
                "0.5",
                "1",
                "5",
                "10")
  
  
  
  ### some parameters for plotting ##############
  
  
  # Create the plot using ggplot2
  p1 <- ggplot(df, aes(x = x, y = y)) +
    # ylim(c(0,y_max))+   # use the same scale
    ylim(c(0,130))+
    geom_polygon(fill = "lightblue", color = "white") +
    geom_point(color = data_point_color)+
    labs(x = "Paraoxon (µM)", 
         y = "Relative Abundance",
         title = title) + 
    theme_minimal() +
    theme(axis.title.x = element_text(size = 12, color = "black"),
          axis.title.y = element_text(size = 12, color = "black"),
          axis.text.x = element_text(size = 10, color = "black"),
          axis.text.y = element_text(size = 10, color = "black"))+
    geom_abline(intercept = c(y_min, 100), slope = 0, color = "blue", linetype = "dashed")+
    scale_x_continuous(breaks = c(0:8),
                       labels = mylables)+
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
          axis.text.x = element_text(angle = 90, 
                                     color = "black"),
          plot.title = element_text(hjust = 0.5,
                                    size = 13, face = 2,
                                    color = "black"))+
    annotate("text", x=6, y= 120, 
             label= paste0("AUC = ", relative_AUC))+
    theme(axis.title.x = element_text(size = 12, color = "black", margin = margin(t = 10)))
  
  return(p1)
}


# a function to just calculate the area 
area_calc <- function(i) {
  x <- 0:8 # x axis 
  y <- as.numeric(m[i,1:9])
  # to close the polygon, need additional three points
  p1.x <- 8
  p1.y <- 0
  
  p2.x <- 0
  p2.y <- 0
  
  p3.x <- 0
  p3.y <- y[1]
  
  x = c(x, p1.x, p2.x, p3.x)
  y = c(y, p1.y, p2.y, p3.y)
  
  n <- length(x)
  
  # Shoelace formula to calculate the area
  area <- abs(sum(x * c(y[-1], y[1]) - y * c(x[-1], x[1])) / 2)
  relative_AUC = sprintf("%.2f", area / 900)
  
  # finish area calculation
  
  return(relative_AUC)
}
