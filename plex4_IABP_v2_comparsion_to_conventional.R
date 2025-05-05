library(dplyr)
m <- read.csv("IABP_IABP_plex4.txt", sep = "\t",
              header = T, stringsAsFactors = F)
dim(m)

table(apply(m, 1, function(x) sum(is.na(x))))
m <- m[complete.cases(m),]

m$mean_NP <- apply(m[,1:3], 1, function(x)  mean(x))
m$mean_ABP <- apply(m[,4:6], 1, function(x)  mean(x))
m$mean_Cmp <- apply(m[,7:10], 1, function(x)  mean(x))

min <- min(c(m$mean_NP,
             m$mean_ABP,
             m$mean_Cmp))

max <- max(c(m$mean_NP,
             m$mean_ABP,
             m$mean_Cmp))



par(mfrow = c(3,1))

hist(m$mean_NP, breaks = seq(min, max, length.out = 100), main = "NP", xlab = "")
abline(v = median(m$mean_NP), col = "red", lty = 2, lwd =2)

hist(m$mean_ABP, breaks = seq(min, max, length.out = 100), main = "ABP", xlab = "")
abline(v = median(m$mean_ABP), col = "red", lty = 2, lwd =2)

hist(m$mean_Cmp, breaks = seq(min, max, length.out = 100), main = "Competitor", xlab = "")
abline(v = median(m$mean_Cmp), col = "red", lty = 2, lwd =2)


m$SN <- m$mean_ABP - m$mean_NP
m$drop <- m$mean_Cmp - m$mean_ABP

m <- filter(m, SN > 1)
m <- filter(m, drop < log(0.8,2))

write.csv(m, "IABP_19_protein_list.csv")

iabp_list <- row.names(m)

m_conventional <- read.csv("all_41_proteins_with_area_information.csv", header = T, stringsAsFactors = F,
                           row.names = 1)
m_conventional <- filter(m_conventional, area < 0.8)
old_list <- row.names(m_conventional)

library(VennDiagram)
venn.diagram(
  x = list(iabp_list, old_list),
  category.names = c("" , "" ),
  filename = '#14_venn_diagramm.png',
  fill = c("yellow", "blue"),
  cex = 3,
  output=TRUE
)


## merge the two datafrmae
intersect(row.names(m), row.names(m_conventional))
setdiff(row.names(m), row.names(m_conventional))
setdiff(row.names(m_conventional), row.names(m))
        
write.csv(data.frame(protein = setdiff(row.names(m), row.names(m_conventional))),
          "Unique_proteins_in_IABP.csv")

m$protein <- row.names(m)
m_conventional$protein <- row.names(m_conventional)

md <- merge(m, m_conventional, by = "protein")
dev.off()
plot(2^md$drop, md$area)

colnames(md)
data <- select(md, protein, drop, area)
data$protein <- gsub(".*\\|", "", data$protein)
data$protein <- gsub("_MOUSE", "", data$protein)
data$drop <- 2^data$drop



# Create scatter plot with labels
ggplot(data, aes(x = area, y = drop)) +
  geom_point(size = 3) +  
  coord_cartesian(xlim = c(0.1, 0.9), ylim = c(0.1, 0.9)) + 
  
  geom_text(aes(label = protein), 
          
            vjust = -0.5, 
            hjust = 1,
            size = 3) +  
  labs(
       x = "AUC from conventional ABBP",
       y = "Competitor / ABP from IABP") +
  theme_minimal() +
  theme(panel.border = element_rect(size = 1, colour = "black", fill = NA))


# Create scatter plot with fitting line

  
  
  ggplot(data, aes(x = area, y = drop)) +
  geom_point(size = 3) +  
  coord_cartesian(xlim = c(0.1, 0.9), ylim = c(0.1, 0.9)) + 
    geom_smooth(method = "lm", se = FALSE, color = "blue") +  
  labs(
    x = "AUC from conventional ABBP",
    y = "Competitor / ABP from IABP") +
  theme_minimal() +
  theme(panel.border = element_rect(size = 1, colour = "black", fill = NA))+
    annotate("text",
             x = 0.25, y = 0.75,
            label = paste0("r = ", round(cor(data$drop, data$area), 2)))


