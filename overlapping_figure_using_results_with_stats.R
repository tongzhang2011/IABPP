library(dplyr)
x1 <- read.csv("Target_list_from_IABP_analysis.csv",
              header = T, stringsAsFactors = F)

x2 <- read.csv("Target_list_from_conventional_analysis.csv",
               header = T, stringsAsFactors = F)

iabp_list <- x1$X
old_list <- x2$X

library(VennDiagram)
venn.diagram(
  x = list(iabp_list, old_list),
  category.names = c("" , "" ),
  filename = 'Overlapping between IABP and conventional.png',
  fill = c("yellow", "blue"),
  cex = 3,
  output=TRUE
)

