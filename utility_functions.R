box_plot <- function(data,my_title) {
  
  dat <- data
  
  boxplot(dat, xaxt = "n",
          ylab = "Relative Intensity",
          main = my_title)
  
  xmin <- par("usr")[1]
  xmax <- par("usr")[2]
  ymin <- par("usr")[3]
  ymax <- par("usr")[4]
  
  y_pos <- ymin*0.95
  
  axis(side = 1, labels = FALSE)
  text(x = 1:length(dat),
       y = y_pos,
       labels = colnames(dat),
       xpd = NA,
       srt = 60,
       adj = 1)
  
  
  text(x = 1:length(dat),
       y = apply(dat, 2, function(x) median(x, na.rm = T)) + 0.5,
       labels = sprintf("%.1f", 
                        apply(dat, 2, function(x) median(x, na.rm = T))),
       col = "red",
       font = 2,
       cex = 0.8)
  
}

global_median <- function(m) {
  # mc <- m[complete.cases(m), ]
  sample.bias <- apply(m, 2, median, na.rm = T)
  m <- sweep(m, 2, sample.bias, "-")
  
  median_of_samples <- median(sample.bias, na.rm = T)
  
  m <- m + median_of_samples
  
  return(m)
}


##################################################################################

####################################################################################calcluate p and FC
# step 1: define t.test and FC for this case
###########################################
#a t.test function that requires at least 2 valid numbers from both groups
at_least_number <- 2

t.test_at_least_2_v2 <- function(num1, num2) {
  if( sum(!is.na(num1)) >= at_least_number &  sum(!is.na(num2)) >= at_least_number ) {
    my_p.value <- t.test(num1,num2,na.rm=TRUE,var.equal=TRUE)$p.value
  } else {  my_p.value <- NA
  }
  return( my_p.value)
}

#apply this funtion to two dataframe.
#df1 contains the first group;
#df2 contains the second group;
t.test_v2_df <- function(df1,df2) {
  n <- dim(df1)[1]
  my_p_value <- vector(length = n)
  
  for (i in 1:n) {
    my_p_value[i] <- t.test_at_least_2_v2(as.numeric(df1[i,]), as.numeric(df2[i,]))
  }
  return(my_p_value)
}
#end of the function
##########################################

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#similarly, FC also requires at least two 
#valid values from each comparing group
fc_at_least_2  <- function(num1, num2) {
  if( sum(!is.na(num1)) >= at_least_number &  sum(!is.na(num2)) >= at_least_number ) {
    fc <-  mean(num2, na.rm=TRUE) - mean(num1, na.rm=TRUE)
  } else { fc <- NA}
  return( fc)
}

#define a function use two dataframes as input
#i.e: do the calculation at the datafrmae level
FC <- function(df1, df2) {
  n <- dim(df1)[1]
  my_fc <- vector(length = n)
  
  for (i in 1:n) {
    my_fc[i] <- fc_at_least_2(as.numeric(df1[i,]), as.numeric(df2[i,]))
  }
  return(my_fc)
}
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


################################################################################################






vplot <- function(d1,region1, region2, p_cut = 0.01, FC_cut = log2(1.5)) {
  ROI1 <- region1
  ROI2 <- region2
  
  my_title <- paste(ROI2, ROI1, sep = " / ")
  
  df1 <-  dplyr::select(d1, starts_with(ROI1))
  df2 <-  dplyr::select(d1, starts_with(ROI2))
  
  df <- cbind(df1, df2)
  
  df$p_value <- t.test_v2_df(df1, df2)
  df$FC <- FC(df1, df2)
  
  df$p_adj <- p.adjust(df$p_value, method = "none") # adjust the p value
  df$p_value <-  df$p_adj  # just replace the pvalue column with the adjusted p value
  
  valid_p <- as.numeric(table(is.na(df$p_value))[1])
  my_title <- paste(my_title,  
                    paste0("(n = ", valid_p, ")"))
  
  
  #  p_cut <- 0.01
  # FC_cut <- 0.5
  
  
  d_up <- filter(df, p_value < p_cut, FC > FC_cut)
  d_dn <- filter(df, p_value < p_cut, FC < -FC_cut)
  
  Increase <- dim(d_up)[1]
  Decrease <- dim(d_dn)[1]
  
  point_size <- 1
  
  # create a plot
  with(df,plot(FC, -log(p_value,10),pch = 20,
               xlab = "",
               ylab = "",
               cex.axis = 1, font.axis = 1,
               cex.lab = 0.8, font.lab = 1,
               cex = point_size))
  
  xmin <- par("usr")[1]
  xmax <- par("usr")[2]
  ymin <- par("usr")[3]
  ymax <- par("usr")[4]
  
  xrange <- xmax - xmin
  yrange <- ymax - ymin
  xlab_pos_1 <- xmin + 0.2*xrange
  xlab_pos_2 <- xmin + 0.8*xrange
  ylab_pos <- ymin + 0.9*yrange
  
  
  with(subset(df,p_value > p_cut),points(FC, -log(p_value,10),pch = 20, col = "gray",
                                         cex = point_size))
  with(subset(df,p_value < p_cut  & FC > FC_cut ),points(FC, -log(p_value,10),pch = 20, col = "red",
                                                         cex = point_size))
  with(subset(df,p_value < p_cut  & FC < -FC_cut ),points(FC, -log(p_value,10),pch = 20, col = "blue",
                                                          cex = point_size))
  abline(h = -log(p_cut,10),lty =2 ,col = "black")
  abline(v = FC_cut, lty =2 ,col = "red")
  abline(v = -FC_cut, lty =2 ,col = "blue")
  box(lwd = 2)
  text(xlab_pos_1, ylab_pos, paste("Decrease: ", Decrease, sep = ""),
       cex = 1.2, col = "blue")
  text(xlab_pos_2, ylab_pos, paste("Increase: ", Increase, sep = ""),
       cex = 1.2, col = "red")
  title(main = my_title)
  
}





sig_list <- function(d1, region1, region2, p_cut = 0.01, FC_cut = 1) {
  ROI1 <- region1
  ROI2 <- region2
  
  my_title <- paste(ROI2, ROI1, sep = "_vs_")
  
  df1 <-  dplyr::select(d1, starts_with(ROI1))
  df2 <-  dplyr::select(d1, starts_with(ROI2))
  
  df <- cbind(df1, df2)
  
  df$p_value <- t.test_v2_df(df1, df2)
  df$FC <- FC(df1, df2)
  
  
  df$p_adj <- p.adjust(df$p_value, method = "none") # adjust the p value
  df$p_value <-  df$p_adj  # just replace the pvalue column with the adjusted p value
  
  valid_p <- as.numeric(table(is.na(df$p_value))[1])
  my_title <- paste(my_title,  
                    paste0("(n = ", valid_p, ")"))
  
  
  df_sig <- filter(df,p_value < p_cut)
  df_sig <- filter(df_sig,abs(FC) > FC_cut)
  df_sig <- filter(df_sig,FC < 0) # decrese ones only
  
  write.csv(df_sig, paste(my_title, "_putative_target.csv", sep = "_"))
  
}









