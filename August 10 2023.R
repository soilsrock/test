## Added this code to make package install cleaner
list.of.packages <- c("caret", "randomForest", "DescTools", "dplyr","raster","epiR","tcltk", "data.table","rgdal","ggpubr","sf","sp", "ggplot2", "rasterVis")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

setwd("D:/IHF/IHF_Research/RandomForest")
#install.packages("caret")
library(caret)
library(tcltk)
#install.packages("randomForest")
library(randomForest)
#install.packages("raster")
library(raster)
library(data.table)
#install.packages("rgdal")
library(rgdal)
#install.packages("ggpubr")
library(ggpubr)
#install.packages("sf")
library(sf)
#install.packages("sp")
library(sp)
library(ggplot2)
library(rasterVis)
library(epiR)
library(dplyr)
library(DescTools)
setwd("D:/IHF/IHF_Research/RandomForest")
###Rasters as covariates
flow.path<-raster("FPL.tif" )
LS<-raster("LSF.tif")
landform<-raster("TPI.tif")
saga.wetness.index<-raster("SWI.tif")
aspect<-raster("Aspect.tif")
em<-raster("EM.tif")
DEM<-raster("DEM.tif")
Multiresolution.Index.of.RidgeTop.Flatness=raster("MRRTF.tif")
Multiresolution.Index.of.Valley.Bottom.Flatness=raster("MRVBF.tif")
RGB.2021<-raster("PCA.tif")

### import shapefile of the soil points on farm
#shp<-shapefile("D:/IHF_Research/Phase 2 texture data/new_seventy.shp")
## Changed path relative to current directory
setwd("D:/IHF/IHF_Research/Phase 2 texture data")
shp<-shapefile("new_seventy.shp")
plot(shp)



#
rasStack=stack(LS,flow.path,aspect,landform,saga.wetness.index,DEM,em,
               RGB.2021,Multiresolution.Index.of.RidgeTop.Flatness,Multiresolution.Index.of.Valley.Bottom.Flatness)
rasStack[[1]]
plot(rasStack[[1]])

points(shp)


#
setwd("D:/IHF/IHF_Research/RandomForest/csv_files")
# ##read point data and convert them into SPDF
Pointcorrdinates = read.csv("70_100cm.csv")
#
# this transforms the data.frame to a SpatialPointsDataFrame (i.e. a points shapefile)
coordinates(Pointcorrdinates)= ~ lat+long
points<-Pointcorrdinates

##Extract raster value by points
point.cov <- raster::extract(rasStack,points)
head(point.cov)

#### look at structure of pnt. coord.
str(Pointcorrdinates)
crs(Pointcorrdinates)

####fresh data with lab/covaruate data w/ 69 pts
lab.covar<- cbind(points,point.cov)
head(lab.covar)

##end of data for soil and covariates going into model
lab.data <- lab.covar
head(lab.data)

# Lists to store results
rmse_list <- c()
ccc_list <- c()
sequence <- c(seq(10, 65, 5), 69) # Sample sizes
all_data <- rbind(lab.data@data)

for (aaa in 1:50) {  # Run 50 iterations
  set.seed(aaa)
  
  for (i in sequence) {  # For each sample size
    sample_index <- sample(seq_len(nrow(all_data)), size = i)
    data_sample <- all_data[sample_index,]
    
    rmse_values <- c()
    ccc_values <- c()
    
    # Set up k-fold cross-validation
    cv <- trainControl(method = "CV", number = 10) # 10-fold CV
    
    # Train Random Forest model with cross-validation
    rf.model <- caret::train(Sand ~ LSF + FPL + Aspect + TPI + SWI + DEM + EM + PCA + MRRTF + MRVBF, data = data_sample, method = "rf", trControl = cv)
    
    
    # Make predictions on the same data (cross-validated predictions)
    prediction <- predict(rf.model, dplyr::select(data_sample, -Sand))
    
    # Calculate RMSE
    test_rmse <- sqrt(mean((data_sample$Sand - prediction)^2))
    
    # Manually calculate CCC
    rho <- cor(data_sample$Sand, prediction)
    mu_x <- mean(data_sample$Sand)
    mu_y <- mean(prediction)
    sigma_x <- sd(data_sample$Sand)
    sigma_y <- sd(prediction)
    test_ccc <- (2 * rho * sigma_x * sigma_y) / (sigma_x^2 + sigma_y^2 + (mu_x - mu_y)^2)
    
    rmse_list <- c(rmse_list, test_rmse)  # Store RMSE for this sample size
    ccc_list <- c(ccc_list, test_ccc)  # Store CCC for this sample size
    
    # Print CCC value for this iteration and sample size
    print(paste("CCC value for iteration", aaa, "and sample size", i, ":", test_ccc))
  }
  
  print(paste("Iteration for aaa =", aaa, "completed"))
}

# Create the data frame
model_summary.sand70100 <- data.frame(sample_size = rep(sequence, times = 50), 
                                     ccc = unlist(ccc_list), 
                                     rmse = unlist(rmse_list))
# Fixed column name

# Export to CSV
write.csv(model_summary.sand70100, "70100updated.sand_ccc_rmse.csv", row.names = FALSE)

# Summarizing data
df.3 <- model_summary.clay010 %>%
  group_by(sample_size) %>%  # Fixed column name
  summarise(
    mean_ccc = mean(ccc),
    mean_rmse = mean(rmse),
    mean_mae = mean(mae),  # Added mean for MAE
    lci = t.test(ccc, conf.level = 0.9)$conf.int[1],
    uci = t.test(ccc, conf.level = 0.9)$conf.int[2]
  )


### this is to create 50 iterations per sample size of 10,15,20...65 at 10-40cm with mean values of confidene level at 90%
sequence <- c(seq(10, 65, 5), 69)
model_summary.0140.clay <- data.frame(sample_size = rep(c(seq(10, 65, 5), 69), each = 50), ccc = unlist(ccc_list), rmse = unlist(rmse_list))

sequence <- c(seq(10, 65, 5), 69)

model_summary.0140.clay <- data.frame(sample_size = rep(sequence, each = 50), ccc = unlist(ccc_list), rmse = unlist(rmse_list))

df.1111=model_summary.0140.clay %>%
  group_by(sample_size) %>%
  summarise(
    mean_ccc=mean(ccc),
    mean_rmse=mean(rmse),
    lci=t.test(ccc,conf.level = 0.9)$conf.int[1],
    uci=t.test(ccc,conf.level = 0.9)$conf.int[2])

df.1111

### this is to create 50 iterations per sample size of 10,15,20...65 at 40-70cm with mean values of confidene level at 90%
sequence <- c(seq(10, 65, 5), 69)
model_summary.4070.clay <- data.frame(sample_size = rep(c(seq(10, 65, 5), 69), each = 50), ccc = unlist(ccc_list), rmse = unlist(rmse_list))

sequence <- c(seq(10, 65, 5), 69)

model_summary.4070.clay <- data.frame(sample_size = rep(sequence, each = 50), ccc = unlist(ccc_list), rmse = unlist(rmse_list))

model_summary.4070.clay
df.2111=model_summary.4070.clay %>%
  group_by(sample_size) %>%
  summarise(
    mean_ccc=mean(ccc),
    mean_rmse=mean(rmse),
    lci=t.test(ccc,conf.level = 0.9)$conf.int[1],
    uci=t.test(ccc,conf.level = 0.9)$conf.int[2])

df.2111

### this is to create 50 iterations per sample size of 10,15,20...65 at 70-100cm with mean values of confidene level at 90%
sequence <- c(seq(10, 65, 5), 69)
model_summary.70100.clay <- data.frame(sample_size = rep(c(seq(10, 65, 5), 69), each = 50), ccc = unlist(ccc_list), rmse = unlist(rmse_list))

sequence <- c(seq(10, 65, 5), 69)

# Your model fitting loop here...

model_summary.70100.clay <- data.frame(sample_size = rep(sequence, each = 50), ccc = unlist(ccc_list), rmse = unlist(rmse_list))

df.3111 = model_summary.70100.clay %>%
  group_by(sample_size) %>%
  summarise(
    mean_ccc = mean(ccc),
    mean_rmse = mean(rmse),
    lci = t.test(ccc, conf.level = 0.9)$conf.int[1],
    uci = t.test(ccc, conf.level = 0.9)$conf.int[2])


df.3111


new1 <- c( 0.68)
new2 <- c(0.77)

# Perform the Mann Whitney U test
test_result <- wilcox.test(new1, new2)

# Print the test result
print(test_result)

# Now test_results is a list of test results for each pair of datasets


# Now test_results is a list of test results for each sample size




###Now time to plot the sample size vs CCC for 0-10cm#############################
pl1 <- ggplot(data = df)
pl1 <- pl1 + geom_bar(aes(x=sample_size, y=mean, fill = sample_size), stat="identity",fill="lightblue")
pl1 <- pl1 + geom_errorbar(aes(x=sample_size, ymin=lci, ymax= uci), width = 0.3, color ="red", size = 1.5)
pl1 <- pl1 + geom_text(aes(x=sample_size, y=lci, label = round(lci,2)), size= 3, vjust = 1.5)
pl1 <- pl1 + geom_text(aes(x=sample_size, y=uci, label = round(uci,2)), size= 3, vjust = -1)
pl1 <- pl1 + theme_classic()
pl1 <- pl1 + labs(title = "Mean with 90% confidence intervals for clay 0-10 cm",color="blue")
pl1 <- pl1 + labs(x= "Sample Size", y = "Mean CCC")



# Add this line to set the breaks for the x-axis
pl1 <- pl1 + scale_x_continuous(breaks = seq(10, 65, by = 5))

pl1


###Now time to plot the sample size vs CCC for 10-40 cm cm#############################
pl2 <- ggplot(data = df.1)
pl2<- pl2 + geom_bar(aes(x=sample_size, y=mean, fill = sample_size), stat="identity",fill="lightblue")
pl2 <- pl2 + geom_errorbar(aes(x=sample_size, ymin=lci, ymax= uci), width = 0.3, color ="red", size = 1.5)
pl2 <- pl2 + geom_text(aes(x=sample_size, y=lci, label = round(lci,2)), size= 3, vjust = 1.5)
pl2 <- pl2 + geom_text(aes(x=sample_size, y=uci, label = round(uci,2)), size= 3, vjust = -1)
pl2 <- pl2 + theme_classic()
pl2 <- pl2 + labs(title = "Mean with 90% confidence intervals for clay 10-40 cm",color="blue")
pl2 <- pl2 + labs(x= "Sample Size", y = "Mean CCC")

# Add this line to set the breaks for the x-axis
pl2 <- pl2 + scale_x_continuous(breaks = seq(10, 65, by = 5))

pl2

###Now time to plot the sample size vs CCC for 40-70 cm ############################
pl3 <- ggplot(data = df.2)
pl3<- pl3 + geom_bar(aes(x=sample_size, y=mean, fill = sample_size), stat="identity",fill="lightblue")
pl3 <- pl3 + geom_errorbar(aes(x=sample_size, ymin=lci, ymax= uci), width = 0.3, color ="red", size = 1.5)
pl3 <- pl3 + geom_text(aes(x=sample_size, y=lci, label = round(lci,2)), size= 3, vjust = 1.5)
pl3 <- pl3 + geom_text(aes(x=sample_size, y=uci, label = round(uci,2)), size= 3, vjust = -1)
pl3 <- pl3 + theme_classic()
pl3 <- pl3 + labs(title = "Mean with 90% confidence intervals for clay 40-70 cm",color="blue")
pl3 <- pl3 + labs(x= "Sample Size", y = "Mean CCC")

# Add this line to set the breaks for the x-axis
pl3 <- pl3 + scale_x_continuous(breaks = seq(10, 65, by = 5))

pl3

# Now time to plot the sample size vs CCC for 70-100 cm
pl4 <- ggplot(data = df.3111)
pl4 <- pl4 + geom_bar(aes(x = sample_size, y = mean_ccc, fill = sample_size), stat = "identity", fill = "lightblue")
pl4 <- pl4 + geom_errorbar(aes(x = sample_size, ymin = lci, ymax = uci), width = 0.3, color = "red", size = 1.5)
pl4 <- pl4 + geom_text(aes(x = sample_size, y = lci, label = round(lci, 2)), size = 3, vjust = 1.5)
pl4 <- pl4 + geom_text(aes(x = sample_size, y = uci, label = round(uci, 2)), size = 3, vjust = -1)
pl4 <- pl4 + theme_classic()
pl4 <- pl4 + labs(title = "Mean with 90% confidence intervals for clay 70-100 cm", color = "blue")
pl4 <- pl4 + labs(x = "Sample Size", y = "Mean CCC")

# Add this line to set the breaks for the x-axis to include 69
pl4 <- pl4 + scale_x_continuous(breaks = c(seq(10, 65, by = 5), 69))

pl4

length(ccc_list)
length(rmse_list)


### for basic boxplot without blue bars
library(ggplot2)
####box plot of sample size vs ccc/rmse 0-10cm
# Box plot for CCC
ggplot(model_summary.sand010, aes(x = factor(sample_size, levels = sequence), y = ccc)) +
  geom_boxplot() +
  labs(x = "Sample Size", y = "CCC", title = "0-10 cm")

# Box plot for RMSE
ggplot(model_summary.sand010, aes(x = factor(sample_size, levels = sequence), y = rmse)) +
  geom_boxplot() +
  labs(x = "Sample Size", y = "RMSE", title = "0-10 cm")

####box plot of sample size vs ccc/rmse 10-40cm
ggplot(model_summary.clay1040 , aes(x = factor(sample_size, levels = sequence), y = ccc)) +
  geom_boxplot() +
  labs(x = "Sample Size", y = "CCC", title = "10-40 cm")

ggplot(model_summary.clay1040, aes(x = factor(sample_size, levels = sequence), y = rmse)) +
  geom_boxplot() +
  labs(x = "Sample Size", y = "RMSE", title = "10-40 cm")

####box plot of sample size vs ccc/rmse 40-70cm
ggplot(model_summary.clay4070, aes(x = factor(sample_size, levels = sequence), y = ccc)) +
  geom_boxplot() +
  labs(x = "Sample Size", y = "CCC", title = "40-70 cm")

ggplot(model_summary.clay4070, aes(x = factor(sample_size, levels = sequence), y = rmse)) +
  geom_boxplot() +
  labs(x = "Sample Size", y = "RMSE", title = "40-70 cm")

####box plot of sample size vs ccc/rmse 70-100cm
ggplot(model_summary.clay70100, aes(x = factor(sample_size, levels = sequence), y = ccc)) +
  geom_boxplot() +
  labs(x = "Sample Size", y = "CCC", title = "Clay 70-100 cm")

ggplot(model_summary.clay70100, aes(x = factor(sample_size, levels = sequence), y = rmse)) +
  geom_boxplot() +
  labs(x = "Sample Size", y = "RMSE", title = "Clay 70-100 cm")

# Load the required package
library(gridExtra)

# Calculate the overall minimum and maximum
y_min <- min(
  min(model_summary.clay010$ccc, na.rm = TRUE),
  min(model_summary.clay1040$ccc, na.rm = TRUE),
  min(model_summary.clay4070$ccc, na.rm = TRUE),
  min(model_summary.clay70100$ccc, na.rm = TRUE)
)
y_max <- max(
  max(model_summary.clay010$ccc, na.rm = TRUE),
  max(model_summary.clay1040$ccc, na.rm = TRUE),
  max(model_summary.clay4070$ccc, na.rm = TRUE),
  max(model_summary.clay70100$ccc, na.rm = TRUE)
) + 0.15

# Create the four plots from the 4 boxplots created for ccc
plot1 <- ggplot(model_summary.clay010, aes(x = factor(sample_size), y = ccc)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(y_min, y_max)) +
  labs(x = "Sample Size", y = "CCC", title = "0-10 cm") +
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(size = 20),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18)) +
  annotate("text", 
           x = 1:13,  # Adjusted to cover 13 sample sizes
           y = y_max - c(0.7, 0.59, 0.29, 0.30, 0.25, 0.25, 0.27, 0.23, 0.23, 0.20, 0.17, 0.17, 0.17),  # example y values, adjust as needed
           label = c("a", "ab", "abc", "abc", "acd", "cde", "cdef", "defg", "efgh", "fgh", "gh", "h", "h"),  # example labels, adjust as needed
           size = 4)



# Repeat for the other plots...
plot2 <- ggplot(model_summary.clay1040, aes(x = factor(sample_size), y = ccc)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(y_min, y_max)) +
  labs(x = "Sample Size", y = "CCC", title = "10-40 cm") +
  theme(axis.title.x = element_text(size = 20),  # change the size as needed
        axis.title.y = element_text(size = 20),  # change the size as needed
        plot.title = element_text(size = 20),    # change the size as needed
        axis.text.x = element_text(size = 18),   # change the size as needed
        axis.text.y = element_text(size = 18)) +
  annotate("text",
           x = 1:13, 
           y = c(y_max - c(0.3, 0.20, 0.18, 0.15, 0.11, 0.1, 0.1, 0.1, 0.1, 0.1, 0.07, 0.08, 0.08)),  # Added a new y value
           label = c("a", "ab", "ac", "cd", "cd", "cd", "cd", "cd", "cd", "cd", "cd", "d", "cd"), size = 4)  # Added a new label


##
plot3 <- ggplot(model_summary.clay4070, aes(x = factor(sample_size), y = ccc)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(y_min, y_max)) +
  labs(x = "Sample Size", y = "CCC", title = "40-70 cm") +
  theme(axis.title.x = element_text(size = 20),  # change the size as needed
        axis.title.y = element_text(size = 20),  # change the size as needed
        plot.title = element_text(size = 20),    # change the size as needed
        axis.text.x = element_text(size = 18),   # change the size as needed
        axis.text.y = element_text(size = 18)) +  # change the size as needed
  annotate("text",
           x = 1:13, 
           y = c(y_max - c(0.10, 0.05, 0.04, 0.04, 0.09, 0.1, 0.1, 0.1, 0.1, 0.1, 0.07, 0.08, 0.08)),  # Added a new y value
           label = c("a", "ab", "abc", "acd", "cde", "def", "defg", "efg", "fgh", "gh", "hi", "ij", "j"), size = 4)  # Added a new label

##
plot4 <- ggplot(model_summary.clay70100, aes(x = factor(sample_size), y = ccc)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(y_min, y_max)) +
  labs(x = "Sample Size", y = "CCC", title = "70-100 cm") +
  theme(axis.title.x = element_text(size = 20),  # change the size as needed
        axis.title.y = element_text(size = 20),  # change the size as needed
        plot.title = element_text(size = 20),    # change the size as needed
        axis.text.x = element_text(size = 18),   # change the size as needed
        axis.text.y = element_text(size = 18)) +  # change the size as needed
  annotate("text",
           x = 1:13, 
           y = c(y_max - c(0.8, 0.6, 0.6, 0.27, 0.32, 0.52, 0.51, 0.52, 0.51, 0.6, 0.57, 0.58, 0.58)),  # Added a new y value
           label = c("a", "ab", "ab", "ac", "ac", "ac", "acd", "cd", "cd", "cd", "cd", "cd", "d"), size = 4)  # Added a new label

# Arrange the four plots into a 2x2 grid
grid.arrange(plot1, plot2, plot3, plot4, nrow = 2)

#better resolution 
# ... (your code to create the plots)

# Arrange the four plots into a 2x2 grid
plotsccc <- gridExtra::arrangeGrob(plot1, plot2, plot3, plot4, nrow = 2)
plotsccc
# Save the plot to a file with high resolution
ggsave("plots.png", plots, dpi = 500, width = 20, height = 15)


####ploting RMSE in a 2x2 matrix from the boxplots
# Calculate the overall minimum and maximum
y_min <- min(
  min(model_summary.clay010$rmse, na.rm = TRUE),
  min(model_summary.clay1040$rmse, na.rm = TRUE),
  min(model_summary.clay4070$rmse, na.rm = TRUE),
  min(model_summary.clay70100$rmse, na.rm = TRUE)
)
y_max <- max(
  max(model_summary.clay010$rmse, na.rm = TRUE),
  max(model_summary.clay1040$rmse, na.rm = TRUE),
  max(model_summary.clay4070$rmse, na.rm = TRUE),
  max(model_summary.clay70100$rmse, na.rm = TRUE)
) + 0.15

# Create the four plots from the 4 boxplots created for ccc
# Calculate the y values for the annotations

plot10 <- ggplot(model_summary.clay010, aes(x = factor(sample_size), y = rmse)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(y_min, y_max)) +
  labs(x = "Sample Size", y = "RMSE", title = "0-10 cm") +
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(size = 20),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18)) +
  annotate("text", 
           x = 1:13, 
           y = y_max - c(5.20, 3.90, 5.7, 7.6, 9.3, 5.7, 7.6, 10.9, 11, 11, 11, 11, 11),  # example y values
           label = c("a", "a", "a", "a", "ab", "abc", "abc", "bcd", "cde", "cdef", "def", "ef", "f"), 
           size = 5)


########
plot20 <- ggplot(model_summary.clay1040, aes(x = factor(sample_size), y = rmse)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(y_min, y_max)) +
  labs(x = "Sample Size", y = "RMSE", title = "10-40 cm") +
  theme(axis.title.x = element_text(size = 20),  # change the size as needed
        axis.title.y = element_text(size = 20),  # change the size as needed
        plot.title = element_text(size = 20),    # change the size as needed
        axis.text.x = element_text(size = 18),   # change the size as needed
        axis.text.y = element_text(size = 18)) +
  annotate("text", 
           x = 1:13, 
           y = y_max - c(2, 4.80, 5.40, 6.60, 7, 7.95, 8.7, 8.30, 8, 8, 10, 11.2, 13.5),  # example y values
           label = c("a", "ab", "abc", "abc", "abc", "acd", "cde", "cde", "def", "efg", "efg", "fg","g"), 
           size = 5)


####################

plot30 <- ggplot(model_summary.clay4070, aes(x = factor(sample_size), y = rmse)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(y_min, y_max)) +
  labs(x = "Sample Size", y = "RMSE", title = "40-70 cm") +
  theme(axis.title.x = element_text(size = 20),  # change the size as needed
        axis.title.y = element_text(size = 20),  # change the size as needed
        plot.title = element_text(size = 20),    # change the size as needed
        axis.text.x = element_text(size = 18),   # change the size as needed
        axis.text.y = element_text(size = 18)) +  # change the size as needed
  annotate("text", 
           x = 1:13, 
           y = y_max - c(10, 12.3, 12.3, 12.6, 12.6, 12.6, 12.6, 12.5 , 12.6 , 12.6 , 12.6, 12.5, 12.6),  # example y values
           label = c("a", "ab", "ab", "ab", "ab", "ab", "ab", "ab", "ab", "a", "ab", "a", "a"), 
           size = 5)

plot30
###################

plot40 <- ggplot(model_summary.clay70100, aes(x = factor(sample_size), y = rmse)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(y_min, y_max)) +
  labs(x = "Sample Size", y = "RMSE", title = "70-100 cm") +
  theme(axis.title.x = element_text(size = 20),  # change the size as needed
        axis.title.y = element_text(size = 20),  # change the size as needed
        plot.title = element_text(size = 20),    # change the size as needed
        axis.text.x = element_text(size = 18),   # change the size as needed
        axis.text.y = element_text(size = 18)) +  # change the size as needed
  annotate("text", 
           x = 1:13, 
           y = y_max - c(0.2, 0.0, 0.0, 2.2, 1.9, 1.6, 2 , 4.3, 3.7 , 3.8 , 3.5, 3.5, 3.5),    # example y values
           label = c("a", "ab", "ab", "ab", "ab", "ab", "ab", "ab", "ab", "ab", "ab", "a", "ab"), 
           size = 5)

plot40
# Arrange the four plots into a 2x2 grid
grid.arrange(plot10, plot20, plot30, plot40, nrow = 2)

#better resolution 
# ... (your code to create the plots)
plot10
plot20
plot30
plot40

# Arrange the four plots into a 2x2 grid
plots.1 <- gridExtra::arrangeGrob(plot10, plot20, plot30, plot40, nrow = 2)

# Save the plot to a file with high resolution
ggsave("rmse plots.png", plots.1, dpi = 500, width = 20, height = 15)


# Add mean values to the error bars
pl1 <- pl1 + geom_text(aes(x=sample_size, y=mean, label = round(mean,2)), size= 3, vjust = 0)

pl1 <- pl1 + theme_classic()
pl1 <- pl1 + labs(title = "Mean with 90% confidence intervals for clay 0-10 cm",color="blue")
pl1 <- pl1 + labs(x= "Sample Size", y = "Mean CCC")


# Create a data frame with your ccc_list and num_points
df <- data.frame(ccc = unlist(ccc_list), sample_size = unlist(num_points))

# Run a one-way ANOVA
res.aov <- aov(ccc ~ sample_size, data = df)

# Print the summary of the ANOVA
summary(res.aov)

# Create a data frame with your mse_list and num_points
df <- data.frame(mse = unlist(mse_list), sample_size = unlist(num_points))

# Plot MSE vs sample size
ggplot(df, aes(x = sample_size, y = mse)) +
  geom_point() + # This will create the scatterplot
  geom_line() + # This will connect the points with a line
  labs(x = "Sample Size", y = "Mean Squared Error", title = "MSE by Sample Size") +
  theme_minimal()

# Assume mse_design1 and mse_design2 are your MSE values for each design
test_result <- wilcox.test(mse_design1, mse_design2)

# Print the result
print(test_result)

mse_list <- c()  # Initialize the list

for (i in 1:100) {  # Loop over your models/data subsets
  # Here, perform your model fitting and prediction.
  # Then compute the MSE. This is just a placeholder example.
  mse <- i^2  # This should be your computed MSE
  
  # Append the MSE to the list
  mse_list <- c(mse_list, mse)
}

# Now you can create a dataframe with mse_list
df <- data.frame(mse = unlist(mse_list))

# Get the CCC values for two different sample sizes
ccc_sample_size_10 <- model_summary.70100.clay[model_summary.70100.clay$sample_size == 10, "ccc"]
ccc_sample_size_15 <- model_summary.4070.clay[model_summary.70100.clay$sample_size == 15, "ccc"]

# Conduct a Mann-Whitney U test
test_result <- wilcox.test(ccc_sample_size_10, ccc_sample_size_15)

# Print the p-value
print(test_result$p.value)

library(ggplot2)

ggplot(model_summary.70100.clay, aes(x = factor(sample_size), y = ccc)) +
  geom_boxplot() +
  labs(x = "Sample Size", y = "CCC", title = "CCC values for different sample sizes")

# Run the Wilcoxon test of predicted and actual values from the model 
wilcox_result <- wilcox.test(dt$RealValue, dt$PredictedValue)

# Print the p-value
print(wilcox_result$p.value)

library(Metrics)
mae_value <- mae(actual_values, predicted_values)
mae_value <- mae(dt$RealValue, dt$PredictedValue)
print(mae_value)






####average rmse values

# Create the plots
install.packages("gridExtra")
library(gridExtra)
plot1111 <- ggplot(model_summary.010.silt, aes(x = factor(sample_size), y = ccc)) +
  geom_boxplot() +
  labs(x = "Sample Size", y = "RMSE", title = "0-10 cm")

plot2111<- ggplot(model_summary.040.silt, aes(x = factor(sample_size), y = ccc)) +
  geom_boxplot() +
  labs(x = "Sample Size", y = "RMSE", title = "10-40 cm")

plot3111 <- ggplot(model_summary.4070.silt, aes(x = factor(sample_size), y = ccc)) +
  geom_boxplot() +
  labs(x = "Sample Size", y = "RMSE", title = "40-70 cm")

plot4111 <- ggplot(model_summary.70100.silt, aes(x = factor(sample_size), y = ccc)) +
  geom_boxplot() +
  labs(x = "Sample Size", y = "RMSE", title = "70-100 cm")

# Arrange the plots in a 2x2 grid with an overall title
grid.arrange(plot1111, plot2111, plot3111, plot4111, ncol = 2, 
             top = textGrob("Silt", gp = gpar(fontsize = 15, font = 1)))


# ggplot of raster importance. Extract variable importance from the random forest model
#var_importance <- importance(rf)

# Convert to data frame for ggplot2
#var_importance_df <- data.frame(Variable = rownames(var_importance), 
# Importance = var_importance[, "IncNodePurity"])

# Order the variables by importance
#var_importance_df <- var_importance_df[order(var_importance_df$Importance, decreasing = TRUE), ]

# Create the plot
# Create the plot with colored bars

#ggplot(var_importance_df, aes(x = reorder(Variable, Importance), y = Importance, fill = Variable)) +
# geom_bar(stat = "identity") +
#coord_flip() +
#labs(x = "Variable", y = "Importance", title = "Sand 70-100 cm Variable Importance") +
#scale_fill_viridis(discrete = TRUE) + # Use the viridis color palette
#theme_minimal() +
# theme(legend.position = "none") # Remove the legend


# Define seeds
seeds <- 1:50

# Prepare for the loop
Data = data.frame()

for(i in seeds) { 
  set.seed(i)
  rf=randomForest(lab.data@data[c("LSF", "FPL","Aspect","TPI","SWI",
                                  "DEM", "PCA_1", "EM","MRRTF","MRVBF")],
                  y=lab.data@data$Clay, mtry=10,importance=TRUE,
                  na.action = na.omit) 
  
  # Extract the variable importance
  data <- rf$importance
  # Convert to a data frame
  data <- data.frame(Variable = rownames(data), Importance = data[, "%IncMSE"])
  # Combine with the overall data
  Data <- rbind(Data, data)
}


# Boxplot for variable importance
ggplot(Data, aes(x = reorder(Variable, Importance), y = Importance, fill = Variable)) +
  geom_boxplot() +
  theme_light() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 11)
  ) +
  ggtitle("Clay 10-40 cm  cm Variable Importance") +
  xlab("") +
  coord_flip()

