##Retrieve sample intensity data
#Set the working directory to where intensity data csv files are found
setwd("C:/Users/Jordan/Documents/NTU/Y4S2/Data Analysis/Sample Stats")

#Create empty list to store sample file names
list_sample_file_names <- list()

#Loop through the folder to pick up sample names and store them in "list_sample_names"
for (sample_name in list.files()){
  list_sample_file_names <- append(list_sample_file_names, sample_name)
}

#Create empty list to store sample names
list_sample_data <- list()

#Load each sample csv into the environment
for (sample_file in list_sample_file_names){
  #Read the CSV file and store it in the list
  trimmed_sample_name <- strsplit(sample_file, split="_")[[1]][1]
  list_sample_data[[trimmed_sample_name]] <- read.csv(sample_file)
}

##Clean up sample intensity data
#Replace NaN values with 0 in sample files
for (sample_name in names(list_sample_data)){
  sample_df <- list_sample_data[[sample_name]]
  #Replace NaN values with 0
  sample_df[is.na(sample_df)] <- 0
  #Update the list element with modified data frame
  list_sample_data[[sample_name]] <- sample_df
}

##Processing of intensity data - normalize and calculate summative statistics
#Function for normalizing Intensity/Volume data
normalize <- function(Int_vol){
  (Int_vol - min(Int_vol, na.rm = TRUE)) / (max(Int_vol, na.rm = TRUE) - min(Int_vol, na.rm = TRUE))
}

#Add Normalized_Int_vol using min-max normalization and Int_vol columns
for (sample_name in names(list_sample_data)){
  sample_df <- list_sample_data[[sample_name]]
  #Calculate Intensity/Volume (skip if Intensity is 0)
  sample_df$Int_vol <- ifelse(sample_df$Intensity != 0, sample_df$Intensity/sample_df$Volume, 0)
  #Calculate Normalized_Int_vol
  sample_df$Normalized_Int_vol <- normalize(sample_df$Int_vol)
  #Update the data frame in the list
  list_sample_data[[sample_name]] <- sample_df
}

#Function to extract RegionAbbr and Normalized_Int_vol columns into a new dataframe
organize_samples <- function(sample_names, list_sample_data) {
  #Create an empty list to store dataframes for the samples
  sample_data <- list()

  #Loop through each sample name and extract the corresponding dataframe
  for (sample_name in sample_names) {
    selected_df <- subset(list_sample_data[[sample_name]], select = c("RegionAbbr", "Normalized_Int_vol"))
    #Rename Normalized_Int_vol column
    colnames(selected_df)[-1] <- paste(sample_name, "Normalized_Int_vol", sep = "_")
    sample_data[[sample_name]] <- selected_df
  }

  #Combine the dataframes into a single dataframe by merging on 'RegionAbbr'
  samples_df <- Reduce(function(x, y) merge(x, y, by = "RegionAbbr", all = TRUE), sample_data)

  return(samples_df)
}

#Specify obese and lean samples
obese_samples <- c("sample4974", "sample4985", "sample4986", "sample7073",
                    "sample9294", "sample9390" , "sample9272", "sample9269")
lean_samples <- c("sample4932", "sample4933","sample8293",
                  "sample8294", "sample8405", "sample8406")

#Organize obese and lean samples into dataframes
obese_samples_df <- organize_samples(obese_samples, list_sample_data)
lean_samples_df <- organize_samples(lean_samples, list_sample_data)

##Plot PCA of samples using all brain regions
#Load ggplot2
library(ggplot2)

#Load ggfortify
library(ggfortify)

#Load ggrepel
library(ggrepel)

#Create a dataframe that only selects the sample columns from lean_samples_df
PCA_lean <- lean_samples_df[,c(2:7)]

#Set row names as brain region abbreviations
rownames(PCA_lean) <- lean_samples_df[,1]

#Create a dataframe that only selects the sample columns from obese_samples_df
PCA_obese <- obese_samples_df[,c(2:9)]

#Combine lean and obese data frames into a single data frame
PCA_combined <- cbind(PCA_lean,PCA_obese)

#Transpose the combined data frame to have samples as rows and brain region abbreviations as columns
PCA_combined <- t(PCA_combined)

#Extract the first part of row names "sampleXXXX" and set them as new row names
for (row_number in 1:nrow(PCA_combined)) {
  row_name_parts <- strsplit(rownames(PCA_combined)[row_number], "_")[[1]][1]
  row_name <- sub("sample", "", row_name_parts)
  rownames(PCA_combined)[row_number] <- row_name
}

#Store sample names
PCA_Sample_names <- rownames(PCA_combined)

#Create group labels
group_labels <- c("Lean", "Lean", "Lean", "Lean", "Lean", "Lean", "Obese", "Obese", "Obese", "Obese", "Obese", "Obese", "Obese", "Obese")

#Combine sample names and group labels
PCA_group_labels <- data.frame(PCA_Sample_names, group_labels)

#Perform Principal Component Analysis (PCA) on PCA_combined
pca <- prcomp(PCA_combined)

#Set up the PNG device for saving the plot
png("C:/Users/Jordan/Documents/NTU/Y4S2/Data Analysis/Plots/PCA of samples using all brain regions.png", units="in", width=8, height=5, res=300)

#Plot PCA plot to see the spread of obese and lean samples
autoplot(pca, data = PCA_group_labels, label = FALSE) +
  geom_point(aes(colour = ifelse(PCA_group_labels[,2] == 'Obese', 'Obese', 'Lean')), size = 2) +
  geom_text_repel(aes(label = PCA_group_labels[,1], colour = ifelse(PCA_group_labels[,2] == 'Obese', 'Obese', 'Lean')),
                  size = 3, box.padding = 0.5, segment.alpha = 0) +
  scale_colour_manual(name = "Groups",
                      values = c("Obese" = "#F8766D", "Lean" = "#619CFF"))+
  theme(legend.key.size = unit(2, "lines"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 16, face = "bold"))

#End the PNG device
dev.off()

##Conduct Mann Whitney U test
#Create a dataframe to store the combined summary
Combined_summary <- data.frame(RegionAbbr = obese_samples_df$RegionAbbr)

#Perform Mann Whitney U test
Combined_summary$p_value <- rep(NA, nrow(Combined_summary))
for (row_number in 1:nrow(Combined_summary)) {
  #Perform unpaired two-sample Mann Whitney U test between corresponding rows in obese_samples_df and lean_samples_df
  test_results <- wilcox.test(
    x = as.numeric(obese_samples_df[row_number, c("sample4974_Normalized_Int_vol", "sample4985_Normalized_Int_vol", "sample4986_Normalized_Int_vol", "sample7073_Normalized_Int_vol",
                                                  "sample9294_Normalized_Int_vol", "sample9390_Normalized_Int_vol", "sample9272_Normalized_Int_vol", "sample9269_Normalized_Int_vol")]),
    y = as.numeric(lean_samples_df[row_number,  c("sample4932_Normalized_Int_vol", "sample4933_Normalized_Int_vol","sample8293_Normalized_Int_vol",
                                                  "sample8294_Normalized_Int_vol", "sample8405_Normalized_Int_vol", "sample8406_Normalized_Int_vol")]),
    paired = FALSE, exact = FALSE
  )

  #Store Mann Whitney U test result in Combined_summary data frame
  Combined_summary$U_value[row_number] <- test_results$statistic
  Combined_summary$p_value[row_number] <- test_results$p.value
}

#Replace NaN value with 0
Combined_summary[is.na(Combined_summary)] <- 0

#Store significant brain regions
significant_regions <- Combined_summary[Combined_summary$p_value < 0.05 & Combined_summary$p_value > 0, ]

#Save p-value significant criteria
significant_criteria <- Combined_summary$p_value < 0.05 & Combined_summary$p_value > 0

#Filter only significant brain regions
obese_samples_df <- obese_samples_df[,c("RegionAbbr", "sample4974_Normalized_Int_vol", "sample4985_Normalized_Int_vol", "sample4986_Normalized_Int_vol", "sample7073_Normalized_Int_vol",
                                       "sample9294_Normalized_Int_vol", "sample9390_Normalized_Int_vol", "sample9272_Normalized_Int_vol", "sample9269_Normalized_Int_vol")]
lean_samples_df <- lean_samples_df[, c("RegionAbbr","sample4932_Normalized_Int_vol", "sample4933_Normalized_Int_vol","sample8293_Normalized_Int_vol",
                                        "sample8294_Normalized_Int_vol", "sample8405_Normalized_Int_vol", "sample8406_Normalized_Int_vol")]

#Merge the two data frames based on the common column "BrainAbbr"
merged_df <- merge(lean_samples_df, obese_samples_df, by = "RegionAbbr")

##Mapping significant brain regions back to parent labels
#Load CSV file containing brain regions parent classification
parent_classification <- read.csv("C:/Users/Jordan/Documents/NTU/Y4S2/Data Analysis/Atlas_Parent_Classification.csv")

#Create empty vector to store parent brain classification
parent_label <- rep(NA,nrow(merged_df))

#Iterate over rows of merged_df and retrieve significant brain region abbreviation
for (row_number_merged in 1:nrow(merged_df)){
  brain_abbr <- merged_df[row_number_merged,1]
  #Iterate over rows of parent_classification
  for (row_number in 1:nrow(parent_classification)){
    abbreviation <- parent_classification[row_number,3]
    #Check for matching abbreviation
    if (brain_abbr == abbreviation){
      parent_label[row_number_merged] <- parent_classification[row_number,5]
      break  #Exit the inner loop once a match is found
    }
  }
}

#Add the parent_label as a new column at index 2 of merged_df
merged_df <- cbind(merged_df[, 1], parent_label, merged_df[, -1])

#Renaming column 1 - containing brain region abbreviations
colnames(merged_df)[1] <- "RegionAbbr"

#Renaming sample columns
colnames(merged_df)[3:ncol(merged_df)] <-c(lean_samples,obese_samples)

#Criteria for filtering rows that do not belong to the parent label "Fiber Tracts"
not_fiber_tracts <- merged_df$parent_label != "Fiber Tracts"

#Replacing NA values with TRUE
not_fiber_tracts[is.na(not_fiber_tracts)] <- TRUE

#Store significant brain regions
significant_regions <- Combined_summary[Combined_summary$p_value < 0.05 & Combined_summary$p_value > 0 &
                                          not_fiber_tracts, ]
#Save p-value significant criteria
significant_criteria <- Combined_summary$p_value < 0.05 & Combined_summary$p_value > 0

#Filter only brain regions that are not "Fiber Tracts
obese_samples_df <- obese_samples_df[not_fiber_tracts & Combined_summary$p_value >0,c("RegionAbbr", "sample4974_Normalized_Int_vol", "sample4985_Normalized_Int_vol", "sample4986_Normalized_Int_vol", "sample7073_Normalized_Int_vol",
                                                                                      "sample9294_Normalized_Int_vol", "sample9390_Normalized_Int_vol", "sample9272_Normalized_Int_vol", "sample9269_Normalized_Int_vol")]
lean_samples_df <- lean_samples_df[not_fiber_tracts & Combined_summary$p_value >0, c("RegionAbbr","sample4932_Normalized_Int_vol", "sample4933_Normalized_Int_vol","sample8293_Normalized_Int_vol",
                                                                                     "sample8294_Normalized_Int_vol", "sample8405_Normalized_Int_vol", "sample8406_Normalized_Int_vol")]

#Filer only brain regions that are not "Fiber Tracts and regions with 0 intensity
Combined_summary <- Combined_summary[not_fiber_tracts & Combined_summary$p_value >0,]

#Plot volcano plot
#Calculate median intensity for each brain region in lean samples
for (row_number in 1:nrow(lean_samples_df)) {
  #Extract intensity data for the current brain region
  intensity_data_for_1_brain_region <- as.numeric(lean_samples_df[row_number, c(2:ncol(lean_samples_df))])
  #Calculate median intensity for the current brain region
  median_intensity <- median(intensity_data_for_1_brain_region)
  #Store the median intensity in the Combined_summary data frame under the Lean_Median column
  Combined_summary$Lean_Median[row_number] <- median_intensity
}

#Calculate median intensity for each brain region in obese samples
for (row_number in 1:nrow(obese_samples_df)) {
  #Extract intensity data for the current brain region
  intensity_data_for_1_brain_region <- as.numeric(obese_samples_df[row_number, c(2:ncol(obese_samples_df))])
  #Calculate median intensity for the current brain region
  median_intensity <- median(intensity_data_for_1_brain_region)
  #Store the median intensity in the Combined_summary data frame under the Obese_Median column
  Combined_summary$Obese_Median[row_number] <- median_intensity
}

#Calculate fold change and log2 fold change
Combined_summary$Fold_change <- Combined_summary$Obese_Median / Combined_summary$Lean_Median
Combined_summary$log2_Fold_Change <- log2(Combined_summary$Fold_change)

#Create a data frame for volcano plot with RegionAbbr, log2 fold change, and p-value
volcano_df <- data.frame(
  RegionAbbr = Combined_summary$RegionAbbr,
  log2_Fold_Change = Combined_summary$log2_Fold_Change,
  p_value = Combined_summary$p_value
)
#Initialize the significance column with "Not significant"
volcano_df$significance <- "Not significant"

#Loop through each row in the volcano_df data frame
for (row_number in 1:nrow(volcano_df)) {
  #Check if p-value is less than 0.05
  if (volcano_df[row_number, "p_value"] < 0.05) {
    #Check if log2 fold change is greater than 0.10
    if (volcano_df[row_number, "log2_Fold_Change"] > 0.10) {
      #Assign "Higher Neuronal Activity" if both conditions are met
      volcano_df$significance[row_number] <- "Higher Neuronal Activity"
    } else {
      #Otherwise, assign "Lower Neuronal Activity"
      volcano_df$significance[row_number] <- "Lower Neuronal Activity"
    }
  }
}

#Set up the PNG device for saving the plot
png("C:/Users/Jordan/Documents/NTU/Y4S2/Data Analysis/Plots/Volcano Plot.png", units="in", width=8, height=5, res=300)

#Create a volcano plot using ggplot
ggplot(data = volcano_df, aes(x = log2_Fold_Change, y = -log10(p_value), color = significance)) +
  geom_vline(xintercept = c(-0.10, 0.10), col = "gray", linetype = 'dashed') +  #Add vertical lines
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +  #Add horizontal line for significance threshold
  geom_point(size = 2) +
  scale_color_manual(values = c("#F8766D", "grey"),
                     labels = c("Higher", "Not significant")) +
  geom_text_repel(data = subset(volcano_df, significance != "Not significant"),
                  aes(label = RegionAbbr), size = 3, vjust = -2, box.padding = 0.4, segment.alpha = 0) +
  labs(x = "log2(Fold Change)", y = "-log10(p-value)", color = "Neuronal Activity") +
  theme_minimal() +
  theme(legend.text = element_text(size = 14),
        legend.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))

#End the PNG device
dev.off()

##Plot box plot
#Subset merged_df based on conditions:
#1. Exclude rows where the "parent_label" column is equal to "Fiber Tracts"
#2. Include rows where the significance criteria are met
merged_df <- merged_df[merged_df$parent_label != "Fiber Tracts" & significant_criteria, ]

#Reshape the merged data into long format
library(tidyr)
long_df <- pivot_longer(merged_df, cols = -c(RegionAbbr, parent_label), names_to = "Sample", values_to = "Value")

#Initialize a vector to store group labels
group_labels <- character(length(long_df$Sample))

#Loop through vector containing sample names
for (sample_name in obese_samples) {
  #Check each row of data that contains the sample name
  match_indices <- grepl(sample_name, long_df$Sample)
  #Assign "obese" to matching samples
  group_labels[match_indices] <- "Obese"
}

#Assign "lean" to non-matching samples
group_labels[group_labels != "Obese"] <- "Lean"

#Add the group column to the data frame
long_df$Groups <- group_labels

#Subset the data to include only the first 10 brain regions
subset_df <- long_df #[long_df$RegionAbbr %in% unique(long_df$RegionAbbr)[1:22], ]

#Create a data frame of unique RegionAbbr and their corresponding p-values
unique_regions <- unique(long_df$RegionAbbr)
unique_p_values <- significant_regions[,2]

#Order the unique RegionAbbr based on their corresponding p-values
ordered_regions <- unique_regions[order(unique_p_values)]

#Reorder RegionAbbr based on increasing p-values
long_df$RegionAbbr <- factor(long_df$RegionAbbr, levels = ordered_regions)

#Define a function to create the plot for a specific parent label
plot_by_parent_label <- function(label) {
  gg <- ggplot(long_df[long_df$parent_label == label, ],
               aes(x = RegionAbbr, y = Value, fill = Groups, color = Groups)) +
    geom_boxplot(position = position_dodge(width = 0.75), width = 0.5, alpha = 0.3) +
    geom_point(position = position_dodge(width = 0.75)) +
    labs(title = paste("Parental Brain Region -", label),
         x = "Brain Regions", y = "Normalized Intensity", fill = "Groups", color = "Groups") +
    scale_fill_manual(values = c("Obese" = "#FF5733", "Lean" = "#619CFF")) +
    scale_color_manual(values = c("Obese" = "#FF5733", "Lean" = "#619CFF")) +
    theme_minimal() +
    theme(
      legend.key.size = unit(2, "lines"),
      legend.text = element_text(size = 18),
      legend.title = element_text(size = 22, face = "bold"),
      axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 18),
      axis.title.x = element_text(size = 22),
      axis.title.y = element_text(size = 22),
      plot.title = element_text(size = 24, face = "bold")
    )
  print(gg)
}

#Function to generate and save box plots for each parent label
generate_and_save_plots <- function(labels) {
  for (label in labels) {
    png(paste0("C:/Users/Jordan/Documents/NTU/Y4S2/Data Analysis/Plots/Boxplot ", label, ".png"),
        units = "in", width = 10, height = 8, res = 300)
    plot_by_parent_label(label)
    dev.off()
  }
}

#Call the function to generate and save plots
generate_and_save_plots(unique(merged_df$parent_label))

##Plot heatmap
#Load pheatmap
library(pheatmap)

#Load grid
library(grid)

#Create matrix data to plot heatmap
heatmap_data <- Combined_summary[Combined_summary$p_value < 0.05,c("Lean_Median", "Obese_Median")]

#Create a new column that combines parent labels and brain region abbreviations
combined_labels <- paste(merged_df[, 2], Combined_summary[Combined_summary$p_value < 0.05,"RegionAbbr"], sep = " ")

#Set row names of heatmap_data using combined_labels
rownames(heatmap_data) <- combined_labels

#Order row names based on the combined parent labels and brain region abbreviations,
#and then by the significance of the brain region
ordered_rownames <- combined_labels[order(merged_df[, 2], significant_regions$p_value)]

#Reorder rows of heatmap_data based on the ordered row names
heatmap_data <- heatmap_data[match(ordered_rownames, rownames(heatmap_data)), ]

#Add brain region location indexes to row names
rownames(heatmap_data) <- paste(1:nrow(heatmap_data), rownames(heatmap_data), sep = "-")

#Set column names of heatmap_data
colnames(heatmap_data) <- c("Lean", "Obese")

#Define the colors for the heatmap
high_values <- "#F8766D"
low_values <- "#619CFF"

#Generate a color palette with a range of colors between the extreme colors
color_palette <- colorRampPalette(c(low_values, high_values))

#Set up the PNG device for saving the plot
png("C:/Users/Jordan/Documents/NTU/Y4S2/Data Analysis/Plots/Heatmap.png", units = "in", width = 10, height = 8, res = 300)

#Add custom viewport settings for the heatmap and plot heatmap
setHook("grid.newpage", function() pushViewport(viewport(x=1, y=1, width=0.9, height=0.9, name="vp", just=c("right", "top"))), action="prepend")
pheatmap(heatmap_data, cluster_cols = FALSE, cluster_rows = FALSE, angle_col = 0,
         fontsize_row = 18, fontsize_col = 18,
          color = color_palette(100))
setHook("grid.newpage", NULL, "replace")

#Add labels for x-axis and y-axis
grid.text("Groups", x = 0.33, y = -0.02, gp = gpar(fontsize = 22))
grid.text("Brain Regions", x = -0.02, y = 0.51, rot = 90, gp = gpar(fontsize = 22))

#End the PNG device
dev.off()

##Plot pie chart
#Initialize a vector of ones, counting the number of brain regions mapped back to parental brain regions
Number_of_Significant_Regions <- rep(1, nrow(significant_regions))

#Create dataframe to store number of significant regions in parental brain regions
Count_significant_regions <- data.frame(parent_label=merged_df[,2],Number_of_Significant_Regions)

#Sum up the counts for each unique parent label
Count <- aggregate(Number_of_Significant_Regions ~ parent_label, data = Count_significant_regions, sum)

#Calculate total count of significant regions
total_count <- sum(Count$Number_of_Significant_Regions)

#Define the two colors
color1 <- "#F8766D"
color2 <- "#619CFF"

#Generate a palette with shades between the two colors
custom_palette <- colorRampPalette(c(color1, color2))(nrow(Count))

#Set up the PNG device for saving the plot
png("C:/Users/Jordan/Documents/NTU/Y4S2/Data Analysis/Plots/Piechart.png", units = "in", width = 10, height = 8, res = 300)

#Plot pie chart with custom color palette
ggplot(Count, aes(x = "", y =Number_of_Significant_Regions, fill = parent_label)) +
  geom_bar(width = 1, stat = "identity") +
  geom_text(aes(label = paste0(Number_of_Significant_Regions, "/", total_count)),
            position = position_stack(vjust = 0.5), size = 8) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = setNames(custom_palette, unique(merged_df[,2]))) +
  labs(fill = "Parental Brain Regions") +
  theme_void() +
  theme(
    legend.key.size = unit(2, "lines"),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 20, face = "bold")
  )

#End the PNG device
dev.off()

##Plot PCA of samples using only significant brain regions
#Filter significant brain regions
PCA_combined <- PCA_combined[,significant_criteria & not_fiber_tracts]

#Perform Principal Component Analysis (PCA) on PCA_combined
pca <- prcomp(PCA_combined)

#Set up the PNG device for saving the plot
png("C:/Users/Jordan/Documents/NTU/Y4S2/Data Analysis/Plots/PCA of samples using only significant brain regions.png", units = "in", width = 8, height = 5, res = 300)

#Plot PCA plot to see the spread of obese and lean samples
autoplot(pca, data = PCA_group_labels, label = FALSE) +
  geom_point(aes(colour = ifelse(PCA_group_labels[,2] == 'Obese', 'Obese', 'Lean')), size = 2) +
  geom_text_repel(aes(label = PCA_group_labels[,1], colour = ifelse(PCA_group_labels[,2] == 'Obese', 'Obese', 'Lean')),
                  size = 3, box.padding = 0.5, segment.alpha = 0) +
  scale_colour_manual(name = "Groups",
                      values = c("Obese" = "#F8766D", "Lean" = "#619CFF")) +
  theme(legend.key.size = unit(2, "lines"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 16, face = "bold"))

#End the PNG device
dev.off()

#Rearrange columns for Combined_summary for export into csv
Combined_summary <-Combined_summary[,c(1,3,2)]

#Export CSV file with Mann Whitney U's U value and p-value
write.csv(Combined_summary, "C:/Users/Jordan/Documents/NTU/Y4S2/Data Analysis/Mann Whitney U Statistics.csv", row.names = FALSE)

#Combining parent labels to significant regions
significant_regions <- cbind(merged_df[,2],significant_regions)

#Rearrange columns for significant_regions for export into csv
significant_regions <- significant_regions[,c(1,2,4,3)]

#Export CSV file with Mann Whitney U's U value and p-value (significant regions only)
write.csv(significant_regions, "C:/Users/Jordan/Documents/NTU/Y4S2/Data Analysis/Mann Whitney U Statistics (significant).csv", row.names = FALSE)
