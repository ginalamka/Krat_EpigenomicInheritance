#!/bin/bash

# this will be used to prep the coverage files from methylseq for DMR analysis
# steps created with help from ChatGPT
# read the convo here: https://chatgpt.com/share/677eea4d-0dcc-8012-b66b-c48da30a54f9

module load bedtools
cd ~/kratepi/work

# Step 1: copy methyldackel files into a working directory in home
# cd /scratch/glamka/kratepi_GCF_02_nfcore_methylseq_rm46/bwameth/results/methyldackel/
# cp *_CpG.bedGraph ~/kratepi/work

# Step 2: create windows.bed by using the reference.fasta.fai after creating genome_sizes.txt

# cut -f1,2 GCF_019054845.1_ASM1905484v1_genomic.fa.fai > genome_sizes.txt

# bedtools makewindows -g genome_sizes.txt -w 500 > windows.bed
# bedtools makewindows -g genome_sizes.txt -w 500 -s 250 > overlapping_windows.bed # == 11380062 windows
# as of 1/30/25, also made overlapping with 
  # 1000 & 500 == 5690540 windows
  # 100 & 50 == 56896318 windows

# Step 3: subset the sorted files so that only a coverage between 5 and 20 are kept
# EXAMPLE for sorting: sort -k1,1 -k2,2n 16_GFL_1568_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov > sorted_16_GFL_1568_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov
# NOTE, bedGraph files are already sorted, so this is now unnecessary

# Loop through all files matching the pattern sorted_*.deduplicated.cov
for input_file in *.markdup.sorted_CpG.bedGraph; do

# Generate the output file name by adding "subset_" at the start
output_file="subset_${input_file}"

# Filter rows where col5 + col6 >= 5 and <= 20
awk '{
        sum = $5 + $6
        if (sum >= 5 && sum <= 20) print $0
    }' "$input_file" > "$output_file"

echo "Subset rows saved to $output_file"
done

# Step 4: aggregate methylation data (for each sample) over windows with bedtools
# EXAMPLE: bedtools map -a windows.bed -b methylation_sample1.bed -c 4 -o mean > aggregated_coverage_sample1.bed

# Loop through all methylation .bed files
for file in subset_*; do

# Generate the output filename by replacing "subset_" with "aggregated"
output_file="aggregated_${file#subset_}"

# Use bedtools map to calculate mean methylation values
bedtools map -a overlapping_windows.bed -b "$file" -c 4 -o mean > "$output_file"

# Print a message for tracking progress
echo "Processed $file into $output_file"
done

echo "All files aggregated successfully!"


# Step 5: combine all samples and use awk to grab columns of interest
# EXAMPLE: paste aggregated_coverage_03.bed aggregated_coverage_12.bed aggregated_coverage_19.bed > all_3samples_aggregated_coverage.bed
# awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8"\t"$12}' all_3samples_aggregated_coverage.bed > formatted_all_3samples_aggregated_coverage.bed

paste aggregated_* > combined_all_samples_aggregated_subset_overlappingwindows.bed
# Extract desired columns; the first 4 columns, then every 4th column after that
awk '{
    output = $1"\t"$2"\t"$3"\t"$4  # Start with columns 1-4
    for (i = 8; i <= NF; i += 4) { # Iterate over every 4th column starting at 8
        output = output"\t"$i
    }
    print output
}' combined_all_samples_aggregated_subset_overlappingwindows.bed > FINAL_formatted_combined_all_samples_aggregated_subset_overlappingwindows.bed

#### Prepare for analysis in R
R
setwd("/home/gfl0003/kratepi/work")
directory = getwd()
data <- read.table("FINAL_formatted_combined_all_samples_aggregated_subset_overlappingwindows.bed", header=FALSE)

# fix column labels
wind <- paste0("w_", seq_len(nrow(data))       # label windows
#rmcontig <- data[,4:51]                        # remove contig start/end columns      
#data <- cbind(wind, rmcontig)                  # merge
data <- cbind(wind, data)
colnames(data) <- c("windows", "contig", "start", "end", paste0("samp_", 1:45), paste0("samp_", 47:49))

# filter the data so that windows with no variation (0 or 100% in all) are removed
# REMEMBER TO CHANGE THE COLUMNS DEPENDING IF FILTERING ABOVE
filtered_data1 <- data[apply(data[, 5:ncol(data)], 1, function(row) all(row != 0)), ]
filtered_data2 <- filtered_data1[apply(filtered_data1[, 5:ncol(filtered_data1)], 1, function(row) all(row != 100)), ]
data <- filtered_data2

write.csv(
  filtered_data2,
  file = paste(directory, "/all_samples_methylationdata_withvariation_andmissingness.csv", sep = ""),
  sep = ",",
  row.names = FALSE,  # Removes the row numbers
  quote = FALSE       # Prevents quotes around numeric or character values
)

# Filter rows where X% or more of columns 5 onwards have a value (not ".")
library(dplyr)
filt_data_100 <- data[apply(data[, 5:ncol(data)], 1, function(row) all(row != ".")), ]
write.csv(filt_data_100, file = paste(directory, "/filt_data_100.csv", sep = ""), sep = ",", row.names = FALSE, quote = FALSE)

filt_data_90 <- filtered_data2 %>%
  filter(rowSums(across(5:ncol(filtered_data2), ~ .x != ".")) / (ncol(filtered_data2) - 4) >= 0.9)
write.csv(filt_data_90, file = paste(directory, "/filt_data_90.csv", sep = ""), sep = ",", row.names = FALSE, quote = FALSE)

filt_data_70 <- filtered_data2 %>%
  filter(rowSums(across(5:ncol(filtered_data2), ~ .x != ".")) / (ncol(filtered_data2) - 4) >= 0.7)
write.csv(filt_data_70, file = paste(directory, "/filt_data_70.csv", sep = ""), sep = ",", row.names = FALSE, quote = FALSE)

filt_data_50 <- filtered_data2 %>%
  filter(rowSums(across(5:ncol(filtered_data2), ~ .x != ".")) / (ncol(filtered_data2) - 4) >= 0.5)
write.csv(filt_data_50, file = paste(directory, "/filt_data_50.csv", sep = ""), sep = ",", row.names = FALSE, quote = FALSE)

save.image("/home/gfl0003/kratepi/work/methylation_analysis_filtering_data.RData")


###########################################################
# Old code from when using the Bismark output

# Step 1: grab the methylation coverage files and unzip them
##gunzip *.deduplicated.bismark.cov.gz

########

# Step 3: sort the coverage files by contig and position
# EXAMPLE: sort -k1,1 -k2,2n 16_GFL_1568_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov > sorted_16_GFL_1568_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov
# NOTE, bedGraph files are already sorted, so this is now unnecessary

# Loop through all files ending with ".deduplicated.bismark.cov"
##for file in *.deduplicated.bismark.cov; do
##  # Generate the output filename by adding "sorted_" prefix
##  sorted_file="sorted_${file}"

# Sort the file by contig (column 1) and position (column 2)
##  sort -k1,1 -k2,2n "$file" > "$sorted_file"

# Print a message for tracking progress
##  echo "Sorted $file and saved as $sorted_file"
##done

##echo "All files sorted successfully!"
