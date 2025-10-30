
rm(list = ls()); gc()

# Get list of all CSV files in the folder
selectfolder = "csv_fraction_supertype_5per"
csv_files <- list.files(paste0("~/Documents/sea-ad/data/", selectfolder), 
                        pattern = "*.csv", full.names = TRUE)



# Initialize list to store all attractor results
all_MILists <- list()

# Loop through each CSV file
for (i in seq_along(csv_files)) {
  cat("Processing file", i, "of", length(csv_files), ":", basename(csv_files[i]), "\n")
  
  # Read the data
  data <- read.csv(csv_files[i], row.names = 1, check.names = FALSE)
  
  metadata <- data[, 1:19]
  
  data <- data[, -c(1:19)]
  
  data <- as.matrix(t(data))
  
  library("dplyr")
  
  data <- rbind(data, 
                "percent.6e10.positive.area" = as.numeric(metadata$`percent 6e10 positive area`))
  
  library("cafr")
  
  
  # Find the mi 
  MI_list <- getAllMIWz(data, data["percent.6e10.positive.area", ], negateMI = T, sorting = T)
  # Store the attractor.final in the list with the filename as the key
  file_name <- tools::file_path_sans_ext(basename(csv_files[i]))
  all_MILists[[file_name]] <- MI_list
  
  cat("Completed processing:", file_name, "\n")
}

# Print summary
cat("Processed", length(all_MILists), "files total\n")
cat("Names of processed files:\n")
print(names(all_MILists))
# save(all_MILists, file = paste0("/Users/lingyi/Documents/sea-ad/models/attractors/MI_full/MIs_", selectfolder, "_Full_6e10.RData"))

