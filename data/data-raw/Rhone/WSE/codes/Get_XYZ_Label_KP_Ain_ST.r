rm(list = ls())
graphics.off()
cat("\014")


setwd("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/data-raw/Rhone/WSE/")

file_bathy_Ain <- "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/Case_studies/Rhone/Mage_models/Rhone_Ain_PCH/default-mage/net/"
bathymetry_Ain <- readLines(file.path(file_bathy_Ain, "Reach_014.ST"))
bathymetry_Ain <- bathymetry_Ain[-1]

# Remove headers and get xyz values

fields <- strsplit(trimws(bathymetry_Ain), "\\s+")
n_fields <- lengths(fields)

is_header <- n_fields >= 5

is_end <- sapply(fields, function(x) {
    length(x) >= 3 && all(x[1:3] == "999.9990")
})

is_xyz <- (n_fields %in% c(3, 4)) & !is_end

block_id <- cumsum(is_header) # each header starts a new block
block_id[is_end] <- NA # remove end marker lines

# keep only headers
header_lines <- fields[is_header]

# get first and before-last values
header_info <- t(sapply(header_lines, function(x) {
    first_val <- x[1]
    before_last <- x[length(x) - 1]
    c(profile_id = first_val, header_val = before_last)
}))

header_info <- as.data.frame(header_info, stringsAsFactors = FALSE)

xyz_list <- lapply(which(is_xyz), function(i) {
    nums <- as.numeric(fields[[i]][1:3])
    lbl <- if (length(fields[[i]]) >= 4) fields[[i]][4] else NA
    blk <- block_id[i] # which header/block this line belongs to
    c(nums, label = lbl, block = blk)
})

xyz_df <- as.data.frame(do.call(rbind, xyz_list), stringsAsFactors = FALSE)
colnames(xyz_df) <- c("X", "Y", "Z", "label", "block")
xyz_df$profile_id <- header_info$profile_id[as.numeric(xyz_df$block)]
xyz_df$header_val <- header_info$header_val[as.numeric(xyz_df$block)]

# optional: remove temporary block column
xyz_df$block <- NULL

# Convert numeric columns
xyz_df$X <- as.numeric(xyz_df$X)
xyz_df$Y <- as.numeric(xyz_df$Y)
xyz_df$Z <- as.numeric(xyz_df$Z)
xyz_df$profile_id <- as.numeric(xyz_df$profile_id)
xyz_df$header_val <- as.numeric(xyz_df$header_val)


write.table(xyz_df, file = "pre_processing/model_ST_PCH_PGA.txt", row.names = FALSE, quote = FALSE)

#############

bathymetry_Ain <- readLines(file.path(file_bathy_Ain, "Reach_005.ST"))
bathymetry_Ain <- bathymetry_Ain[-1]

# Remove headers and get xyz values

fields <- strsplit(trimws(bathymetry_Ain), "\\s+")
n_fields <- lengths(fields)

is_header <- n_fields >= 5

is_end <- sapply(fields, function(x) {
    length(x) >= 3 && all(x[1:3] == "999.9990")
})

is_xyz <- (n_fields %in% c(3, 4)) & !is_end

block_id <- cumsum(is_header) # each header starts a new block
block_id[is_end] <- NA # remove end marker lines

# keep only headers
header_lines <- fields[is_header]

# get first and before-last values
header_info <- t(sapply(header_lines, function(x) {
    first_val <- x[1]
    before_last <- x[length(x) - 1]
    c(profile_id = first_val, header_val = before_last)
}))

header_info <- as.data.frame(header_info, stringsAsFactors = FALSE)

xyz_list <- lapply(which(is_xyz), function(i) {
    nums <- as.numeric(fields[[i]][1:3])
    lbl <- if (length(fields[[i]]) >= 4) fields[[i]][4] else NA
    blk <- block_id[i] # which header/block this line belongs to
    c(nums, label = lbl, block = blk)
})

xyz_df <- as.data.frame(do.call(rbind, xyz_list), stringsAsFactors = FALSE)
colnames(xyz_df) <- c("X", "Y", "Z", "label", "block")
xyz_df$profile_id <- header_info$profile_id[as.numeric(xyz_df$block)]
xyz_df$header_val <- header_info$header_val[as.numeric(xyz_df$block)]

# optional: remove temporary block column
xyz_df$block <- NULL

# Convert numeric columns
xyz_df$X <- as.numeric(xyz_df$X)
xyz_df$Y <- as.numeric(xyz_df$Y)
xyz_df$Z <- as.numeric(xyz_df$Z)
xyz_df$profile_id <- as.numeric(xyz_df$profile_id)
xyz_df$header_val <- as.numeric(xyz_df$header_val)


write.table(xyz_df, file = "pre_processing/model_ST_PGA_CAIN.txt", row.names = FALSE, quote = FALSE)

#############

bathymetry_Ain <- readLines(file.path(file_bathy_Ain, "Reach_006.ST"))
bathymetry_Ain <- bathymetry_Ain[-1]

# Remove headers and get xyz values

fields <- strsplit(trimws(bathymetry_Ain), "\\s+")
n_fields <- lengths(fields)

is_header <- n_fields >= 5

is_end <- sapply(fields, function(x) {
    length(x) >= 3 && all(x[1:3] == "999.9990")
})

is_xyz <- (n_fields %in% c(3, 4)) & !is_end

block_id <- cumsum(is_header) # each header starts a new block
block_id[is_end] <- NA # remove end marker lines

# keep only headers
header_lines <- fields[is_header]

# get first and before-last values
header_info <- t(sapply(header_lines, function(x) {
    first_val <- x[1]
    before_last <- x[length(x) - 1]
    c(profile_id = first_val, header_val = before_last)
}))

header_info <- as.data.frame(header_info, stringsAsFactors = FALSE)

xyz_list <- lapply(which(is_xyz), function(i) {
    nums <- as.numeric(fields[[i]][1:3])
    lbl <- if (length(fields[[i]]) >= 4) fields[[i]][4] else NA
    blk <- block_id[i] # which header/block this line belongs to
    c(nums, label = lbl, block = blk)
})

xyz_df <- as.data.frame(do.call(rbind, xyz_list), stringsAsFactors = FALSE)
colnames(xyz_df) <- c("X", "Y", "Z", "label", "block")
xyz_df$profile_id <- header_info$profile_id[as.numeric(xyz_df$block)]
xyz_df$header_val <- header_info$header_val[as.numeric(xyz_df$block)]

# optional: remove temporary block column
xyz_df$block <- NULL

# Convert numeric columns
xyz_df$X <- as.numeric(xyz_df$X)
xyz_df$Y <- as.numeric(xyz_df$Y)
xyz_df$Z <- as.numeric(xyz_df$Z)
xyz_df$profile_id <- as.numeric(xyz_df$profile_id)
xyz_df$header_val <- as.numeric(xyz_df$header_val)


write.table(xyz_df, file = "pre_processing/model_ST_CAIN_Confluence.txt", row.names = FALSE, quote = FALSE)

#############
