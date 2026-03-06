rm(list = ls())
graphics.off()

ST <- readLines("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/data-raw/Synthetic_case/Rhone/PamHyr/_PAMHYR_/Rhone_PCH_Hautes_eaux/default-mage/net/Reach_004.ST")

# ST <- c(ST1, ST2)

# Skip the header (starts with '#')
lines <- ST[!grepl("^#", ST)] # remove comment lines

# Split each line by whitespace
split_lines <- strsplit(lines, "\\s+")
# Filter: keep lines with 3–4 elements and no 999.9990
split_lines <- Filter(function(x) {
  length(x) <= 4 && !any(x == "999.9990")
}, split_lines)
# Pad shorter lines with NAs
maxlen <- max(sapply(split_lines, length))
split_lines <- lapply(split_lines, function(x) {
  length(x) <- maxlen
  x
})

# Convert to data frame
df <- as.data.frame(do.call(rbind, split_lines), stringsAsFactors = FALSE)

# Convert numeric columns where possible
df[] <- lapply(df, type.convert, as.is = TRUE)

# Name columns (you can adjust)
colnames(df) <- c("X", "Y", "Z", "Label")[1:ncol(df)]

head(df)

write.table(df, file = "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/data-raw/Synthetic_case/Rhone/QGIS/XYZ_Rhone_last_reach.txt", row.names = FALSE, col.names = TRUE, sep = ";")


# ####################################################
# # After creating the interpolated ST


# ST <- readLines("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/data-raw/Synthetic_case/Rhone/PamHYR_synthetic/ST_created/rhone_reach_LAG_ANT.ST")

# # ST <- c(ST1, ST2)

# # Skip the header (starts with '#')
# lines <- ST[!grepl("^#", ST)] # remove comment lines

# # Split each line by whitespace
# split_lines <- strsplit(lines, "\\s+")
# # Filter: keep lines with 3–4 elements and no 999.9990
# split_lines <- Filter(function(x) {
#   length(x) <= 5 && !any(x == "999.9990")
# }, split_lines)
# # Pad shorter lines with NAs
# maxlen <- max(sapply(split_lines, length))
# split_lines <- lapply(split_lines, function(x) {
#   length(x) <- maxlen
#   x
# })

# # Convert to data frame
# df <- as.data.frame(do.call(rbind, split_lines), stringsAsFactors = FALSE)

# # Convert numeric columns where possible
# df[] <- lapply(df, type.convert, as.is = TRUE)

# if (all(is.na(df$V1))) {
#   df <- df[, -which(names(df) == "V1")]
# } else {
#   stop("something is wrong")
# }

# # Name columns (you can adjust)
# colnames(df) <- c("X", "Y", "Z", "Label")[1:ncol(df)]

# head(df)

# write.table(df, file = "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/data-raw/Synthetic_case/Rhone/PamHYR_synthetic/rhone_reach_LAG_ANT.csv", row.names = FALSE, col.names = TRUE, sep = ";")
