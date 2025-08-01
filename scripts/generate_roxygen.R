library(stringr)

generate_roxygen_for_function <- function(filepath, fun_name) {
  file_content <- readLines(filepath)

  # Locate function line
  pattern <- paste0("^", fun_name, "\\s*<-\\s*function")
  func_line_index <- grep(pattern, file_content)

  if (length(func_line_index) == 0) {
    stop("Function not found in file.")
  }

  # Load function into environment
  env <- new.env()
  source(filepath, local = env)
  fun_obj <- env[[fun_name]]

  if (!is.function(fun_obj)) {
    stop(paste(fun_name, "is not a valid function."))
  }

  # Extract arguments
  args <- names(formals(fun_obj))

  # Create Roxygen block
  roxygen_lines <- c(
    paste0("#' Description of the function `", fun_name, "`"),
    "#'",
    paste0("#' @param ", args, " Description of the parameter"),
    "#'",
    "#' @return Result",
    "#'",
    "#' @examples",
    paste0("#' ", fun_name, "(", paste(rep("...", length(args)), collapse = ", "), ")"),
    "#'",
    "#' @export"
  )

  # Insert block above function
  new_content <- append(file_content, roxygen_lines, after = func_line_index[1] - 1)

  # Overwrite file
  writeLines(new_content, filepath)

  message("✅ Roxygen added for function: ", fun_name)
}

# ---- Get command line arguments ----
args <- commandArgs(trailingOnly = TRUE)
filepath <- args[1]
# Read all the lines of the file
file_lines <- readLines(filepath)

# Use regex to detect all function definitions
fun_lines <- grep("^[a-zA-Z0-9_.]+\\s*<-\\s*function", file_lines, value = TRUE)

if (length(fun_lines) == 0) {
  stop("No function found in the file.")
}

# Loop over each detected function
for (fun_line in fun_lines) {
  # Extract function name before "<-"
  fun_name <- sub("\\s*<-.*", "", fun_line)
  fun_name <- trimws(fun_name)

  cat("Generating Roxygen for function:", fun_name, "\n")

  # Call your function to write the roxygen comment
  generate_roxygen_for_function(filepath, fun_name)
}
