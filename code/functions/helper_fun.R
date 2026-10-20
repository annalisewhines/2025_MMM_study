################################################################################
###############---------------- HELPER FUNCTIONS ----------------###############
################################################################################

### Helper functions used in subsequent functions and code

## Population standard deviation -----------------------------------------------
# @ x: A numerical vector of values

sdP <- function(x){
  sqrt(mean(x^2) - (mean(x))^2) 
}

## Standard error --------------------------------------------------------------
# @ x: A numerical vector of values

se <- function(x) var(x) / sqrt(length(x))

## Mode ------------------------------------------------------------------------
# @ x: A numerical vector of values

getmode <- function(x) {
  uniqX <- unique(x)
  uniqX[which.max(tabulate(match(x, uniqX)))]
  
}

## Task data standardization ---------------------------------------------------
#@ x: A numerical vector of values
#@ ref: A single numerical values to be used as the reference
#@ levels: Possible values the data to be standardized can take

standardization <- function(x, ref = 0, levels = 1:4){
  (x - ref) / sdP(levels)
}

## Shuffle subject IDs ---------------------------------------------------------
#@ x: Character string

id_shuffle <- function(x) { return(x) }
#new_id <- unlist(strsplit(x, ""))
#new_id <- new_id[1:10]
#set.seed(1)
#new_id <- sample(new_id, length(new_id), FALSE)
#new_id <- paste(new_id, collapse = "")
#new_id


## Scale variables to 0-1 ---------------------------------------------------------
#@ x: Variable vector

rescale <- function(x){
  (x-min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE)-min(x, na.rm = TRUE))
}


## Normalize negatively skewed data ---------------------------------------------------------
#@ x: Variable vector

norm_neg_skew <- function(x){
  (1 - sqrt(1 - x))
}


## Removal of orphan prolificID's ---------------------------------------------------------
json_parsing <- function(input_file, output_file) {
  # Read file line-by-line
  lines <- readLines(input_file, warn = FALSE)
  
  keep <- logical(length(lines))
  
  i <- 1
  while (i <= length(lines)) {
    line <- lines[i]
    
    # If this line contains prolificID
    if (grepl('"prolificID"', line)) {
      
      # Check if there is a next line and it contains "start"
      if (i < length(lines) && grepl('"start"', lines[i + 1])) {
        keep[i] <- TRUE
        keep[i + 1] <- TRUE
        i <- i + 2
      } else {
        # Drop this prolificID (no following start)
        i <- i + 1
      }
      
    } else {
      # Non-prolificID lines are kept by default
      keep[i] <- TRUE
      i <- i + 1
    }
  }
  
  # Write cleaned file
  writeLines(lines[keep], output_file)
}

# How to call the function:
# json_parsing("pill_game_raw.json", "pill_game_clean.json")