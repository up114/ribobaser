# data-raw/convert_csv_to_rda.R

csvs <- list.files("data", pattern = "[.]csv$", full.names = TRUE)
stopifnot(length(csvs) > 0)

# helper: if the first column looks like ids, promote to rownames
promote_rownames <- function(df) {
  first <- names(df)[1]
  # heuristics: non-numeric and unique
  if (!is.numeric(df[[first]]) && anyDuplicated(df[[first]]) == 0) {
    rn <- df[[first]]
    df[[first]] <- NULL
    rn <- as.character(rn)
    rownames(df) <- rn
  }
  df
}

for (f in csvs) {
  name <- tools::file_path_sans_ext(basename(f))
  obj_name <- make.names(name)

  df <- readr::read_csv(f, show_col_types = FALSE) |> as.data.frame(check.names = FALSE)
  df <- promote_rownames(df)

  assign(obj_name, df)                  # creates an object with that name
  save(list = obj_name,
       file = file.path("data", paste0(obj_name, ".rda")),
       compress = "xz")
}
