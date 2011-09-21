

parse_rosetta_database_path <- function(opt) {

  if(is.null(opt$database)) {
    opt$database <- Sys.getenv("ROSETTA3_DB")
    if(opt$database == ""){
      stop("Unable to locate path to the rosetta_database. Either pass in the path via the --database flag, or set the ROSETTA3_DB evironment variable.")
    }
  }
  if(!file.exists(opt$database)){
    stop(paste("The rosetta_database path '", opt$database, "' does not exist."))
  }
  opt
}
