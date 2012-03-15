
all_output_formats <- read.csv(
	file.path(base_dir, "scripts/methods/output_formats.csv"), header=T, sep="\t")

make_output_formats_options_list <- function(output_formats){
	alply(output_formats, 1, function(output_format){
		make_option(
			c(paste("--", output_format$id, sep="")),
			action="store_true",
			type="logical",
			default=output_format$use_by_default,
			dest=as.character(output_format$id),
			help=paste("Generate output plots using the", output_format$id, "format.  [Default \"%default\"]"))
	})
}    

get_output_formats <- function(options, output_formats){
	adply(output_formats, 1, function(output_format){
	       	if(options[as.character(output_format$id)] == TRUE){
			return(output_format)
		} else {
			return(data.frame())
		}
	})
}
