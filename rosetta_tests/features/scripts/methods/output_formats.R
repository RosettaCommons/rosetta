

##############################################
#     Output Formats for Feature Plots
#
#USAGE:
#
#In the setup script define a data.frame output_formats:
#
#   output_formats <- rbind(
#     <output_format_1>,
#     <output_format_2>,
#     ...)
#
#Then at at a feature plot script:
# 
#   save_plots()
#


output_web_raster <- list(
  id="output_web_raster",
  extensions = ".png",
  height=8,
  width=10.5,
  dpi=300,
  scale=1)

output_web_vector <- list(
  id="output_web_vector",
  extensions = ".svg",
  height=8,
  width=10.5,
  dpi=300,
  scale=.1)

output_slide_raster <- list(
  id="output_slide_raster",
  extension = ".png",
  height=6,
  width=8,
  dpi=300,
  scale=1)

output_slide_vector <- list(
  id="output_slide_vector",
  extension = ".svg",
  height=6,
  width=8,
  dpi=300,
  scale=1)

output_slide_pdf <- list(
  id="output_slide_pdf",
  extension = ".pdf",
  height=6,
  width=8,
  dpi=300,
  scale=1)

output_print_raster <- list(
  id="output_print_raster",
  extensions = ".png",
  height=8,
  width=10.5,
  dpi=300,
  scale=1)

output_print_vector <- list(
  id="output_print_vector",
  extensions = ".svg",
  height=8,
  width=10.5,
  dpi=300,
  scale=1)


output_print_pdf <- list(
  id="output_print_pdf",
  extensions = ".pdf",
  height=15,
  width=18,
  dpi=300,
  scale=1)

output_huge_raster <- list(
  id="output_huge_raster",
  extensions = ".png",
  height=15,
  width=18,
  dpi=300,
  scale=1)

output_huge_vector <- list(
  id="output_huge_vector",
  extensions = ".svg",
  height=15,
  width=18,
  dpi=300,
  scale=1)


output_huge_pdf <- list(
  id="output_huge_pdf",
  extensions = ".pdf",
  height=8,
  width=10.5,
  dpi=300,
  scale=1)



#INPUT:
# opt should be a list containing values like
#
# $output_web_vector
# [1] TRUE
#OUTPUT:
# matrix with output_formats suitable for passing to save_plots(...)
get_output_formats <- function(opt){
  output_formats <- NULL
  if(opt$output_web_raster) output_formats <- rbind(output_formats, output_web_raster)
  if(opt$output_web_vector) output_formats <- rbind(output_formats, output_web_vector)
  if(opt$output_slide_raster) output_formats <- rbind(output_formats, output_slide_raster)
  if(opt$output_slide_vector) output_formats <- rbind(output_formats, output_slide_vector)
  if(opt$output_slide_pdf) output_formats <- rbind(output_formats, output_slide_pdf)
  if(opt$output_print_raster) output_formats <- rbind(output_formats, output_print_raster)
  if(opt$output_print_vector) output_formats <- rbind(output_formats, output_print_vector)
  if(opt$output_print_pdf) output_formats <- rbind(output_formats, output_print_pdf)
  if(opt$output_huge_raster) output_formats <- rbind(output_formats, output_huge_raster)
  if(opt$output_huge_vector) output_formats <- rbind(output_formats, output_huge_vector)
  if(opt$output_huge_pdf) output_formats <- rbind(output_formats, output_huge_pdf)
  return(output_formats)
}

