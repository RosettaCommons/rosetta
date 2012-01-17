

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


output_web_raster <- data.frame(
  id="output_web_raster",
  extension = ".png",
  height=8,
  width=10.5,
  dpi=300,
  scale=1)

output_web_vector <- data.frame(
  id="output_web_vector",
  extension = ".svg",
  height=8,
  width=10.5,
  dpi=300,
  scale=.1)

output_slide_raster <- data.frame(
  id="output_slide_raster",
  extension = ".png",
  height=6,
  width=8,
  dpi=300,
  scale=1)

output_slide_vector <- data.frame(
  id="output_slide_vector",
  extension = ".svg",
  height=6,
  width=8,
  dpi=300,
  scale=1)

output_slide_pdf <- data.frame(
  id="output_slide_pdf",
  extension = ".pdf",
  height=6,
  width=8,
  dpi=300,
  scale=1)

output_print_raster <- data.frame(
  id="output_print_raster",
  extension = ".png",
  height=8,
  width=10.5,
  dpi=300,
  scale=1)

output_print_vector <- data.frame(
  id="output_print_vector",
  extension = ".svg",
  height=8,
  width=10.5,
  dpi=300,
  scale=1)


output_print_pdf <- data.frame(
  id="output_print_pdf",
  extension = ".pdf",
  height=8,
  width=10.5,
  dpi=300,
  scale=1)

output_huge_raster <- data.frame(
  id="output_huge_raster",
  extension = ".png",
  height=15,
  width=18,
  dpi=300,
  scale=1)

output_huge_vector <- data.frame(
  id="output_huge_vector",
  extension = ".svg",
  height=15,
  width=18,
  dpi=300,
  scale=1)


output_huge_pdf <- data.frame(
  id="output_huge_pdf",
  extension = ".pdf",
  height=15,
  width=18,
  dpi=300,
  scale=1)



#INPUT:
# opt should be a list containing values like
#
# $output_web_vector
# [1] TRUE
#OUTPUT:
# data.frame with output_formats suitable for passing to save_plots(...)
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

