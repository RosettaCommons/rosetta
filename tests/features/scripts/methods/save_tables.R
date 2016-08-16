# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


date_code <- function(d=NA){
	# reference http://www.r-cookbook.com/node/17
	if(is.na(d)) d <- Sys.Date()
	pattern <- '20([[:digit:]]{2})-([[:digit:]]{2})-([[:digit:]]{2})'
	paste(
		sub(pattern, '\\1', d), sub(pattern, '\\2', d), sub(pattern, '\\3', d),
		sep="")
}


# Save a data.frame as a table. For each output format,
# generate a table and put in the output directory
save_tables <- function(
	features_analysis,
	table,
	table_id,
	sample_sources,
	output_dir,
	output_formats,
	table_title=NULL,
	caption=caption,
  quote_strings=F,
	...
) {
	extra_args <- list(...)
	tryCatch(table_id, error=function(e){
		stop(paste(
			"ERROR: Unable to save the table because ",
			"the 'table_id' is not specified.\n", e, sep=""))
	})

	tryCatch(features_analysis, error=function(e){
		stop(paste(
			"ERROR: Unable to save the table '", table_id,"' ",
			"because the specified 'features_analysis' is not valid.\n",
			e, sep=""))
	})

	tryCatch(sample_sources, error=function(e){
		stop(paste(
			"ERROR: Unable to save the table '", table_id, "' ",
			"because the specified 'sample_sources' is not valid.\n",
			e, sep=""))
	})

	if(nrow(sample_sources)==0){
		stop(paste(
			"ERROR: Unable to save the table '", table_id, "' ",
			"because no sample_sources were specified.\n", e, sep=""))
	}

	tryCatch(output_dir, error=function(e){
		stop(paste(
			"ERROR: Unable to save the table '", table_id, "' ",
			"because the specified 'output_dir' ",
			"is not a valid variable.\n",
			e, sep=""))
	})

	tryCatch(output_formats, error=function(e){
		stop(paste(
			"ERROR: Unable to save the table '", table_id, "' ",
			"because the 'output_formats' parameter is not valid.\n",
			e, sep=""))
	})
	table_formats <- output_formats[output_formats$type == "table",]

	if(nrow(table_formats)==0){
		stop(paste(
			"ERROR: Unable to save the table '", table_id, "' ",
			"because no output formats were specified.", sep=""))
	}

	a_ply(table_formats, 1, function(fmt){
		full_output_dir <- file.path(output_dir, features_analysis@id, fmt$id)
		if(!file.exists(full_output_dir)){
			dir.create(full_output_dir, recursive=TRUE)
		}
		date <- date_code()
		fname <- paste(table_id, date, sep="_")
		full_path <- file.path(full_output_dir, paste(fname, fmt$extension, sep=""))
		extension <- substring(fmt$extension, 2)
		cat("Saving Table with extension ", extension, ": ", full_path, sep="")
		timing <- system.time({
			tryCatch({
				if(as.character(fmt$id) == "output_latex_sideways_table"){
					cat("\\usepackage{booktabs}\n\\usepackage{rotating}\n", file=full_path)
					print(xtable(
						table, caption=gsub("\n", "\\\\\n", caption), ...),
						file=full_path,
						type="latex",
						append=TRUE,
						booktabs=TRUE,
						floating.environment="sidewaystable*",
						...)
				} else if(as.character(fmt$id) == "output_latex_table"){
					cat("\\usepackage{booktabs}\n", file=full_path)
					print(
						xtable(table, caption=gsub("\n", "\\\\\n", caption), ...),
						file=full_path,
						type="latex",
						append=TRUE,
						booktabs=TRUE,
						...)
				} else if(as.character(fmt$id) == "output_html"){
					cat(table_css_header(), file=full_path)
					print(
						xtable(table, gsub("\n", "<br>\n", caption), ...),
						file=full_path,
						type=extension,
						append=TRUE,
						html.table.attributes="border=0",
						...)
				} else if(as.character(fmt$id) == "output_csv"){
					write.csv(
						x=table,
						file=full_path,
						row.names=FALSE,
            quote=quote_strings)
				} else if(as.character(fmt$id) == "output_tab_delimited_table") {
          write.table(
            table,
            full_path,
            sep="\t",
            row.names=F,
            quote=quote_strings)
				}
			}, error=function(e){
				cat("\n")
				cat(paste(
					"ERROR: Generating and saving the table:\n",
					e, sep=""))
			})
		})
		cat(" ... ", as.character(round(timing[3],2)), "s\n", sep="")
	})
}


table_css_header <- function()
	"
<head>
<style type=\"text/css\">
{
	font-family: \"Helvetica\", \"Lucida Sans Unicode\", \"Lucida Grande\", Sans-Serif;
	font-size: 10px;
	background: #fff;
	margin: 45px;
	width: 480px;
	border-collapse: collapse;
	text-align: left;
}


th {
	font-family: \"Helvetica\", \"Lucida Sans Unicode\", \"Lucida Grande\", Sans-Serif;
	font-size: 12px;
	font-weight: normal;
	color: #036;
	padding-top:10px;
	padding-bottom:2px;
	padding-right:8px;
	padding-left:8px;
	border-bottom: 2px solid #6678b1;
}


td {
	font-size: 12px;
	border-bottom: 1px solid #ccc;
	color: #336;
	padding: 4px 6px;
}

tbody tr:hover td {
	color: #009;
}
</style>
<head>
"
