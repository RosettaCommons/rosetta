#! /usr/bin/Rscript --vanilla
# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.



initialize_analysis_manager_db <- function(
	analysis_manager_db_fname) {

  sele <- "
CREATE TABLE IF NOT EXISTS features_reporters (
	features_reporter_type_name TEXT,
	PRIMARY KEY(features_reporter_type_name));

CREATE TABLE IF NOT EXISTS sample_sources (
	sample_source_id TEXT,
	filename TEXT,
	description TEXT,
	PRIMARY KEY(sample_source_id));

CREATE TABLE IF NOT EXISTS sample_source_features_reporters (
	sample_source_id TEXT,
	features_reporter_type_name TEXT,
	FOREIGN KEY(sample_source_id)
		REFERENCES sample_sources(sample_source_id)
		DEFERRABLE INITIALLY DEFERRED,
	FOREIGN KEY(features_reporter_type_name)
		REFERENCES features_reporters(features_reporter_type_name)
		DEFERRABLE INITIALLY DEFERRED);

CREATE TABLE IF NOT EXISTS features_analyses (
	features_analysis_id TEXT,
	filename TEXT,
	author TEXT,
	brief_description TEXT,
	long_description TEXT,
	script TEXT,
	PRIMARY KEY(features_analysis_id));

CREATE TABLE IF NOT EXISTS features_analysis_keywords (
	features_analysis_id TEXT,
	keyword TEXT,
	FOREIGN KEY(features_analysis_id)
		REFERENCES features_analyses(features_analysis_id)
		DEFERRABLE INITIALLY DEFERRED);

CREATE TABLE IF NOT EXISTS features_analysis_features_reporter_dependencies (
	features_analysis_id TEXT,
	features_reporter_type_name TEXT,
	FOREIGN KEY(features_analysis_id)
		REFERENCES feature_analyses(features_analysis_id)
		DEFERRABLE INITIALLY DEFERRED,
	FOREIGN KEY(features_reporter_type_name)
		REFERENCES features_reporters(features_reporter_type_name)
		DEFERRABLE INITIALLY DEFERRED);

CREATE TABLE IF NOT EXISTS features_analysis_plot_formats (
	format_id TEXT,
	extension TEXT,
	height REAL,
	width REAL,
	dpi REAL,
	scale REAL,
	PRIMARY KEY(format_id));

CREATE TABLE IF NOT EXISTS features_analysis_plots (
	plot_id TEXT,
	features_analysis_id TEXT,
	date_code TEXT,
	filename TEXT,
	format_id TEXT,
	PRIMARY KEY(plot_id, date_code, filename, format_id),
	FOREIGN KEY(features_analysis_id)
		REFERENCES feature_analyses(features_analysis_id)
		DEFERRABLE INITIALLY DEFERRED,
	FOREIGN KEY(format_id)
		REFERENCES features_analysis_plot_formats(format_id)
		DEFERRABLE INITIALLY DEFERRED);

CREATE TABLE IF NOT EXISTS features_analysis_plot_sample_sources (
	plot_id TEXT,
	sample_source_rank INTEGER,
	sample_source_id TEXT,
	PRIMARY KEY(plot_id, sample_source_id),
	UNIQUE(plot_id, sample_source_rank),
	FOREIGN KEY(plot_id)
		REFERENCES features_analysis_plots(plot_id)
		DEFERRABLE INITIALLY DEFERRED,
	FOREIGN KEY(sample_source_id)
		REFERENCES sample_sources(sample_source_id)
		DEFERRABLE INITIALLY DEFERRED);

CREATE TABLE IF NOT EXISTS features_analysis_plot_data(
	plot_id TEXT,
	plot_data BLOB,
	PRIMARY KEY(plot_id),
	FOREIGN KEY(plot_id)
		REFERENCES features_analysis_plot(plot_id)
		DEFERRABLE INITIALLY DEFERRED);
"

  con <- dbConnect(engine, as.character(analysis_manager_db_fname))
	sele_split <- paste(strsplit(sele, ";\\W*", perl=TRUE)[[1]], ";", sep="")
	l_ply(sele_split, function(sele) dbGetQuery(con, sele))
  con
}


add_features_analysis_to_analysis_manager <- function(
	analysis_manager_con,
	features_analysis) {

	sql <- "INSERT OR IGNORE INTO features_analyses VALUES (?,?,?,?,?,?);"
	bind <- data.frame(
			as.character(features_analysis@id),
			as.character(features_analysis@filename),
			as.character(features_analysis@author),
			as.character(features_analysis@brief_description),
			as.character(features_analysis@long_description),
			paste(deparse(features_analysis@run), collapse=""))
	dbGetPreparedQuery(analysis_manager_con, sql, bind)

	if(!is.na(features_analysis@keywords)){
		sql <- "INSERT OR IGNORE INTO features_analysis_keywords VALUES (?,?);"
		bind <- data.frame(
			features_analysis_id=features_analysis@id,
			keyword=features_analysis@keywords)
		dbGetPreparedQuery(analysis_manager_con, sql, bind)
	}
}


add_fetures_analysis_plot_formats_to_analysis_manger <- function(
	analysis_manager_con,
	output_formats) {

	sql <- "INSERT OR IGNORE INTO features_analysis_plot_formts VALUES (?,?,?,?,?,?);"
	dbGetPreparedQuery(analysis_manager_con, sql, output_formats);

}


add_features_analysis_plot <- function(
	features_analysis,
	plot_id,
	sample_sources,
	date_code,
	filename,
	plot_format) {

	sql <- "INSERT OR IGNORE INTO features_analysis_plots VALUES (?,?,?,?,?);"
	bind <- data.frame(
		plot_id=plot_id,
		features_analysis_id=features_analysis@id,
		date_code=date_code,
		filename,
		plot_format$id)
	dbGetPreparedQuery(analysis_manager_con, sql, bind)

	sql <- "INSERT OR IGNORE INTO features_analysis_plot_sample_sources VALUES (?,?,?);"
	for(i in 1:nrow(sample_sources)){
		bind <- data.frame(
			plot_id=plot_id,
			sample_source_rank=i,
			as.character(sample_sources[i,"sample_source"]))
		dbGetPreparedQuery(analysis_manager_con, sql, bind)
	}
}
