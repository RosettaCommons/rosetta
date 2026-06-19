// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/gemmi_util.cc
/// @brief  Utilities for working with Gemmi CIF file data
///
/// @author Rocco Moretti (rmorettiase@gmail.com)


// Unit headers
#include <utility/gemmi_util.hh>

#include <utility/exit.hh>

#include <numeric/types.hh>

#include <gemmi/cif.hpp>

#include <json.hpp>

namespace utility {

std::string
as_cif_value(json const & obj) {
	if ( obj.is_null() ) { return "?"; }
	if ( obj.is_boolean() ) { return obj.get<bool>() ? "YES" : "NO"; }
	if ( obj.is_string() ) { return gemmi::cif::quote(obj.get<std::string>()); }
	if ( obj.is_number() ) { return std::to_string(obj.get<double>()); }
	if ( obj.is_array() ) {
		std::string s = obj[0].get<std::string>();
		for ( numeric::Size ii(1); ii < obj.size(); ++ii ) {
			s += ' ';
			s += obj[ii].get<std::string>();
		}
		return gemmi::cif::quote(s);
	}
	utility_exit_with_message("In ready mmJSON, cannot covert entry `" + obj.dump() + "` to cif-style entry.");
}

/// @brief Read an mmJSON file into a Gemmi Document object.
/// @details Implementation based on Gemmi's read_mmjson_insitu()/fill_document_from_sajson()
/// but reimplemented to keep the Gemmi dependency header-only.
gemmi::cif::Document
gemmi_load_mmjson(std::string const & contents_of_file, std::string const & filename) {
	gemmi::cif::Document doc;
	doc.source = filename;

	json obj = json::parse(contents_of_file);

	if ( ! obj.is_object() ) { utility_exit_with_message("mmJSON file " + filename + " is not a JSON (dict) object at top level." ); }
	if ( obj.empty() ) { utility_exit_with_message("mmJSON file " + filename + " does not contain valid contents" ); }
	// We don't have tracers at utility level ...
	//if ( obj.size() != 1 ) { TR.Warning << "mmJSON file " << filename << " contains more than one structure -- only loading the first." << std::endl; }

	std::string const & block_name = obj.begin().key();
	if ( block_name.substr(0,5) != "data_" ) { utility_exit_with_message("mmJSON file " + filename + " entry " + block_name + " does not start with `data_`, as required."); }
	doc.blocks.emplace_back( block_name.substr(5) );
	std::vector<gemmi::cif::Item>& items = doc.blocks[0].items;

	auto & top = obj.begin().value();
	if ( ! top.is_object() ) { utility_exit_with_message("mmJSON file " + filename + " entry " + block_name + " is not a JSON (dict) object." ); }

	for ( auto it = top.begin(); it != top.end(); ++it) {
		std::string category_name = "_" + it.key() + ".";
		auto const & category = it.value();
		if ( !category.is_object() || category.empty() ) {
			utility_exit_with_message("mmJSON file " + filename + " entry " + it.key() + " is malformed." );
		}
		numeric::Size cif_cols = category.size();
		numeric::Size cif_rows = category.begin().value().size();
		if ( cif_rows > 1 ) {
			items.emplace_back(gemmi::cif::LoopArg{});
			gemmi::cif::Loop & loop = items.back().loop;
			loop.tags.reserve(cif_cols);
			loop.values.resize(cif_cols * cif_rows);
		}
		numeric::Size jj = 0;
		for ( auto it2 = category.begin(); it2 != category.end(); ++it2, ++jj ) {
			std::string tag = category_name + it2.key();
			auto const & arr = it2.value();
			if ( !arr.is_array() ) { utility_exit_with_message("mmJSON file " + filename + " entry " + tag + " is malformed"); }
			if ( arr.size() != cif_rows ) { utility_exit_with_message("mmJSON file " + filename + " entry " + tag + " does not have the expected number of entries."); }
			if ( cif_rows == 1 ) {
				items.emplace_back(tag, as_cif_value(arr[0]));
			} else if ( cif_rows > 1 ) {
				gemmi::cif::Loop & loop = items.back().loop;
				loop.tags.emplace_back(tag);
				for ( numeric::Size kk = 0; kk < cif_rows; ++kk ) {
					loop.values[jj + kk*cif_cols] = as_cif_value(arr[kk]);
				}
			}
		}
	}

	return doc;
}


} // namespace utility
