// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/antibody/grafting/util.cc
/// @brief Helpers for antibody grafting code
/// @author Rocco Moretti (rmorettiase@gmail.com)
/// @author Sergey Lyskov, Jeliazko Jeliazkov

#include <protocols/antibody/grafting/util.hh>

#ifdef __ANTIBODY_GRAFTING__
#include <regex>
#include <string>
#endif

namespace protocols {
namespace antibody {
namespace grafting {


using std::string;


// It's only usable if the regex works.
bool antibody_grafting_usable() {
#ifndef __ANTIBODY_GRAFTING__
	return false;
#else
	string s("is the regex code working?");
	std::regex reg("regex");
	return std::regex_search(s,reg); // Problematic GCC versions return false for everything
#endif
}

} // grafting
} // antibody
} // protocols


#ifdef __ANTIBODY_GRAFTING__

#include <protocols/antibody/grafting/exception.hh>

#include <basic/Tracer.hh>

#include <utility/vector0.hh>
#include <utility/string_util.hh>
#include <utility/stream_util.hh>

#include <map>
#include <fstream>

static THREAD_LOCAL basic::Tracer TR("protocols.antibody.grafting.util");

namespace protocols {
namespace antibody {
namespace grafting {

/// @details Parse plain text file with following format, Replace some of the charactes in legend
/// <legend-prefix> field-1 field-2 field-3 ...
/// f11 f12 f13 ...
/// f21 f22 f23 ...
/// @return utility::vector0< std::map<string field, string value> >
/// @trows _AE_scs_failed_ on unexpeceted formating
utility::vector0< std::map<string, string> >
parse_plain_text_with_columns(
	string file_name,
	string legend_prefix/*="# "*/,
	char legend_separator/*=' '*/,
	string data_prefix/*=""*/,
	char data_separator/*=' '*/
)
{
	utility::vector0< std::map<string, string> > result;
	utility::vector0<string> legend;

	std::ifstream f(file_name);

	string line;
	while( std::getline(f, line) ) {
		if( utility::startswith(line, legend_prefix) ) {
			if( legend_separator == ' ') legend = utility::split_whitespace( line.substr( legend_prefix.size() ) );
			else legend = utility::string_split( line.substr( legend_prefix.size() ), legend_separator);

			for(auto & l : legend) {
				l = utility::strip(l, ' ');
				l = utility::replace_in(l, ". ", ".");
				l = utility::replace_in(l, " ", "-");
			}

			TR.Trace << "Found legend: " << line << std::endl << "Legend: " << legend << std::endl;
		}
		else {
			if( utility::startswith(line, data_prefix) ) {
				if( !legend.size() ) throw _AE_scs_failed_("File: " + file_name + " is missing legend!");
				else {
					//TR.Trace << "Got line:" << line << std::endl;

					utility::vector0< string > parts;

					//auto parts( utility::string_split(line, '\t') );
					if( data_separator == ' ') parts = utility::split_whitespace(line);
					else parts = utility::string_split(line, data_separator);

					//TR.Trace << "Line parts.size: " << parts.size() << " legend size:" << legend.size() << std::endl;

					if( parts.size() == legend.size() ) {
						std::map<string, string> fields;  int i=0;
						for(auto e : parts) fields[ legend[i++] ] = e;
						result.push_back(fields);
					} else throw _AE_scs_failed_("Number of fileds does not match legend!\nFile: "+file_name+"\nLine: "+line);
				}
			}
		}
	}

	return result;
}

} // namespace grafting
} // namespace antibody
} // namespace protocols

#endif // __ANTIBODY_GRAFTING__
