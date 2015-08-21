// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/io/raw_data/RawStruct.cc
///
/// @brief Struct base class
/// @author James Thompson, Monica Berrondo

// C++ Headers
#include <string>
#include <map>

// mini headers
#include <core/io/raw_data/RawStruct.hh>

#include <ObjexxFCL/format.hh>

#include <utility/vector1.hh>

#include <basic/Tracer.hh>

namespace core {
namespace io {
namespace raw_data {

using namespace ObjexxFCL::format;

static int precision (3);

RawStruct::~RawStruct() {
}

// Print the energy line/header into silent files
// You must supply a score_map (from the CachedData in the pose)
// which is customized in the protocol and returned to the job distributor
// (or to your main function)
void
RawStruct::print_header(
	std::ostream& out,
	std::map < std::string, core::Real > const & score_map,
	std::map < std::string, std::string > const & string_map,
	bool print_sequence // = true
) const {
	using namespace ObjexxFCL;

	std::map< std::string, Real >::const_iterator pair;
	std::map< std::string, std::string>::const_iterator string_pair;
	Size width (8);

	if ( print_sequence ) out << "SEQUENCE: " << sequence_ << std::endl;
	out << "SCORE:";

	// show score first
	pair = score_map.find( "total_score" );
	if ( pair != score_map.end() ) {
		if ( pair->first.length() > 8 ) width = pair->first.length();
		out << ' ' << A( width, pair->first );
	}
	pair = score_map.find( "score" );
	if ( pair != score_map.end() ) {
		if ( pair->first.length() > 8 ) width = pair->first.length();
		out << ' ' << A( width, pair->first );
	}
	pair = score_map.find( "rms" );
	if ( pair != score_map.end() ) {
		if ( pair->first.length() > 8 ) width = pair->first.length();
		out << ' ' << A( width, pair->first );
	}

	for ( pair=score_map.begin(); pair!=score_map.end(); ++pair ) {
		if ( pair->first == "total_score" || pair->first == "score" || pair->first == "rms" ) continue;
		if ( pair->first.length() > 8 ) width = pair->first.length();
		out << ' ' << A( width, pair->first );
	}
	for ( string_pair=string_map.begin(); string_pair!=string_map.end(); ++string_pair ) {
		width = std::max( Size( std::max(string_pair->first.length(), string_pair->second.length() ) ), Size(8) );
		out << ' ' << A( width, string_pair->first );
	}
	out << " description " << std::endl;
}

void
RawStruct::print_scores(
	std::ostream& out,
	std::map < std::string, core::Real > const & score_map,
	std::map < std::string, std::string > const & string_map
) const {
	std::map< std::string, Real >::const_iterator pair;
	std::map< std::string, std::string>::const_iterator string_pair;
	Size width (8);

	out << "SCORE:";

	// show score first
	pair = score_map.find( "total_score" );
	if ( pair != score_map.end() ) {
		if ( pair->first.length() > 8 ) width = pair->first.length();
		out << ' ' << F( width, precision, pair->second );
	}
	pair = score_map.find( "score" );
	if ( pair != score_map.end() ) {
		if ( pair->first.length() > 8 ) width = pair->first.length();
		out << ' ' << F( width, precision, pair->second );
	}
	pair = score_map.find( "rms" );
	if ( pair != score_map.end() ) {
		if ( pair->first.length() > 8 ) width = pair->first.length();
		out << ' ' << F( width, precision, pair->second );
	}

	for ( pair=score_map.begin(); pair!=score_map.end(); ++pair ) {
		if ( pair->first == "total_score" || pair->first == "score" || pair->first == "rms" ) continue;
		if ( pair->first.length() > 8 ) width = pair->first.length();
		out << ' ' << F( width, precision, pair->second );
	}
	for ( string_pair=string_map.begin(); string_pair!=string_map.end(); ++string_pair ) {
		width = std::max( Size( std::max(string_pair->first.length(), string_pair->second.length() ) ), Size(8) );
		out << ' ' << A( width, string_pair->second );
	}
	out << " " << decoy_tag_ << std::endl;
}

void
RawStruct::print_conformation( std::ostream& out ) const {
	out << "Don't know how to print_conformation from RawStruct! Use a derived class!" << decoy_tag_ << std::endl;
}

Real RawStruct::get_debug_rmsd() {
	basic::T("core.io.silent.RawStruct")
		<< "Warning: calling get_debug_rmsd() method of RawStruct!";
	return 0;
}


} // namespace silent
} // namespace io
} // namespace core
