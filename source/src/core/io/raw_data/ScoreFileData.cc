// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/io/raw_data/ScoreFileData.cc
///
/// @brief a way to write out just the score using the SilentFile stuff
/// @author Monica Berrondo, Luki Goldschmidt

// C++ Headers
#include <iostream>
#include <string>
#include <map>

// mini headers
#include <utility/io/ozstream.hh>

#include <basic/Tracer.hh>
#include <core/io/raw_data/Raw.fwd.hh>
#include <core/io/raw_data/ScoreFileData.hh>
#include <core/io/raw_data/ScoreStructText.hh>
#include <core/io/raw_data/ScoreStructJSON.hh>

#include <utility/file/file_sys_util.hh>

///Basic headers
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>


static basic::Tracer TR( "core.io.raw_data.ScoreFileData" );

namespace core {
namespace io {
namespace raw_data {

/// @brief write the given structure to the supplied filename.
bool ScoreFileData::write_struct(
	const RawStructOP s,
	std::map < std::string, core::Real > const & score_map,
	std::map < std::string, std::string > const & string_map
) {

	bool success = false;
	utility::io::ozstream output;

	runtime_assert( s != nullptr );

	if ( ! utility::file::file_exists( filename_ ) ) {
		output.open( filename_ );
		s->print_header( output, score_map, string_map );
	}

	output.open_append( filename_ );
	s->print_scores( output, score_map, string_map );
	success = output.good();

	return success;
}


bool ScoreFileData::write_scorefile(
	std::string tag,
	std::map < std::string, core::Real > const & score_map,
	std::map < std::string, std::string > const & string_map,
	bool use_json
) {

	RawStructOP outputter(nullptr);

	if ( use_json ) {
		// JSON
		outputter = utility::pointer::make_shared< ScoreStructJSON >( tag );
	} else {
		// Old plain-text format
		outputter = utility::pointer::make_shared< ScoreStructText >( tag );
	}

	if ( !outputter ) {
		return false;
	}

	return write_struct( outputter, score_map, string_map );
}

bool
ScoreFileData::write_scorefile(
	std::string tag,
	std::map < std::string, core::Real > const & score_map,
	std::map < std::string, std::string > const & string_map
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	bool use_json = false;
	std::string format = option[ out::file::scorefile_format ].value();
	if ( format == "text" ) {
	} else if ( format == "json" || format == "JSON" ) {
		use_json = true;
	} else {
		TR << "Invalid score file format specified: \"" << format << "\". No output generated!" << std::endl;
		return false;
	}

	return write_scorefile(tag, score_map, string_map, use_json);
}

} // namespace silent
} // namespace io
} // namespace core
