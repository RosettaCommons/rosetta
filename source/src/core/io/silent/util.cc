// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/io/silent/util.cc
/// @brief utility functions for silent-file classes.
/// @author James Thompson

// C++ Headers
#include <string>
#include <map>

// mini headers
#include <basic/Tracer.hh>
#include <core/io/silent/SilentFileData.hh>

#include <utility/io/izstream.hh>
#include <utility/file/gzip_util.hh>
#include <utility/file/file_sys_util.hh>

// option key includes

#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <utility/vector1.hh>


namespace core {
namespace io {
namespace silent {

static THREAD_LOCAL basic::Tracer tr( "core.io.silent" );

void
gzip() {
	using basic::options::option;
	using namespace basic::options::OptionKeys;
	// gzip output file if desired
	std::string const basename = option[ out::file::silent ]();
	if ( ! option[ out::silent_gz ]() )  return;

	utility::vector1< std::string > file_list;
	file_list.push_back( basename );

	for ( utility::vector1< std::string >::const_iterator
			fn = file_list.begin(), end = file_list.end(); fn != end; ++fn
			) {
		utility::io::izstream in_stream( *fn );
		if ( ! in_stream ) continue;

		in_stream.close();
		in_stream.clear();
		tr.Info << "GZIP SILENT FILE: " << *fn << std::endl;
		utility::file::gzip( *fn, true );
	} // loop over each file
} // gzip


/////////////////////////////////////////////////////////////////
std::map< std::string, bool >
initialize_tag_is_done( std::string const & silent_file ){

	std::map< std::string, bool > tag_is_done;
	utility::vector1< std::string > tags_done;
	SilentFileData silent_file_data;

	if ( utility::file::file_exists( silent_file ) ) {
		tags_done = silent_file_data.read_tags_fast( silent_file );
		for ( utility::vector1< std::string >::const_iterator iter = tags_done.begin(), end = tags_done.end(); iter != end; ++iter ) {
			std::cout << "Already done: " << *iter << std::endl;
			tag_is_done[ *iter ] = true;
		}
	}

	return tag_is_done;
}


} // namespace silent
} // namespace io
} // namespace core
