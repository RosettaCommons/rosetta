// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
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

// mini headers
#include <basic/Tracer.hh>

#include <utility/io/izstream.hh>
// AUTO-REMOVED #include <utility/io/ozstream.hh>
#include <utility/file/gzip_util.hh>

// option key includes

#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

//Auto Headers
#include <utility/io/mpistream.hh>


namespace core {
namespace io {
namespace silent {

static basic::Tracer tr("core.io.silent");

void
gzip() {
	using basic::options::option;
	using namespace basic::options::OptionKeys;
	// gzip output file if desired
	std::string const basename = option[ out::file::silent ]();
	if ( option[ out::silent_gz ]() ) {
		utility::vector1< std::string > file_list;
		file_list.push_back( basename );

		for ( utility::vector1< std::string >::const_iterator
				fn = file_list.begin(), end = file_list.end(); fn != end; ++fn
		) {
			utility::io::izstream in_stream( *fn );
			if ( in_stream ) {
				in_stream.close();
				in_stream.clear();
				tr.Info << "GZIP SILENT FILE: " << *fn << std::endl;
				utility::file::gzip( *fn, true );
			} // file exists
		} // loop over each file
	} // if ( option[ out::silent_gz ]() )
} // gzip

} // namespace silent
} // namespace io
} // namespace core
