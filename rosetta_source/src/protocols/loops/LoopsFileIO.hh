// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loops/LoopsFileIO.hh
/// @brief This class exists to handle the reading and writing of loops files.
/// @author Brian D. Weitzner

#ifndef INCLUDED_protocols_loops_LoopsFileIO_HH
#define INCLUDED_protocols_loops_LoopsFileIO_HH

// Unit header
#include <protocols/loops/LoopsFileIO.fwd.hh>

// Package headers
#include <protocols/loops/Loop.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C/C++ headers
#include <fstream>
#include <string>

namespace protocols {
namespace loops {

///////////////////////////////////////////////////////////////////////////
// a list of loops
class LoopsFileIO : public utility::pointer::ReferenceCount {

public:

	//constructor
	LoopsFileIO();

		//copy constructor
	LoopsFileIO( const LoopsFileIO & src );


	// assignment operator
	LoopsFileIO & operator =( LoopsFileIO const & src );

		// destructor
		~LoopsFileIO();

	friend std::ostream & operator<<( std::ostream & os, const LoopsFileIO & loops_file_io );


	SerializedLoopList read_loop_file( std::string filename );

	SerializedLoopList use_custom_legacy_file_format(
		std::istream & is,
		std::string filename,
		bool strict_looprelax_checks,
		std::string token
	);

private:

	void read_stream_to_END(
		std::istream & is,
		std::string filename /*for error msg */,
		bool strict_looprelax_checks = true,
		std::string token = "LOOP" );

private:
	SerializedLoopList loops_;

}; // LoopsFileIO

} //namespace loops
} //namespace protocols

#endif //INCLUDED_protocols_loops_LoopsFileIO_HH
