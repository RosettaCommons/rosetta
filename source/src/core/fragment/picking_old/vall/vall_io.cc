// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragment/picking_old/vall/vall_io.hh
/// @brief  reading/writing of Vall libraries
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// unit headers
#include <core/fragment/picking_old/vall/vall_io.hh>

// type headers
#include <core/types.hh>

// package headers
#include <core/fragment/picking_old/vall/VallLibrary.hh>

// project headers
#include <basic/Tracer.hh>

// utility headers
#include <utility/exit.hh>
#include <utility/io/izstream.hh>

// C++ headers
#include <sstream>
#include <string>

#include <utility/vector1.hh>


#ifdef WIN32
#include <ctime>
#endif

namespace core {
namespace fragment {
namespace picking_old {
namespace vall {


// Tracer instance for this file
// Named after the original location of this code
static THREAD_LOCAL basic::Tracer TR( "core.fragment.picking_old.vall.vall_io" );


/// @brief load standard Vall library from file
/// @param[in] filename
/// @param[out] vall
/// @remarks Currently a difference in id (pdb name) or resi (residue sequence id)
///  will start a new stretch of lines defining a VallSection
void vall_library_from_file( std::string const & filename, VallLibrary & library, core::Size const preallocate ) {
	using core::Size;

	utility::io::izstream stream( filename );
	if ( !stream ) {
		utility_exit_with_message( "can't open file: " + filename );
	}

	// statistics
	Size n_lines = 0;
	Size prior_library_size = library.size();
	Size prior_n_residues = library.n_residues();

	TR << "Reading Vall library from " << filename << " ... " << std::endl;

	time_t time_start = time( NULL );

	std::string prior_id;
	Size prior_resi;
	VallResidue current_residue;
	VallSection section;

	// try to reduce number of allocation events if requested
	if ( preallocate > 0 ) {
		library.reserve( preallocate );
	}

	// parse Vall from file
	std::string line;
	while ( getline( stream, line ) ) {
		if ( line[0] != '#' ) {
			++n_lines;
			prior_id = current_residue.id();
			prior_resi = current_residue.resi();

			current_residue.fill_from_string( line );

			// check for start of new continuous stretch
			if ( ( current_residue.resi() != prior_resi + 1 ) || ( current_residue.id() != prior_id ) ) {
				if ( section.size() > 0 ) {
					library.add_section( section );
				}
				section.clear(); // new section
			}

			section.append_residue( current_residue );

			if ( n_lines % 100000 == 0 ) {
				TR << "   " << n_lines << std::endl;
			}
		}
	}
	if ( section.size() > 0 ) { // handle final set of entries
		library.add_section( section );
	}

	// free unused allocated memory
	if ( library.size() <= ( preallocate / 2 ) ) { // heuristic
		library.tighten_memory( false ); // entire library
	} else {
		library.tighten_memory(); // books only
	}

	time_t time_end = time( NULL );

	TR << "... done.  Read " << n_lines << " lines.  Time elapsed: " << ( time_end - time_start ) << " seconds." << std::endl;
	TR << "Prior library contained " << prior_library_size << " sections totaling " << prior_n_residues << " residues." << std::endl;
	TR << "Added " << library.size() << " sections to library totaling " << library.n_residues() << " residues." << std::endl;
}


} // vall
} // picking_old
} // fragment
} // core

