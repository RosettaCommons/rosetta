// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief helper app to convert two-residue pdbs into a MotifLibrary text file
/// @author sthyme
///

#include <devel/init.hh>
#include <core/options/option.hh>

#include <protocols/motifs/Motif.hh>
#include <protocols/motifs/MotifLibrary.hh>
#include <protocols/motifs/motif_utils.hh>

// c++ headers
#include <utility/io/ozstream.hh>

int
main( int argc, char * argv [] )
{

	try {

	devel::init( argc, argv );
	protocols::motifs::MotifLibrary motifs( protocols::motifs::get_MotifLibrary_user() );
	protocols::motifs::MotifCOPs motifcops = motifs.library();

	std::string filename( "MotifLibrary.motifs" );
	utility::io::ozstream motif_output_file( filename );
	for( protocols::motifs::MotifCOPs::const_iterator motifcop_itr = motifcops.begin(), end_itr = motifcops.end();
			motifcop_itr != end_itr; ++motifcop_itr ) {
			protocols::motifs::MotifCOP motifcop( *motifcop_itr );
			motif_output_file << *motifcop;
	}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
