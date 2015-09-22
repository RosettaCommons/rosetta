// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief app for motif-biased protein-DNA design
/// @author sthyme


#include <devel/init.hh>

// Project Headers (protocols)
#include <protocols/dna/DnaInterfaceFinder.hh>
#include <protocols/motifs/motif_utils.hh>
#include <protocols/motifs/MotifDnaPacker.hh>

// Project Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/dna/setup.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>
static THREAD_LOCAL basic::Tracer TR( "apps.pilot.motif_dna_packer_design" );

// Utility Headers

// C++ Headers
#include <string>

// Option Key Includes
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/dna.OptionKeys.gen.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/string_util.hh>

#include <utility/vector0.hh>


////////////////////////////////////////////////////////////////////////////////

void
motif_dna_packer_design()
{

	utility::vector1< std::string > pdb_files( basic::options::start_files() );

	// Set up scoring
	core::scoring::ScoreFunctionOP scorefxn( core::scoring::get_score_function() );
	bool minimize( false );
	if ( basic::options::option[basic::options::OptionKeys::dna::design::minimize].user() ) minimize = basic::options::option[ basic::options::OptionKeys::dna::design::minimize ]();

	for ( utility::vector1< std::string >::const_iterator pdb_file( pdb_files.begin() );
			pdb_file != pdb_files.end(); ++pdb_file ) {
		std::string pdb_name( *pdb_file );

		TR << "Working on file: " << pdb_name << std::endl;
		std::string pdb_prefix( utility::string_split( utility::string_split( pdb_name, '/' ).back(), '.' ).front() );

		if ( basic::options::option[ basic::options::OptionKeys::out::prefix ].user() ) {
			pdb_prefix = basic::options::option[ basic::options::OptionKeys::out::prefix ]();
		}

		core::pose::PoseOP pose( new core::pose::Pose );
		core::import_pose::pose_from_pdb( *pose, pdb_name );
		protocols::motifs::make_dna_mutations( *pose );
		// This sets BasePartner in pose cacheable data
		core::scoring::dna::set_base_partner( *pose ); // For some reason I am resetting this after making a DNA mutation, I think it's a bug fix

		protocols::motifs::MotifDnaPacker motif_dna_packer( scorefxn, minimize, pdb_prefix );
		motif_dna_packer.apply( *pose );
	}

	std::cout << "SUCCESSFUL COMPLETION" << std::endl;

}

////////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{

	try {

		devel::init( argc, argv );
		motif_dna_packer_design();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
