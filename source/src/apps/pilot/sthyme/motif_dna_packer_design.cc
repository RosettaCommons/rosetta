// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @file
/// @brief app for motif-biased protein-DNA design
/// @author sthyme
///

#include <devel/init.hh>

// Project Headers (protocols)
#include <protocols/dna/DnaInterfaceFinder.hh>
// AUTO-REMOVED #include <protocols/dna/DnaInterfacePacker.hh>
// AUTO-REMOVED #include <protocols/dna/PDBOutput.hh>
// AUTO-REMOVED #include <protocols/dna/RestrictDesignToProteinDNAInterface.hh>
// AUTO-REMOVED #include <protocols/dna/util.hh>
#include <protocols/motifs/motif_utils.hh>
#include <protocols/motifs/MotifDnaPacker.hh>

// Project Headers
// AUTO-REMOVED #include <core/chemical/ResidueType.hh>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>
#include <core/import_pose/import_pose.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
// AUTO-REMOVED #include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/dna/setup.hh>
// AUTO-REMOVED #include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>
static basic::Tracer TR("apps.pilot.motif_dna_packer_design");

// Utility Headers
// AUTO-REMOVED #include <utility/io/ozstream.hh>

// C++ Headers
#include <string>

// Option Key Includes
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/in.OptionKeys.gen.hh>
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
	core::scoring::ScoreFunctionOP scorefxn( core::scoring::getScoreFunction() );
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

		core::pose::PoseOP pose = new core::pose::Pose;
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
