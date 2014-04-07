// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @brief
/// @author jk

#include <iostream>
#include <iomanip>

#include <core/pose/Pose.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/excn/Exceptions.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


using namespace core;
using namespace core::scoring;
using namespace basic::options;
using namespace basic::options::OptionKeys;

OPT_KEY( String, desired_sequence )
OPT_KEY( Boolean, repack_preserved )

static basic::Tracer TR( "apps.pilot.johnk_gapless_threading.main" );

/// General testing code
int
main( int argc, char * argv [] )
{
    try {
	NEW_OPT( desired_sequence, "the sequence to thread onto the current backbone", "" );
	NEW_OPT( repack_preserved, "if true, repack residues which match in both sequences", false );

	devel::init(argc, argv);

	TR << "Starting gapless threading" << std::endl;

	std::string const sequence_to_build = option[ desired_sequence ];

	// scoring function
	scoring::ScoreFunctionOP scorefxn( getScoreFunction() );

	// create pose from pdb
	pose::Pose pose;
	std::string const input_pdb_name ( basic::options::start_file() );
	core::import_pose::pose_from_pdb( pose, input_pdb_name );

	if ( sequence_to_build.length() != pose.total_residue() ) {
		TR << "Error!! Desired sequence does not have the same number of residues as input PDB" << std::endl;
		exit(1);
	}

	// Setup packer task
  pack::task::PackerTaskOP packer_task( pack::task::TaskFactory::create_packer_task( pose ));
	packer_task->set_bump_check( false );
	packer_task->initialize_from_command_line();
	packer_task->or_include_current( true );

	int num_preserved(0);
	int num_mut(0);

	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		chemical::AA const current_aa( pose.residue(ii).aa());
		char desired_aa = sequence_to_build[ii-1];
		if ( oneletter_code_from_aa(current_aa) == desired_aa ) {
			if ( option[ repack_preserved ] ) {
				packer_task->nonconst_residue_task(ii).restrict_to_repacking();
			} else {
				packer_task->nonconst_residue_task(ii).prevent_repacking();
			}
			++num_preserved;
		} else {
			utility::vector1< bool > allowed_aas( chemical::num_canonical_aas, false );
			chemical::AA aa( chemical::aa_from_oneletter_code( desired_aa ) );
			allowed_aas[ aa ] = true;
			packer_task->nonconst_residue_task(ii).restrict_absent_canonical_aas( allowed_aas );
			++num_mut;
		}
	}

	// Redesign
	pack::pack_rotamers( pose, *scorefxn, packer_task );

	// Report Scores
	(*scorefxn)(pose);
	pose.dump_scored_pdb( "requested_sequence.pdb", *scorefxn );

	TR << "Successfully finished replacing sequence." << std::endl;
	TR << "Made " << num_mut << " mutations, preserved sequence at " << num_preserved << " positions." << std::endl;

    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
        return -1;
    }
    return 0;
}



