// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief
/// @author jk

#include <iostream>
#include <iomanip>

#include <devel/init.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID_Map.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/TwelveANeighborGraph.hh>
#include <core/scoring/Energies.hh>
#include <protocols/moves/SuperimposeMover.hh>

#include <basic/Tracer.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/after_opts.hh>

#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/rms_util.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>


using namespace core;
using namespace basic::options;
using namespace core::scoring;
using namespace core::optimization;
using namespace basic::options::OptionKeys;

static thread_local basic::Tracer TR( "apps.pilot.chet_trp_to_gly.main" );


OPT_KEY( Integer, rms_resnum )
OPT_KEY( String, native_pdb )

core::Size rosetta_resnum;

bool
is_desired_resnum(
	core::pose::Pose const & pose1,
	core::pose::Pose const & ,//pose2,
	core::Size resno,
	core::Size atomno
)
{
	core::conformation::Residue const & rsd = pose1.residue(resno);
	if ( resno != rosetta_resnum ) return false;
	return rsd.is_protein() && !rsd.atom_is_hydrogen(atomno);
}


/// General testing code
int
main( int argc, char * argv [] )
{

	NEW_OPT( rms_resnum, "which residue to report the rmsd for", 1 );
	NEW_OPT( native_pdb, "compare rmsd to this pdb structure", "native.pdb" );

	devel::init(argc, argv);

	TR << "Starting chet_report_rmsd" << std::endl;

	// create pose for native pose
	pose::Pose native_pose;
	std::string const native_pdb_name ( option[ native_pdb ] );
	core::import_pose::pose_from_pdb( native_pose, native_pdb_name );

	// create pose for curr_pose
	pose::Pose curr_pose;
	std::string const input_pdb_name ( basic::options::start_file() );
	core::import_pose::pose_from_pdb( curr_pose, input_pdb_name );

	if ( curr_pose.total_residue() != native_pose.total_residue() ) {
		TR << "ERROR!! Native and comparison pose do not have the same number of residues!!" << std::endl;
		exit(1);
	}

	// find rosetta resid for the residue of interest
	int const monitor_rmsd_pdb_number = option[ rms_resnum ];
	rosetta_resnum = 0;
	for ( core::Size j = 1; j <= curr_pose.total_residue(); ++j ) {
		if ( ( curr_pose.pdb_info()->number(j) == monitor_rmsd_pdb_number ) ) {
					 rosetta_resnum = j;
		}
	}
	if ( rosetta_resnum == 0 ) {
		TR << "ERROR!! Could not find residue" << std::endl;
		exit(1);
	}

	// align curr pose onto native_pose
	protocols::moves::SuperimposeMoverOP sp_mover = new protocols::moves::SuperimposeMover( native_pose );
	sp_mover->apply( curr_pose );

	// report rmsd just for one residue
	core::Real const rms = rmsd_no_super( native_pose, curr_pose, is_desired_resnum );

	TR << "Rms for this single residue on " <<input_pdb_name<<" is: " << rms << std::endl;

	return 0;

}


