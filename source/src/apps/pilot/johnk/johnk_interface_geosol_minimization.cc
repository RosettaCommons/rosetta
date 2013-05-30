// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @brief
/// @author jk

// Project Headers
#include <devel/init.hh>
#include <core/types.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/pose/Pose.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Atom.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/etable/EtableOptions.hh>
#include <core/scoring/hbonds/hbonds_geom.hh>
#include <core/scoring/hbonds/types.hh>
#include <core/scoring/hbonds/constants.hh>
#include <basic/Tracer.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <protocols/toolbox/task_operations/RestrictToInterface.hh>
#include <core/pack/pack_rotamers.hh>

#include <numeric/constants.hh>

//#include <core/scoring/ScoreFunction.hh>
//#include <core/scoring/ScoreFunctionFactory.hh>

#include <basic/options/option_macros.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/excn/Exceptions.hh>

// C++ Headers
#include <cmath>
#include <iostream>
#include <iomanip>
#include <map>

//Auto Headers
#include <core/import_pose/import_pose.hh>


using namespace core;
using namespace core::scoring;
using namespace core::scoring::hbonds;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::optimization;

OPT_KEY( Boolean, minimize_chi )
OPT_KEY( Boolean, minimize_bb )
OPT_KEY( Boolean, minimize_jump )
OPT_KEY( Boolean, do_repack )
OPT_KEY( Boolean, use_geosol )

static basic::Tracer TR( "apps.pilot.johnk_test_interface_geosol_minimization.main" );

// NOTE: GLOBAL VARIABLE - this would better be done with a static (and could be used to ensure it's filled, too...)
utility::vector1 <bool> interface;

void define_interface( core::pose::Pose const & pose ) {
	interface.resize( pose.total_residue(), false );
	core::Real interface_dist = 8.0;
	core::Size rb_jump = 1;
	scoring::ScoreFunctionOP scorefxn( getScoreFunction() );
	pack::task::TaskFactory tf;
	tf.push_back( new protocols::toolbox::task_operations::RestrictToInterface( rb_jump, interface_dist ) );
	pack::task::PackerTaskOP task = tf.create_task_and_apply_taskoperations( pose );
	for ( core::Size i=1; i <= pose.total_residue(); ++i ) {
		if ( task->pack_residue(i) ) interface.at(i)=true;
	}
}


bool
is_interface_sc(
	core::pose::Pose const & pose1,
	core::pose::Pose const & ,//pose2,
	core::Size resno,
	core::Size atomno
)
{
	core::conformation::Residue const & rsd = pose1.residue(resno);
	if ( ! interface.at(resno) ) return false;
	return rsd.is_protein() && !rsd.atom_is_backbone(atomno) && !rsd.atom_is_hydrogen(atomno);
}


/// General testing code
int
main( int argc, char * argv [] )
{
    try {
	NEW_OPT( minimize_chi, "include chi angles in minimization", false );
	NEW_OPT( minimize_bb, "include backbone in minimization", false );
	NEW_OPT( minimize_jump, "include jumps in minimization", false );
	NEW_OPT( do_repack, "repack the interface", false );
	NEW_OPT( use_geosol, "use geosol instead of LK", false );

	devel::init(argc, argv);

	TR << "jk doing geosol interface minimiations" << std::endl;

	// scoring function
	scoring::ScoreFunctionOP scorefxn( getScoreFunction() );

	if ( option[ use_geosol ] ) {
		//	scorefxn->reset();
		//	scorefxn->set_weight( core::scoring::fa_sol, 0.65 );
		scorefxn->set_weight( core::scoring::fa_sol, 0.0 );
		scorefxn->set_weight( core::scoring::occ_sol_fitted, 0.65 );
	}

	// Read input pose, with support for waters
	pose::Pose input_pose, pose;
	std::string const input_pdb_name ( basic::options::start_file() );
	core::import_pose::pose_from_pdb( input_pose, input_pdb_name );
	(*scorefxn)(input_pose);
	pose = input_pose;

	define_interface( pose );

	// setup a packertask
	utility::vector1 <bool>	allow_moving( pose.total_residue(), false );
	pack::task::PackerTaskOP packer_task( pack::task::TaskFactory::create_packer_task( pose) );

	// setting degrees of freedom which can move during minimization
	kinematics::MoveMap mm_all;
	mm_all.set_chi( false );
	mm_all.set_bb( false );
	mm_all.set_jump( false );

	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		if ( ! interface.at(ii) ) continue;
		mm_all.set_chi( ii, option[ minimize_chi ] );
		mm_all.set_bb( ii, option[ minimize_bb ] );
		mm_all.set_jump( ii, option[ minimize_jump ] );
		allow_moving.at(ii) = true;
		packer_task->nonconst_residue_task( ii ).restrict_to_repacking();
	}
	packer_task->restrict_to_residues( allow_moving );

	// do the repack
	if ( option[ do_repack ] ) {
		TR << "Starting repack...." << std::endl;
		pack::pack_rotamers( pose, *scorefxn, packer_task );
	}

	// minimize protein
	if ( option[ minimize_chi ] || option[ minimize_bb ] || option[ minimize_jump ] ) {
		TR << "Starting minimization...." << std::endl;
		AtomTreeMinimizer minimizer;
		MinimizerOptions min_options( "dfpmin", 0.001, true );
		minimizer.run( pose, mm_all, *scorefxn, min_options );
	}

	(*scorefxn)(pose);

	TR << "jk done minimizing" << std::endl;

	// compute sidechain interface rmsd
	//	core::Real interface_sc_rms = rmsd_no_super_subset( input_pose, pose, interface, is_interface_sc );
	core::Real interface_sc_rms = rmsd_no_super( input_pose, pose, is_interface_sc );
	TR << "jk interface sidechain rmsd is " << interface_sc_rms << std::endl;

	// compute per-residue interface rmsd
	// JK FILL THIS IN IF NEEDED....

	if ( option[ use_geosol ] ) {
		pose.dump_pdb("geosol_min.pdb");
	} else {
		pose.dump_pdb("LK_min.pdb");
	}

	TR << "jk done analysis" << std::endl;

    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
    }
    return 0;
}

