// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file calibrate_pdb.cc
/// @brief Calibrates the input structure for the Rosetta forcefield, using repacking, rotamer trials and minimization.
/// @author Roland A. Pache, PhD

//Unit Headers
#include <protocols/loops/Loops.hh>
#include <protocols/viewer/viewers.hh>


//Rosetta Headers
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/id/TorsionID.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <basic/options/option.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <devel/init.hh>

#include <basic/Tracer.hh> //tracer output

#include <utility/io/izstream.hh>
#include <protocols/jobdist/Jobs.hh>
#include <protocols/jobdist/standard_mains.hh>

//C++ Headers
#include <iostream>
#include <map>
#include <string>

//option key includes

#include <basic/options/option_macros.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


using namespace core;

static basic::Tracer TR( "pilot_apps.calibrate_pdb_via_sidechain_optimization" );

OPT_1GRP_KEY(Boolean, calibrate_pdb, no_repacking)
OPT_1GRP_KEY(Boolean, calibrate_pdb, no_rottrials)
OPT_1GRP_KEY(Boolean, calibrate_pdb, no_minimization)

/////////
/////////main function
/////////

void*
my_main( void* )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	//load input PDB into pose
	utility::vector1< protocols::jobdist::BasicJobOP > input_jobs = protocols::jobdist::load_s_and_l();
	if ( input_jobs.size() != 1 ) {
		utility_exit_with_message( "Expected exactly one pdb to be specified from the -s or -l flags" );
	}
	pose::Pose pose;
	import_pose::pose_from_file( pose, input_jobs[ 1 ]->input_tag() , core::import_pose::PDB_file);
	int const nres( pose.size() );

	//define score function
	scoring::ScoreFunctionOP scorefxn=core::scoring::get_score_function();

	//get starting energy
	Real starting_total_energy = (*scorefxn)(pose);

	//create packer task
	pack::task::PackerTaskOP base_packer_task( pack::task::TaskFactory::create_packer_task( pose ));
	base_packer_task->initialize_from_command_line().restrict_to_repacking().or_include_current( false );
	base_packer_task->set_bump_check( true );
	pack::task::PackerTaskOP this_packer_task( base_packer_task->clone() );
	utility::vector1<bool> allow_repacked( nres, true );

	//create minimizer
	optimization::AtomTreeMinimizer minimizer;
	optimization::MinimizerOptions options( "lbfgs_armijo_nonmonotone", 0.001, true /*use_nblist*/, false /*deriv_check*/ );

	//prepare move map for repacking, rotamer trials and minimization, changing non-disulfide bonded residues only
	kinematics::MoveMap mm_all_sc;
	mm_all_sc.set_bb( false );
	mm_all_sc.set_chi( false );
	for ( Size i=1; i<= pose.size(); i++ ) {
		if ( pose.residue(i).type().is_disulfide_bonded() || pose.residue(i).has_variant_type( chemical::SIDECHAIN_CONJUGATION ) ) {
			allow_repacked[i] = false;
			TR << "Disabling side-chain optimization of disulfide bonded residue " << i << std::endl;
		} else {
			mm_all_sc.set_chi( true );
		}
	}
	this_packer_task->restrict_to_residues( allow_repacked );

	//repack
	if ( !option[calibrate_pdb::no_repacking].user() ) {
		TR << std::endl << "repacking pose..." << std::endl;
		pack::pack_rotamers( pose, *scorefxn, this_packer_task );
	}

	//perform rotamer trials
	if ( !option[calibrate_pdb::no_rottrials].user() ) {
		TR << std::endl << "performing one round of rotamer trials..." << std::endl;
		pack::rotamer_trials( pose, *scorefxn, this_packer_task );
	}

	//minimize sidechains
	if ( !option[calibrate_pdb::no_minimization].user() ) {
		TR << std::endl << "minimizing pose..." << std::endl;
		minimizer.run( pose, mm_all_sc, *scorefxn, options );
	}

	//get final energy
	Real last_total_energy = (*scorefxn)(pose);
	TR << std::endl << "total_energy_before_calibration: " << starting_total_energy << std::endl;
	TR << "total_energy_after_calibration: " << last_total_energy << std::endl;

	//output the minimized pose
	std::string pdb_prefix = option[out::prefix];
	std::string outname=option[out::path::path]().name()+pdb_prefix+"_calibrated.pdb";
	std::ofstream out(outname.c_str(), std::ios::out | std::ios::binary);
	io::pdb::dump_pdb( pose, out );
	out << "total_energy_before_calibration: " << starting_total_energy << std::endl;
	out << "total_energy_after_calibration: " << last_total_energy << std::endl;

	return 0;
}

int
main( int argc, char * argv [] )
{

	try {

		// define viable options
		NEW_OPT(calibrate_pdb::no_repacking, "skip repacking", false);
		NEW_OPT(calibrate_pdb::no_rottrials, "skip rotamer trials", false);
		NEW_OPT(calibrate_pdb::no_minimization, "skip minimization", false);

		// initialize option and random number system
		devel::init( argc, argv );

		protocols::viewer::viewer_main( my_main );

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
