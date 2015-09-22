// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief


// libRosetta headers

#include <protocols/viewer/viewers.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/dna/setup.hh>
#include <core/scoring/dna/base_geometry.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
//#include <core/scoring/elec/FA_ElecEnergy.hh>
//#include <core/scoring/etable/EtableEnergy.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
//#include <core/scoring/etable/count_pair/CountPairAll.hh>
//#include <core/scoring/etable/count_pair/CountPairFunction.hh>
//#include <core/scoring/etable/count_pair/CountPairFactory.hh>
//#include <core/scoring/etable/count_pair/CountPair1BC4.hh>

#include <core/types.hh>

#include <core/chemical/AA.hh>

#include <core/conformation/Residue.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>


#include <core/pose/Pose.hh>


//#include <basic/options/after_opts.hh>


#include <devel/init.hh>

#include <core/io/pdb/pose_io.hh>

#include <utility/vector1.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>


#include <ObjexxFCL/format.hh>


// C++ headers
//#include <cstdlib>
#include <iostream>
#include <string>
#include <set>
#include <cstdlib>
#include <sstream>

//silly using/typedef


#include <basic/Tracer.hh>

#include <core/import_pose/import_pose.hh>
#include <core/kinematics/Jump.hh>
#include <utility/vector0.hh>
#include <basic/options/keys/OptionKeys.hh>


using namespace core;
using namespace protocols;
using namespace ObjexxFCL;
using namespace ObjexxFCL::format;

////////////////////////////////////////////////
// danger USING ////////////////////////////////
using namespace core;
using namespace protocols;
using namespace pose;
using namespace conformation;
using namespace chemical;
using namespace scoring;
using namespace basic::options;
using utility::vector1;
using std::string;
using core::import_pose::pose_from_pdb;
////////////////////////////////////////////////


static THREAD_LOCAL basic::Tracer TR( "apps.pilot.chen.simple_dna_regression_test" );

/// @details  Show details of the internal DNA geometry
void
show_dna_geometry( pose::Pose const & pose )
{
	scoring::dna::show_dna_geometry( pose, TR.Info );
	TR.Info << std::endl;
}


// /// @details  Test static scoring of the pose
// void
// score_test( pose::Pose & pose )
// {

//  { // intra-dna scoring
//   ScoreFunction scorefxn;

//  }


// }


void
repack_test( pose::Pose & pose )
{

	// Parameters:
	std::string const weights_file( "my_dna.wts" );

	ScoreFunction scorefxn;

	//scorefxn.energy_method_options().exclude_DNA_DNA( false );
	// Safer:
	methods::EnergyMethodOptions options( scorefxn.energy_method_options() );
	options.exclude_DNA_DNA( false );
	scorefxn.set_energy_method_options( options );

	// scorefxn.add_weights_from_file( basic::database::full_name( "scoring/weights/"+weights_file ) );
	scorefxn.add_weights_from_file( weights_file );
	scorefxn.set_weight( atom_pair_constraint, 1.0 );
	scorefxn.set_weight(     angle_constraint, 1.0 );

	// score starting pose
	TR << "Score before packing: " << scorefxn( pose ) << std::endl;
	pose.energies().total_energies().show_nonzero( TR );
	TR << std::endl;

	pack::task::PackerTaskOP task
		( pack::task::TaskFactory::create_packer_task( pose ));
	task->set_bump_check( true );

	Size const nres( pose.total_residue() );

	task->initialize_from_command_line();
	for ( Size ii = 1; ii <= nres; ++ii ) {
		if ( pose.residue(ii).is_protein() ) {
			task->nonconst_residue_task( ii ).restrict_to_repacking();
			task->nonconst_residue_task( ii ).restrict_to_repacking();
			assert( task->pack_residue(ii) );
		} else {   // repack_DNA
			task->nonconst_residue_task( ii ).restrict_to_repacking();
		}

		assert( !task->design_residue(ii) );
	}

	// dont include current
	task->or_include_current( false );

	pack::pack_rotamers( pose, scorefxn, task );

	TR << "Pack score: " << scorefxn( pose ) << std::endl;
	pose.energies().total_energies().show_nonzero( TR );
	TR << std::endl;

	io::pdb::dump_pdb( pose, "1aay_repacked.pdb" );
}


/// @details  The main routine, calls smaller tests.
void
run_tests()
{

	// reads a single dna-containing pose
	Pose pose;
	core::import_pose::pose_from_pdb( pose, "input/1aay.pdb" );

	scoring::dna::set_base_partner( pose );

	// geometry of dna
	show_dna_geometry( pose );

	// static scoring
	//  score_test( pose );

	// repack test
	repack_test( pose );

	// rotamer trials test

	// correlate basepair dna design test

	// minimization test


}


void*
my_main( void* )
{
	run_tests();
	exit(0);
}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {

		// initialize option and random number system
		devel::init( argc, argv );

		protocols::viewer::viewer_main( my_main );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
