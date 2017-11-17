// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
#include <devel/init.hh>

#include <core/io/pdb/build_pose_as_is.hh>
#include <core/io/pdb/pdb_writer.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/jobdist/Jobs.hh>
#include <protocols/jobdist/standard_mains.hh>

#include <protocols/loops/loop_closure/kinematic_closure/KinematicMover.hh>

#include <utility/string_util.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


using namespace core;

int main( int argc, char * argv [] )
{
  try {
	devel::init( argc, argv );
	utility::vector1< protocols::jobdist::BasicJobOP > input_jobs = protocols::jobdist::load_s_and_l();
	pose::Pose pose;
	core::import_pose::pose_from_file( pose, input_jobs[ 1 ]->input_tag() , core::import_pose::PDB_file);

	protocols::loops::loop_closure::kinematic_closure::KinematicMover myKinematicMover( 0.6 );
	myKinematicMover.set_rama_check( true );
	// DJM: debug seeing if this is causing problems with the loopmodel extend loops
	//myKinematicMover.set_idealize_loop_first( true ); // start without any native angles or lengths

	Size alc_start, alc_middle, alc_end; // three pivot residues for kinematic loop closure
	alc_start = 2;
	alc_end = 10;
	Size middle_offset = (alc_end - alc_start) / 2; // need to ensure this isn't a proline
	alc_middle = alc_start + middle_offset;

	myKinematicMover.set_pivots(alc_start, alc_middle, alc_end);
	myKinematicMover.set_temperature(0.6);

	myKinematicMover.set_sweep_nonpivot_torsions( true );

	utility::vector1< Size > res_to_sweep( 4 );
	utility::vector1< Size > tors_to_sweep( 4 );
	utility::vector1< Real > start_angles( 4 );
	utility::vector1< Real > step_sizes( 4 );
	utility::vector1< Size > nsteps( 4 );

	res_to_sweep[ 1 ] = 4; tors_to_sweep[ 1 ] = 1; start_angles[ 1 ] = pose.phi( 4 ) - 5; step_sizes[ 1 ] = 2.5; nsteps[ 1 ] = 5;
	res_to_sweep[ 2 ] = 4; tors_to_sweep[ 2 ] = 2; start_angles[ 2 ] = pose.psi( 4 ) - 5; step_sizes[ 2 ] = 2.5; nsteps[ 2 ] = 5;
	res_to_sweep[ 3 ] = 8; tors_to_sweep[ 3 ] = 1; start_angles[ 3 ] = pose.phi( 8 ) - 5; step_sizes[ 3 ] = 2.5; nsteps[ 3 ] = 5;
	res_to_sweep[ 4 ] = 8; tors_to_sweep[ 4 ] = 2; start_angles[ 4 ] = pose.psi( 8 ) - 5; step_sizes[ 4 ] = 2.5; nsteps[ 4 ] = 5;

	myKinematicMover.set_nonpivot_res_to_sweep( res_to_sweep );
	myKinematicMover.set_nonpivot_bb_torsion_id( tors_to_sweep );
	myKinematicMover.set_sweep_start_angle( start_angles );
	myKinematicMover.set_sweep_step_size( step_sizes );
	myKinematicMover.set_sweep_nsteps( nsteps );

	//myKinematicMover.set_sample_nonpivot_torsions( true );

	Size ii = 1;
	scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function();
	while ( myKinematicMover.sweep_incomplete() ) {
		pose::Pose loop_pose( pose );
		myKinematicMover.apply( loop_pose );
		std::cout << "Found one sc: " << (*sfxn)( loop_pose ) << std::endl;
		loop_pose.dump_pdb( "test_sweep_" + utility::to_string( ii ) + ".pdb" );

		++ii;
	}

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
