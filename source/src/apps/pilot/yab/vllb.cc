// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   tests for variable length loop building
/// @brief  apps/pilot/yab/vlb.cc
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// project headers
#include <devel/init.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/conformation/ResidueFactory.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <basic/Tracer.hh>
#include <protocols/forge/build/BuildManager.hh>
#include <protocols/forge/build/ConnectRight.hh>
#include <protocols/forge/build/SegmentRebuild.hh>
#include <protocols/forge/build/Interval.hh>
#include <protocols/forge/components/VarLengthBuild.hh>
#include <protocols/forge/methods/pose_mod.hh>
#include <protocols/forge/methods/util.hh>
#include <core/fragment/picking/FragmentLibraryManager.hh>
#include <core/fragment/picking/vall/VallLibrarian.hh>
#include <core/fragment/picking/vall/eval/EnergyEval.hh>
#include <core/fragment/picking/vall/gen/LengthGen.hh>
#include <core/fragment/picking/vall/scores/VallFragmentScore.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/viewer/viewers.hh>
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <cstdio>
#include <iostream>
#include <string>

//Auto Headers
#include <core/import_pose/import_pose.hh>


static THREAD_LOCAL basic::Tracer TR( "apps.pilot.yab.vllb" );


void * ligand_test( void * ) {
	using core::pose::Pose;
	using core::import_pose::pose_from_file;
	using core::scoring::dssp::Dssp;
	using protocols::moves::MS_SUCCESS;

	using protocols::forge::build::BuildManager;
	using protocols::forge::build::SegmentRebuild;
	using protocols::forge::components::Interval;
	using protocols::forge::components::VarLengthBuild;

	typedef std::string String;

	Pose pose;
	core::import_pose::pose_from_file( pose, "1mau_CHC_d1.pdb" , core::import_pose::PDB_file);

#if defined GL_GRAPHICS
	protocols::viewer::add_conformation_viewer( pose.conformation(), "Ligand Test" );
#endif

	// assign sec.struct
	Dssp dssp( pose );
	dssp.insert_ss_into_pose( pose );

	// Feed instructions to build manager.
	// Here we try to build an 11-mer loop in place of the original 6-mer.
	BuildManager manager;
	manager.add( new SegmentRebuild( Interval( 138, 149 ), "LLEEEEEELLLHHH" ) );

	VarLengthBuild vlb( manager );
	vlb.apply( pose );

	if ( vlb.get_last_move_status() == MS_SUCCESS ) {
		TR << "Rebuild test successfully closed loops." << std::endl;
	} else {
		TR << "Rebuild test failed to close loops." << std::endl;
	}

	pose.dump_pdb( "1mau.after_vlb.pdb" );

	return 0;
}


void * vlb_test( void * ) {
	using core::pose::Pose;
	using core::import_pose::pose_from_file;
	using core::scoring::dssp::Dssp;
	using protocols::moves::MS_SUCCESS;

	using protocols::forge::build::BuildManager;
	using protocols::forge::build::SegmentRebuild;
	using protocols::forge::components::Interval;
	using protocols::forge::components::VarLengthBuild;
	typedef std::string String;

	// input pose
	Pose pose;
	core::import_pose::pose_from_file( pose, "2bodx.pdb" , core::import_pose::PDB_file);

#if defined GL_GRAPHICS
	protocols::viewer::add_conformation_viewer( pose.conformation(), "VLB Test" );
#endif

	// assign sec.struct
	Dssp dssp( pose );
	dssp.insert_ss_into_pose( pose );

	// Feed instructions to build manager.
	// Here we try to build an 11-mer loop in place of the original 6-mer.
	BuildManager manager;
	manager.add( new SegmentRebuild( Interval( 18, 23 ), String( 11, 'H' ) ) );

	// Here we try to grow an n-terminal extension.
	manager.add( new SegmentRebuild( Interval( 1, 1 ), String( 7, 'H' ) ) );

	// Here we try to grow a c-terminal extension.
	manager.add( new SegmentRebuild( Interval( pose.n_residue(), pose.n_residue() ), String( 6, 'H' ) ) );

	// Init VLB.  Be aware this is a bootstrap implementation, even
	// remotely sane results are not guaranteed. To get things pinned
	// down with the proper implementation and benchmarked is going to
	// take some time.
	VarLengthBuild vlb( manager );
	vlb.apply( pose );

	if ( vlb.get_last_move_status() == MS_SUCCESS ) {
		TR << "Rebuild test successfully closed loops." << std::endl;
	} else {
		TR << "Rebuild test failed to close loops." << std::endl;
	}

	pose.dump_pdb( "2bodx.after_vlb.pdb" );

	return 0;
}


void * connect_test( void * ) {
	using core::chemical::ResidueTypeSetCAP;
	using core::pose::Pose;
	using core::import_pose::pose_from_file;
	using core::scoring::dssp::Dssp;
	using protocols::moves::MS_SUCCESS;

	using protocols::forge::build::BuildManager;
	using protocols::forge::build::ConnectRight;
	using protocols::forge::build::SegmentRebuild;
	using protocols::forge::components::Interval;
	using protocols::forge::components::VarLengthBuild;

	typedef std::string String;

	// load structures
	Pose brsA, brsD;
	core::import_pose::pose_from_file( brsA, "1brsA.pdb" , core::import_pose::PDB_file);
	core::import_pose::pose_from_file( brsD, "1brsD.pdb" , core::import_pose::PDB_file);

#if defined GL_GRAPHICS
	protocols::viewer::add_conformation_viewer( brsA.conformation(), "Connect Test" );
#endif

	// assign sec.struct
	Dssp dsspA( brsA ), dsspD( brsD );
	dsspA.insert_ss_into_pose( brsA );
	dsspD.insert_ss_into_pose( brsD );

	// residue type set for manager
	ResidueTypeSetCAP fa_rsd_type_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	// build a new loop, connect 1brsD, n-term ext, c-term ext
	BuildManager manager;
	manager.add( new SegmentRebuild( Interval( 66, 69 ), String( 20, 'H' ), fa_rsd_type_set, false ) );
	manager.add( new ConnectRight( 73, 39, brsD ) );
	manager.add( new SegmentRebuild( Interval( 1, 1 ), String( 7, 'H' ), fa_rsd_type_set, false ) );
	manager.add( new SegmentRebuild( Interval( brsA.n_residue(), brsA.n_residue() ), String( 6, 'H' ), fa_rsd_type_set, false ) );

	// do actual building
	VarLengthBuild vlb( manager );
	vlb.apply( brsA );

	if ( vlb.get_last_move_status() == MS_SUCCESS ) {
		TR << "Rebuild test successfully closed loops." << std::endl;
	} else {
		TR << "Rebuild test failed to close loops." << std::endl;
	}

	// output
	brsA.dump_pdb( "1brs_complex.after_vlb.pdb" );

	return 0;
}

void * picking_test( void * ) {
	using core::Size;
	using core::fragment::picking::FragmentLibraryManager;
	using core::fragment::picking::vall::VallLibrarian;
	using core::fragment::picking::vall::eval::EnergyEval;
	using core::fragment::picking::vall::gen::LengthGen;
	using core::fragment::picking::vall::scores::VallFragmentScore;
	using core::pose::Pose;
	using core::import_pose::pose_from_file;
	using core::scoring::ScoreFunction;
	using protocols::forge::build::SegmentRebuild;
	using protocols::forge::components::Interval;

	typedef std::string String;

	// input pose
	Pose pose;
	core::import_pose::pose_from_file( pose, "2bodx.pdb" , core::import_pose::PDB_file);

	// fake a 6-mer section
	SegmentRebuild rebuild( Interval( 18, 23 ), String( 6, 'H' ), &pose.residue( 1 ).residue_type_set(), false );
	rebuild.modify( pose );

	// add cutpoint variants
	Size const cutpoint = protocols::forge::methods::find_cutpoint( pose, 18, 23 );
	protocols::forge::methods::add_cutpoint_variants( pose, cutpoint );

	// define score function with quadratic chainbreak
	ScoreFunction fx;
	fx.set_weight( core::scoring::chainbreak, 1.0 );

	// pick
	VallLibrarian librarian;
	librarian.add_fragment_gen( new LengthGen( 6 ) ); // 6-mers
	librarian.add_fragment_eval( new EnergyEval( pose, 18, fx ) );

	// catalog fragments
	librarian.catalog( FragmentLibraryManager::get_instance()->get_Vall() );

	return 0;
}


int main( int argc, char * argv [] ) {
	try {
	// initialize rosetta
	devel::init( argc, argv );

	protocols::viewer::viewer_main( ligand_test );

	 } catch ( utility::excn::EXCN_Base const & e ) {
		 std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}

