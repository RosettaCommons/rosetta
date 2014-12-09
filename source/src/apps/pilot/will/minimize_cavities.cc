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


#include <core/types.hh>
#include <devel/init.hh>


#include <core/chemical/AA.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueSelector.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/residue_io.hh>
#include <core/chemical/VariantType.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/Constraints.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/ScoreFunction.hh>

#include <core/io/pdb/pose_io.hh>

#include <core/pose/Pose.hh>

#include <basic/options/option.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/ScoreMover.hh>

#include <core/scoring/packstat/types.hh>
#include <core/scoring/packstat/compute_sasa.hh>
#include <core/scoring/packstat/packing_score_params.hh>
#include <core/scoring/packstat/AtomRadiusMap.hh>
#include <core/scoring/packstat/SimplePDB.hh>
#include <protocols/simple_moves/AddCavitiesMover.hh>

#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <utility/io/izstream.hh>

#include <numeric/xyzVector.hh>
#include <numeric/io.hh>

#include <ObjexxFCL/format.hh>

#include <fstream>

// C++ headers
#include <iostream>
#include <string>
#include <vector>

// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>




using namespace std;
using namespace utility;
using namespace numeric;
using namespace core;
using core::pose::Pose;
using core::import_pose::pose_from_pdb;
using core::io::pdb::dump_pdb;
using core::io::pdb::dump_bfactor_pdb;
using ObjexxFCL::string_of;
using core::id::AtomID;
using namespace core::scoring::packstat;
using scoring::constraints::ConstraintOP;
using scoring::constraints::FuncOP;
using scoring::constraints::ConstraintSetOP;
using core::Real;
// using namespace ObjexxFCL::format;

void
test_suck_res( std::string fname ) {
	using namespace core;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace basic;
	using namespace core::scoring::packstat;

	Pose pose;
	core::import_pose::pose_from_pdb( pose, fname );

	kinematics::MoveMapOP mm = new kinematics::MoveMap;
	mm->set_bb ( true ); mm->set_chi( true );	mm->set_jump( true );

  ScoreFunctionOP sfstd( get_score_function_legacy( PRE_TALARIS_2013_STANDARD_WTS ) );
  ScoreFunctionOP sf( get_score_function_legacy( PRE_TALARIS_2013_STANDARD_WTS ) );
	// ScoreFunctionOP sf( new ScoreFunction );
	sf->set_weight( suck , 1.0 );
	sf->set_weight( coordinate_constraint, 1.0 );
	sfstd->set_weight( suck , 0.0 );
	sfstd->set_weight( coordinate_constraint, 0.0 );
	// std::cerr << *sf << std::endl;
	simple_moves::AddCavitiesMover add_cav(50,1.0,150,3.0);
	simple_moves::AddCavitiesMover rem_cav(0);

	core::Real ps = compute_packing_score(pose,0);
	std::cerr << "SC " << I(2,-1) << " "
		 								 << I(2,-1) << " "
		 								 << F(7,3,0.0) << " "
										 << F(7,3,0.0) << " "
									   << F(8,2,(*sfstd)(pose) ) << " "
										 << F(8,2,(*sf)(pose) ) << " "
										 << F(8,2,pose.energies().total_energies()[suck]) << " "
										 << F(8,2,pose.energies().total_energies()[coordinate_constraint]) << " "
										 << F(8,4,ps ) << " "
										 << std::endl;

	protocols::simple_moves::MinMover minstd( mm, sfstd , "dfpmin", 0.001, true, false, false );
	minstd.min_options()->nblist_auto_update(true);
	pack::task::PackerTaskOP taskstd = pack::task::TaskFactory::create_packer_task( pose );
	taskstd->restrict_to_repacking();
	protocols::simple_moves::RotamerTrialsMover packstd( sfstd, *taskstd );
	minstd .apply( pose );
	packstd.apply( pose );
	// core::Real score_before = (*sfstd)(pose), score_mid=0,score_after=0;
	// core::Real    ps_before = get_packing_score(pose), ps_mid=0,ps_after=0;
	// std::cerr << "packing score before: " << ps_before << std::endl;

	ps = compute_packing_score(pose,0);
	std::cerr << "SC " << I(2, 0) << " "
		 								 << I(2, 0) << " "
		 								 << F(7,3,0.0) << " "
										 << F(7,3,0.0) << " "
									   << F(8,2,(*sfstd)(pose) ) << " "
										 << F(8,2,(*sf)(pose) ) << " "
										 << F(8,2,pose.energies().total_energies()[suck]) << " "
										 << F(8,2,pose.energies().total_energies()[coordinate_constraint]) << " "
										 << F(8,4,ps ) << " "
										 << std::endl;

	protocols::simple_moves::MinMover min( mm, sf , "dfpmin", 0.001, true, false, false );
	min.min_options()->nblist_auto_update(true);
	pack::task::PackerTaskOP task = pack::task::TaskFactory::create_packer_task( pose );
	task->restrict_to_repacking();
	protocols::simple_moves::RotamerTrialsMover pack( sf, *task );


	// vector1<Real> pscores;
	// pscores.push_back(ps_before);
	for( int sucker_iter = 1; sucker_iter <= 3; ++sucker_iter ) {

		// std::cerr << "add cavs" << std::endl;
		rem_cav.apply( pose );
		add_cav.apply( pose );	// will remove if already there
		dump_pdb( pose, fname + "_with_suckers_before" + to_string<int>( sucker_iter ) + ".pdb" );

		taskstd = pack::task::TaskFactory::create_packer_task( pose );
		taskstd->restrict_to_repacking();
		packstd = protocols::simple_moves::RotamerTrialsMover( sfstd, *taskstd );
		task = pack::task::TaskFactory::create_packer_task( pose );
		task->restrict_to_repacking();
		pack = protocols::simple_moves::RotamerTrialsMover( sf, *task );

		core::Real sckwts[] =	{ 0, 1.0, 5, 12.0, 35.5,  3.0, 1.5, 0.0 };
		core::Real cc_wts[] = { 0, 1.0, 2,  5.0,  9.0, 16.0, 3.0, 0.0 };
		for( int i = 1; i <= 7; ++i ) {
			sf->set_weight( suck                 , sckwts[i] );
			sf->set_weight( coordinate_constraint, cc_wts[i] );
			// std::cerr << "pack" << std::endl;
			pack.apply(pose);
			// std::cerr << "min" << std::endl;
			min.apply(pose);
			// std::cerr << "score" << std::endl;
			ps = compute_packing_score(pose,0);
			std::cerr << "SC " << I(2,sucker_iter) << " "
												 << I(2, i) << " "
												 << F(7,3,sckwts[i]) << " "
												 << F(7,3,cc_wts[i]) << " "
											   << F(8,2,(*sfstd)(pose) ) << " "
												 << F(8,2,(*sf)(pose) ) << " "
												 << F(8,2,pose.energies().total_energies()[suck]) << " "
												 << F(8,2,pose.energies().total_energies()[coordinate_constraint]) << " "
												 << F(8,4,ps ) << " "
												 << std::endl;
			// std::cerr << pose.energies() << std::endl;
		}
		// std::cerr << "rm cav" << std::endl;
		rem_cav.apply(pose);
		taskstd = pack::task::TaskFactory::create_packer_task( pose );
		taskstd->restrict_to_repacking();
		packstd = protocols::simple_moves::RotamerTrialsMover( sfstd, *taskstd );
		task = pack::task::TaskFactory::create_packer_task( pose );
		task->restrict_to_repacking();
		pack = protocols::simple_moves::RotamerTrialsMover( sf, *task );

		// std::cerr << "pack" << std::endl;
		packstd.apply(pose);
		// std::cerr << "min" << std::endl;
		minstd .apply(pose);
		// std::cerr << "score" << std::endl;
		ps = compute_packing_score(pose,0);
		std::cerr << "SC " << I(2,99) << " "
											 << I(2,99) << " "
											 << F(7,3,0.0) << " "
											 << F(7,3,0.0) << " "
										   << F(8,2,(*sfstd)(pose) ) << " "
											 << F(8,2,(*sf)(pose) ) << " "
											 << F(8,2,pose.energies().total_energies()[suck]) << " "
											 << F(8,2,pose.energies().total_energies()[coordinate_constraint]) << " "
											 << F(8,4,ps ) << " "
											 << std::endl;

		dump_pdb( pose, fname + "_with_suckers_after" + to_string<int>( sucker_iter ) + ".pdb" );

	}
	// std::cerr << "packing scores: " << pscores.size() << endl;
	// for( int i = 1; 1 <= pscores.size(); ++i ) {
	// 	std::cerr << pscores[i] << " ";
	// }
	// std::cerr << std::endl;
}

// void
// test_min_cav( string fname ){
// 	using namespace core::scoring;
// 	using namespace constraints;
//
// 	Pose pose;
// 	core::import_pose::pose_from_pdb( pose, fname );
//
// 	kinematics::MoveMapOP mm = new kinematics::MoveMap;
// 	mm->set_bb ( true );
// 	mm->set_chi( true );
//
//   ScoreFunctionOP sf( get_score_function() );
// 	sf->set_weight( atom_pair_constraint, 1.0 );
//
//   ScoreFunctionOP sfnoc( get_score_function() );
// 	sfnoc->set_weight( atom_pair_constraint, 0.0 );
//
// 	std::cerr << "get cavities" << std::endl;
// 	CavBalls cbs = get_cavities( pose, 10.0, 150, 3.0 );
// 	int Ncb = 5;
//
// 	std::cerr << "make constraints" << std::endl;
// 	ConstraintSetOP cstset = new ConstraintSet();
// 	FuncOP h08 = new HarmonicFunc( 4.0, 6.0 );
// 	for( int i = 1; i <= Ncb; ++i ) {
// 		if( cbs[i].radius() < 1.8 ) {
// 			Ncb = i-1;
// 			break;
// 		}
// 		std::cerr << "adding cb" << cbs[i].str() << std::endl;
// 		add_surrounding_constraints( pose, cbs[i].xyz(), cbs[i].radius()+3.2, h08, cstset );
// 	}
// 	pose.constraint_set( cstset );
//
// 	protocols::simple_moves::MinMover mincst( mm, sf   , "dfpmin", 0.1, true );
// 	mincst.min_options()->nblist_auto_update(true);
// 	protocols::simple_moves::MinMover minnoc( mm, sfnoc, "dfpmin", 0.1, true );
// 	minnoc.min_options()->nblist_auto_update(true);
//
// 	pack::task::PackerTaskOP task = pack::task::TaskFactory::create_packer_task( pose );
// 	task->restrict_to_repacking();
// 	// protocols::simple_moves::PackRotamersMover packcst( sf, *task );
// 	// protocols::simple_moves::PackRotamersMover packnoc( sfnoc, *task );
// 	protocols::simple_moves::RotamerTrialsMover packcst( sf, *task );
// 	protocols::simple_moves::RotamerTrialsMover packnoc( sfnoc, *task );
//
// 	ofstream out( (fname.substr(0,fname.size()-4)+"_before.pdb").c_str(), ios_base::app );
// 	dump_pdb( pose, out );
// 	for( int i = 1; i <= Ncb; ++i )
// 		out << cbs[i].hetero_atom_line() << std::endl;
// 	out.close();
//
// 	Real CST_WEIGHT = 50.0;
// 	Real FA_REP_WT  = sf->get_weight( fa_rep );
//
// 	for( int i = 0; i <= 9; ++i ) {
// 		std::cerr << "applying... " << i << " " << (*sf)(pose) << " " << (*sfnoc)(pose)
//       	      << " SCORE FUNC: " << sf->get_weight(fa_rep) << " " << sf->get_weight(atom_pair_constraint)
//     					<< " SCORE FUNC2: " << sfnoc->get_weight(fa_rep) << " " << sfnoc->get_weight(atom_pair_constraint) << std::endl;
// 		sf->set_weight( atom_pair_constraint, ((Real)i) * 0.1 * CST_WEIGHT );
// 		sf->set_weight( fa_rep, ((Real)(10-i)) * 0.10 * FA_REP_WT );
// 		// std::cerr << "========================= about to pack =======================" << std::endl;
// 		packcst.apply( pose );
// 		// std::cerr << "========================== about to min =======================" << std::endl;
// 		mincst.apply( pose );
// 		// std::cerr << "============================ done min =========================" << std::endl;
// 		dump_pdb( pose, fname.substr(0,fname.size()-4)+"_after"+string_of(i)+".pdb" );
// 	}
// 	for( int i = 10; i <= 19; ++i ) {
// 		std::cerr << "applying... " << i << " " << (*sf)(pose) << " " << (*sfnoc)(pose)
// 			        << " SCORE FUNC: " << sf->get_weight(fa_rep) << " " << sf->get_weight(atom_pair_constraint) << std::endl;
// 		sf->set_weight( fa_rep, ((Real)(i-10)+1) * 0.1 * FA_REP_WT );
// 		sf->set_weight( atom_pair_constraint, ((Real)(19-i)) * 0.1 * CST_WEIGHT );
// 		// std::cerr << "========================= about to pack =======================" << std::endl;
// 		packcst.apply( pose );
// 		// std::cerr << "========================== about to min =======================" << std::endl;
// 		mincst.apply( pose );
// 		// std::cerr << "============================ done min =========================" << std::endl;
// 		dump_pdb( pose, fname.substr(0,fname.size()-4)+"_after"+string_of(i)+".pdb" );
// 	}
// 	std::cerr << "applying... " << 20 << " " << (*sf)(pose) << " " << (*sfnoc)(pose) << std::endl;
// 	packnoc.apply( pose );
// 	minnoc.apply( pose );
// 	dump_pdb( pose, fname.substr(0,fname.size()-4)+"_after20.pdb" );
// 	std::cerr << "done... " << (*sf)(pose) << " " << (*sfnoc)(pose) << std::endl;
//
// }
//

int
main (int argc, char *argv[])
{

	try {



	devel::init( argc, argv );

  using namespace core::scoring::packstat;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using namespace utility;

  // test_io();

	// test_sasa_dots();

	if( option[ in::file::s ].user() ) {
  	vector1<file::FileName> files( option[ in::file::s ]() );
  	for( size_t i = 1; i <= files.size(); ++i ) {
    	test_suck_res( files[i] );
  	}
	} else if( option[ in::file::l ].user() ) {
  	vector1<file::FileName> files( option[ in::file::l ]() );
  	for( size_t i = 1; i <= files.size(); ++i ) {
			utility::io::izstream list( files[i] );
			std::string fname;
			while( list >> fname ) {
				// std::cerr << "'" << fname << "'" << std::endl;
    		test_suck_res( fname );
			}
  	}
	}
	return 0;



	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
