// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/denovo_design/DisulfidizeMover.cxxtest.hh
/// @brief  test suite for protocols::denovo_design::components::DisulfidizeMover
/// @author Tom Linsky (tlinsky@uw.edu)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Unit headers
#include <protocols/denovo_design/movers/DisulfidizeMover.hh>

// Protocol headers
#include <protocols/denovo_design/util.hh>
#include <protocols/simple_moves/MutateResidue.hh>

// Project headers

// Core headers
#include <core/io/pdb/build_pose_as_is.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/disulfides/DisulfideMatchingPotential.hh>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.DisulfidizeMover.cxxtest" );

//Uncomment the following line to dump PDBs, to debug:
//#define test_protocols_denovo_design_DisulfidizeMoverTests_dump_pdbs

// --------------- Test Class --------------- //
class DisulfidizeMoverTests : public CxxTest::TestSuite {
public:

	// Shared initialization goes here.
	void setUp() {
		// load params for ligand
		protocols_init();

		// set preserve header always
		basic::options::option[basic::options::OptionKeys::run::preserve_header].value(true);
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_multiplex() {
		using namespace protocols::denovo_design;
		core::pose::Pose input_pose;
		core::io::pdb::build_pose_from_pdb_as_is( input_pose, "protocols/denovo_design/disulf_test.pdb" );

		DisulfidizeMover disulf;
		disulf.set_min_loop( 6 );

		core::pose::PoseOP runpose = input_pose.clone();
		disulf.apply( *runpose );
		TS_ASSERT_EQUALS( disulf.get_last_move_status(), protocols::moves::MS_SUCCESS );

		utility::vector1< core::pose::PoseOP > poses;
		poses.push_back( runpose );

#ifdef test_protocols_denovo_design_DisulfidizeMoverTests_dump_pdbs
		runpose->dump_pdb("temp_0000.pdb"); //For debugging only
#endif

		core::pose::PoseOP additional = disulf.get_additional_output();

#ifdef test_protocols_denovo_design_DisulfidizeMoverTests_dump_pdbs
		core::Size temp(0);
#endif

		while ( additional ) {
#ifdef test_protocols_denovo_design_DisulfidizeMoverTests_dump_pdbs
			++temp;
			char name[256];
			sprintf(name, "temp_%04lu.pdb", temp);
			std::string namestr(name);
			additional->dump_pdb(namestr);
#endif
			poses.push_back( additional );
			additional = disulf.get_additional_output();
		}

		// there should be five results
		// hpark: beta_nov15 makes only 3; this may not of issue
		TS_ASSERT_EQUALS( poses.size(), /*5*/ 3 );

		// each should have disulfides
		std::set< core::Size > num_disulf;
		BOOST_FOREACH ( core::pose::PoseOP p, poses ) {
			core::Size cyd_count = 0;
			for ( core::Size i=1; i<=p->size(); ++i ) {
				TS_ASSERT( p );
				if ( p->residue(i).type().is_disulfide_bonded() ) {
					++cyd_count;
				}
			}
			TS_ASSERT( cyd_count );
			num_disulf.insert( cyd_count );

			// should be DISULFIDIZE pdb residue info in the pose
			TS_ASSERT_DIFFERS( p->pdb_info(), core::pose::PDBInfoOP() );
			utility::vector1< core::Size > tagged;
			for ( core::Size resid=1; resid<=p->size(); ++resid ) {
				if ( p->pdb_info()->res_haslabel( resid, "DISULFIDIZE" ) ) {
					tagged.push_back( resid );
				}
			}
			TS_ASSERT_EQUALS( tagged.size(), cyd_count );
			TR.Debug << "Found the following residues in disulfides: " << tagged << std::endl;

			// records should be dumped with the pose
			std::stringstream pdbstr;
			//p->dump_pdb( pdbstr, "1" ); JAB - tag unused
			p->dump_pdb( pdbstr);
			TS_ASSERT_DIFFERS( pdbstr.str().find( "DISULFIDIZE" ), std::string::npos );
		}
		TS_ASSERT( num_disulf.find(1) == num_disulf.end() );
		TS_ASSERT( num_disulf.find(3) == num_disulf.end() );
		TS_ASSERT( num_disulf.find(2) != num_disulf.end() );
		TS_ASSERT( num_disulf.find(4) != num_disulf.end() );
	}

	void test_disulfidize() {
		using namespace protocols::denovo_design;
		core::pose::Pose input_pose;
		core::io::pdb::build_pose_from_pdb_as_is( input_pose, "protocols/denovo_design/disulf_test.pdb" );

		// test util.cc: convert_to_poly_ala
		core::pose::PoseOP posecopy = input_pose.clone();
		core::select::residue_selector::ResidueSubset subset( posecopy->size(), true );
		construct_poly_ala_pose( *posecopy, true, subset, subset );
		TS_ASSERT_EQUALS( posecopy->size(), input_pose.size() );

		DisulfidizeMover disulf;
		// check conversion to poly-ala and residue type check function
		for ( core::Size i=1, endi=posecopy->size(); i<=endi; ++i ) {
			TR << "Res " << i << " name " << posecopy->residue(i).name() << std::endl;
			if ( posecopy->residue(i).name3() != "GLY" ) {
				if ( ! posecopy->residue(i).type().is_disulfide_bonded() ) {
					TS_ASSERT_EQUALS( posecopy->residue(i).name3(), "ALA" );
					TS_ASSERT( disulf.check_residue_type( input_pose, i ) );
				}
			} else {
				TS_ASSERT( !disulf.check_residue_type( input_pose, i ) );
			}
		}

		// test seqpos check
		for ( core::Size i=1, endi=input_pose.size(); i<=endi; ++i ) {
			for ( core::Size j=i, endj=input_pose.size(); j<=endj; ++j ) {
				if ( j - i + 1  < 8 ) {
					TS_ASSERT( !disulf.check_disulfide_seqpos( i, j ) );
					TS_ASSERT( !disulf.check_disulfide_seqpos( j, i ) );
				} else {
					TS_ASSERT( disulf.check_disulfide_seqpos( i, j ) );
					TS_ASSERT( disulf.check_disulfide_seqpos( j, i ) );
				}
			}
		}

		// test distance check
		for ( core::Size i=1, endi=input_pose.size(); i<=endi; ++i ) {
			for ( core::Size j=i, endj=input_pose.size(); j<=endj; ++j ) {
				core::Real const dist = input_pose.residue(j).nbr_atom_xyz().distance(
					input_pose.residue(i).nbr_atom_xyz() );
				if ( dist > 5 ) {
					TS_ASSERT( !disulf.check_disulfide_cb_distance( input_pose, i, j ) );
					TS_ASSERT( !disulf.check_disulfide_cb_distance( input_pose, j, i ) );
				} else {
					TS_ASSERT( disulf.check_disulfide_cb_distance( input_pose, i, j ) );
					TS_ASSERT( disulf.check_disulfide_cb_distance( input_pose, j, i ) );
				}
			}
		}

		// test make disulfide
		// 30->47 is a valid disulfide
		protocols::simple_moves::MutateResidue m1(30, "ALA");
		protocols::simple_moves::MutateResidue m2(47, "ALA");
		m1.set_preserve_atom_coords( true );
		m2.set_preserve_atom_coords( true );
		m1.apply( *posecopy );
		m2.apply( *posecopy );
		TS_ASSERT_EQUALS( posecopy->residue(30).name3(), "ALA" );
		TS_ASSERT_EQUALS( posecopy->residue(47).name3(), "ALA" );
		core::scoring::ScoreFunctionOP sfxn_full( core::scoring::get_score_function() );
		disulf.make_disulfide( *posecopy, 30, 47, false, sfxn_full );
		TS_ASSERT_EQUALS( posecopy->residue(30).name(), "CYS:disulfide" );
		TS_ASSERT_EQUALS( posecopy->residue(47).name(), "CYS:disulfide" );

		/// checks disulfide rosetta score
		core::scoring::ScoreFunctionOP sfxn( new core::scoring::ScoreFunction() );
		sfxn->set_weight( core::scoring::dslf_fa13, 1.0 );
		m1.apply( *posecopy );
		m2.apply( *posecopy );
		TS_ASSERT( disulf.check_disulfide_score( *posecopy, 47, 30, sfxn, sfxn_full ) );
		m1.apply( *posecopy );
		m2.apply( *posecopy );
		TS_ASSERT( !disulf.check_disulfide_score( *posecopy, 46, 30, sfxn, sfxn_full ) );

		// check match rt
		core::scoring::disulfides::DisulfideMatchingPotential disulfPot;
		m1.apply( *posecopy );
		m2.apply( *posecopy );
		TS_ASSERT( !disulf.check_disulfide_match_rt( *posecopy, 30, 46, disulfPot, false ) );
		m1.apply( *posecopy );
		m2.apply( *posecopy );
		TS_ASSERT( disulf.check_disulfide_match_rt( *posecopy, 47, 30, disulfPot, false ) );

		core::select::residue_selector::ResidueSubset residueset1( input_pose.size(), true );
		core::select::residue_selector::ResidueSubset residueset2( input_pose.size(), true );
		TS_ASSERT( residueset1[1] );
		TR << residueset2[1] << std::endl;
		DisulfidizeMover::DisulfideList disulfs =
			disulf.find_possible_disulfides( input_pose, residueset1, residueset2, sfxn_full );

		TR << "Pair\tRes1\tRes2" << std::endl;
		for ( core::Size i=1, imax=disulfs.size(); i<=imax; ++i ) {
			TR << i << "\t" << disulfs[i].first << "\t" << disulfs[i].second << std::endl;
		}

#ifdef test_protocols_denovo_design_DisulfidizeMoverTests_dump_pdbs
		input_pose.dump_pdb("disulfidize_temp.pdb");
#endif

		TS_ASSERT_EQUALS( disulfs.size(), 2 );
		TS_ASSERT_EQUALS( disulfs[1].first, 5 );
		TS_ASSERT_EQUALS( disulfs[1].second, 26 );
		TS_ASSERT_EQUALS( disulfs[2].first, 30 );
		TS_ASSERT_EQUALS( disulfs[2].second, 47 );

		DisulfidizeMover::DisulfideList empty_disulfide_list;
		utility::vector1< DisulfidizeMover::DisulfideList > disulfide_configurations =
			disulf.recursive_multiple_disulfide_former( empty_disulfide_list, disulfs );

		TS_ASSERT_EQUALS( disulfide_configurations.size(), 3 );


		// test recursive function
		DisulfidizeMover::DisulfideList mylist, tmplist;
		mylist.push_back( std::make_pair( 1, 2 ) );
		mylist.push_back( std::make_pair( 2, 3 ) );
		mylist.push_back( std::make_pair( 5, 20 ) );
		mylist.push_back( std::make_pair( 21, 30 ) );

		utility::vector1< DisulfidizeMover::DisulfideList > all_combinations=
			disulf.recursive_multiple_disulfide_former( tmplist, mylist );

		TS_ASSERT_EQUALS( all_combinations.size(), 11 );
		TR << "All combinations: " << all_combinations << std::endl;
	}

};
