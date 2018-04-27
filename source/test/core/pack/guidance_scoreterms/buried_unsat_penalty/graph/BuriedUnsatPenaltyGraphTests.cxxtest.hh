// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/scoring/buried_unsat_penalty/graph/BuriedUnsatPenaltyGraphTests.cxxtest.hh
/// @brief  Unit tests for the BuriedUnsatPenaltyGraph, a graph structure used by the BuriedUnsatPenalty energy.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers
#include <core/pack/guidance_scoreterms/buried_unsat_penalty/graph/BuriedUnsatPenaltyGraph.hh>
#include <core/pack/guidance_scoreterms/buried_unsat_penalty/graph/BuriedUnsatPenaltyGraphOptions.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/chemical/ResidueType.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>
#include <utility/string_util.hh>
#include <utility/graph/Graph.hh>

static basic::Tracer TR("BuriedUnsatPenaltyGraphTests");


class BuriedUnsatPenaltyGraphTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init();

	}

	void tearDown(){

	}

	void test_graph_setup() {
		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, "core/pack/guidance_scoreterms/buried_unsat_penalty/graph/hbnet_testcase1.pdb"  );

		core::pack::guidance_scoreterms::buried_unsat_penalty::graph::BuriedUnsatPenaltyGraphOptionsOP options( new core::pack::guidance_scoreterms::buried_unsat_penalty::graph::BuriedUnsatPenaltyGraphOptions(2.0, 0.25, 1.0, 9.0, 5.0, -0.25) );
		core::scoring::hbonds::HBondOptionsOP hboptions( new core::scoring::hbonds::HBondOptions );
		hboptions->initialize_from_options();
		core::pack::guidance_scoreterms::buried_unsat_penalty::graph::BuriedUnsatPenaltyGraph graph( pose.total_residue(), options, hboptions );
		graph.initialize_graph_for_scoring(pose);

		for ( core::Size i(1), imax(graph.num_nodes()); i<=imax; ++i ) {
			core::pack::guidance_scoreterms::buried_unsat_penalty::graph::BuriedUnsatPenaltyNode const & curnode( *(dynamic_cast<core::pack::guidance_scoreterms::buried_unsat_penalty::graph::BuriedUnsatPenaltyNode*>(graph.get_node(i))) );
			TS_ASSERT_EQUALS( curnode.stored_data_->donor_acceptor_groups_.size(), curnode.stored_data_->donor_acceptor_groups_intrares_hbonds_accepted_.size() );
			TS_ASSERT_EQUALS( curnode.stored_data_->donor_acceptor_groups_.size(), curnode.stored_data_->donor_acceptor_groups_intrares_hbonds_donated_.size() );

			if ( pose.residue_type(i).aa() != core::chemical::aa_pro ) {
				TS_ASSERT_LESS_THAN( 1, curnode.stored_data_->donor_acceptor_groups_.size() ); //Each residue should have at least 1 donor and 1 acceptor group (backbone groups).
			} else {
				TS_ASSERT_EQUALS( 1, curnode.stored_data_->donor_acceptor_groups_.size() ); //Proline has only one acceptor (the O).
				core::pack::guidance_scoreterms::buried_unsat_penalty::graph::BuriedUnsatPenaltyGraphHbondDonorAcceptorGroup const &group1( curnode.stored_data_->donor_acceptor_groups_[1] );
				TS_ASSERT( !group1.is_donor() );
				TS_ASSERT( group1.is_acceptor() );
				TS_ASSERT_EQUALS( utility::strip( pose.residue_type(i).atom_name( group1.heavyatom_index() ) ), "O" );
			}
			if ( pose.residue_type(i).is_polar() ) {
				TS_ASSERT_LESS_THAN(2, curnode.stored_data_->donor_acceptor_groups_.size() ); //Polar residues should have at least 2.
			}
			if ( pose.residue_type(i).aa() == core::chemical::aa_glu || pose.residue_type(i).aa() == core::chemical::aa_asp || pose.residue_type(i).aa() == core::chemical::aa_gln || pose.residue_type(i).aa() == core::chemical::aa_asn || pose.residue_type(i).aa() == core::chemical::aa_his ) {
				TS_ASSERT_EQUALS( 4, curnode.stored_data_->donor_acceptor_groups_.size() );
			} else if ( pose.residue_type(i).aa() == core::chemical::aa_lys || pose.residue_type(i).aa() == core::chemical::aa_thr || pose.residue_type(i).aa() == core::chemical::aa_ser ) {
				TS_ASSERT_EQUALS( 3, curnode.stored_data_->donor_acceptor_groups_.size() );
			} else if ( pose.residue_type(i).aa() == core::chemical::aa_arg ) {
				TS_ASSERT_EQUALS( 5, curnode.stored_data_->donor_acceptor_groups_.size() );
			}
			TR << "\nResidue\tGroup\tHeavyatom\tIs_acceptor\tIs_counted\tN_protons\n";
			for ( core::Size j(1), jmax( curnode.stored_data_->donor_acceptor_groups_.size()); j<=jmax; ++j ) {
				core::pack::guidance_scoreterms::buried_unsat_penalty::graph::BuriedUnsatPenaltyGraphHbondDonorAcceptorGroup const &curgroup( curnode.stored_data_->donor_acceptor_groups_[j] );
				TR << i << "\t" << j << "\t" << utility::strip(pose.residue_type(i).atom_name( curgroup.heavyatom_index() )) << "\t" << (curgroup.is_acceptor() ? "TRUE" : "FALSE") << "\t" << (curgroup.is_counted() ? "TRUE" : "FALSE") << "\t" << curgroup.n_protons() << "\n";
				if ( i == 54 ) {
					if ( j == 1 ) {
						TS_ASSERT_EQUALS( utility::strip(pose.residue_type(i).atom_name( curgroup.heavyatom_index() )), "O");
						TS_ASSERT( curgroup.is_acceptor() );
						TS_ASSERT( curgroup.is_counted() );
						TS_ASSERT_EQUALS(curgroup.n_protons(), 0);
					} else if ( j == 2 ) {
						TS_ASSERT_EQUALS( utility::strip(pose.residue_type(i).atom_name( curgroup.heavyatom_index() )), "OD1");
						TS_ASSERT( curgroup.is_acceptor() );
						TS_ASSERT( curgroup.is_counted() );
						TS_ASSERT_EQUALS(curgroup.n_protons(), 0);
					} else if ( j == 3 ) {
						TS_ASSERT_EQUALS( utility::strip(pose.residue_type(i).atom_name( curgroup.heavyatom_index() )), "N");
						TS_ASSERT( !curgroup.is_acceptor() );
						TS_ASSERT( curgroup.is_counted() );
						TS_ASSERT_EQUALS(curgroup.n_protons(), 1);
					} else if ( j == 3 ) {
						TS_ASSERT_EQUALS( utility::strip(pose.residue_type(i).atom_name( curgroup.heavyatom_index() )), "ND1");
						TS_ASSERT( !curgroup.is_acceptor() );
						TS_ASSERT( curgroup.is_counted() );
						TS_ASSERT_EQUALS(curgroup.n_protons(), 2);
					}
				}
			}
			TR << std::endl;
		}

		//Additional edge tests:
		{ //Tests for residue ASN54
			TR << "\nTESTING NODE 54's CONNECTIONS:\n";
			using namespace core::pack::guidance_scoreterms::buried_unsat_penalty::graph;
			BuriedUnsatPenaltyNode const & curnode( *(dynamic_cast<BuriedUnsatPenaltyNode*>(graph.get_node(54))) );
			TS_ASSERT_EQUALS( curnode.num_edges(), 5 );
			TR << "Edge\tNode1\tNode2\tHBondList\n";
			core::Size count(0);
			for ( utility::graph::EdgeListConstIterator it( curnode.const_edge_list_begin() ); it != curnode.const_edge_list_end(); ++it ) {
				BuriedUnsatPenaltyEdge const &curedge( static_cast< BuriedUnsatPenaltyEdge const & >( **it ) );
				TR << ++count << "\t" << curedge.get_first_node_ind() << "\t" << curedge.get_second_node_ind() << "\t";
				for ( core::Size j(1), jmax(curedge.edge_data_->hbonds_list_.size()); j<=jmax; ++j ) {
					BuriedUnsatPenaltyGraphHbond const & curhbond( curedge.edge_data_->hbonds_list_[j] );
					TR << curhbond.acceptor_group() << "[" << (curhbond.first_node_is_the_acceptor_ ? curedge.get_first_node_ind() : curedge.get_second_node_ind()) << "]-" << curhbond.donor_group() << "[" << (curhbond.first_node_is_the_acceptor_ ? curedge.get_second_node_ind() : curedge.get_first_node_ind()) << "]";
					if ( j<jmax ) TR << ",";
				}
				TR << "\n";

				//In this case, each hydrogen bond is to a unique residue:
				TS_ASSERT( curedge.edge_data_->hbonds_list_.size() == 1 );
				BuriedUnsatPenaltyGraphHbond const &curhbond( curedge.edge_data_->hbonds_list_[1] );
				core::Size const &acc_node_index( curhbond.first_node_is_the_acceptor() ? curedge.get_first_node_ind() : curedge.get_second_node_ind() );
				core::Size const &don_node_index( curhbond.first_node_is_the_acceptor() ? curedge.get_second_node_ind() : curedge.get_first_node_ind() );
				BuriedUnsatPenaltyNode const & acc_node( *(dynamic_cast<BuriedUnsatPenaltyNode const *>( graph.get_node(acc_node_index) ) ) );
				BuriedUnsatPenaltyNode const & don_node( *(dynamic_cast<BuriedUnsatPenaltyNode const *>( graph.get_node(don_node_index) ) ) );
				BuriedUnsatPenaltyGraphHbondDonorAcceptorGroup const & acc_group( acc_node.donor_acceptor_group( curhbond.acceptor_group() ) );
				BuriedUnsatPenaltyGraphHbondDonorAcceptorGroup const & don_group( don_node.donor_acceptor_group( curhbond.donor_group() ) );

				if ( count == 1 ) {
					TS_ASSERT( don_node_index == 54 );
					TS_ASSERT( acc_node_index == 50 );
					TS_ASSERT_EQUALS( utility::strip( pose.residue_type(54).atom_name( don_group.heavyatom_index() ) ), "N" );
					TS_ASSERT_EQUALS( utility::strip( pose.residue_type(50).atom_name( acc_group.heavyatom_index() ) ), "O" );
				} else if ( count == 2 ) {
					TS_ASSERT( don_node_index == 58 );
					TS_ASSERT( acc_node_index == 54 );
					TS_ASSERT_EQUALS( utility::strip( pose.residue_type(58).atom_name( don_group.heavyatom_index() ) ), "N" );
					TS_ASSERT_EQUALS( utility::strip( pose.residue_type(54).atom_name( acc_group.heavyatom_index() ) ), "O" );
				} else if ( count == 3 ) {
					TS_ASSERT( don_node_index == 54 );
					TS_ASSERT( acc_node_index == 128 );
					TS_ASSERT_EQUALS( utility::strip( pose.residue_type(54 ).atom_name( don_group.heavyatom_index() ) ), "ND2" );
					TS_ASSERT_EQUALS( utility::strip( pose.residue_type(128).atom_name( acc_group.heavyatom_index() ) ), "ND1" );
				} else if ( count == 4 ) {
					TS_ASSERT( don_node_index == 54 );
					TS_ASSERT( acc_node_index == 129 );
					TS_ASSERT_EQUALS( utility::strip( pose.residue_type(54 ).atom_name( don_group.heavyatom_index() ) ), "ND2" );
					TS_ASSERT_EQUALS( utility::strip( pose.residue_type(129).atom_name( acc_group.heavyatom_index() ) ), "OD1" );
				} else if ( count == 5 ) {
					TS_ASSERT( don_node_index == 204 );
					TS_ASSERT( acc_node_index == 54 );
					TS_ASSERT_EQUALS( utility::strip( pose.residue_type(204).atom_name( don_group.heavyatom_index() ) ), "ND2" );
					TS_ASSERT_EQUALS( utility::strip( pose.residue_type(54 ).atom_name( acc_group.heavyatom_index() ) ), "OD1" );
				}
			}
			TR << std::endl;
		}
	}

	void test_smallpose_graph_setup() {
		//core::pose::Pose bigpose;
		//core::import_pose::pose_from_file( bigpose, "core/pack/guidance_scoreterms/buried_unsat_penalty/graph/hbnet_testcase1.pdb"  );

		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, "core/pack/guidance_scoreterms/buried_unsat_penalty/graph/hbnet_testcase2.pdb"  );

		core::pack::guidance_scoreterms::buried_unsat_penalty::graph::BuriedUnsatPenaltyGraphOptionsOP options( new core::pack::guidance_scoreterms::buried_unsat_penalty::graph::BuriedUnsatPenaltyGraphOptions(2.0, 0.25, 1.0, 9.0, 5.0, -0.25) );
		core::scoring::hbonds::HBondOptionsOP hboptions( new core::scoring::hbonds::HBondOptions );
		hboptions->initialize_from_options();
		core::scoring::hbonds::HBondDatabase const &hbdatabase( *(core::scoring::hbonds::HBondDatabase::get_database( hboptions->params_database_tag() )) );

		//core::scoring::hbonds::HBondSet hbondset1( *hboptions, bigpose, false );
		//core::scoring::hbonds::HBondSet hbondset2( *hboptions, bigpose.residue(204), bigpose.residue(54), hbdatabase );
		core::scoring::hbonds::HBondSet hbondset3( *hboptions, pose, false );
		core::scoring::hbonds::HBondSet hbondset4( *hboptions, pose.residue(2), pose.residue(5), hbdatabase );

		/*TR << "\nCONTENTS OF HBONDSET1:\n";
		hbondset1.show( TR );
		TR << std::endl;
		TR << "\nCONTENTS OF HBONDSET2:\n";
		hbondset2.show( TR );
		TR << std::endl;
		*/
		TR << "\nCONTENTS OF HBONDSET3:\n";
		hbondset3.show( TR );
		TR << std::endl;
		TR << "\nCONTENTS OF HBONDSET4:\n";
		hbondset4.show( TR );
		TR << std::endl;

		core::pack::guidance_scoreterms::buried_unsat_penalty::graph::BuriedUnsatPenaltyGraph graph( pose.total_residue(), options, hboptions );
		graph.initialize_graph_for_scoring(pose);

		{ //Tests for residue ASN5
			TR << "\nTESTING NODE 5's CONNECTIONS:\n";
			using namespace core::pack::guidance_scoreterms::buried_unsat_penalty::graph;
			BuriedUnsatPenaltyNode const & curnode( *(dynamic_cast<BuriedUnsatPenaltyNode*>(graph.get_node(5))) );
			TS_ASSERT_EQUALS( curnode.num_edges(), 2 );
			TR << "Edge\tNode1\tNode2\tHBondList\n";
			core::Size count(0);
			for ( utility::graph::EdgeListConstIterator it( curnode.const_edge_list_begin() ); it != curnode.const_edge_list_end(); ++it ) {
				BuriedUnsatPenaltyEdge const &curedge( static_cast< BuriedUnsatPenaltyEdge const & >( **it ) );
				TR << ++count << "\t" << curedge.get_first_node_ind() << "\t" << curedge.get_second_node_ind() << "\t";
				for ( core::Size j(1), jmax(curedge.edge_data_->hbonds_list_.size()); j<=jmax; ++j ) {
					BuriedUnsatPenaltyGraphHbond const & curhbond( curedge.edge_data_->hbonds_list_[j] );
					TR << curhbond.acceptor_group() << "[" << (curhbond.first_node_is_the_acceptor_ ? curedge.get_first_node_ind() : curedge.get_second_node_ind()) << "]-" << curhbond.donor_group() << "[" << (curhbond.first_node_is_the_acceptor_ ? curedge.get_second_node_ind() : curedge.get_first_node_ind()) << "]";
					if ( j<jmax ) TR << ",";
				}
				TR << "\n";

				//In this case, each hydrogen bond is to a unique residue:
				TS_ASSERT( curedge.edge_data_->hbonds_list_.size() == 1 );
				BuriedUnsatPenaltyGraphHbond const &curhbond( curedge.edge_data_->hbonds_list_[1] );
				core::Size const &acc_node_index( curhbond.first_node_is_the_acceptor() ? curedge.get_first_node_ind() : curedge.get_second_node_ind() );
				core::Size const &don_node_index( curhbond.first_node_is_the_acceptor() ? curedge.get_second_node_ind() : curedge.get_first_node_ind() );
				BuriedUnsatPenaltyNode const & acc_node( *(dynamic_cast<BuriedUnsatPenaltyNode const *>( graph.get_node(acc_node_index) ) ) );
				BuriedUnsatPenaltyNode const & don_node( *(dynamic_cast<BuriedUnsatPenaltyNode const *>( graph.get_node(don_node_index) ) ) );
				BuriedUnsatPenaltyGraphHbondDonorAcceptorGroup const & acc_group( acc_node.donor_acceptor_group( curhbond.acceptor_group() ) );
				BuriedUnsatPenaltyGraphHbondDonorAcceptorGroup const & don_group( don_node.donor_acceptor_group( curhbond.donor_group() ) );

				if ( count == 1 ) {
					TS_ASSERT( don_node_index == 2 );
					TS_ASSERT( acc_node_index == 5 );
					TS_ASSERT_EQUALS( utility::strip( pose.residue_type(2).atom_name( don_group.heavyatom_index() ) ), "ND2" );
					TS_ASSERT_EQUALS( utility::strip( pose.residue_type(5).atom_name( acc_group.heavyatom_index() ) ), "OD1" );
				} else if ( count == 2 ) {
					TS_ASSERT( don_node_index == 5 );
					TS_ASSERT( acc_node_index == 8 );
					TS_ASSERT_EQUALS( utility::strip( pose.residue_type(5).atom_name( don_group.heavyatom_index() ) ), "ND2" );
					TS_ASSERT_EQUALS( utility::strip( pose.residue_type(8).atom_name( acc_group.heavyatom_index() ) ), "OD1" );
				}
			}
			TR << std::endl;
		}

	}

};
