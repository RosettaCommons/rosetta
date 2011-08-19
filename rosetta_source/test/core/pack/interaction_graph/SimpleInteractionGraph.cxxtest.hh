// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/interaction_graph/FASTERInteractionGraph.cxxtest.hh
/// @brief  test suite for the FASTER interaction graph
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Test framework headers
#include <cxxtest/TestSuite.h>

// Core Headers
#include <core/pack/interaction_graph/SimpleInteractionGraph.hh>

#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pack/packer_neighbors.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>


// Test headers
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

// C++ headers
#include <iostream>

class SimpleInteractionGraphTests : public CxxTest::TestSuite {
public:

	void setUp() {
		core_init();
	}

	void test_instantiate_simple_ig() {

		using namespace core::chemical;
		using namespace core::graph;
		using namespace core::pack;
		using namespace core::pack::interaction_graph;
		using namespace core::pack::rotamer_set;
		using namespace core::pack::task;
		using namespace core::pose;
		using namespace core::scoring;
		using core::Size;

		PoseOP trpcage = create_trpcage_ideal_poseop();
		PackerTaskOP task = TaskFactory::create_packer_task( *trpcage );
		//task->or_preserve_c_beta( true );

		utility::vector1< bool > allowed_aas( num_canonical_aas, false );
		allowed_aas[ aa_ala ] = allowed_aas[ aa_phe ] = allowed_aas[ aa_arg ] = true;

		for ( Size ii = 1; ii <= 20; ++ii ) {
			if ( ii == 11 || ii == 12 || ii == 13 ) {
				task->nonconst_residue_task( ii ).restrict_absent_canonical_aas( allowed_aas );
			} else {
				task->nonconst_residue_task( ii ).prevent_repacking();
			}
		}

		ScoreFunctionOP sfxn = ScoreFunctionFactory::create_score_function( "standard" );
		sfxn->set_weight( fa_pair, 0.0 );
		sfxn->set_weight( hbond_sc, 0.0 );
		sfxn->set_weight( hbond_bb_sc, 0.0 );

		core::Real startscore = (*sfxn)( *trpcage ); // score the pose first;
		sfxn->setup_for_packing( *trpcage, task->repacking_residues(), task->designing_residues() );
		GraphOP packer_neighbor_graph = create_packer_graph( *trpcage, *sfxn, task );

		RotamerSetsOP rotsets( new RotamerSets() );
		rotsets->set_task( task );
		rotsets->build_rotamers( *trpcage, *sfxn, packer_neighbor_graph );
		rotsets->prepare_sets_for_packing( *trpcage, *sfxn );

		//trpcage->dump_pdb( "test_trpcage_1.pdb" );

		SimpleInteractionGraphOP simple_ig = new SimpleInteractionGraph;
		simple_ig->set_scorefunction( *sfxn );
		simple_ig->initialize( *trpcage );

		/*std::cout << "Initial pose one-body energy: " << sfxn->weights().dot( trpcage->energies().onebody_energies( 12 )) << std::endl;
		for ( Node::EdgeListConstIter iter = trpcage->energies().energy_graph().get_node( 12 )->edge_list_begin(),
				iter_end = trpcage->energies().energy_graph().get_node( 12 )->edge_list_end(); iter != iter_end; ++iter ) {
			std::cout << "Energy Edge to " << (*iter)->get_other_ind( 12 ) << " ";
			EnergyEdge const * eedge = static_cast< EnergyEdge const * > (*iter);
			EnergyMap emap = eedge->fill_energy_map();
			//emap.show_weighted( std::cout, sfxn->weights() );
			
			std::cout << sfxn->weights().dot( emap ) << std::endl;
		}
		std::cout << "Initial one-body energy: " << static_cast< SimpleNode const * > (simple_ig->get_node( 12 ))->current_one_body_energy() << std::endl;
		for ( Node::EdgeListConstIter iter = simple_ig->get_node( 12 )->edge_list_begin(),
				iter_end = simple_ig->get_node( 12 )->edge_list_end(); iter != iter_end; ++iter ) {
			SimpleEdge const * simple_edge = static_cast< SimpleEdge const * > ( *iter );
			std::cout << "Curr energy: " << simple_edge->get_other_ind( 12 ) << " " << simple_edge->get_current_energy() << std::endl;
		}*/


		core::conformation::ResidueCOP res12rot4 = rotsets->rotamer_set_for_residue( 12 )->rotamer( 4 );

		core::Real negdelta = simple_ig->consider_substitution( 12, res12rot4 );
		trpcage->replace_residue( 12, *res12rot4, false );
		core::Real score_after_res12rot4 = (*sfxn)(*trpcage);
		//trpcage->dump_pdb( "test_trpcage_2.pdb" );

		TS_ASSERT_DELTA( negdelta, startscore - score_after_res12rot4, 1e-6 );

		/*std::cout << "Afterwards pose one-body energy: " << sfxn->weights().dot( trpcage->energies().onebody_energies( 12 )) << std::endl;
		for ( Node::EdgeListConstIter iter = trpcage->energies().energy_graph().get_node( 12 )->edge_list_begin(),
				iter_end = trpcage->energies().energy_graph().get_node( 12 )->edge_list_end(); iter != iter_end; ++iter ) {
			std::cout << "Energy Edge to " << (*iter)->get_other_ind( 12 ) << " ";
			EnergyEdge const * eedge = static_cast< EnergyEdge const * > (*iter);
			EnergyMap emap = eedge->fill_energy_map();
			//emap.show_weighted( std::cout, sfxn->weights() );
			
			std::cout << sfxn->weights().dot( emap ) << std::endl;
		}
		std::cout << "Alternate one-body energy: " << static_cast< SimpleNode const * > (simple_ig->get_node( 12 ))->proposed_one_body_energy() << std::endl;
		for ( Node::EdgeListConstIter iter = simple_ig->get_node( 12 )->edge_list_begin(),
				iter_end = simple_ig->get_node( 12 )->edge_list_end(); iter != iter_end; ++iter ) {
			SimpleEdge const * simple_edge = static_cast< SimpleEdge const * > ( *iter );
			std::cout << "Alt energy: " << simple_edge->get_other_ind( 12 ) << " " << simple_edge->get_proposed_energy() << std::endl;
		}*/
	}

	void test_simple_ig_neutral_compare() {
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::conformation;
		using namespace core::pack::interaction_graph;

		core::pose::Pose  oneten(create_1ten_pdb_pose());
		ScoreFunctionOP sfxn = ScoreFunctionFactory::create_score_function( "standard" );
		sfxn->score(oneten);

		SimpleInteractionGraphOP simple_ig = new SimpleInteractionGraph;
		simple_ig->set_scorefunction( *sfxn );
		simple_ig->initialize( oneten );

		for (core::Size ii(1); ii <= oneten.n_residue(); ++ii) {
			ResidueCOP res(new Residue(oneten.residue(ii)));
			TS_ASSERT_DELTA( simple_ig->consider_substitution(ii, res), 0, 0.0001);
			simple_ig->reject_change( ii );
		}
	}
};
