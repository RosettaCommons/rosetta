// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/hbonds/graph/HBondGraph.cxxtest.hh
/// @brief  test suite for HBondGraph, AtomLevelHBondGraph, AtomInfo
/// @author Jack Maguire, jackmaguire1444@gmail.com

// Test headers
#include <cxxtest/TestSuite.h>

#include <platform/types.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <utility/graph/Graph.hh>

#include <core/conformation/Residue.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/hbonds/HBondGraph_util.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/TenANeighborGraph.hh>

#include <core/scoring/hbonds/graph/HBondGraph.hh>
#include <core/scoring/hbonds/graph/AtomLevelHBondGraph.hh>
#include <core/scoring/hbonds/graph/AtomInfo.hh>

#include <core/types.hh>

#include <numeric/angle.functions.hh>
#include <numeric/xyzVector.hh>

#include <test/UTracer.hh>
#include <test/core/init_util.hh>

//Auto Headers
#include <utility/vector1.hh>


static basic::Tracer TR("core.pack.hbonds.HBondGraph.cxxtest");

class HBondGraphTests : public CxxTest::TestSuite
{
public:
	void setUp()
	{
		core_init();
	}

	void test_hbond_graph(){
		core::pose::Pose testPose;
		core::import_pose::pose_from_file( testPose, "core/pack/dunbrack/1UBQ_repack.pdb" , core::import_pose::PDB_file);

		core::pack::task::PackerTaskOP task = core::pack::task::TaskFactory::create_packer_task( testPose );
		core::scoring::ScoreFunctionOP sfxn = core::scoring::ScoreFunctionFactory::create_score_function( "ref2015.wts" );
		sfxn->score( testPose );

		utility::graph::GraphOP packer_neighbor_graph = core::pack::create_packer_graph( testPose, * sfxn, task );

		core::pack::rotamer_set::RotamerSets rotsets;
		rotsets.set_task( task );
		rotsets.build_rotamers( testPose, *sfxn, packer_neighbor_graph );
		core::Size const nrot = rotsets.nrotamers();

		core::scoring::hbonds::graph::HBondGraph hbg( nrot );
		core::pack::hbonds::init_node_info( hbg, rotsets );
		TS_ASSERT_EQUALS( hbg.num_nodes(), nrot );

		core::scoring::hbonds::graph::AtomLevelHBondGraph alhbg( nrot );
		core::pack::hbonds::init_node_info( alhbg, rotsets );
		TS_ASSERT_EQUALS( alhbg.num_nodes(), nrot );

		for ( core::Size rot = 1; rot <= nrot; ++rot ) {
			core::scoring::hbonds::graph::HBondNode const * hbg_node = dynamic_cast< core::scoring::hbonds::graph::HBondNode const * > ( hbg.get_node( rot ) );
			TS_ASSERT( hbg_node );
			TS_ASSERT_EQUALS( hbg_node->global_rotamer_id(), rot );
			TS_ASSERT_EQUALS( hbg_node->moltenres(), rotsets.moltenres_for_rotamer( rot ) );
			TS_ASSERT_EQUALS( hbg_node->local_rotamer_id(), rot - rotsets.nrotamer_offset_for_moltenres( rotsets.moltenres_for_rotamer( rot ) ) );

			core::scoring::hbonds::graph::HBondNode const * alhbg_node = dynamic_cast< core::scoring::hbonds::graph::HBondNode const * > ( alhbg.get_node( rot ) );
			TS_ASSERT( hbg_node );
			TS_ASSERT_EQUALS( alhbg_node->global_rotamer_id(), rot );
			TS_ASSERT_EQUALS( alhbg_node->moltenres(), rotsets.moltenres_for_rotamer( rot ) );
			TS_ASSERT_EQUALS( alhbg_node->local_rotamer_id(), rot - rotsets.nrotamer_offset_for_moltenres( rotsets.moltenres_for_rotamer( rot ) ) );
		}

		//Starting position for second moltenres
		core::Size const offset = rotsets.nrotamer_offset_for_moltenres( 2 ) + 1;

		hbg.register_hbond( 1, offset, -0.5 );
		hbg.register_hbond( 2, offset, -0.5 );
		hbg.register_hbond( 3, offset, -0.5 );
		hbg.register_hbond( 4, offset, -0.5 );
		hbg.register_hbond( 5, offset + 1, -0.5 );

		alhbg.register_hbond( 1, offset, -0.5 );
		alhbg.register_hbond( 2, offset, -0.5 );
		alhbg.register_hbond( 3, offset, -0.5 );
		alhbg.register_hbond( 4, offset, -0.5 );
		alhbg.register_hbond( 5, offset + 1, -0.5 );

		TS_ASSERT( hbg.find_hbondedge( 1, offset ) );
		TS_ASSERT( hbg.find_hbondedge( 2, offset ) );
		TS_ASSERT( hbg.find_hbondedge( 3, offset ) );
		TS_ASSERT( hbg.find_hbondedge( 4, offset ) );
		TS_ASSERT( hbg.find_hbondedge( 5, offset + 1 ) );

		TS_ASSERT( alhbg.find_hbondedge( 1, offset ) );
		TS_ASSERT( alhbg.find_hbondedge( 2, offset ) );
		TS_ASSERT( alhbg.find_hbondedge( 3, offset ) );
		TS_ASSERT( alhbg.find_hbondedge( 4, offset ) );
		TS_ASSERT( alhbg.find_hbondedge( 5, offset + 1 ) );

		core::pack::hbonds::delete_edges_with_degree_zero( hbg );
		core::pack::hbonds::delete_edges_with_degree_zero( alhbg );

		TS_ASSERT( hbg.find_hbondedge( 1, offset ) );
		TS_ASSERT( hbg.find_hbondedge( 2, offset ) );
		TS_ASSERT( hbg.find_hbondedge( 3, offset ) );
		TS_ASSERT( hbg.find_hbondedge( 4, offset ) );
		TS_ASSERT( ! hbg.find_hbondedge( 5, offset + 1 ) );

		TS_ASSERT( alhbg.find_hbondedge( 1, offset ) );
		TS_ASSERT( alhbg.find_hbondedge( 2, offset ) );
		TS_ASSERT( alhbg.find_hbondedge( 3, offset ) );
		TS_ASSERT( alhbg.find_hbondedge( 4, offset ) );
		TS_ASSERT( ! alhbg.find_hbondedge( 5, offset + 1 ) );

		utility::vector1< bool > buried( testPose.size(), true );
		core::pack::hbonds::determine_atom_level_node_info_for_all_nodes( alhbg, rotsets, buried );

		core::pack::hbonds::determine_atom_level_edge_info_for_all_edges( alhbg, rotsets, *core::scoring::hbonds::HBondDatabase::get_database(), testPose.energies().tenA_neighbor_graph(), testPose );
	}

	void test_AtomInfo(){
		numeric::xyzVector< float > xyz( 1, 2, 3 );
		core::scoring::hbonds::graph::AtomInfo ai(
			17,
			xyz,
			true,
			false,
			true,
			false
		);

		TS_ASSERT_EQUALS( ai.local_atom_id(), 17 );
		TS_ASSERT( ai.is_hydrogen() );
		TS_ASSERT( ! ai.is_donor() );
		TS_ASSERT( ai.is_acceptor() );
		TS_ASSERT( ! ai.is_hydroxyl() );

		ai.is_hydrogen( false );
		TS_ASSERT( ! ai.is_hydrogen() );

		ai.is_donor( true );
		TS_ASSERT( ai.is_donor() );

		ai.is_acceptor( false );
		TS_ASSERT( ! ai.is_acceptor() );

		ai.is_hydroxyl( true );
		TS_ASSERT( ai.is_hydroxyl() );
	}

};
