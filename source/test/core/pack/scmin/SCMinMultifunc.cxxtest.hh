// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/scmin/SCMinMultifunc.cxxtest.hh
/// @brief  Sidechain minimization multifunc class tests
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <core/pack/scmin/SCMinMultifunc.hh>

// Package headers
#include <core/pack/scmin/AtomTreeSCMinMinimizerMap.hh>
#include <core/pack/scmin/AtomTreeCollection.hh>

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

#include <platform/types.hh>

#include <core/graph/Graph.hh>
// AUTO-REMOVED #include <core/kinematics/AtomTree.hh>
// AUTO-REMOVED #include <core/optimization/DOF_Node.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
// AUTO-REMOVED #include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/packer_neighbors.hh>

//#include <core/chemical/ResidueTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Residue.hh>

// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

// AUTO-REMOVED #include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/pose/PDBInfo.hh>

#include <core/scoring/MinimizationGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/types.hh>

// AUTO-REMOVED #include <basic/Tracer.hh>
// AUTO-REMOVED #include <basic/basic.hh>

// AUTO-REMOVED #include <numeric/constants.hh>
// AUTO-REMOVED #include <numeric/random/random.hh>
// AUTO-REMOVED #include <numeric/angle.functions.hh>

#include <test/UTracer.hh>

//Auto Headers
#include <core/kinematics/DomainMap.hh>
#include <utility/vector1.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.pack.scmin.SCMinMultifunc.cxxtest");

using namespace core;

class SCMinMultifuncTests : public CxxTest::TestSuite
{

public:
	SCMinMultifuncTests() {}

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	void test_construct_SCMinMultifunc()
	{
		using namespace conformation;
		using namespace chemical;
		using namespace pack::rotamer_set;
		using namespace pack::scmin;
		using namespace pose;
		using namespace scoring;
		using namespace optimization;
		using namespace graph;

		//typedef utility::vector1< core::conformation::ResidueCOP > ResidueCOPs;
		scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();

		// read in pose
		Pose pose = create_trpcage_ideal_pose();
		(*scorefxn)( pose );
		pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));

		for ( Size ii = 1; ii <= 20; ++ii ) {
			if ( ii == 7 || ii == 11 ) continue;
			task->nonconst_residue_task( ii ).restrict_to_repacking();
		}
		RotamerSetsOP rotsets( new RotamerSets() );
		rotsets->set_task( task );
		graph::GraphOP packer_neighbor_graph = pack::create_packer_graph( pose, *scorefxn, task );

		rotsets->build_rotamers( pose, *scorefxn, packer_neighbor_graph );

		AtomTreeCollectionOP collection( new AtomTreeCollection( pose, *task ) );
		collection->residue_atomtree_collection( 7 ).set_active_restype_index( 21 );
		AtomTreeSCMinMinimizerMap scminmap;
		scminmap.set_total_residue( 20 );
		scminmap.activate_residue_dofs( 7 );

		scminmap.setup( collection );

		utility::vector1< ResidueCOP > bgres( 20 );
		for ( Size ii = 1; ii <= 20; ++ii ) bgres[ ii ] = core::conformation::ResidueCOP( core::conformation::ResidueOP( new Residue( pose.residue( ii ) ) ) );

		MinimizationGraph g( 20 );
		g.copy_connectivity( * packer_neighbor_graph );

		EnergyMap emap; // dummy
		scorefxn->setup_for_minimizing_for_node(
			* g.get_minimization_node( 7 ),
			collection->residue_atomtree_collection( 7 ).active_residue(),
			scminmap, pose, false, emap );
		g.get_minimization_node( 7 )->setup_for_minimizing(
			collection->residue_atomtree_collection( 7 ).active_residue(),
			pose, *scorefxn, scminmap );
		scminmap.set_natoms_for_residue( 7, collection->residue_atomtree_collection( 7 ).active_residue().natoms()  );
		for ( Graph::EdgeListIter
				eiter = g.get_node( 7 )->edge_list_begin(),
				eiter_end = g.get_node( 7 )->edge_list_end();
				eiter != eiter_end; ++eiter ) {
			Size other_ind = (*eiter)->get_other_ind( 7 );
			scorefxn->setup_for_minimizing_for_node(
				* g.get_minimization_node( other_ind ), * bgres[ other_ind ],
				scminmap, pose, false, emap );
			MinimizationEdge & minedge = static_cast< MinimizationEdge & > ( **eiter );
			if ( other_ind < 7 ) {
				scorefxn->setup_for_minimizing_sr2b_enmeths_for_minedge(
					*bgres[ other_ind ], collection->residue_atomtree_collection( 7 ).active_residue(),
					minedge, scminmap, pose, true, false, ( EnergyEdge * ) 0, emap );
				minedge.setup_for_minimizing( *bgres[ other_ind ],
					collection->residue_atomtree_collection( 7 ).active_residue(),
					pose, *scorefxn, scminmap );
			} else {
				scorefxn->setup_for_minimizing_sr2b_enmeths_for_minedge(
					collection->residue_atomtree_collection( 7 ).active_residue(), *bgres[ other_ind ],
					minedge, scminmap, pose, true, false, ( EnergyEdge * ) 0, emap );
				minedge.setup_for_minimizing( collection->residue_atomtree_collection( 7 ).active_residue(),
					*bgres[ other_ind ], pose, *scorefxn, scminmap );
			}
			scminmap.set_natoms_for_residue( other_ind, bgres[ other_ind ]->natoms()  );
		}

		SCMinMultifunc scmf( pose, bgres, *scorefxn, g, scminmap );
		optimization::Multivec chi( 3 );
		//std::cout << "chi: ";
		for ( Size ii = 1; ii <= 3; ++ii ) {
			chi[ ii ] = collection->residue_atomtree_collection( 7 ).active_residue().chi( ii );
		//	std::cout << chi[ ii ] << " ";
		}
		//std::cout << std::endl;
		Real score = scmf( chi );
		//std::cout << "score: " << score << std::endl;
		optimization::Multivec dEdchi( 3 );
		scmf.dfunc( chi, dEdchi );
		//std::cout << "dEdchi: ";
		//for ( Size ii = 1; ii <= 3; ++ii ) {
		//	std::cout << dEdchi[ ii ] << " ";
		//}
		//std::cout << std::endl;
		Real step = 0.1;
		optimization::Multivec chi_step( 3 );
		//std::cout << "chi_stepped: ";
		for ( Size ii = 1; ii <= 3; ++ii ) {
			chi_step[ ii ] = chi[ ii ] - step * dEdchi[ ii ];
		//	std::cout << chi_step[ ii ] << " ";
		}
		//std::cout << std::endl;
		Real score2 = scmf( chi_step );
		TS_ASSERT( score2 < score );
		//std::cout << "score2: " << score2 << std::endl;
		//std::cout << "chi2: ";
		//for ( Size ii = 1; ii <= 3; ++ii ) {
		//	chi[ ii ] = collection->residue_atomtree_collection( 7 ).active_residue().chi( ii );
		//	std::cout << chi[ ii ] << " ";
		//}
		//std::cout << std::endl;
	}

};
