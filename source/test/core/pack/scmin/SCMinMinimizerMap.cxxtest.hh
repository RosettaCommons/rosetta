// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/scmin/SCMinMinimizerMap.cxxtest.hh
/// @brief  Tests for the AtomTreeCollection classes
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <core/pack/scmin/AtomTreeSCMinMinimizerMap.hh>
#include <core/pack/scmin/AtomTreeCollection.hh>

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

#include <platform/types.hh>

#include <utility/graph/Graph.hh>
#include <core/optimization/DOF_Node.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/packer_neighbors.hh>

//#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>


#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>


#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/types.hh>


#include <test/UTracer.hh>

//Auto Headers
#include <core/kinematics/DomainMap.hh>
#include <utility/vector1.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.pack.rotamer_set.RotamerSet.cxxtest");

using namespace core;

class SCMinMinimizerMapTests : public CxxTest::TestSuite
{

public:
	SCMinMinimizerMapTests() {}

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	void test_construct_SCMinMinimizerMap()
	{
		using namespace conformation;
		using namespace chemical;
		using namespace pack::rotamer_set;
		using namespace pack::scmin;
		using namespace pose;
		using namespace optimization;

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
		utility::graph::GraphOP packer_neighbor_graph = pack::create_packer_graph( pose, *scorefxn, task );

		rotsets->build_rotamers( pose, *scorefxn, packer_neighbor_graph );

		AtomTreeCollectionOP collection( new AtomTreeCollection( pose, *task ) );
		collection->residue_atomtree_collection( 7 ).set_active_restype_index( 20 );
		AtomTreeSCMinMinimizerMap scminmap;
		scminmap.set_total_residue( 20 );
		scminmap.activate_residue_dofs( 7 );

		scminmap.setup( collection );

		/// Having gotten here, we have successfully traversed the atom tree and registered atoms
		/// of the 1-residue atom tree that are controlled by chi-dihedrals and assigned them to
		/// the chi DOF nodes.

		Residue const & res7( collection->residue_atomtree_collection( 7 ).active_residue() );
		for ( Size ii = 1; ii <= scminmap.n_dof_nodes(); ++ii ) {
			DOF_Node const & iinode( scminmap.dof_node( ii ));
			TS_ASSERT( iinode.atoms().size() > 0 );
			for ( Size jj = 1; jj <= iinode.atoms().size(); ++jj ) {
				TS_ASSERT( res7.type().last_controlling_chi( iinode.atoms()[ jj ].atomno() ) == ii );
				TS_ASSERT( iinode.atoms()[ jj ].rsd() == 7 );
			}
			TS_ASSERT( res7.chi_atoms( ii )[ 4 ] == iinode.atoms()[ 1 ].atomno() );
		}

		collection->residue_atomtree_collection( 11 ).set_active_restype_index( 8 );

		scminmap.clear_active_dofs();
		scminmap.activate_residue_dofs( 11 );
		scminmap.setup( collection );

		Residue const & res11( collection->residue_atomtree_collection( 11 ).active_residue() );
		for ( Size ii = 1; ii <= scminmap.n_dof_nodes(); ++ii ) {
			DOF_Node const & iinode( scminmap.dof_node( ii ));
			TS_ASSERT( iinode.atoms().size() > 0 );
			for ( Size jj = 1; jj <= iinode.atoms().size(); ++jj ) {
				TS_ASSERT( res11.type().last_controlling_chi( iinode.atoms()[ jj ].atomno() ) == ii );
				TS_ASSERT( iinode.atoms()[ jj ].rsd() == 11 );
			}
			TS_ASSERT( res11.chi_atoms( ii )[ 4 ] == iinode.atoms()[ 1 ].atomno() );
		}


	}

};
