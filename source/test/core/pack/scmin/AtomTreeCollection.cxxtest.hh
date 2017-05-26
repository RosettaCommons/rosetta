// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/scmin/AtomTreeCollection.cxxtest.hh
/// @brief  Tests for the AtomTreeCollection classes
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <core/pack/scmin/AtomTreeCollection.hh>

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

#include <platform/types.hh>

#include <utility/graph/Graph.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/packer_neighbors.hh>

//#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>


#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>


#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/types.hh>


#include <numeric/constants.hh>
#include <numeric/angle.functions.hh>

#include <test/UTracer.hh>

//Auto Headers
#include <utility/vector1.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.pack.rotamer_set.RotamerSet.cxxtest");

using namespace core;

class AtomTreeCollectionTests : public CxxTest::TestSuite
{

public:
	AtomTreeCollectionTests() {}

	// Shared initialization goes here.
	void setUp() {
		core_init();
		//core_init_with_additional_options();
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	void test_construct_residue_atom_tree_collection()
	{
		using namespace conformation;
		using namespace chemical;
		using namespace pack::rotamer_set;
		using namespace pack::scmin;
		using namespace pose;

		//typedef utility::vector1< core::conformation::ResidueCOP > ResidueCOPs;


		scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();

		// read in pose
		Pose pose = create_trpcage_ideal_pose();
		(*scorefxn)( pose );
		pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));

		for ( Size ii = 1; ii <= 20; ++ii ) {
			if ( ii == 7 ) continue;
			task->nonconst_residue_task( ii ).restrict_to_repacking();
		}

		RotamerSetsOP rotsets( new RotamerSets() );
		rotsets->set_task( task );
		utility::graph::GraphOP packer_neighbor_graph = pack::create_packer_graph( pose, *scorefxn, task );

		rotsets->build_rotamers( pose, *scorefxn, packer_neighbor_graph );

		// </stolen code>
		// OK -- lets build an AtomTreeCollection
		RotamerSetCOP r7rotset = rotsets->rotamer_set_for_residue( 7 );
		TS_ASSERT( r7rotset->get_n_residue_types() == 21 ); // 20 + 1 extra his restype.
		//std::cout << "MINE N tyr rots; " << r7rotset->get_n_rotamers_for_residue_type( 20 ) << std::endl;
		TS_ASSERT( r7rotset->get_n_rotamers_for_residue_type( 20 ) == 9 ); // this test will crash below if nrots < 4
		ResidueAtomTreeCollection ratc( task->residue_task( 7 ), pose.conformation(), pose.residue( 7 ) );
		ratc.set_active_restype_index( 20 );
		{ // scope
			kinematics::AtomTree const & tree( ratc.active_atom_tree() );
			Residue const & res( ratc.active_residue() );
			for ( Size ii = 1; ii <= res.nchi(); ++ii ) {
				TS_ASSERT_DELTA(
					numeric::principal_angle_radians( numeric::constants::d::deg2rad * res.chi( ii ) ),
					numeric::principal_angle_radians( tree.dof( id::DOF_ID( id::AtomID( res.chi_atoms(ii)[4], 1 ), id::PHI ))),
					1e-12 );
			}
		}

		// now, lets change shit up
		ratc.set_rescoords( *r7rotset->rotamer( r7rotset->get_residue_type_begin(20) + 3 ) );
		ratc.update_atom_tree();

		{ // scope
			kinematics::AtomTree const & tree( ratc.active_atom_tree() );
			Residue const & res( ratc.active_residue() );
			for ( Size ii = 1; ii <= res.nchi(); ++ii ) {
				TS_ASSERT_DELTA(
					numeric::principal_angle_radians( numeric::constants::d::deg2rad * res.chi( ii ) ),
					numeric::principal_angle_radians( tree.dof( id::DOF_ID( id::AtomID( res.chi_atoms(ii)[4], 1 ), id::PHI ))),
					1e-12 );
			}
		}

		ratc.set_chi( 2, 55.0 );
		ratc.update_residue();
		{ // scope
			kinematics::AtomTree const & tree( ratc.active_atom_tree() );
			Residue const & res( ratc.active_residue() );
			for ( Size ii = 1; ii <= res.nchi(); ++ii ) {
				TS_ASSERT_DELTA(
					numeric::principal_angle_radians( numeric::constants::d::deg2rad * res.chi( ii ) ),
					numeric::principal_angle_radians( tree.dof( id::DOF_ID( id::AtomID( res.chi_atoms(ii)[4], 1 ), id::PHI ))),
					1e-12 );
			}
		}
	}

	void test_atom_tree_collection_from_rotamer_sets_ctor()
	{
		using namespace conformation;
		using namespace chemical;
		using namespace pack::rotamer_set;
		using namespace pack::scmin;
		using namespace pose;

		scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();

		// read in pose
		Pose pose = create_trpcage_ideal_pose();
		(*scorefxn)( pose );
		pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));

		for ( Size ii = 1; ii <= 20; ++ii ) {
			if ( ii == 7 || ii == 13 ) continue;
			task->nonconst_residue_task( ii ).restrict_to_repacking();
		}

		RotamerSetsOP rotsets( new RotamerSets() );
		rotsets->set_task( task );
		utility::graph::GraphOP packer_neighbor_graph = pack::create_packer_graph( pose, *scorefxn, task );

		rotsets->build_rotamers( pose, *scorefxn, packer_neighbor_graph );

		AtomTreeCollection collection( pose, *task );
		TS_ASSERT( collection.residue_atomtree_collection_op( 7 ) );
		TS_ASSERT( collection.residue_atomtree_collection_op( 13 ) );
		TS_ASSERT( collection.residue_atomtree_collection_op( 7 ) != collection.residue_atomtree_collection_op( 13 ) );

	}

	void test_atom_tree_collection_from_rotamer_set_ctor()
	{
		using namespace conformation;
		using namespace chemical;
		using namespace pack::rotamer_set;
		using namespace pack::scmin;
		using namespace pose;

		scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();

		// read in pose
		Pose pose = create_trpcage_ideal_pose();
		(*scorefxn)( pose );
		pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));

		for ( Size ii = 1; ii <= 20; ++ii ) {
			if ( ii == 7 || ii == 13 ) continue;
			task->nonconst_residue_task( ii ).restrict_to_repacking();
		}

		RotamerSetFactory rsf;
		RotamerSetOP rotset = rsf.create_rotamer_set( pose.residue( 7 ) );
		rotset->set_resid( 7 );
		utility::graph::GraphOP packer_neighbor_graph = pack::create_packer_graph( pose, *scorefxn, task );
		rotset->build_rotamers( pose, *scorefxn, *task, packer_neighbor_graph, true );

		AtomTreeCollection collection( pose, task->residue_task( 7 ), 7 );
		TS_ASSERT( collection.residue_atomtree_collection_op( 7 ) );

	}

};
