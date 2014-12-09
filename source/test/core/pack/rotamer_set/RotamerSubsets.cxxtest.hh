// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/rotamer_set/RotamerSubsets.cxxtest.hh
/// @brief  test suite for RotamerSubsets
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>


#include <platform/types.hh>

#include <core/graph/Graph.hh>

#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSubsets.hh>
#include <core/pack/packer_neighbors.hh>

#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

// Core headers
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

//Auto Headers
#include <utility/vector0.hh>
#include <utility/vector1.hh>


using namespace core;

class RotamerSubsetsTest : public CxxTest::TestSuite
{
	//chemical::ResidueTypeSetCAP residue_set;

public:

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	///////////////////////////////////////////////////////////////////////////////
	// ------------------------------------------ //
	/// @brief test for rotamer_trials

	void test_rotamer_subset_constructor()
	{
		using namespace conformation;
		using namespace chemical;
		using namespace pack::rotamer_set;
		using namespace pack::task;
		using namespace pose;
		using namespace core::scoring;

		PoseOP trpcage = create_trpcage_ideal_poseop();
		PackerTaskOP task = TaskFactory::create_packer_task( *trpcage );

		utility::vector1< bool > allowed_aas( num_canonical_aas, false );
		allowed_aas[ aa_ala ] = allowed_aas[ aa_phe ] = allowed_aas[ aa_arg ] = true;

		for ( Size ii = 1; ii <= 20; ++ii ) {
			if ( ii == 11 || ii == 12 || ii == 13 ) {
				task->nonconst_residue_task( ii ).restrict_absent_canonical_aas( allowed_aas );
			} else {
				task->nonconst_residue_task( ii ).prevent_repacking();
			}
		}

		ScoreFunctionOP sfxn = get_score_function();
		(*sfxn)( *trpcage ); // score the pose first;

		RotamerSetsOP rotsets( new RotamerSets() );
		rotsets->set_task( task );
		graph::GraphOP packer_neighbor_graph = pack::create_packer_graph( *trpcage, *sfxn, task );
		rotsets->build_rotamers( *trpcage, *sfxn, packer_neighbor_graph );

		//std::cout << "Rotsets built: " << rotsets->nrotamers() << " rotamers.";
		//for ( Size ii = 1; ii <= rotsets->nrotamers(); ++ii ) {
		//	std::cout << "Rotamer " << ii << " mr#: " << rotsets->moltenres_for_rotamer( ii ) << " rot#: "
		//		<< rotsets->rotid_on_moltenresidue( ii ) << std::endl;
		//}

		utility::vector0< Size > rs;
		rs.push_back( 5 );
		rs.push_back( 15 );
		rs.push_back( 25 );
		rs.push_back( 35 );
		rs.push_back( 45 );
		rs.push_back( 55 );

		RotamerSubsetsOP rsubset( new RotamerSubsets( *rotsets, rs ) );
		/*std::cout << "Rotsubset with " << rsubset->nrotamers() << " rotamers" << std::endl;
		for ( Size ii = 1; ii <= rsubset->nrotamers(); ++ii ) {
			std::cout << "Rotamer " << ii << " mr#: " << rsubset->moltenres_for_rotamer( ii ) << " rot#: "
				<< rsubset->rotid_on_moltenresidue( ii ) << std::endl;
			std::cout << "TS_ASSERT( rsubset->moltenres_for_rotamer( " << ii << " ) = " << rsubset->moltenres_for_rotamer( ii ) << ");" << std::endl;
			std::cout << "TS_ASSERT( rsubset->rotid_on_moltenresidue( " << ii << " ) = " << rsubset->rotid_on_moltenresidue( ii ) << ");" << std::endl;
		}*/

		TS_ASSERT( rsubset->moltenres_for_rotamer( 1 ) == 2);
		TS_ASSERT( rsubset->rotid_on_moltenresidue( 1 ) == 1);
		TS_ASSERT( rsubset->moltenres_for_rotamer( 2 ) == 2);
		TS_ASSERT( rsubset->rotid_on_moltenresidue( 2 ) == 2);
		TS_ASSERT( rsubset->moltenres_for_rotamer( 3 ) == 2);
		TS_ASSERT( rsubset->rotid_on_moltenresidue( 3 ) == 3);
		TS_ASSERT( rsubset->moltenres_for_rotamer( 4 ) == 2);
		TS_ASSERT( rsubset->rotid_on_moltenresidue( 4 ) == 4);
		TS_ASSERT( rsubset->moltenres_for_rotamer( 5 ) == 3);
		TS_ASSERT( rsubset->rotid_on_moltenresidue( 5 ) == 1);
		TS_ASSERT( rsubset->moltenres_for_rotamer( 6 ) == 3);
		TS_ASSERT( rsubset->rotid_on_moltenresidue( 6 ) == 2);

		TS_ASSERT( rsubset->rotamer( 1 ) == rotsets->rotamer( 5 ) );
		TS_ASSERT( rsubset->rotamer( 2 ) == rotsets->rotamer( 15 ) );
		TS_ASSERT( rsubset->rotamer( 3 ) == rotsets->rotamer( 25 ) );
		TS_ASSERT( rsubset->rotamer( 4 ) == rotsets->rotamer( 35 ) );
		TS_ASSERT( rsubset->rotamer( 5 ) == rotsets->rotamer( 45 ) );
		TS_ASSERT( rsubset->rotamer( 6 ) == rotsets->rotamer( 55 ) );

	}

};
