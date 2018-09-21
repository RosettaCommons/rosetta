// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/task_operations/PruneBuriedUnsatsOperation.cxxtest.hh
/// @brief  test suite for PruneBuriedUnsatsOperation
/// @author Brian Coventry (bcov@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>

#include <platform/types.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <core/conformation/Residue.hh>
#include <core/id/AtomID_Map.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/interaction_graph/AnnealableGraphBase.hh>
#include <core/pack/task/operation/OperateOnResidueSubset.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/scoring/atomic_depth/AtomicDepth.hh>
#include <core/scoring/atomic_depth/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>

#include <protocols/task_operations/PruneBuriedUnsatsOperation.hh>

#include <core/types.hh>

#include <test/UTracer.hh>
#include <test/core/init_util.hh>

//Auto Headers
#include <utility/vector1.hh>
#include <utility/pointer/memory.hh>

using namespace core;


static basic::Tracer TR("protocols.task_operations.PruneBuriedUnsatsOperation.cxxtest");

class PruneBuriedUnsatsOperationTests : public CxxTest::TestSuite
{
public:
	void setUp()
	{
		core_init();
	}


	const float burial_probe_radius = 0.1f;
	const float burial_depth = 0.1f;
	const float burial_resolution = 0.25f;

	// If this test fails, we don't expect anything else to pass. But the issue is that the
	//  atoms aren't being considered buried.
	void test_burial_sanity_check() {
		{
			pose::Pose pose;
			import_pose::pose_from_file( pose, "protocols/task_operations/arg_and_two_ser.pdb" );

			core::id::AtomID_Map< bool > buried = scoring::atomic_depth::atoms_deeper_than(
				pose, burial_depth, false, burial_probe_radius, true, burial_resolution );

			TS_ASSERT( buried( 1, pose.residue(1).atom_index("NH1") ) );
			TS_ASSERT( buried( 1, pose.residue(1).atom_index("NE") ) );
			TS_ASSERT( buried( 1, pose.residue(1).atom_index("NH2") ) );
			TS_ASSERT( buried( 2, pose.residue(2).atom_index("OG") ) );
			TS_ASSERT( buried( 3, pose.residue(3).atom_index("OG") ) );
		}
		{
			pose::Pose pose;
			import_pose::pose_from_file( pose, "protocols/task_operations/asn_to_bb.pdb" );

			core::id::AtomID_Map< bool > buried = scoring::atomic_depth::atoms_deeper_than(
				pose, burial_depth, false, burial_probe_radius, true, burial_resolution );

			TS_ASSERT( buried( 1, pose.residue(1).atom_index("ND2") ) );
			TS_ASSERT( buried( 1, pose.residue(1).atom_index("OD1") ) );
			TS_ASSERT( buried( 3, pose.residue(3).atom_index("O") ) );
		}
		{
			pose::Pose pose;
			import_pose::pose_from_file( pose, "protocols/task_operations/asn_to_satisfied_bb.pdb" );

			core::id::AtomID_Map< bool > buried = scoring::atomic_depth::atoms_deeper_than(
				pose, burial_depth, false, burial_probe_radius, true, burial_resolution );

			TS_ASSERT( buried( 1, pose.residue(1).atom_index("ND2") ) );
			TS_ASSERT( buried( 1, pose.residue(1).atom_index("OD1") ) );
			TS_ASSERT( buried( 3, pose.residue(3).atom_index("O") ) );
			TS_ASSERT( buried( 6, pose.residue(6).atom_index("N") ) );
		}
		{
			pose::Pose pose;
			import_pose::pose_from_file( pose, "protocols/task_operations/asn_to_satisfied_ser.pdb" );

			core::id::AtomID_Map< bool > buried = scoring::atomic_depth::atoms_deeper_than(
				pose, burial_depth, false, burial_probe_radius, true, burial_resolution );

			TS_ASSERT( buried( 1, pose.residue(1).atom_index("ND2") ) );
			TS_ASSERT( buried( 1, pose.residue(1).atom_index("OD1") ) );
			TS_ASSERT( buried( 2, pose.residue(2).atom_index("OG") ) );
			TS_ASSERT( buried( 3, pose.residue(3).atom_index("OG1") ) );
		}
		{
			pose::Pose pose;
			import_pose::pose_from_file( pose, "protocols/task_operations/asn_to_ser.pdb" );

			core::id::AtomID_Map< bool > buried = scoring::atomic_depth::atoms_deeper_than(
				pose, burial_depth, false, burial_probe_radius, true, burial_resolution );

			TS_ASSERT( buried( 1, pose.residue(1).atom_index("ND2") ) );
			TS_ASSERT( buried( 1, pose.residue(1).atom_index("OD1") ) );
			TS_ASSERT( buried( 2, pose.residue(2).atom_index("OG") ) );
		}

	}


	pack::rotamer_set::RotamerSetsOP
	make_rotsets( core::pose::Pose & pose,
		scoring::ScoreFunctionOP const & sfxn,
		std::string const & allow_repacking_res,
		utility::vector1<std::pair<std::string, std::string> > const & designables,
		bool allow_even_trades,
		bool add_taskop ) {

		using utility::pointer::make_shared;
		using namespace core::pack::task::operation;
		using namespace core::select::residue_selector;

		core::pack::task::TaskFactoryOP tf( new core::pack::task::TaskFactory() );
		tf->push_back( TaskOperationCOP( new IncludeCurrent() ) );
		tf->push_back( make_shared< OperateOnResidueSubset >(
			make_shared< PreventRepackingRLT >(), make_shared< ResidueIndexSelector >( allow_repacking_res ), true ) );

		if ( add_taskop ) {
			protocols::task_operations::PruneBuriedUnsatsOperationOP prune( new protocols::task_operations::PruneBuriedUnsatsOperation );
			prune->allow_even_trades( allow_even_trades );
			prune->atomic_depth_probe_radius( burial_probe_radius );
			prune->atomic_depth_resolution( burial_resolution );
			prune->atomic_depth_cutoff( burial_depth );

			tf->push_back( prune );
		}

		for ( std::pair<std::string, std::string> const & pair : designables ) {
			std::string const & res_string = pair.first;
			std::string const & aa_string = pair.second;

			RestrictAbsentCanonicalAASRLTOP only_XX( new RestrictAbsentCanonicalAASRLT() );
			only_XX->aas_to_keep(aa_string);

			tf->push_back( make_shared< OperateOnResidueSubset >(
				only_XX, make_shared< ResidueIndexSelector >( res_string ), false ) );

		}

		pack::task::PackerTaskOP task = tf->create_task_and_apply_taskoperations( pose );

		pack::rotamer_set::RotamerSetsOP rotsets( new pack::rotamer_set::RotamerSets() );
		pack::interaction_graph::AnnealableGraphBaseOP ig = nullptr;

		pack::pack_rotamers_setup( pose, *sfxn, task, rotsets, ig );

		return rotsets;
	}


	// Testing function to start a packing trajectory and then makes sure the rotamer sets
	// looks how we expect it to
	void do_test_pruning(
		pose::Pose & pose,
		std::string const & allow_repacking_res,
		utility::vector1<std::pair<std::string, std::string> > const & designables,
		bool allow_even_trades,
		utility::vector1< Size > const & our_res,
		utility::vector1< bool > const & expect_current,
		utility::vector1< std::set<char> > const & expect_aa
	) {



		scoring::ScoreFunctionOP sfxn ( new scoring::ScoreFunction() );
		sfxn->set_weight( scoring::hbond_sc, 1.0f );
		sfxn->set_weight( scoring::hbond_bb_sc, 1.0f );
		sfxn->set_weight( scoring::hbond_lr_bb, 1.0f );
		sfxn->set_weight( scoring::hbond_sr_bb, 1.0f );

		pack::rotamer_set::RotamerSetsOP rotsets_no_task = make_rotsets( pose, sfxn, allow_repacking_res, designables, allow_even_trades, false );
		pack::rotamer_set::RotamerSetsOP rotsets = make_rotsets( pose, sfxn, allow_repacking_res, designables, allow_even_trades, true );

		TS_ASSERT( rotsets->nrotamers() < rotsets_no_task->nrotamers() );

		// Now lets make sure that the RotamerSets look the way we expect them to

		for ( Size iour = 1; iour <= our_res.size(); iour++ ) {
			Size this_resnum = our_res[ iour ];
			bool this_expect_current = expect_current[ iour ];
			std::set< char > this_expect_aa = expect_aa[ iour ];

			pack::rotamer_set::RotamerSetCOP rotset = rotsets->rotamer_set_for_residue( this_resnum );

			conformation::Residue const & match_res = pose.residue( this_resnum );

			bool found_current = false;
			std::set< char > found_aa;

			for ( Size irot = 1; irot <= rotset->num_rotamers(); irot++ ) {
				conformation::ResidueCOP rotamer = rotset->rotamer( irot );
				found_aa.insert( rotamer->name1() );

				if ( match_res.name1() != rotamer->name1() ) continue;

				bool chi_match = true;
				for ( Size chi = 1; chi < match_res.nchi(); chi++ ) {
					if ( std::abs( match_res.chi( chi ) - rotamer->chi( chi ) ) > 1 ) {
						chi_match = false;
						break;
					}
				}

				if ( chi_match ) found_current = true;

			}

			TS_ASSERT( found_current == this_expect_current );
			TS_ASSERT( this_expect_aa == found_aa );
		}

	}

	// Testing whether the math functions work. 2 good atoms and 1 bad atom = keep
	void test_pruning_arg_and_two_ser(){

		pose::Pose pose;
		import_pose::pose_from_file( pose, "protocols/task_operations/arg_and_two_ser.pdb" );



		std::string allow_repacking_res = "1";

		utility::vector1<std::pair<std::string, std::string> > designables { { "1", "RA" } };

		utility::vector1<Size> our_res { 1 };
		utility::vector1<bool> expect_current { true };
		utility::vector1< std::set<char> > expect_aa { { 'A', 'R' } };



		do_test_pruning( pose, allow_repacking_res, designables, false, our_res, expect_current, expect_aa);

	}

	// Testing even trades on sc atoms
	void test_pruning_asn_to_ser(){

		pose::Pose pose;
		import_pose::pose_from_file( pose, "protocols/task_operations/asn_to_ser.pdb" );

		{
			std::string allow_repacking_res = "1";

			utility::vector1<std::pair<std::string, std::string> > designables { { "1", "NA" } };

			utility::vector1<Size> our_res { 1 };
			utility::vector1<bool> expect_current { true };
			utility::vector1< std::set<char> > expect_aa { { 'A', 'N' } };



			do_test_pruning( pose, allow_repacking_res, designables, true, our_res, expect_current, expect_aa);
		}
		{
			std::string allow_repacking_res = "1";

			utility::vector1<std::pair<std::string, std::string> > designables { { "1", "NA" } };

			utility::vector1<Size> our_res { 1 };
			utility::vector1<bool> expect_current { false };
			utility::vector1< std::set<char> > expect_aa { { 'A' } };



			do_test_pruning( pose, allow_repacking_res, designables, false, our_res, expect_current, expect_aa);
		}

	}

	// Testing presatisfaction on sc atoms
	void test_pruning_asn_to_satisfied_ser(){

		pose::Pose pose;
		import_pose::pose_from_file( pose, "protocols/task_operations/asn_to_satisfied_ser.pdb" );

		{
			std::string allow_repacking_res = "1";

			utility::vector1<std::pair<std::string, std::string> > designables { { "1", "NA" } };

			utility::vector1<Size> our_res { 1 };
			utility::vector1<bool> expect_current { false };
			utility::vector1< std::set<char> > expect_aa { { 'A' } };



			do_test_pruning( pose, allow_repacking_res, designables, true, our_res, expect_current, expect_aa);
		}
		{
			std::string allow_repacking_res = "1";

			utility::vector1<std::pair<std::string, std::string> > designables { { "1", "NA" } };

			utility::vector1<Size> our_res { 1 };
			utility::vector1<bool> expect_current { false };
			utility::vector1< std::set<char> > expect_aa { { 'A' } };



			do_test_pruning( pose, allow_repacking_res, designables, false, our_res, expect_current, expect_aa);
		}

	}


	// Testing even trades on bb atoms
	void test_pruning_asn_to_bb(){

		pose::Pose pose;
		import_pose::pose_from_file( pose, "protocols/task_operations/asn_to_bb.pdb" );

		{
			std::string allow_repacking_res = "1";

			utility::vector1<std::pair<std::string, std::string> > designables { { "1", "NA" } };

			utility::vector1<Size> our_res { 1 };
			utility::vector1<bool> expect_current { true };
			utility::vector1< std::set<char> > expect_aa { { 'A', 'N' } };



			do_test_pruning( pose, allow_repacking_res, designables, true, our_res, expect_current, expect_aa);
		}
		{
			std::string allow_repacking_res = "1";

			utility::vector1<std::pair<std::string, std::string> > designables { { "1", "NA" } };

			utility::vector1<Size> our_res { 1 };
			utility::vector1<bool> expect_current { false };
			utility::vector1< std::set<char> > expect_aa { { 'A' } };



			do_test_pruning( pose, allow_repacking_res, designables, false, our_res, expect_current, expect_aa);
		}

	}

	// Testing presatisfaction on bb atoms
	void test_pruning_asn_to_satisfied_bb(){

		pose::Pose pose;
		import_pose::pose_from_file( pose, "protocols/task_operations/asn_to_satisfied_bb.pdb" );

		{
			std::string allow_repacking_res = "1";

			utility::vector1<std::pair<std::string, std::string> > designables { { "1", "NA" } };

			utility::vector1<Size> our_res { 1 };
			utility::vector1<bool> expect_current { false };
			utility::vector1< std::set<char> > expect_aa { { 'A' } };



			do_test_pruning( pose, allow_repacking_res, designables, true, our_res, expect_current, expect_aa);
		}
		{
			std::string allow_repacking_res = "1";

			utility::vector1<std::pair<std::string, std::string> > designables { { "1", "NA" } };

			utility::vector1<Size> our_res { 1 };
			utility::vector1<bool> expect_current { false };
			utility::vector1< std::set<char> > expect_aa { { 'A' } };



			do_test_pruning( pose, allow_repacking_res, designables, false, our_res, expect_current, expect_aa);
		}

	}


};
