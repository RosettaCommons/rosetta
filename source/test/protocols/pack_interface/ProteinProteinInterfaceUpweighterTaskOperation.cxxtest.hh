// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/pack_interface/ProteinProteinInterfaceUpweighterTaskOperation.cxxtest.hh
/// @brief  test suite for ApproximateBuriedUnsatPenalty
/// @author Longxing Cao (longxing@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>

#include <platform/types.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <core/conformation/Residue.hh>
#include <core/id/AtomID_Map.hh>
#include <core/pack/interaction_graph/AnnealableGraphBase.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/task/operation/OperateOnResidueSubset.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/pack_interface/ProteinProteinInterfaceUpweighterTaskOperation.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>

#include <core/types.hh>

#include <test/UTracer.hh>
#include <test/core/init_util.hh>

//Auto Headers
#include <utility/vector1.hh>
#include <utility/pointer/memory.hh>

using namespace core;


static basic::Tracer TR("protocols.pack_interface.ProteinProteinInterfaceUpweighterTests.cxxtest");

class ProteinProteinInterfaceUpweighterTests : public CxxTest::TestSuite
{
public:
	void setUp()
	{
		core_init();
		import_pose::pose_from_file( _4SER_, "protocols/pack_interface/4SER_buriedASP.pdb.gz" );
	}



	// Testing function to start a packing trajectory and then use ala scanning
	//  to ensure that the energies were calculated correctly
	core::Real do_test_packing(
		pose::Pose & pose,
		std::string const & allow_repacking_res,
		utility::vector1<std::pair<std::string, std::string> > const & designables,
		Size const & res,
		utility::vector1<core::pack::task::operation::TaskOperationOP> task_ops
	) {


		using utility::pointer::make_shared;
		using namespace core::pack::task::operation;
		using namespace core::select::residue_selector;

		// We also need hbond_sc this time so we know what rotamer is what

		scoring::ScoreFunctionOP sfxn ( new scoring::ScoreFunction() );
		sfxn->set_weight( scoring::hbond_sc, 1.0 );



		core::pack::task::TaskFactoryOP tf( new core::pack::task::TaskFactory() );
		tf->push_back( TaskOperationCOP( new IncludeCurrent() ) );
		tf->push_back( make_shared< OperateOnResidueSubset >(
			make_shared< PreventRepackingRLT >(), make_shared< ResidueIndexSelector >( allow_repacking_res ), true ) );

		for ( TaskOperationOP & task : task_ops ) {
			tf->push_back( task );
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
		ig->prepare_for_simulated_annealing();
		ig->blanket_assign_state_0();

		// First lets make sure that the ig and the rotsets are mirrors of each other
		TS_ASSERT( rotsets->nmoltenres() == (Size)ig->get_num_nodes() );

		for ( Size imolt = 1; imolt <= rotsets->nmoltenres(); imolt++ ) {
			TS_ASSERT( rotsets->nrotamers_for_moltenres( imolt ) == (Size)ig->get_num_states_for_node( imolt ) );
		}

		// Now lets find the original rotamers and tell the interaction graph that we're placing them
		//  Also take note of the alanine irots

		Size moltres = rotsets->resid_2_moltenres( res );
		pack::rotamer_set::RotamerSetCOP rotset = rotsets->rotamer_set_for_residue( res );

		conformation::Residue const & match_res = pose.residue( res );
		Size the_irot = 0;
		Size ala_irot = 0;
		for ( Size irot = 1; irot <= rotset->num_rotamers(); irot++ ) {
			conformation::ResidueCOP rotamer = rotset->rotamer( irot );
			if ( rotamer->name1() == 'A' ) {
				ala_irot = irot;
				continue;
			}
			if ( match_res.name1() != rotamer->name1() ) continue;

			bool chi_match = true;
			for ( Size chi = 1; chi < match_res.nchi(); chi++ ) {
				if ( std::abs( match_res.chi( chi ) - rotamer->chi( chi ) ) > 1 ) {
					chi_match = false;
					break;
				}
			}
			if ( ! chi_match ) continue;

			the_irot = irot;
		}
		TS_ASSERT( ala_irot != 0 );
		TS_ASSERT( the_irot != 0 );

		ig->set_state_for_node( moltres, the_irot );

		// Ok, now do alanine scanning and see if the energies make sense


		core::PackerEnergy delta_energy = 0;
		core::PackerEnergy prev_energy = 0;
		ig->consider_substitution( moltres, ala_irot, delta_energy, prev_energy );

		return -delta_energy;

	}

	// In this test, residue 2 has a buried polar that is being satisfied by 4 residues
	// Tests the sc-sc case for ERROR FIX 3
	void test_packing_4SER(){
		pose::Pose pose = _4SER_;


		std::string allow_repacking_res = "2,4";

		utility::vector1<std::pair<std::string, std::string> > designables { { "2", "DA" }, { "4", "SA" } };

		Size our_res = 4;


		utility::vector1<core::pack::task::operation::TaskOperationOP> task_ops;
		Real score_no_upweight = do_test_packing( pose, allow_repacking_res, designables, our_res, task_ops );
		core::pack::task::operation::TaskOperationOP upweighter( utility::pointer::make_shared<protocols::pack_interface::ProteinProteinInterfaceUpweighter>(3.0) );
		task_ops.push_back(upweighter);
		Real score_upweight = do_test_packing( pose, allow_repacking_res, designables, our_res, task_ops );

		TS_ASSERT_DELTA(3 * score_no_upweight, score_upweight, 0.1);


	}


private:
	pose::Pose _4SER_;


};
