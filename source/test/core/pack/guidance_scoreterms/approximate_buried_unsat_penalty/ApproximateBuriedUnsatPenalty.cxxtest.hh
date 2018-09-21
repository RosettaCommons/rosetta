// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/guidance_scoreterms/approximate_buried_unsat_penalty/ApproximateBuriedUnsatPenalty.cxxtest.hh
/// @brief  test suite for ApproximateBuriedUnsatPenalty
/// @author Brian Coventry (bcov@uw.edu)

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
#include <core/scoring/atomic_depth/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>

#include <protocols/pack_interface/ProteinProteinInterfaceUpweighterTaskOperation.hh>

#include <core/types.hh>

#include <test/UTracer.hh>
#include <test/core/init_util.hh>

//Auto Headers
#include <utility/vector1.hh>
#include <utility/pointer/memory.hh>

using namespace core;


static basic::Tracer TR("core.pack.guidance_scoreterms.approximate_buried_unsat_penalty.ApproximateBuriedUnsatPenaltyTests.cxxtest");

class ApproximateBuriedUnsatPenaltyTests : public CxxTest::TestSuite
{
public:
	void setUp()
	{
		core_init();
		import_pose::pose_from_file( _4SER_, "core/pack/guidance_scoreterms/approximate_buried_unsat_penalty/4SER_buriedASP.pdb.gz" );
		import_pose::pose_from_file( _1ARG_, "core/pack/guidance_scoreterms/approximate_buried_unsat_penalty/1ARG_buriedASP.pdb.gz" );
		import_pose::pose_from_file( _1ARG_bbMET_, "core/pack/guidance_scoreterms/approximate_buried_unsat_penalty/1ARG_bbMET_buriedASP.pdb.gz" );
		import_pose::pose_from_file( _1ARG_bbMET_O_, "core/pack/guidance_scoreterms/approximate_buried_unsat_penalty/1ARG_bbMET_buriedO.pdb.gz" );
		import_pose::pose_from_file( _APA_, "core/pack/guidance_scoreterms/approximate_buried_unsat_penalty/APA.pdb.gz" );
		import_pose::pose_from_file( _proline_nh_test_, "core/pack/guidance_scoreterms/approximate_buried_unsat_penalty/proline_nh_test.pdb.gz" );
	}

	void assert_edge( scoring::EnergyGraph const & energy_graph, Size res1, Size res2, float expect ) {
		scoring::EnergyEdge const * edge = energy_graph.find_energy_edge( res1, res2 );
		TS_ASSERT( edge );
		if ( ! edge ) return;

		float actual = edge->fill_energy_map()[ scoring::approximate_buried_unsat_penalty ];
		TS_ASSERT_DELTA( actual, expect, 0.01 );
	}

	const float burial_probe_radius = 2.3f;
	const float burial_depth = 2.5f;
	const float burial_resolution = 0.25f;

	// If this test fails, we don't expect anything else to pass. But the issue is that the
	//  buried polar atom being tested isn't being considered buried.
	//  If you have to fix this test (probably because you changed atomic radii),
	//  your best bet is to play with the burial metric numbers until everything works.
	//  ( Or maybe hardcode the atomic weights back if things get rough )
	void test_burial_sanity_check() {
		{
			pose::Pose pose = _4SER_;

			core::id::AtomID_Map< bool > buried = scoring::atomic_depth::atoms_deeper_than(
				pose, burial_depth, false, burial_probe_radius, true, burial_resolution );


			TS_ASSERT( buried( 2, pose.residue(2).atom_index("OD2") ) );
			TS_ASSERT( !buried( 2, pose.residue(2).atom_index("OD1") ) );
			TS_ASSERT( !buried( 2, pose.residue(2).atom_index("O") ) );
			TS_ASSERT( !buried( 2, pose.residue(2).atom_index("N") ) );
			TS_ASSERT( !buried( 4, pose.residue(4).atom_index("OG") ) );
			TS_ASSERT( !buried( 5, pose.residue(5).atom_index("OG") ) );
			TS_ASSERT( !buried( 6, pose.residue(6).atom_index("OG") ) );
			TS_ASSERT( !buried( 7, pose.residue(7).atom_index("OG") ) );
		}
		{
			pose::Pose pose = _1ARG_;

			core::id::AtomID_Map< bool > buried = scoring::atomic_depth::atoms_deeper_than(
				pose, burial_depth, false, burial_probe_radius, true, burial_resolution );


			TS_ASSERT( buried( 3, pose.residue(3).atom_index("OD2") ) );
			TS_ASSERT( !buried( 3, pose.residue(3).atom_index("OD1") ) );
			TS_ASSERT( !buried( 3, pose.residue(3).atom_index("O") ) );
			TS_ASSERT( !buried( 3, pose.residue(3).atom_index("N") ) );
			TS_ASSERT( !buried( 1, pose.residue(1).atom_index("NH1") ) );
			TS_ASSERT( !buried( 1, pose.residue(1).atom_index("NH2") ) );
			TS_ASSERT( !buried( 1, pose.residue(1).atom_index("NE") ) );
			TS_ASSERT( !buried( 1, pose.residue(1).atom_index("N") ) );
			TS_ASSERT( !buried( 1, pose.residue(1).atom_index("OXT") ) );
			TS_ASSERT( !buried( 1, pose.residue(1).atom_index("O") ) );
		}
		{
			pose::Pose pose = _1ARG_bbMET_;

			core::id::AtomID_Map< bool > buried = scoring::atomic_depth::atoms_deeper_than(
				pose, burial_depth, false, burial_probe_radius, burial_resolution );


			TS_ASSERT( buried( 3, pose.residue(3).atom_index("OD2") ) );
			TS_ASSERT( !buried( 3, pose.residue(3).atom_index("OD1") ) );
			TS_ASSERT( !buried( 3, pose.residue(3).atom_index("O") ) );
			TS_ASSERT( !buried( 3, pose.residue(3).atom_index("N") ) );
			TS_ASSERT( !buried( 1, pose.residue(1).atom_index("NH1") ) );
			TS_ASSERT( !buried( 1, pose.residue(1).atom_index("NH2") ) );
			TS_ASSERT( !buried( 1, pose.residue(1).atom_index("NE") ) );
			TS_ASSERT( !buried( 1, pose.residue(1).atom_index("N") ) );
			TS_ASSERT( !buried( 1, pose.residue(1).atom_index("OXT") ) );
			TS_ASSERT( !buried( 1, pose.residue(1).atom_index("O") ) );
			TS_ASSERT( !buried( 8, pose.residue(8).atom_index("N") ) );
			TS_ASSERT( !buried( 8, pose.residue(8).atom_index("OXT") ) );
			TS_ASSERT( !buried( 8, pose.residue(8).atom_index("O") ) );
		}
		{
			pose::Pose pose = _1ARG_bbMET_O_;

			core::id::AtomID_Map< bool > buried = scoring::atomic_depth::atoms_deeper_than(
				pose, burial_depth, false, burial_probe_radius, burial_resolution );


			TS_ASSERT( buried( 7, pose.residue(7).atom_index("O") ) );
			TS_ASSERT( !buried( 7, pose.residue(7).atom_index("N") ) );
			TS_ASSERT( !buried( 5, pose.residue(5).atom_index("NH1") ) );
			TS_ASSERT( !buried( 5, pose.residue(5).atom_index("NH2") ) );
			TS_ASSERT( !buried( 5, pose.residue(5).atom_index("NE") ) );
			TS_ASSERT( !buried( 5, pose.residue(5).atom_index("N") ) );
			TS_ASSERT( !buried( 5, pose.residue(5).atom_index("OXT") ) );
			TS_ASSERT( !buried( 5, pose.residue(5).atom_index("O") ) );
			TS_ASSERT( !buried( 3, pose.residue(3).atom_index("N") ) );
			TS_ASSERT( !buried( 3, pose.residue(3).atom_index("OXT") ) );
			TS_ASSERT( !buried( 3, pose.residue(3).atom_index("O") ) );
		}
		// Everything is buried for these next two
		{
			pose::Pose pose = _APA_;

			core::id::AtomID_Map< bool > buried = scoring::atomic_depth::atoms_deeper_than(
				pose, 0, false, 0.1, burial_resolution );


			TS_ASSERT( buried( 2, pose.residue(2).atom_index("N") ) );
			TS_ASSERT( buried( 2, pose.residue(2).atom_index("O") ) );
		}
		{
			pose::Pose pose = _proline_nh_test_;

			core::id::AtomID_Map< bool > buried = scoring::atomic_depth::atoms_deeper_than(
				pose, 0, false, 0.1, burial_resolution );


			TS_ASSERT( buried( 2, pose.residue(2).atom_index("N") ) );
			TS_ASSERT( buried( 2, pose.residue(2).atom_index("O") ) );
			TS_ASSERT( buried( 4, pose.residue(4).atom_index("N") ) );
			TS_ASSERT( buried( 4, pose.residue(4).atom_index("O") ) );
			TS_ASSERT( buried( 4, pose.residue(4).atom_index("OXT") ) );
			TS_ASSERT( buried( 4, pose.residue(4).atom_index("OG") ) );
		}

	}


	// Testing function to start a packing trajectory and then use ala scanning
	//  to ensure that the energies were calculated correctly
	void do_test_packing(
		pose::Pose & pose,
		std::string const & allow_repacking_res,
		utility::vector1<std::pair<std::string, std::string> > const & designables,
		utility::vector1<Size> const & our_res,
		utility::vector1<Real> expected_scores,
		bool add_ig_edge_reweight = false,
		bool assume_const_backbone = false,
		bool bury_everything = false
	) {


		using utility::pointer::make_shared;
		using namespace core::pack::task::operation;
		using namespace core::select::residue_selector;

		// We also need hbond_sc this time so we know what rotamer is what

		scoring::ScoreFunctionOP sfxn ( new scoring::ScoreFunction() );
		scoring::methods::EnergyMethodOptions options = sfxn->energy_method_options();
		options.approximate_buried_unsat_penalty_hbond_energy_threshold( 0 );
		options.approximate_buried_unsat_penalty_burial_atomic_depth( bury_everything ?   0 : burial_depth );
		options.approximate_buried_unsat_penalty_burial_probe_radius( bury_everything ? 0.1 : burial_probe_radius );
		options.approximate_buried_unsat_penalty_burial_resolution( burial_resolution );
		options.approximate_buried_unsat_penalty_assume_const_backbone( assume_const_backbone );
		sfxn->set_energy_method_options( options );

		sfxn->set_weight( scoring::approximate_buried_unsat_penalty, 10.0 );



		core::pack::task::TaskFactoryOP tf( new core::pack::task::TaskFactory() );
		tf->push_back( TaskOperationCOP( new IncludeCurrent() ) );
		tf->push_back( make_shared< OperateOnResidueSubset >(
			make_shared< PreventRepackingRLT >(), make_shared< ResidueIndexSelector >( allow_repacking_res ), true ) );

		if ( add_ig_edge_reweight ) {            // Make it something wacky so we can't accidentally pass
			tf->push_back( make_shared<protocols::pack_interface::ProteinProteinInterfaceUpweighter>( 7.3 ) );
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

		utility::vector1<Size> our_irots;
		utility::vector1<Size> ala_irots;

		for ( Size resnum : our_res ) {
			Size moltres = rotsets->resid_2_moltenres( resnum );
			pack::rotamer_set::RotamerSetCOP rotset = rotsets->rotamer_set_for_residue( resnum );

			conformation::Residue const & match_res = pose.residue( resnum );
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

			ala_irots.push_back( ala_irot );
			// our_irots.push_back( the_irot );

			ig->set_state_for_node( moltres, the_irot );
		}

		// Ok, now do alanine scanning and see if the energies make sense


		for ( Size ipos = 1; ipos <= our_res.size(); ipos++ ) {

			Size resnum = our_res[ ipos ];
			Size moltres = rotsets->resid_2_moltenres( resnum );

			// Size the_irot = our_irots[ ipos ];
			Size ala_irot = ala_irots[ ipos ];

			core::PackerEnergy delta_energy = 0;
			core::PackerEnergy prev_energy = 0;
			ig->consider_substitution( moltres, ala_irot, delta_energy, prev_energy );

			// std::cout << ipos << " " << resnum << " " << -delta_energy << " " << expected_scores[ipos] << std::endl;
			TS_ASSERT_DELTA( -delta_energy, expected_scores[ ipos ], 0.1 );

		}


	}

	// In this test, residue 2 has a buried polar that is being satisfied by 4 residues
	void test_scoring_4SER(){
		pose::Pose pose = _4SER_;

		// We only need our scoreterm

		scoring::ScoreFunctionOP sfxn ( new scoring::ScoreFunction() );
		scoring::methods::EnergyMethodOptions options = sfxn->energy_method_options();
		options.approximate_buried_unsat_penalty_hbond_energy_threshold( 0 );
		options.approximate_buried_unsat_penalty_burial_atomic_depth( burial_depth );
		options.approximate_buried_unsat_penalty_burial_probe_radius( burial_probe_radius );
		options.approximate_buried_unsat_penalty_burial_resolution( burial_resolution );
		sfxn->set_energy_method_options( options );

		sfxn->set_weight( scoring::approximate_buried_unsat_penalty, 1.0 );

		sfxn->score( pose );

		scoring::Energies const & energies = pose.energies();

		// Remember that pair energies get halved
		TS_ASSERT_DELTA( energies.residue_total_energy( 2 ), -1, 0.01 ); // 1 - 4 * 0.5
		TS_ASSERT_DELTA( energies.residue_total_energy( 4 ), 1, 0.01 );  // -0.5 + 3 * 0.5
		TS_ASSERT_DELTA( energies.residue_total_energy( 5 ), 1, 0.01 );  // -0.5 + 3 * 0.5
		TS_ASSERT_DELTA( energies.residue_total_energy( 6 ), 1, 0.01 );  // -0.5 + 3 * 0.5
		TS_ASSERT_DELTA( energies.residue_total_energy( 7 ), 1, 0.01 );  // -0.5 + 3 * 0.5

		scoring::EnergyGraph const & energy_graph = energies.energy_graph();

		// hbonds
		assert_edge( energy_graph, 2, 4, -1 );
		assert_edge( energy_graph, 2, 5, -1 );
		assert_edge( energy_graph, 2, 6, -1 );
		assert_edge( energy_graph, 2, 7, -1 );

		// oversat penalties
		assert_edge( energy_graph, 4, 5, 1 );
		assert_edge( energy_graph, 4, 6, 1 );
		assert_edge( energy_graph, 4, 7, 1 );
		assert_edge( energy_graph, 5, 6, 1 );
		assert_edge( energy_graph, 5, 7, 1 );
		assert_edge( energy_graph, 6, 7, 1 );



	}

	// In this test, residue 2 has a buried polar that is being satisfied by 4 residues
	// Tests the sc-sc case for ERROR FIX 3
	void test_packing_4SER(){
		pose::Pose pose = _4SER_;


		std::string allow_repacking_res = "2,4-7";

		utility::vector1<std::pair<std::string, std::string> > designables { { "2", "DA" }, { "4-7", "SA" } };

		utility::vector1<Size> our_res { 2, 4, 5, 6, 7 };
		utility::vector1<Real> expected_scores { -30, 20, 20, 20, 20 };

		do_test_packing( pose, allow_repacking_res, designables, our_res, expected_scores, false, false );
		do_test_packing( pose, allow_repacking_res, designables, our_res, expected_scores, false, true );

	}

	// In this test, residue 2 has a buried polar that is being satisfied by 4 residues
	// Tests ability of scoreterm to avoid being affected by IGEdgeReweights
	//     Residue 2 is on a different chain from residues 4-7
	// Tests the sc-sc case for ERROR FIX 3
	void test_packing_IGReweight(){
		pose::Pose pose = _4SER_;


		std::string allow_repacking_res = "2,4-7";

		utility::vector1<std::pair<std::string, std::string> > designables { { "2", "DA" }, { "4-7", "SA" } };

		utility::vector1<Size> our_res { 2, 4, 5, 6, 7 };
		utility::vector1<Real> expected_scores { -30, 20, 20, 20, 20 };

		do_test_packing( pose, allow_repacking_res, designables, our_res, expected_scores, true, false );
		do_test_packing( pose, allow_repacking_res, designables, our_res, expected_scores, true, true );

	}


	// In this test, residue 3 has a buried polar that is being satisfied by 2 atoms of residue 1
	void test_scoring_1ARG(){
		pose::Pose pose = _1ARG_;

		// We only need our scoreterm

		scoring::ScoreFunctionOP sfxn ( new scoring::ScoreFunction() );
		scoring::methods::EnergyMethodOptions options = sfxn->energy_method_options();
		options.approximate_buried_unsat_penalty_hbond_energy_threshold( 0 );
		options.approximate_buried_unsat_penalty_burial_atomic_depth( burial_depth );
		options.approximate_buried_unsat_penalty_burial_probe_radius( burial_probe_radius );
		options.approximate_buried_unsat_penalty_burial_resolution( burial_resolution );
		sfxn->set_energy_method_options( options );

		sfxn->set_weight( scoring::approximate_buried_unsat_penalty, 1.0 );

		sfxn->score( pose );

		scoring::Energies const & energies = pose.energies();

		// Remember that pair energies get halved
		TS_ASSERT_DELTA( energies.residue_total_energy( 3 ), 0, 0.01 ); // 1 - 2 * 0.5
		TS_ASSERT_DELTA( energies.residue_total_energy( 1 ), 0, 0.01 );  // 2 * -0.5 + 1

		scoring::EnergyGraph const & energy_graph = energies.energy_graph();

		// hbonds
		assert_edge( energy_graph, 1, 3, -2 );

		// oversat penalties
		//  There is one, but it's internal...

	}

	// In this test, residue 3 has a buried polar that is being satisfied by 2 atoms of residue 1
	//  Tests ERROR FIX 3 for the case when the two atoms are on the same rotamer
	//  Also tests conversion of a 2-body edge to a 1-body edge
	void test_packing_1ARG(){
		pose::Pose pose = _1ARG_;


		std::string allow_repacking_res = "1,3";

		utility::vector1<std::pair<std::string, std::string> > designables { { "1", "RA" }, { "3", "DA" } };

		utility::vector1<Size> our_res { 1, 3 };
		utility::vector1<Real> expected_scores { -10, -10 };

		do_test_packing( pose, allow_repacking_res, designables, our_res, expected_scores, false, true );
		do_test_packing( pose, allow_repacking_res, designables, our_res, expected_scores, false, false );

	}



	// In this test, residue 3 has a buried polar that is being satisfied by 2 atoms of residue 1 and a bb of residue 8
	void test_scoring_1ARG_bbMET(){
		pose::Pose pose = _1ARG_bbMET_;

		// We only need our scoreterm

		scoring::ScoreFunctionOP sfxn ( new scoring::ScoreFunction() );
		scoring::methods::EnergyMethodOptions options = sfxn->energy_method_options();
		options.approximate_buried_unsat_penalty_hbond_energy_threshold( 0 );
		options.approximate_buried_unsat_penalty_burial_atomic_depth( burial_depth );
		options.approximate_buried_unsat_penalty_burial_probe_radius( burial_probe_radius );
		options.approximate_buried_unsat_penalty_burial_resolution( burial_resolution );
		sfxn->set_energy_method_options( options );

		sfxn->set_weight( scoring::approximate_buried_unsat_penalty, 1.0 );

		sfxn->score( pose );

		scoring::Energies const & energies = pose.energies();

		// Remember that pair energies get halved
		TS_ASSERT_DELTA( energies.residue_total_energy( 3 ), -0.5, 0.01 ); // 1 - 3 * 0.5
		TS_ASSERT_DELTA( energies.residue_total_energy( 1 ), 1, 0.01 );  // 2 * -0.5 + 1 + 2 * 0.5
		TS_ASSERT_DELTA( energies.residue_total_energy( 8 ), 0.5, 0.01 );  // -0.5 + 2 * 0.5

		scoring::EnergyGraph const & energy_graph = energies.energy_graph();

		// hbonds
		assert_edge( energy_graph, 1, 3, -2 );
		assert_edge( energy_graph, 3, 8, -1 );

		// oversat penalties
		assert_edge( energy_graph, 1, 8, 2 );

	}

	// In this test, residue 3 has a buried polar that is being satisfied by 2 atoms of residue 1
	//  Tests ERROR FIX 3 for the case when the two atoms are on the same rotamer
	//  Tests ERROR FIX 3 for the case when one atom is a backbone
	//  Also tests conversion of a 2-body edge to a 1-body edge
	void test_packing_1ARG_bbMET(){
		pose::Pose pose = _1ARG_bbMET_;


		std::string allow_repacking_res = "1,3,8";

		utility::vector1<std::pair<std::string, std::string> > designables { { "1", "RA" }, { "3", "DA" }, { "8", "MA" } };

		utility::vector1<Size> our_res { 1, 3, 8 };

		utility::vector1<Real> expected_scores { 10,  // -10 * 2 + 10 * 3
			-20,  // 10 - 10 * 3
			0  // 0 // bb is part of background for packing
			};

		do_test_packing( pose, allow_repacking_res, designables, our_res, expected_scores, false, true );
		do_test_packing( pose, allow_repacking_res, designables, our_res, expected_scores, false, false );

	}


	// In this test, residue 2 has a buried polar bb that is being satisfied by 2 atoms of residue 4 and a bb of residue 8
	void test_scoring_1ARG_bbMET_O(){
		pose::Pose pose = _1ARG_bbMET_O_;

		// We only need our scoreterm

		scoring::ScoreFunctionOP sfxn ( new scoring::ScoreFunction() );
		scoring::methods::EnergyMethodOptions options = sfxn->energy_method_options();
		options.approximate_buried_unsat_penalty_hbond_energy_threshold( 0 );
		options.approximate_buried_unsat_penalty_burial_atomic_depth( burial_depth );
		options.approximate_buried_unsat_penalty_burial_probe_radius( burial_probe_radius );
		options.approximate_buried_unsat_penalty_burial_resolution( burial_resolution );
		sfxn->set_energy_method_options( options );

		sfxn->set_weight( scoring::approximate_buried_unsat_penalty, 1.0 );

		sfxn->score( pose );

		scoring::Energies const & energies = pose.energies();

		// Remember that pair energies get halved
		TS_ASSERT_DELTA( energies.residue_total_energy( 7 ), -0.5, 0.01 ); // 1 - 3 * 0.5
		TS_ASSERT_DELTA( energies.residue_total_energy( 5 ), 1, 0.01 );  // 2 * -0.5 + 1 + 2 * 0.5
		TS_ASSERT_DELTA( energies.residue_total_energy( 3 ), 0.5, 0.01 );  // -0.5 + 2 * 0.5

		scoring::EnergyGraph const & energy_graph = energies.energy_graph();

		// hbonds
		assert_edge( energy_graph, 7, 5, -2 );
		assert_edge( energy_graph, 7, 3, -1 );

		// oversat penalties
		assert_edge( energy_graph, 3, 5, 2 );

	}

	// In this test, residue 7 has a buried polar bb that is being satisfied by 2 atoms of residue 5 and a bb of residue 3
	//  Tests ERROR FIX 1
	//  Tests ERROR FIX 2
	//  Also tests conversion of a 2-body edge to a 1-body edge
	void test_packing_1ARG_bbMET_O(){
		pose::Pose pose = _1ARG_bbMET_O_;


		std::string allow_repacking_res = "7,5,3";

		utility::vector1<std::pair<std::string, std::string> > designables { { "7", "MA" }, { "5", "RA" }, { "3", "MA" } };

		utility::vector1<Size> our_res { 7, 5, 3 };

		utility::vector1<Real> expected_scores {  0,  // 0  // bb is part of background for packing
			10,  // -10 * 2 + 3 * 10
			0  // 0 // bb is part of background for packing
			};

		do_test_packing( pose, allow_repacking_res, designables, our_res, expected_scores, false, false );
		do_test_packing( pose, allow_repacking_res, designables, our_res, expected_scores, false, true );

	}



	// In this test, residue 2 is proline and every atom is buried
	void test_scoring_APA(){
		pose::Pose pose = _APA_;

		// We only need our scoreterm

		scoring::ScoreFunctionOP sfxn ( new scoring::ScoreFunction() );
		scoring::methods::EnergyMethodOptions options = sfxn->energy_method_options();
		options.approximate_buried_unsat_penalty_hbond_energy_threshold( 0 );
		options.approximate_buried_unsat_penalty_burial_atomic_depth( 0 );
		options.approximate_buried_unsat_penalty_burial_probe_radius( 0.1 );
		options.approximate_buried_unsat_penalty_burial_resolution( burial_resolution );
		sfxn->set_energy_method_options( options );

		sfxn->set_weight( scoring::approximate_buried_unsat_penalty, 1.0 );

		sfxn->score( pose );

		scoring::Energies const & energies = pose.energies();

		// Remember that pair energies get halved
		TS_ASSERT_DELTA( energies.residue_total_energy( 2 ), 1, 0.01 ); // 1

	}

	// In this test, residue 2 is proline and every atom is buried
	//  Tests assume_const_backbone
	void test_packing_APA(){
		pose::Pose pose = _APA_;

		std::string allow_repacking_res = "2";

		utility::vector1<std::pair<std::string, std::string> > designables { { "2", "AP" } };

		utility::vector1<Size> our_res { 2 };

		utility::vector1<Real> expected_scores_const_bb {  0  // no change on mutation to ala because bb didn't change
			};
		utility::vector1<Real> expected_scores_nonconst_bb {  -10  // ala has a buried N
			};

		do_test_packing( pose, allow_repacking_res, designables, our_res, expected_scores_nonconst_bb, false, false, true );
		do_test_packing( pose, allow_repacking_res, designables, our_res, expected_scores_const_bb,    false, true, true );

	}

	// In this test, residue 2 is proline and residue 4 is trying to h-bond to it
	//   all atoms are buried
	void test_scoring_proline_nh_test(){
		pose::Pose pose = _proline_nh_test_;

		// We only need our scoreterm

		scoring::ScoreFunctionOP sfxn ( new scoring::ScoreFunction() );
		scoring::methods::EnergyMethodOptions options = sfxn->energy_method_options();
		options.approximate_buried_unsat_penalty_hbond_energy_threshold( 0 );
		options.approximate_buried_unsat_penalty_burial_atomic_depth( 0 );
		options.approximate_buried_unsat_penalty_burial_probe_radius( 0.1 );
		options.approximate_buried_unsat_penalty_burial_resolution( burial_resolution );
		sfxn->set_energy_method_options( options );

		sfxn->set_weight( scoring::approximate_buried_unsat_penalty, 1.0 );

		sfxn->score( pose );

		scoring::Energies const & energies = pose.energies();

		// Remember that pair energies get halved
		TS_ASSERT_DELTA( energies.residue_total_energy( 2 ), 1, 0.01 );  // 1      // just the O
		TS_ASSERT_DELTA( energies.residue_total_energy( 4 ), 6, 0.01 );  // 1 + 1 +1 + 3   // O, OXT, OG, NH3

		scoring::EnergyGraph const & energy_graph = energies.energy_graph();

		// hbonds
		assert_edge( energy_graph, 2, 4, 0 );

	}

	// In this test, residue 2 is proline and every atom is buried
	//  Tests assume_const_backbone
	void test_packing_proline_nh_test(){
		pose::Pose pose = _proline_nh_test_;

		std::string allow_repacking_res = "2,4";

		utility::vector1<std::pair<std::string, std::string> > designables { { "2", "AP" }, { "4", "AS" } };

		utility::vector1<Size> our_res { 2, 4 };

		utility::vector1<Real> expected_scores_const_bb {  0,  // no change on mutation to ala because bb didn't change
			-10  // now the proline N is considered buried
			};
		utility::vector1<Real> expected_scores_nonconst_bb {  10,  // alanine is able to satisfy OG
			10 // alanine doesn't have an unsatisfied OG
			};

		do_test_packing( pose, allow_repacking_res, designables, our_res, expected_scores_nonconst_bb, false, false, true );
		do_test_packing( pose, allow_repacking_res, designables, our_res, expected_scores_const_bb,    false, true, true );

	}


private:
	pose::Pose _4SER_;
	pose::Pose _1ARG_;
	pose::Pose _1ARG_bbMET_;
	pose::Pose _1ARG_bbMET_O_;
	pose::Pose _APA_;
	pose::Pose _proline_nh_test_;


};
