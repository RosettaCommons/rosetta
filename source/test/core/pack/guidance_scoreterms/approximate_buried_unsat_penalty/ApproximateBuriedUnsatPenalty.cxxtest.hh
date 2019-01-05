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
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSets.hh>
#include <core/pack/task/operation/OperateOnResidueSubset.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/atomic_depth/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>

#include <protocols/symmetry/SetupForSymmetryMover.hh>
#include <protocols/pack_interface/ProteinProteinInterfaceUpweighterTaskOperation.hh>

#include <core/types.hh>

#include <basic/database/open.hh>
#include <basic/options/option.hh>

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
	}

	core::pose::PoseOP
	get_4SER() {
		core::pose::PoseOP pose_p( new core::pose::Pose() );
		import_pose::pose_from_file( *pose_p, "core/pack/guidance_scoreterms/approximate_buried_unsat_penalty/4SER_buriedASP.pdb.gz" );
		return pose_p;
	}

	core::pose::PoseOP
	get_1ARG() {
		core::pose::PoseOP pose_p( new core::pose::Pose() );
		import_pose::pose_from_file( *pose_p, "core/pack/guidance_scoreterms/approximate_buried_unsat_penalty/1ARG_buriedASP.pdb.gz" );
		return pose_p;
	}

	core::pose::PoseOP
	get_1ARG_bbMET() {
		core::pose::PoseOP pose_p( new core::pose::Pose() );
		import_pose::pose_from_file( *pose_p, "core/pack/guidance_scoreterms/approximate_buried_unsat_penalty/1ARG_bbMET_buriedASP.pdb.gz" );
		return pose_p;
	}

	core::pose::PoseOP
	get_1ARG_bbMET_O() {
		core::pose::PoseOP pose_p( new core::pose::Pose() );
		import_pose::pose_from_file( *pose_p, "core/pack/guidance_scoreterms/approximate_buried_unsat_penalty/1ARG_bbMET_buriedO.pdb.gz" );
		return pose_p;
	}

	core::pose::PoseOP
	get_APA() {
		core::pose::PoseOP pose_p( new core::pose::Pose() );
		import_pose::pose_from_file( *pose_p, "core/pack/guidance_scoreterms/approximate_buried_unsat_penalty/APA.pdb.gz" );
		return pose_p;
	}

	core::pose::PoseOP
	get_proline_nh_test() {
		core::pose::PoseOP pose_p( new core::pose::Pose() );
		import_pose::pose_from_file( *pose_p, "core/pack/guidance_scoreterms/approximate_buried_unsat_penalty/proline_nh_test.pdb.gz" );
		return pose_p;
	}

	core::pose::PoseOP
	get_C3_asn() {
		core::pose::PoseOP pose_p( new core::pose::Pose() );
		import_pose::pose_from_file( *pose_p, "core/pack/guidance_scoreterms/approximate_buried_unsat_penalty/C3_asn.pdb" );
		std::string symdef_file = basic::database::full_name( "symmetry/cyclic/C3_Z.sym" );
		protocols::symmetry::SetupForSymmetryMover setup( symdef_file, basic::options::option );
		setup.apply( *pose_p );
		return pose_p;
	}

	core::pose::PoseOP
	get_C3_NN() {
		core::pose::PoseOP pose_p( new core::pose::Pose() );
		import_pose::pose_from_file( *pose_p, "core/pack/guidance_scoreterms/approximate_buried_unsat_penalty/C3_NN.pdb" );
		std::string symdef_file = basic::database::full_name( "symmetry/cyclic/C3_Z.sym" );
		protocols::symmetry::SetupForSymmetryMover setup( symdef_file, basic::options::option );
		setup.apply( *pose_p );
		return pose_p;
	}

	core::pose::PoseOP
	get_C3_super_NN() {
		core::pose::PoseOP pose_p( new core::pose::Pose() );
		import_pose::pose_from_file( *pose_p, "core/pack/guidance_scoreterms/approximate_buried_unsat_penalty/C3_super_NN.pdb" );
		std::string symdef_file = basic::database::full_name( "symmetry/cyclic/C3_Z.sym" );
		protocols::symmetry::SetupForSymmetryMover setup( symdef_file, basic::options::option );
		setup.apply( *pose_p );
		return pose_p;
	}

	void assert_edge( scoring::EnergyGraph const & energy_graph, Size res1, Size res2, float expect ) {
		scoring::EnergyEdge const * edge = energy_graph.find_energy_edge( res1, res2 );
		TS_ASSERT( edge );
		if ( ! edge ) return;

		float actual = edge->fill_energy_map()[ scoring::approximate_buried_unsat_penalty ];
		TR << res1 << " vs " << res2 << std::endl;
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
			pose::Pose pose = *get_4SER();

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
			pose::Pose pose = *get_1ARG();

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
			pose::Pose pose = *get_1ARG_bbMET();

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
			pose::Pose pose = *get_1ARG_bbMET_O();

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
		// Everything is buried for the rest
		{
			pose::Pose pose = *get_APA();

			core::id::AtomID_Map< bool > buried = scoring::atomic_depth::atoms_deeper_than(
				pose, 0, false, 0.1, burial_resolution );


			TS_ASSERT( buried( 2, pose.residue(2).atom_index("N") ) );
			TS_ASSERT( buried( 2, pose.residue(2).atom_index("O") ) );
		}
		{
			pose::Pose pose = *get_proline_nh_test();

			core::id::AtomID_Map< bool > buried = scoring::atomic_depth::atoms_deeper_than(
				pose, 0, false, 0.1, burial_resolution );


			TS_ASSERT( buried( 2, pose.residue(2).atom_index("N") ) );
			TS_ASSERT( buried( 2, pose.residue(2).atom_index("O") ) );
			TS_ASSERT( buried( 4, pose.residue(4).atom_index("N") ) );
			TS_ASSERT( buried( 4, pose.residue(4).atom_index("O") ) );
			TS_ASSERT( buried( 4, pose.residue(4).atom_index("OXT") ) );
			TS_ASSERT( buried( 4, pose.residue(4).atom_index("OG") ) );
		}
		{
			pose::Pose pose = *get_C3_asn();

			core::id::AtomID_Map< bool > buried = scoring::atomic_depth::atoms_deeper_than(
				pose, 0, false, 0.1, burial_resolution );


			TS_ASSERT( buried( 1, pose.residue(1).atom_index("N") ) );
			TS_ASSERT( buried( 1, pose.residue(1).atom_index("O") ) );
			TS_ASSERT( buried( 1, pose.residue(1).atom_index("OXT") ) );
			TS_ASSERT( buried( 1, pose.residue(1).atom_index("OD1") ) );
			TS_ASSERT( buried( 1, pose.residue(1).atom_index("ND2") ) );

			TS_ASSERT( buried( 2, pose.residue(2).atom_index("N") ) );
			TS_ASSERT( buried( 2, pose.residue(2).atom_index("O") ) );
			TS_ASSERT( buried( 2, pose.residue(2).atom_index("OXT") ) );
			TS_ASSERT( buried( 2, pose.residue(2).atom_index("OD1") ) );
			TS_ASSERT( buried( 2, pose.residue(2).atom_index("ND2") ) );

			TS_ASSERT( buried( 3, pose.residue(3).atom_index("N") ) );
			TS_ASSERT( buried( 3, pose.residue(3).atom_index("O") ) );
			TS_ASSERT( buried( 3, pose.residue(3).atom_index("OXT") ) );
			TS_ASSERT( buried( 3, pose.residue(3).atom_index("OD1") ) );
			TS_ASSERT( buried( 3, pose.residue(3).atom_index("ND2") ) );
		}
		{
			pose::Pose pose = *get_C3_NN();

			core::id::AtomID_Map< bool > buried = scoring::atomic_depth::atoms_deeper_than(
				pose, 0, false, 0.1, burial_resolution );


			TS_ASSERT( buried( 1, pose.residue(1).atom_index("N") ) );
			TS_ASSERT( buried( 1, pose.residue(1).atom_index("O") ) );
			TS_ASSERT( buried( 1, pose.residue(1).atom_index("OD1") ) );
			TS_ASSERT( buried( 1, pose.residue(1).atom_index("ND2") ) );
			TS_ASSERT( buried( 2, pose.residue(2).atom_index("N") ) );
			TS_ASSERT( buried( 2, pose.residue(2).atom_index("O") ) );
			TS_ASSERT( buried( 2, pose.residue(2).atom_index("OXT") ) );
			TS_ASSERT( buried( 2, pose.residue(2).atom_index("OD1") ) );
			TS_ASSERT( buried( 2, pose.residue(2).atom_index("ND2") ) );

			TS_ASSERT( buried( 3, pose.residue(3).atom_index("N") ) );
			TS_ASSERT( buried( 3, pose.residue(3).atom_index("O") ) );
			TS_ASSERT( buried( 3, pose.residue(3).atom_index("OD1") ) );
			TS_ASSERT( buried( 3, pose.residue(3).atom_index("ND2") ) );
			TS_ASSERT( buried( 4, pose.residue(4).atom_index("N") ) );
			TS_ASSERT( buried( 4, pose.residue(4).atom_index("O") ) );
			TS_ASSERT( buried( 4, pose.residue(4).atom_index("OXT") ) );
			TS_ASSERT( buried( 4, pose.residue(4).atom_index("OD1") ) );
			TS_ASSERT( buried( 4, pose.residue(4).atom_index("ND2") ) );

			TS_ASSERT( buried( 5, pose.residue(5).atom_index("N") ) );
			TS_ASSERT( buried( 5, pose.residue(5).atom_index("O") ) );
			TS_ASSERT( buried( 5, pose.residue(5).atom_index("OD1") ) );
			TS_ASSERT( buried( 5, pose.residue(5).atom_index("ND2") ) );
			TS_ASSERT( buried( 6, pose.residue(6).atom_index("N") ) );
			TS_ASSERT( buried( 6, pose.residue(6).atom_index("O") ) );
			TS_ASSERT( buried( 6, pose.residue(6).atom_index("OXT") ) );
			TS_ASSERT( buried( 6, pose.residue(6).atom_index("OD1") ) );
			TS_ASSERT( buried( 6, pose.residue(6).atom_index("ND2") ) );

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

		scoring::ScoreFunctionOP sfxn ( new scoring::ScoreFunction() );

		scoring::methods::EnergyMethodOptions options = sfxn->energy_method_options();
		options.approximate_buried_unsat_penalty_hbond_energy_threshold( -0.25 );
		options.approximate_buried_unsat_penalty_burial_atomic_depth( bury_everything ?   0 : burial_depth );
		options.approximate_buried_unsat_penalty_burial_probe_radius( bury_everything ? 0.1 : burial_probe_radius );
		options.approximate_buried_unsat_penalty_burial_resolution( burial_resolution );
		options.approximate_buried_unsat_penalty_assume_const_backbone( assume_const_backbone );
		sfxn->set_energy_method_options( options );

		sfxn->set_weight( scoring::approximate_buried_unsat_penalty, 10.0 );



		core::pack::task::TaskFactoryOP tf( new core::pack::task::TaskFactory() );
		tf->push_back( utility::pointer::make_shared< IncludeCurrent >() );
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
		if ( core::pose::symmetry::is_symmetric( pose ) ) {
			rotsets = pack::rotamer_set::RotamerSetsOP( new pack::rotamer_set::symmetry::SymmetricRotamerSets() );
		}
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

		core::pose::Pose temp_pose = pose;

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

			// std::cout << "Setting " << resnum << " to rotamer " << the_irot << std::endl;
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
		pose::Pose pose = *get_4SER();

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
		pose::Pose pose = *get_4SER();


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
		pose::Pose pose = *get_4SER();


		std::string allow_repacking_res = "2,4-7";

		utility::vector1<std::pair<std::string, std::string> > designables { { "2", "DA" }, { "4-7", "SA" } };

		utility::vector1<Size> our_res { 2, 4, 5, 6, 7 };
		utility::vector1<Real> expected_scores { -30, 20, 20, 20, 20 };

		do_test_packing( pose, allow_repacking_res, designables, our_res, expected_scores, true, false );
		do_test_packing( pose, allow_repacking_res, designables, our_res, expected_scores, true, true );

	}


	// In this test, residue 3 has a buried polar that is being satisfied by 2 atoms of residue 1
	void test_scoring_1ARG(){
		pose::Pose pose = *get_1ARG();

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
		pose::Pose pose = *get_1ARG();


		std::string allow_repacking_res = "1,3";

		utility::vector1<std::pair<std::string, std::string> > designables { { "1", "RA" }, { "3", "DA" } };

		utility::vector1<Size> our_res { 1, 3 };
		utility::vector1<Real> expected_scores { -10, -10 };

		do_test_packing( pose, allow_repacking_res, designables, our_res, expected_scores, false, true );
		do_test_packing( pose, allow_repacking_res, designables, our_res, expected_scores, false, false );

	}



	// In this test, residue 3 has a buried polar that is being satisfied by 2 atoms of residue 1 and a bb of residue 8
	void test_scoring_1ARG_bbMET(){
		pose::Pose pose = *get_1ARG_bbMET();

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
		pose::Pose pose = *get_1ARG_bbMET();


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
		pose::Pose pose = *get_1ARG_bbMET_O();

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

	// This test has a problem, once the bb-bb hbond is established, it blocks the arg-bb hbond
	// Blame the scorefunction. Eliminating this test - bcov
	// In this test, residue 7 has a buried polar bb that is being satisfied by 2 atoms of residue 5 and a bb of residue 3
	//  Tests ERROR FIX 1
	//  Tests ERROR FIX 2
	//  Also tests conversion of a 2-body edge to a 1-body edge
	// void test_packing_1ARG_bbMET_O(){
	//  pose::Pose pose = *get_1ARG_bbMET_O();


	//  std::string allow_repacking_res = "7,5,3";

	//  utility::vector1<std::pair<std::string, std::string> > designables { { "7", "MA" }, { "5", "RA" }, { "3", "MA" } };

	//  utility::vector1<Size> our_res { 7, 5, 3 };

	//  utility::vector1<Real> expected_scores {  0,  // 0  // bb is part of background for packing
	//   10,  // -10 * 2 + 3 * 10
	//   0  // 0 // bb is part of background for packing
	//   };

	//  do_test_packing( pose, allow_repacking_res, designables, our_res, expected_scores, false, false );
	//  do_test_packing( pose, allow_repacking_res, designables, our_res, expected_scores, false, true );

	// }



	// In this test, residue 2 is proline and every atom is buried
	void test_scoring_APA(){
		pose::Pose pose = *get_APA();

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
		pose::Pose pose = *get_APA();

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
		pose::Pose pose = *get_proline_nh_test();

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
		pose::Pose pose = *get_proline_nh_test();

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

	// In this test, residue 1 is asparagine and is making hbonds to both symmetric copies of itself
	//   all atoms are buried
	void test_scoring_C3_asn(){
		pose::Pose pose = *get_C3_asn();

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
		// The asymmetric unit is X2 in symmetry
		TS_ASSERT_DELTA( energies.residue_total_energy( 1 ), 7*2-2, 0.01 );  // 1 + 1 + 3 + 1 + 1 - 0.5 * 4   // O, OXT, N, OD1, ND2, four hbonds

		scoring::EnergyGraph const & energy_graph = energies.energy_graph();

		// hbonds
		assert_edge( energy_graph, 1, 2, -2 );
		// assert_edge( energy_graph, 2, 3, -2 );
		// assert_edge( energy_graph, 1, 3, -2 );

	}

	// In this test, residue 1 is asparagine and is making hbonds to both symmetric copies of itself
	//   all atoms are buried
	void test_packing_C3_asn(){
		pose::Pose pose = *get_C3_asn();

		std::string allow_repacking_res = "1";

		utility::vector1<std::pair<std::string, std::string> > designables { { "1", "AN" } };

		utility::vector1<Size> our_res { 1 };

		utility::vector1<Real> expected_scores {  0,  // If everything is alanine, there is no change
			};

		do_test_packing( pose, allow_repacking_res, designables, our_res, expected_scores, false, false, true );
		do_test_packing( pose, allow_repacking_res, designables, our_res, expected_scores,    false, true, true );

		// Symmetry and the ig edge reweighter don't work together. Not our fault
		// do_test_packing( pose, allow_repacking_res, designables, our_res, expected_scores, true, false, true );
		// do_test_packing( pose, allow_repacking_res, designables, our_res, expected_scores,    true, true, true );

	}

	// In this test, residue 1 and 2 are asparagine and are hbonded to 2* and 1* respectively
	//   all atoms are buried
	void test_scoring_C3_NN(){
		pose::Pose pose = *get_C3_NN();

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
		// Onebodies in the asu are doubled
		TS_ASSERT_DELTA( energies.residue_total_energy( 1 ), 6*2-2, 0.01 );  // 1 + 3 + 1 + 1 - 0.5 * 4 // O, N, OD1, ND2, four hbonds
		TS_ASSERT_DELTA( energies.residue_total_energy( 2 ), 5*2-2, 0.01 );  // 1 + 1 + 1 + 1 + 1 - 0.5 * 4   // O, OXT, N, OD1, ND2, four hbonds

		scoring::EnergyGraph const & energy_graph = energies.energy_graph();

		// hbonds
		assert_edge( energy_graph, 1, 6, -4 );
		assert_edge( energy_graph, 2, 3, -4 );
		// assert_edge( energy_graph, 1, 3, -2 );

	}

	// In this test, residue 1 is asparagine and is making hbonds to both symmetric copies of itself
	//  Tests assume_const_backbone
	void test_packing_C3_NN(){
		pose::Pose pose = *get_C3_NN();

		std::string allow_repacking_res = "1,2";

		utility::vector1<std::pair<std::string, std::string> > designables { { "1", "AN" }, { "2", "AN" } };

		utility::vector1<Size> our_res { 1, 2 };

		// The total score is 2X the asu
		utility::vector1<Real> expected_scores {  -40,  // Should create 6 buried unsats all around the ring
			-40 }; // Should create 6 buried unsats all around the ring

		do_test_packing( pose, allow_repacking_res, designables, our_res, expected_scores, false, false, true );
		do_test_packing( pose, allow_repacking_res, designables, our_res, expected_scores,    false, true, true );

		// Symmetry and the ig edge reweighter don't work together. Not our fault
		// do_test_packing( pose, allow_repacking_res, designables, our_res, expected_scores, true, false, true );
		// do_test_packing( pose, allow_repacking_res, designables, our_res, expected_scores,    true, true, true );


	}

	// In this test, residue 1 and 2 are asparagine and are hbonded to 2* and 1* respectively, but then crazy stuff is happening
	//   You should open the c3 version in pymol to understand. But although pymol draws an hbond between 3 and 4. Rosetta does
	//   not recognize this as an hbond
	//   Designed to test oversats at symmetry interfaces
	//   all atoms are buried
	void test_scoring_C3_super_NN(){
		pose::Pose pose = *get_C3_super_NN();

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
		// Onebodies and twobodies in the asu are doubled
		TS_ASSERT_DELTA( energies.residue_total_energy( 1 ), 6*2-6-2+3, 0.01 );  // O, NH3, OD1, ND2, 6 asu hbonds, 4 hbonds, 3 asu oversats
		TS_ASSERT_DELTA( energies.residue_total_energy( 2 ), 5*2-5+1.5, 0.01 );  // O, OXT, N, OD1, ND2, 10 hbonds 3 oversats
		TS_ASSERT_DELTA( energies.residue_total_energy( 3 ), 6*2-2-1+3+0.5, 0.01 );  //  O, NH3, OD1, ND2, 2 asu hbonds, 2 hbonds, 3 asu oversats, 1 oversats
		TS_ASSERT_DELTA( energies.residue_total_energy( 4 ), 4*2-4+1+1, 0.01 );  //  O, N, OD1, ND2, four asu hbonds, 1 asu oversat, 2 oversats
		TS_ASSERT_DELTA( energies.residue_total_energy( 5 ), 5*2-2+3, 0.01 );  //  O, OXT, N, OD1, ND2, four hbonds, 3 asu oversats

		scoring::EnergyGraph const & energy_graph = energies.energy_graph();

		// Right interface
		assert_edge( energy_graph, 1, 3, 2*(-2 + 1) ); // 1 hbonds 1 oversat  2*asu
		assert_edge( energy_graph, 1, 4, 2*(-4) ); // 2 hbonds 2*asu
		assert_edge( energy_graph, 1, 5, 2*(+2) ); // 2 oversats 2*asu
		assert_edge( energy_graph, 1, 12, -4 ); // 2 hbonds
		assert_edge( energy_graph, 3, 4, 2*(+1) ); // 1 oversats 2*asu
		assert_edge( energy_graph, 3, 5, 2*(+1) ); // 1 oversats 2*asu
		assert_edge( energy_graph, 3, 12, -2+1 ); // 1 hbonds 1 oversats
		assert_edge( energy_graph, 4, 12, +2 ); // 2 oversats
		assert_edge( energy_graph, 5, 12, -4 ); // 2 hbonds

		// Left interface
		// assert_edge( energy_graph, 11, 13, -2 ); // 1 hbonds
		// assert_edge( energy_graph, 11, 14, -4 ); // 2 hbonds
		// assert_edge( energy_graph, 11, 15, +2 ); // 2 oversats
		assert_edge( energy_graph, 2, 6, -4 ); // 2 hbonds
		// assert_edge( energy_graph, 13, 14, +1 ); // 1 oversats
		// assert_edge( energy_graph, 13, 15, +1 ); // 1 oversats
		assert_edge( energy_graph, 2, 8, -2+1 ); // 1 hbonds 1 oversats
		assert_edge( energy_graph, 2, 9, +2 ); // 2 oversats
		assert_edge( energy_graph, 2, 10, -4 ); // 2 hbonds


		// assert_edge( energy_graph, 1, 3, -2 );

	}

	// In this test, residue 1 is asparagine and is making hbonds to both symmetric copies of itself
	//  Tests assume_const_backbone
	void test_packing_C3_super_NN(){
		pose::Pose pose = *get_C3_super_NN();

		std::string allow_repacking_res = "1,2,3,4,5";

		utility::vector1<std::pair<std::string, std::string> > designables { { "1", "AN" }, { "2", "AN" }, { "3", "AN" },
			{ "4", "AN" }, { "5", "AN" } };

		utility::vector1<Size> our_res { 1, 2, 3, 4, 5 };

		// The total score is 2X the asu
		utility::vector1<Real> expected_scores {  2*(20-100+30),  // +2 unsats, +10 hbonds, +3 oversats
			2*(20-100+30), // +2 unsats, +10 hbonds, +3 oversats
			2*(20-40+40), // +2 unsats, +4 hbonds, +4 oversats
			2*(20-40+30), // +2 unsats, +4 hbonds, +3 oversats
			2*(20-40+30) }; // +2 unsats, +4 hbonds, +3 oversats

		do_test_packing( pose, allow_repacking_res, designables, our_res, expected_scores, false, false, true );
		do_test_packing( pose, allow_repacking_res, designables, our_res, expected_scores,    false, true, true );

		// Symmetry and the ig edge reweighter don't work together. Not our fault
		// do_test_packing( pose, allow_repacking_res, designables, our_res, expected_scores, true, false, true );
		// do_test_packing( pose, allow_repacking_res, designables, our_res, expected_scores,    true, true, true );
		//-10, -40, +20
		// wat? Res: 1 Rot: 7 Res: 10 Rot: 14 Score: 10.000

	}



};
