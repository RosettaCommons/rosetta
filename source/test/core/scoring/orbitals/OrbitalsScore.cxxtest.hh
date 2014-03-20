// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


/// @author Steven Combs (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <core/scoring/orbitals/OrbitalsScore.hh>
#include <core/scoring/orbitals/OrbitalsLookup.hh>

#include <platform/types.hh>

// Package Headers
#include <test/util/pose_funcs.hh>
#include <test/util/deriv_funcs.hh>
#include <test/core/init_util.hh>

#include <core/id/DOF_ID.hh>
#include <core/id/TorsionID.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/io/pdb/pose_io.hh>
#include <basic/options/option.hh>

#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <numeric/conversions.hh>
#include <core/import_pose/import_pose.hh>
//Auto Headers
#include <core/conformation/Atom.fwd.hh>
#include <core/id/NamedStubID.hh>
#include <core/io/pdb/file_data.hh>
#include <core/pose/signals/ConformationEvent.hh>
#include <core/pose/signals/DestructionEvent.hh>
#include <core/pose/signals/EnergyEvent.hh>
#include <core/scoring/types.hh>
#include <core/scoring/constraints/Constraints.fwd.hh>
#include <ObjexxFCL/FArray.fwd.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/database/open.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "core.scoring.orbitals.OrbitalsScore_cxxtest_hh" );

// --------------- Test Class --------------- //

// using declarations
using namespace core;
using namespace core::pose;
using namespace core::scoring;
using namespace core::scoring::methods;

class OrbitalsEnergyTests : public CxxTest::TestSuite {

public:
	core::chemical::ResidueTypeSetOP orbitals_residue_type_set_;

	void setUp() {
		core_init_with_additional_options( "-in:add_orbitals -chemical:exclude_patches LowerDNA  UpperDNA Cterm_amidation SpecialRotamer protein_cutpoint_upper protein_cutpoint_lower VirtualBB ShoveBB VirtualDNAPhosphate VirtualNTerm CTermConnect sc_orbitals pro_hydroxylated_case1 pro_hydroxylated_case2 ser_phosphorylated thr_phosphorylated  tyr_phosphorylated tyr_sulfated lys_dimethylated lys_monomethylated  lys_trimethylated lys_acetylated glu_carboxylated cys_acetylated tyr_diiodinated N_acetylated C_methylamidated MethylatedProteinCterm" );
		if ( ! orbitals_residue_type_set_ ) {
			orbitals_residue_type_set_ = new core::chemical::ResidueTypeSet( "ORBITALS_FA_STANDARD", basic::database::full_name("chemical/residue_type_sets/fa_standard/") );
		}
	}
	void tearDown() {
		//basic::options::option[ basic::options::OptionKeys::in::add_orbitals ](false);
		//orbitals_residue_type_set_ = 0;
	}


	void test_orbital_scoring_function_values()
	{
		core::pose::Pose pose;
		core::import_pose::pose_from_pdbstring(pose, trp_cage_ideal(), *orbitals_residue_type_set_, "trp_cage");
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( pci_cation_pi, 1 );
		sfxn.set_weight( pci_pi_pi, 1 );
		sfxn.set_weight( orbitals_hpol_bb, 1);
		sfxn.set_weight( pci_salt_bridge, 1);
		sfxn.set_weight( pci_hbond, 1);
		Real start_score = sfxn(pose);
		TS_ASSERT_DELTA( -6.9185, start_score, 0.0001 );

	}

	void test_orbital_start_score_start_func_match_w_total_flexibility()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( pci_cation_pi, 1 );
		sfxn.set_weight( pci_pi_pi, 1 );
		sfxn.set_weight( orbitals_hpol_bb, 1);
		sfxn.set_weight( pci_salt_bridge, 1);
		sfxn.set_weight( pci_hbond, 1);
		kinematics::MoveMap movemap( create_movemap_to_allow_all_torsions() );
		AtomDerivValidator adv( pose, sfxn, movemap );
		adv.validate_start_func_matches_start_score( -6.918, false, 1e-3 );
	}


	void test_orbital_deriv_check_w_partial_flexibility()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( pci_cation_pi, 1 );
		sfxn.set_weight( pci_pi_pi, 1 );
		sfxn.set_weight( orbitals_hpol_bb, 1);
		sfxn.set_weight( pci_salt_bridge, 1);
		sfxn.set_weight( pci_hbond, 1);
		kinematics::MoveMap movemap( create_trpcage_movemap_to_allow_bb10_freedom() );
		AtomDerivValidator adv;
		adv.set_pose( pose );
		adv.set_score_function( sfxn );
		adv.set_movemap( movemap );
		adv.simple_deriv_check( false, 1e-4 );
	}


	void test_orbital_deriv_check_w_total_flexibility_all()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( pci_cation_pi, 1 );
		sfxn.set_weight( pci_pi_pi, 1 );
		sfxn.set_weight( orbitals_hpol_bb, 1);
		sfxn.set_weight( pci_salt_bridge, 1);
		sfxn.set_weight( pci_hbond, 1);
		TR << "test_orbital_deriv_check_w_total_flexibility_all: " << sfxn(pose) << std::endl;
		kinematics::MoveMap movemap( create_movemap_to_allow_all_torsions() );
		AtomDerivValidator adv( pose, sfxn, movemap );
		adv.simple_deriv_check( false, 1e-3 );
	}

	void test_orbital_deriv_check_w_total_flexibility_pci_cation_pi()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( pci_cation_pi, 1 );
		TR << "test_orbital_deriv_check_w_total_flexibility_pci_cation_pi: " << sfxn(pose) << std::endl;
		kinematics::MoveMap movemap( create_movemap_to_allow_all_torsions() );
		AtomDerivValidator adv( pose, sfxn, movemap );
		adv.simple_deriv_check( false, 1e-3 );
	}

	void test_orbital_deriv_check_w_total_flexibility_pci_pi_pi()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( pci_pi_pi, 1 );
		TR << "test_orbital_deriv_check_w_total_flexibility_pci_pi_pi: " << sfxn(pose) << std::endl;
		kinematics::MoveMap movemap( create_movemap_to_allow_all_torsions() );
		AtomDerivValidator adv( pose, sfxn, movemap );
		adv.simple_deriv_check( false, 1e-3 );
	}

	void test_orbital_deriv_check_w_total_flexibility_orbitals_hpol_bb()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( orbitals_hpol_bb, 1);
		TR << "test_orbital_deriv_check_w_total_flexibility_orbitals_hpol_bb: " << sfxn(pose) << std::endl;
		kinematics::MoveMap movemap( create_movemap_to_allow_all_torsions() );
		AtomDerivValidator adv( pose, sfxn, movemap );
		adv.simple_deriv_check( false, 1e-3 );
	}

	void test_orbital_deriv_check_w_total_flexibility_pci_salt_bridge()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( pci_salt_bridge, 1);
		TR << "test_orbital_deriv_check_w_total_flexibility_pci_salt_bridge: " << sfxn(pose) << std::endl;
		kinematics::MoveMap movemap( create_movemap_to_allow_all_torsions() );
		AtomDerivValidator adv( pose, sfxn, movemap );
		adv.simple_deriv_check( false, 1e-3 );
	}

	void test_orbital_deriv_check_w_total_flexibility_pci_hbond()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( pci_hbond, 1);
		TR << "test_orbital_deriv_check_w_total_flexibility_pci_hbond: " << sfxn(pose) << std::endl;
		kinematics::MoveMap movemap( create_movemap_to_allow_all_torsions() );
		AtomDerivValidator adv( pose, sfxn, movemap );
		adv.simple_deriv_check( false, 1e-3 );
	}

	void test_make_compiler_happy()
	{
		// without a single active test_xxx(), the compilar quits.
	}




};

