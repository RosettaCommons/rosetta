// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/MMBondAngleScore.cxxtest.hh
/// @brief  test suite for core:scoring::methods::MMBondAngleEnergy
/// @author Colin A. Smith (colin.smith@ucsf.edu)


// Test headers
#include <cxxtest/TestSuite.h>

#include <test/util/pose_funcs.hh>
#include <test/util/deriv_funcs.hh>
#include <test/UTracer.hh>

// Unit headers
#include <core/scoring/mm/MMBondAngleResidueTypeParam.hh>
#include <core/scoring/mm/MMBondAngleResidueTypeParamSet.hh>
#include <core/scoring/methods/MMBondAngleEnergy.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/types.hh>

// Project headers
#include <test/core/init_util.hh>

#include <core/kinematics/DomainMap.hh>
// Auto-header: duplicate removed #include <core/io/pdb/pose_io.hh>

// Utility headers

// C++ headers
#include <iostream>

//Auto Headers
#include <core/pose/annotated_sequence.hh>
#include <utility/vector1.hh>

using namespace core;
using namespace core::pose;
using namespace core::scoring;
using namespace core::scoring::mm;
using namespace core::scoring::methods;

Real
dump_energy_method_energies(
	Pose const & pose,
	ShortRangeTwoBodyEnergy const & emethod,
	EnergyMap const & weightmap,
	std::ostream * os = NULL
)
{
	ScoreFunction sfxn; // unused
	EnergyMap emap;
	EnergyMap twobodyemap;
	core::kinematics::DomainMap domainmap; // unused
	core::Vector F1;
	core::Vector F2;

	Real total_energy(0);

	if ( os ) *os << "Intraresidue Energies:" << "\n";
	for ( Size i = 1; i <= pose.total_residue(); ++i ) {
		emap.zero();
		emethod.eval_intrares_energy( pose.residue(i), pose, sfxn, emap);
		//std::cout <<   std::setprecision(4) << std::fixed << std::setw(8) << emap[ mm_bend ] << ",";
		if ( os ) *os << i << "\t" << emap.dot(weightmap) << "\n";
		total_energy += emap.dot(weightmap);
	}

	if ( os ) *os << "Residue Pair Energies:" << "\n";
	for ( Size i = 1; i <= pose.total_residue(); ++i ) {
		if ( (i+1) <= pose.total_residue() ) {
			twobodyemap.zero();
			emethod.residue_pair_energy( pose.residue(i), pose.residue(i+1), pose, sfxn, twobodyemap);
			//std::cout <<   std::setprecision(4) << std::fixed << std::setw(8) << emap[ mm_twist ] << ",";
			if ( os ) *os << i << "-" << i+1 << "\t" << twobodyemap.dot(weightmap) << "\n";
			total_energy += twobodyemap.dot(weightmap);
		}
	}

	if ( os ) *os << "Atom Derivatives:" << "\n";
	for ( Size i = 1; i <= pose.total_residue(); ++i ) {
		for ( Size j = 1; j <= pose.residue(i).natoms(); ++j ) {
			F1.zero();
			F2.zero();
			emethod.eval_atom_derivative( core::id::AtomID(j, i), pose, domainmap, sfxn, weightmap, F1, F2 );
			if ( os ) {
				*os << i << "\t" << j
					<< "\tF1: " << F1.x() << " " << F1.x() << " " << F1.z()
					<< "\tF2: " << F2.x() << " " << F2.x() << " " << F2.z() << "\n";
			}
		}
	}

	return total_energy;
}


// --------------- Test Class --------------- //

class MMBondAngleEnergyTests : public CxxTest::TestSuite {

public:

	PoseOP pose;
	PoseOP pose_ideal;
	Real delta;

	// --------------- Suite-level Fixture --------------- //

	MMBondAngleEnergyTests() {
		core_init_with_additional_options( "-no_optH" );
	}

	virtual ~MMBondAngleEnergyTests() {}

	static MMBondAngleEnergyTests* createSuite() {
		return new MMBondAngleEnergyTests();
	}

	static void destroySuite( MMBondAngleEnergyTests *suite ) {
		delete suite;
	}

	// --------------- Fixtures --------------- //

	void setUp() {
		// init pose
		pose = create_test_in_pdb_poseop();
		//core::import_pose::pose_from_pdb( *pose, "core/scoring/methods/test_in.pdb" );
		pose_ideal = PoseOP( new Pose() );
		core::pose::make_pose_from_sequence(*pose_ideal, "ACDEFGHIKLMNQRSTVWY", "fa_standard");

		// init delta
		delta = 0.0001;

	}

	void tearDown() {
		pose.reset();
		pose_ideal.reset();
	}

	// --------------- Test Cases --------------- //

	void test_energy() {

		test::UTracer UT("core/scoring/methods/MMBondAngleEnergyTests.u");
		//std::ofstream UT("core/scoring/methods/MMBondAngleEnergyTests.u");

		EnergyMap weightmap;
		weightmap[ mm_bend ] = 1;

		EnergyMethodOptionsOP energy_method_options( new EnergyMethodOptions );

		energy_method_options->bond_angle_residue_type_param_set(scoring::mm::MMBondAngleResidueTypeParamSetOP( new core::scoring::mm::MMBondAngleResidueTypeParamSet() ));

		MMBondAngleEnergyOP mmbondangleenergy( new MMBondAngleEnergy(*energy_method_options) );

		dump_energy_method_energies(*pose, *mmbondangleenergy, weightmap, &UT);

		UT << "\nN-CA-C Central Atoms:" << "\n\n";

		utility::vector1<std::string> central_atoms;
		central_atoms.push_back("N");
		central_atoms.push_back("CA");

		EnergyMethodOptionsOP energy_method_options_n_ca( new EnergyMethodOptions );

		core::scoring::mm::MMBondAngleResidueTypeParamSetOP param_set_n_ca( new core::scoring::mm::MMBondAngleResidueTypeParamSet() );
		param_set_n_ca->central_atoms_to_score(central_atoms);
		energy_method_options_n_ca->bond_angle_residue_type_param_set(param_set_n_ca);

		MMBondAngleEnergyOP mmbondangleenergy_n_ca( new MMBondAngleEnergy(*energy_method_options_n_ca) );

		dump_energy_method_energies(*pose, *mmbondangleenergy_n_ca, weightmap, &UT);

		// test to make sure use ResidueType theta0 gives a nearly zero score

		EnergyMethodOptionsOP energy_method_options_rt_theta0( new EnergyMethodOptions );

		core::scoring::mm::MMBondAngleResidueTypeParamSetOP param_set_rt_theta0( new core::scoring::mm::MMBondAngleResidueTypeParamSet() );
		param_set_rt_theta0->use_residue_type_theta0(true);
		energy_method_options_rt_theta0->bond_angle_residue_type_param_set(param_set_rt_theta0);

		MMBondAngleEnergyOP mmbondangleenergy_rt_theta0( new MMBondAngleEnergy(*energy_method_options_rt_theta0) );

		TS_ASSERT(dump_energy_method_energies(*pose_ideal, *mmbondangleenergy_rt_theta0, weightmap) < 1e-20);
	}

	// output without MMBondAngleResidueTypeParamSet should be identical

	void test_energy_original() {

		test::UTracer UT("core/scoring/methods/MMBondAngleEnergyTests.u");
		//std::ofstream UT("core/scoring/methods/MMBondAngleEnergyTests.u");

		EnergyMap weightmap;
		weightmap[ mm_bend ] = 1;

		EnergyMethodOptionsOP energy_method_options( new EnergyMethodOptions );

		MMBondAngleEnergyOP mmbondangleenergy( new MMBondAngleEnergy(*energy_method_options) );

		dump_energy_method_energies(*pose, *mmbondangleenergy, weightmap, &UT);

		UT << "\nN-CA-C Central Atoms:" << "\n\n";

		utility::vector1<std::string> central_atoms;
		central_atoms.push_back("N");
		central_atoms.push_back("CA");

		EnergyMethodOptionsOP energy_method_options_n_ca( new EnergyMethodOptions );
		energy_method_options_n_ca->bond_angle_central_atoms_to_score(central_atoms);

		MMBondAngleEnergyOP mmbondangleenergy_n_ca( new MMBondAngleEnergy(*energy_method_options_n_ca) );

		dump_energy_method_energies(*pose, *mmbondangleenergy_n_ca, weightmap, &UT);
	}

	void test_mmbond_angle_energy_start_score_start_func_match_w_total_flexibility()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( mm_bend, 0.5 );
		kinematics::MoveMap movemap( create_movemap_to_allow_all_torsions() );
		AtomDerivValidator adv( pose, sfxn, movemap );
		adv.validate_start_func_matches_start_score( 42.85876726805235, false, 1e-6 );
	}


};
