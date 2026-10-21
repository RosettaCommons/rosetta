// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/beta_barrel/MakeBarrelTests.cxxtest.hh
/// @brief  Unit tests for the MakeBarrel mover.
/// @author Andy Watkins

#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/rosettascripts.hh>

#include <protocols/beta_barrel/MakeBarrel.hh>
#include <protocols/beta_barrel/BarrelParametrizationCalculator.hh>
#include <protocols/beta_barrel/parameters/BarrelParameters.hh>
#include <protocols/beta_barrel/parameters/BarrelParametersSet.hh>
#include <core/conformation/parametric/RealValuedParameter.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>

#include <utility/tag/Tag.fwd.hh>
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>

#include <cmath>

static basic::Tracer TR("MakeBarrelTests");

class MakeBarrelTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void tearDown() {
	}

	void test_make_8_strand_parallel_barrel() {
		TR << "Testing 8-strand all-parallel barrel construction." << std::endl;

		utility::tag::TagCOP tag( tagptr_from_string(
			"<MakeBarrel name=\"make\" n_strands=\"8\" shear_number=\"8\" "
			"r0=\"8.0\" antiparallel=\"false\" use_degrees=\"false\" "
			"crick_params_file=\"beta_strand\" residue_name=\"ALA\" "
			"strand_length=\"6\" reset=\"true\" />"
		) );

		basic::datacache::DataMap dummy_data;
		core::pose::Pose testpose;

		protocols::beta_barrel::MakeBarrel makebarrel;
		makebarrel.parse_my_tag( tag, dummy_data );
		makebarrel.apply( testpose );

		// Should have 8 strands * 6 residues = 48 residues
		TS_ASSERT_EQUALS( testpose.total_residue(), 48 );

		// Should have a BarrelParametersSet
		TS_ASSERT( testpose.conformation().n_parameters_sets() >= 1 );

		auto paramset = utility::pointer::dynamic_pointer_cast<
			protocols::beta_barrel::parameters::BarrelParametersSet const >(
			testpose.conformation().parameters_set(1) );
		TS_ASSERT( paramset != nullptr );

		if ( paramset ) {
			TS_ASSERT_EQUALS( paramset->n_strands(), 8 );
			TS_ASSERT_EQUALS( paramset->shear_number(), 8 );
			TS_ASSERT_EQUALS( paramset->antiparallel(), false );
			TS_ASSERT_DELTA( paramset->barrel_radius(), 8.0, 1e-6 );
		}
	}

	void test_make_barrel_has_cylindrical_geometry() {
		TR << "Testing that barrel has approximately cylindrical geometry." << std::endl;

		utility::tag::TagCOP tag( tagptr_from_string(
			"<MakeBarrel name=\"make\" n_strands=\"8\" shear_number=\"8\" "
			"r0=\"10.0\" antiparallel=\"false\" use_degrees=\"false\" "
			"crick_params_file=\"beta_strand\" residue_name=\"ALA\" "
			"strand_length=\"6\" reset=\"true\" />"
		) );

		basic::datacache::DataMap dummy_data;
		core::pose::Pose testpose;

		protocols::beta_barrel::MakeBarrel makebarrel;
		makebarrel.parse_my_tag( tag, dummy_data );
		makebarrel.apply( testpose );

		// Check that CA atoms of residue 1 in each strand are approximately at barrel radius
		// from the z-axis. Strand i starts at residue (i-1)*6 + 1.
		for ( core::Size strand = 1; strand <= 8; ++strand ) {
			core::Size resi = ( strand - 1 ) * 6 + 1;
			numeric::xyzVector< core::Real > ca_pos = testpose.residue( resi ).xyz( "CA" );
			core::Real dist_from_axis = std::sqrt( ca_pos.x() * ca_pos.x() + ca_pos.y() * ca_pos.y() );
			// Should be approximately at r0=10.0, within a few Angstroms (minor helix radius)
			TS_ASSERT_LESS_THAN( std::abs( dist_from_axis - 10.0 ), 4.0 );
		}
	}

};
