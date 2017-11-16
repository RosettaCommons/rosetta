// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/energy_based_clustering/EnergyBasedClusteringTests_oligourea.cxxtest.hh
/// @brief  Unit tests for the EnergyBasedClusteringProtocol, clustering poses containing oligoureas.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers

// Core Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/import_pose/import_pose.hh>

// Protocols Headers
#include <protocols/energy_based_clustering/EnergyBasedClusteringProtocol.hh>
#include <protocols/energy_based_clustering/EnergyBasedClusteringOptions.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("EnergyBasedClusteringTests");

class EnergyBasedClusteringTests_oligourea : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init_with_additional_options( "-energy_based_clustering:use_CB false -in:file:silent protocols/energy_based_clustering/oligourea.silent -in:file:fullatom -force_silent_bitflip_on_read -symmetric_gly_tables" );

	}

	void tearDown(){

	}

	/// @brief Test a Cartesian clustering run with alpha-amino acids, oligoureas and cyclic geometry.
	void test_cartesian_cluster_cyclic_oligourea() {
		protocols::energy_based_clustering::EnergyBasedClusteringOptions options(true);
		TS_ASSERT( !options.use_CB_ );
		options.cyclic_ = true;
		options.cluster_by_ = protocols::energy_based_clustering::EBC_bb_cartesian;
		options.cluster_radius_ = 1.0;
		options.use_CB_ = true;
		options.silent_output_ = true;
		options.output_prefix_ = "cyc_oligourea_cart";

		protocols::energy_based_clustering::EnergyBasedClusteringProtocol clusterer(options); //Load from options object, rather than from global options system.
		clusterer.go();
		TS_ASSERT_EQUALS(clusterer.n_clusters_from_last_run(), 283);
	}

	void test_cartesian_cluster_cyclic_oligourea_s2_symm() {
		protocols::energy_based_clustering::EnergyBasedClusteringOptions options(true);
		TS_ASSERT( !options.use_CB_ );
		options.cyclic_ = true;
		options.cluster_by_ = protocols::energy_based_clustering::EBC_bb_cartesian;
		options.cluster_radius_ = 1.0;
		options.use_CB_ = true;
		options.silent_output_ = true;
		options.output_prefix_ = "cyc_oligourea_cart_symm";
		options.cyclic_symmetry_ = 2;
		options.cyclic_symmetry_mirroring_ = true;
		options.cyclic_symmetry_threshold_ = 50;
		options.cyclic_symmetry_threshold_specified_ = true;

		protocols::energy_based_clustering::EnergyBasedClusteringProtocol clusterer(options); //Load from options object, rather than from global options system.
		clusterer.go();
		TS_ASSERT_EQUALS(clusterer.n_clusters_from_last_run(), 283);
	}

	/// @brief Test a Cartesian clustering run with alpha-amino acids, oligoureas and cyclic geometry.
	void test_dihedral_cluster_cyclic_oligourea() {
		protocols::energy_based_clustering::EnergyBasedClusteringOptions options(true);
		TS_ASSERT( !options.use_CB_ );
		options.cyclic_ = true;
		options.cluster_by_ = protocols::energy_based_clustering::EBC_bb_dihedral;
		options.cluster_radius_ = 250.0;
		options.use_CB_ = true;
		options.silent_output_ = true;
		options.output_prefix_ = "cyc_oligourea_dihed";

		protocols::energy_based_clustering::EnergyBasedClusteringProtocol clusterer(options); //Load from options object, rather than from global options system.
		clusterer.go();
		TS_ASSERT_EQUALS(clusterer.n_clusters_from_last_run(), 250);
	}

};
