// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/helical_bundle_predict/HBPHelixAssignmentsTests.cxxtest.hh
/// @brief  Unit tests for the class that stores helix/coil transition assignments during helical bundle structure prediction.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers
#include <protocols/helical_bundle_predict/HBPHelixAssignments.hh>
#include <protocols/helical_bundle/BundleParametrizationCalculator.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>
#include <utility/pointer/memory.hh>
#include <utility/string_util.hh>

static basic::Tracer TR("HBPHelixAssignmentsTests");


class HBPHelixAssignmentsTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init();

	}

	void tearDown(){

	}


	/// @brief Test the initialization of a HBPHelixAssignments object from a file.
	void test_init_from_file(){
		using namespace protocols::helical_bundle_predict;
		HBPHelixAssignmentsOP hass ( utility::pointer::make_shared< HBPHelixAssignments >() );

		hass->initialize_from_file_contents( utility::file_contents( "protocols/helical_bundle_predict/three_helix.helix_assignments" ) );

		TS_ASSERT_EQUALS(hass->common_r0_, false);
		TS_ASSERT_EQUALS(hass->common_omega0_, true);
		TS_ASSERT_EQUALS(hass->common_delta_omega1_, false);
		TS_ASSERT_EQUALS(hass->confine_to_user_defined_helices_, true);
		TS_ASSERT_DELTA( hass->global_r0_min_, 5.0, 1e-6 );
		TS_ASSERT_DELTA( hass->global_r0_max_, 8.0, 1e-6 );
		TS_ASSERT_DELTA( hass->global_omega0_min_, -3.0/180.0*3.141592654, 1e-6 );
		TS_ASSERT_DELTA( hass->global_omega0_max_, 3.0/180.0*3.141592654, 1e-6 );
		TS_ASSERT_DELTA( hass->global_delta_omega1_min_, 0, 1e-6 );
		TS_ASSERT_DELTA( hass->global_delta_omega1_max_, 3.141592654, 1e-6 );
		TS_ASSERT_DELTA( hass->nucleation_prob_, 0.01, 1e-6 );
		TS_ASSERT_DELTA( hass->extension_prob_, 0.05, 1e-6 );
		TS_ASSERT_DELTA( hass->retraction_prob_, 0.03, 1e-6 );
		TS_ASSERT( hass->global_bundle_calculator_ != nullptr );
		TS_ASSERT_EQUALS( hass->global_bundle_calculator_->residues_per_repeat(), 1 );
		TS_ASSERT_EQUALS( hass->global_bundle_calculator_->residues_per_turn(), 4 );
		TS_ASSERT_EQUALS( hass->global_bundle_calculator_->repeats_per_turn(), 4 );

		TS_ASSERT_EQUALS( hass->helices_.size(), 3 );
		HBPHelixCOP helix1( hass->helices_[1] );
		HBPHelixCOP helix2( hass->helices_[2] );
		HBPHelixCOP helix3( hass->helices_[3] );
		TS_ASSERT_DIFFERS( helix1, nullptr );
		TS_ASSERT_DIFFERS( helix2, nullptr );
		TS_ASSERT_DIFFERS( helix3, nullptr );

		TS_ASSERT_EQUALS( helix1->start_position_, 3 );
		TS_ASSERT_EQUALS( helix1->end_position_, 13 );
		TS_ASSERT_DELTA( helix1->nucleation_prob_, 0.01, 1e-6 );
		TS_ASSERT_DELTA( helix1->extension_prob_, 0.05, 1e-6 );
		TS_ASSERT_DELTA( helix1->retraction_prob_, 0.03, 1e-6 );

		TS_ASSERT_EQUALS( helix2->start_position_, 17 );
		TS_ASSERT_EQUALS( helix2->end_position_, 30 );
		TS_ASSERT_DELTA( helix2->nucleation_prob_, 0.01, 1e-6 );
		TS_ASSERT_DELTA( helix2->extension_prob_, 0.05, 1e-6 );
		TS_ASSERT_DELTA( helix2->retraction_prob_, 0.03, 1e-6 );

		TS_ASSERT_EQUALS( helix3->start_position_, 32 );
		TS_ASSERT_EQUALS( helix3->end_position_, 43 );
		TS_ASSERT_DELTA( helix3->nucleation_prob_, 0.02, 1e-6 );
		TS_ASSERT_DELTA( helix3->extension_prob_, 0.03, 1e-6 );
		TS_ASSERT_DELTA( helix3->retraction_prob_, 0.01, 1e-6 );

		HBPHelixParametersCOP params1( helix1->parameters() );
		HBPHelixParametersCOP params2( helix2->parameters() );
		HBPHelixParametersCOP params3( helix3->parameters() );
		TS_ASSERT( params1 != nullptr );
		TS_ASSERT( params2 != nullptr );
		TS_ASSERT( params3 != nullptr );

		TS_ASSERT_DELTA( params1->delta_omega1_range_.first, 0, 1e-6 );
		TS_ASSERT_DELTA( params1->delta_omega1_range_.second, 3.141592654, 1e-6 );
		TS_ASSERT_DELTA( params1->r0_range_.first, 5.0, 1e-6 );
		TS_ASSERT_DELTA( params1->r0_range_.second, 8.0, 1e-6 );
		TS_ASSERT_DELTA( params1->omega0_range_.first, -3.0/180.0*3.141592654, 1e-6 );
		TS_ASSERT_DELTA( params1->omega0_range_.second, 3.0/180.0*3.141592654, 1e-6 );
		TS_ASSERT_EQUALS( params1->bundle_calculator_->residues_per_repeat(), 1 );
		TS_ASSERT_EQUALS( params1->bundle_calculator_->residues_per_turn(), 4 );
		TS_ASSERT_EQUALS( params1->bundle_calculator_->repeats_per_turn(), 4 );

		TS_ASSERT_DELTA( params2->delta_omega1_range_.first, 45/180.0*3.141592654, 1e-6 );
		TS_ASSERT_DELTA( params2->delta_omega1_range_.second, 135/180.0*3.141592654, 1e-6 );
		TS_ASSERT_DELTA( params2->r0_range_.first, 6.0, 1e-6 );
		TS_ASSERT_DELTA( params2->r0_range_.second, 7.0, 1e-6 );
		TS_ASSERT_DELTA( params2->omega0_range_.first, -3.0/180.0*3.141592654, 1e-6 );
		TS_ASSERT_DELTA( params2->omega0_range_.second, 3.0/180.0*3.141592654, 1e-6 );
		TS_ASSERT_EQUALS( params2->bundle_calculator_->residues_per_repeat(), 1 );
		TS_ASSERT_EQUALS( params2->bundle_calculator_->residues_per_turn(), 4 );
		TS_ASSERT_EQUALS( params2->bundle_calculator_->repeats_per_turn(), 4 );

		TS_ASSERT_DELTA( params3->delta_omega1_range_.first, 0, 1e-6 );
		TS_ASSERT_DELTA( params3->delta_omega1_range_.second, 3.141592654, 1e-6 );
		TS_ASSERT_DELTA( params3->r0_range_.first, 5.0, 1e-6 );
		TS_ASSERT_DELTA( params3->r0_range_.second, 8.0, 1e-6 );
		TS_ASSERT_DELTA( params3->omega0_range_.first, -3.0/180.0*3.141592654, 1e-6 );
		TS_ASSERT_DELTA( params3->omega0_range_.second, 3.0/180.0*3.141592654, 1e-6 );
		TS_ASSERT_EQUALS( params3->bundle_calculator_->residues_per_repeat(), 1 );
		TS_ASSERT_EQUALS( params3->bundle_calculator_->residues_per_turn(), 4 );
		TS_ASSERT_EQUALS( params3->bundle_calculator_->repeats_per_turn(), 4 );

		/*TR << "\nGLOBAL CRICK PARAMS FILE:\n";
		TR << hass->global_crick_params_file_contents_ << "\n";

		TR << "\nHELIX1 CRICK PARAMS FILE:\n";
		TR << params1->crick_params_file_contents_ << "\n";

		TR << "\nHELIX2 CRICK PARAMS FILE:\n";
		TR << params2->crick_params_file_contents_ << "\n";

		TR << "\nHELIX3 CRICK PARAMS FILE:\n";
		TR << params3->crick_params_file_contents_ << "\n";
		*/
	}



};
