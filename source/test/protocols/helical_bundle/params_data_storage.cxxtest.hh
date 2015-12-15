// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/helical_bundle/params_data_storage.cxxtest.hh
/// @brief  Unit tests for the MakeBundle mover that check that the mover is properly storing Crick parameter information in the output Pose object's Conformation object.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Test headers:
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
//#include <test/util/pdb1ubq.hh>

// helical_bundle headers:
#include <protocols/helical_bundle/MakeBundle.hh>
#include <protocols/helical_bundle/parameters/BundleParameters.fwd.hh>
#include <protocols/helical_bundle/parameters/BundleParametersSet.fwd.hh>
#include <protocols/helical_bundle/parameters/BundleParameters.hh>
#include <protocols/helical_bundle/parameters/BundleParametersSet.hh>

// Other Rosetta libraries:
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// --------------- Test Class --------------- //

class ParamsDataStorageTests : public CxxTest::TestSuite {

private:
	core::pose::PoseOP testpose_;
	core::scoring::ScoreFunctionOP scorefxn_;

public:

	void setUp() {
		core_init();
		testpose_ = core::pose::PoseOP( new core::pose::Pose );
		scorefxn_ = core::scoring::get_score_function();
	}

	void tearDown() {
	}

	/// @brief Test storage of Crick parameters in MakeBundle output.
	/// @details Failure of this test indicates that the Crick parameters are not being
	/// properly stored in the output Conformation object.
	void test_params_data_storage()
	{
		using namespace protocols::helical_bundle;
		using namespace protocols::helical_bundle::parameters;

		TS_TRACE( "Starting test_params_data_storage unit test." );
		TS_TRACE( "Written by Vikram K. Mulligan, Baker Laboratory." );
		TS_TRACE( "If this test should fail, it means that Crick parameter data are not being properly stored in the Conformation object in the Pose that results from the MakeBundle mover." );

		MakeBundleOP makebundle ( new MakeBundle ); //Create the mover

		makebundle->set_reset_pose(true); //destroy the original pose.
		makebundle->set_symmetry(2); //Twofold symmetry.

		//Add two helices:
		//Note -- with twofold symmetry, a total of four will be created.
		makebundle->add_helix();
		makebundle->add_helix();

		//Set parameters for the two helices:
		TS_TRACE( "Defining two helices.  (There will be four total, with two-fold symmetry)." );
		TS_TRACE( "Helix1 r0=5.0 omega0=0.05 delta_omega0=0.1 invert=false(default)" );
		makebundle->helix(1)->set_r0(5.0);
		makebundle->helix(1)->set_omega0(0.05);
		makebundle->helix(1)->set_delta_omega0(0.1);
		TS_TRACE( "Helix2 r0=6.5 omega0=0.03 delta_omega0=1.5 invert=true" );
		makebundle->helix(2)->set_r0(6.5);
		makebundle->helix(2)->set_omega0(0.03);
		makebundle->helix(2)->set_delta_omega0(1.5);
		makebundle->helix(2)->set_invert_helix(true);

		//Apply the mover:
		makebundle->apply(*testpose_);

		char outbuffer [1024];
		TS_TRACE( "Reading Crick parameters from the output pose." );
		sprintf(outbuffer, "Number of Parameters objects in the output pose ParametersSet: %lu", testpose_->conformation().parameters_set(1)->n_parameters());
		TS_TRACE( outbuffer );
		for ( core::Size i=1; i<=4; ++i ) {
			BundleParametersOP h1params( utility::pointer::dynamic_pointer_cast<BundleParameters>( testpose_->conformation().parameters_set(1)->parameters(i) ) );
			TS_ASSERT(h1params);
			core::Real r0_1(h1params->r0());
			core::Real omega0_1(h1params->omega0());
			core::Real delta_omega0_1(h1params->delta_omega0());
			bool invert_1(h1params->invert_helix());
			sprintf(outbuffer, "helix %lu: r0=%.4f omega0=%.4f delta_omega0=%.4f invert=%s", i, r0_1, omega0_1, delta_omega0_1, (std::string( invert_1 ? "true" : "false" )).c_str() );
			TS_TRACE( outbuffer );

			if ( i==1 || i==3 ) {
				TS_ASSERT_DELTA( r0_1, 5.0, 1e-5  );
				TS_ASSERT_DELTA( omega0_1, 0.05, 1e-5  );
				TS_ASSERT_DELTA( delta_omega0_1, 0.1 + (i==3 ? 3.141592654 : 0), 1e-5  );
				TS_ASSERT( invert_1 == false );
			} else {
				TS_ASSERT_DELTA( r0_1, 6.5, 1e-5  );
				TS_ASSERT_DELTA( omega0_1, 0.03, 1e-5  );
				TS_ASSERT_DELTA( delta_omega0_1, 1.5 + (i==4 ? 3.141592654 : 0), 1e-5  );
				TS_ASSERT( invert_1 == true );
			}

		}

		//testpose_->dump_pdb("vtemp.pdb"); //DELETE ME

		TS_TRACE( "Finished test_params_data_storage unit test." );

		return;
	}

	/// @brief Test storage of residue linkages in MakeBundle output.
	/// @details Failure of this test indicates that the residue owning pointers are not being
	/// properly stored in the output Conformation object.
	void test_params_residue_storage()
	{
		using namespace protocols::helical_bundle;
		using namespace protocols::helical_bundle::parameters;

		TS_TRACE( "Starting test_params_residue_storage unit test." );
		TS_TRACE( "Written by Vikram K. Mulligan, Baker Laboratory." );
		TS_TRACE( "If this test should fail, it means that residue data are not being properly stored in the Conformation object in the Pose that results from the MakeBundle mover." );

		MakeBundleOP makebundle ( new MakeBundle ); //Create the mover

		makebundle->set_reset_pose(true); //destroy the original pose.
		makebundle->set_symmetry(2); //Twofold symmetry.

		makebundle->set_default_helix_length(10); //10 residues per helix.

		//Add two helices:
		//Note -- with twofold symmetry, a total of four will be created.
		makebundle->add_helix();
		makebundle->add_helix();

		//Set parameters for the two helices:
		TS_TRACE( "Defining two helices.  (There will be four total, with two-fold symmetry)." );
		TS_TRACE( "Helix1 r0=5.0 omega0=0.05 delta_omega0=0.1 invert=false(default)" );
		makebundle->helix(1)->set_r0(5.0);
		makebundle->helix(1)->set_omega0(0.05);
		makebundle->helix(1)->set_delta_omega0(0.1);
		TS_TRACE( "Helix2 r0=6.5 omega0=0.03 delta_omega0=1.5 invert=true" );
		makebundle->helix(2)->set_r0(6.5);
		makebundle->helix(2)->set_omega0(0.03);
		makebundle->helix(2)->set_delta_omega0(1.5);
		makebundle->helix(2)->set_invert_helix(true);
		makebundle->helix(2)->set_helix_length(12);

		//Apply the mover:
		makebundle->apply(*testpose_);

		char outbuffer [1024];
		TS_TRACE( "Reading residue numbers from the output pose Parameters objects." );
		sprintf(outbuffer, "Number of Parameters objects in the output pose ParametersSet: %lu", testpose_->conformation().parameters_set(1)->n_parameters());
		TS_TRACE( outbuffer );

		core::Size resind(1);

		for ( core::Size i=1; i<=4; ++i ) {
			BundleParametersOP h1params( utility::pointer::dynamic_pointer_cast<BundleParameters>( testpose_->conformation().parameters_set(1)->parameters(i) ) );
			TS_ASSERT(h1params);

			sprintf(outbuffer, "Helix %lu residues:", i);

			utility::vector1 < core::Size > residue_indices;

			for ( core::Size ir=1, irmax=h1params->n_residue(); ir<=irmax; ++ir ) {
				residue_indices.push_back( h1params->residue(ir)->seqpos() );
				sprintf(outbuffer, "%s %lu", outbuffer, residue_indices[residue_indices.size()]);
			}

			TS_TRACE(outbuffer);

			core::Size expected_helix_size(10);
			if ( i==2 || i==4 ) expected_helix_size=12;

			TS_ASSERT(residue_indices.size()==expected_helix_size);
			for ( core::Size ii=1; ii<=expected_helix_size; ++ii ) {
				TS_ASSERT(residue_indices[ii]==resind);
				resind++;
			}
		}

		//testpose_->dump_pdb("vtemp.pdb"); //DELETE ME

		TS_TRACE( "Finished test_params_residue_storage unit test." );

		return;
	}

	/// @brief Test storage of Crick parameters when copying poses.
	/// @details Failure of this test indicates that the Crick parameters are not being
	/// properly copied when poses are duplicated.
	void test_params_data_copying()
	{
		using namespace protocols::helical_bundle;
		using namespace protocols::helical_bundle::parameters;

		TS_TRACE( "Starting test_params_data_copying unit test." );
		TS_TRACE( "Written by Vikram K. Mulligan, Baker Laboratory." );
		TS_TRACE( "If this test should fail, it means that Crick parameter data are not being properly copied when a pose is copied." );

		MakeBundleOP makebundle ( new MakeBundle ); //Create the mover

		makebundle->set_reset_pose(true); //destroy the original pose.
		makebundle->set_symmetry(2); //Twofold symmetry.

		//Add two helices:
		//Note -- with twofold symmetry, a total of four will be created.
		makebundle->add_helix();
		makebundle->add_helix();

		//Set parameters for the two helices:
		TS_TRACE( "Defining two helices.  (There will be four total, with two-fold symmetry)." );
		TS_TRACE( "Helix1 r0=5.0 omega0=0.05 delta_omega0=0.1 z1_offset=0.2 invert=false(default)" );
		makebundle->helix(1)->set_r0(5.0);
		makebundle->helix(1)->set_omega0(0.05);
		makebundle->helix(1)->set_delta_omega0(0.1);
		makebundle->helix(1)->set_z1_offset(0.2);
		TS_TRACE( "Helix2 r0=6.5 omega0=0.03 delta_omega0=1.5 z0_offset=-0.2 invert=true" );
		makebundle->helix(2)->set_r0(6.5);
		makebundle->helix(2)->set_omega0(0.03);
		makebundle->helix(2)->set_delta_omega0(1.5);
		makebundle->helix(2)->set_invert_helix(true);
		makebundle->helix(2)->set_z0_offset(-0.2);

		//Apply the mover:
		makebundle->apply(*testpose_);

		//Make a clone of the pose:
		TS_TRACE( "Attempting to clone the pose." );
		core::pose::PoseOP poseclone(testpose_->clone());

		char outbuffer [1024];
		TS_TRACE( "Reading Crick parameters from the cloned pose." );
		sprintf(outbuffer, "Number of Parameters objects in the output pose ParametersSet: %lu", poseclone->conformation().parameters_set(1)->n_parameters());
		TS_TRACE( outbuffer );
		for ( core::Size i=1; i<=4; ++i ) {
			BundleParametersOP h1params( utility::pointer::dynamic_pointer_cast<BundleParameters>( poseclone->conformation().parameters_set(1)->parameters(i) ) );
			TS_ASSERT(h1params);
			core::Real r0_1(h1params->r0());
			core::Real omega0_1(h1params->omega0());
			core::Real delta_omega0_1(h1params->delta_omega0());
			core::Real z1_off( h1params->z1_offset() );
			core::Real z0_off( h1params->z0_offset() );
			bool invert_1(h1params->invert_helix());
			sprintf(outbuffer, "helix %lu: r0=%.4f omega0=%.4f delta_omega0=%.4f invert=%s z1_offset=%.4f z0_offset=%.4f", i, r0_1, omega0_1, delta_omega0_1, (std::string( invert_1 ? "true" : "false" )).c_str(), z1_off, z0_off );
			TS_TRACE( outbuffer );

			if ( i==1 || i==3 ) {
				TS_ASSERT_DELTA( r0_1, 5.0, 1e-5  );
				TS_ASSERT_DELTA( omega0_1, 0.05, 1e-5  );
				TS_ASSERT_DELTA( delta_omega0_1, 0.1 + (i==3 ? 3.141592654 : 0), 1e-5  );
				TS_ASSERT( invert_1 == false );
				TS_ASSERT_DELTA( z1_off, 0.2, 1e-5 );
				TS_ASSERT_DELTA( z0_off, 0, 1e-5 );
			} else {
				TS_ASSERT_DELTA( r0_1, 6.5, 1e-5  );
				TS_ASSERT_DELTA( omega0_1, 0.03, 1e-5  );
				TS_ASSERT_DELTA( delta_omega0_1, 1.5 + (i==4 ? 3.141592654 : 0), 1e-5  );
				TS_ASSERT( invert_1 == true );
				TS_ASSERT_DELTA( z1_off, 0, 1e-5 );
				TS_ASSERT_DELTA( z0_off, -0.2, 1e-5 );
			}

		}

		//Make a copy of the pose:
		TS_TRACE( "Attempting to copy the pose." );
		core::pose::Pose posecopy(*testpose_);

		TS_TRACE( "Reading Crick parameters from the copied pose." );
		sprintf(outbuffer, "Number of Parameters objects in the output pose ParametersSet: %lu", posecopy.conformation().parameters_set(1)->n_parameters());
		TS_TRACE( outbuffer );
		for ( core::Size i=1; i<=4; ++i ) {
			BundleParametersOP h1params( utility::pointer::dynamic_pointer_cast<BundleParameters>( posecopy.conformation().parameters_set(1)->parameters(i) ) );
			TS_ASSERT(h1params);
			core::Real r0_1(h1params->r0());
			core::Real omega0_1(h1params->omega0());
			core::Real delta_omega0_1(h1params->delta_omega0());
			core::Real z1_off( h1params->z1_offset() );
			core::Real z0_off( h1params->z0_offset() );
			bool invert_1(h1params->invert_helix());
			sprintf(outbuffer, "helix %lu: r0=%.4f omega0=%.4f delta_omega0=%.4f invert=%s z1_offset=%.4f z0_offset=%.4f", i, r0_1, omega0_1, delta_omega0_1, (std::string( invert_1 ? "true" : "false" )).c_str(), z1_off, z0_off );
			TS_TRACE( outbuffer );

			if ( i==1 || i==3 ) {
				TS_ASSERT_DELTA( r0_1, 5.0, 1e-5  );
				TS_ASSERT_DELTA( omega0_1, 0.05, 1e-5  );
				TS_ASSERT_DELTA( delta_omega0_1, 0.1 + (i==3 ? 3.141592654 : 0), 1e-5  );
				TS_ASSERT( invert_1 == false );
				TS_ASSERT_DELTA( z1_off, 0.2, 1e-5 );
				TS_ASSERT_DELTA( z0_off, 0, 1e-5 );
			} else {
				TS_ASSERT_DELTA( r0_1, 6.5, 1e-5  );
				TS_ASSERT_DELTA( omega0_1, 0.03, 1e-5  );
				TS_ASSERT_DELTA( delta_omega0_1, 1.5 + (i==4 ? 3.141592654 : 0), 1e-5  );
				TS_ASSERT( invert_1 == true );
				TS_ASSERT_DELTA( z1_off, 0, 1e-5 );
				TS_ASSERT_DELTA( z0_off, -0.2, 1e-5 );
			}

		}

		//Make a copy of the pose:
		TS_TRACE( "Attempting to copy the pose using the assignment operator." );
		core::pose::Pose posecopy2;
		posecopy2 = *testpose_;

		TS_TRACE( "Reading Crick parameters from the assignment-copied pose." );
		sprintf(outbuffer, "Number of Parameters objects in the output pose ParametersSet: %lu", posecopy2.conformation().parameters_set(1)->n_parameters());
		TS_TRACE( outbuffer );
		for ( core::Size i=1; i<=4; ++i ) {
			BundleParametersOP h1params( utility::pointer::dynamic_pointer_cast<BundleParameters>( posecopy2.conformation().parameters_set(1)->parameters(i) ) );
			TS_ASSERT(h1params);
			core::Real r0_1(h1params->r0());
			core::Real omega0_1(h1params->omega0());
			core::Real delta_omega0_1(h1params->delta_omega0());
			bool invert_1(h1params->invert_helix());
			core::Real z1_off( h1params->z1_offset() );
			core::Real z0_off( h1params->z0_offset() );
			sprintf(outbuffer, "helix %lu: r0=%.4f omega0=%.4f delta_omega0=%.4f invert=%s z1_offset=%.4f z0_offset=%.4f", i, r0_1, omega0_1, delta_omega0_1, (std::string( invert_1 ? "true" : "false" )).c_str(), z1_off, z0_off );
			TS_TRACE( outbuffer );

			if ( i==1 || i==3 ) {
				TS_ASSERT_DELTA( r0_1, 5.0, 1e-5  );
				TS_ASSERT_DELTA( omega0_1, 0.05, 1e-5  );
				TS_ASSERT_DELTA( delta_omega0_1, 0.1 + (i==3 ? 3.141592654 : 0), 1e-5  );
				TS_ASSERT( invert_1 == false );
				TS_ASSERT_DELTA( z1_off, 0.2, 1e-5 );
				TS_ASSERT_DELTA( z0_off, 0, 1e-5 );
			} else {
				TS_ASSERT_DELTA( r0_1, 6.5, 1e-5  );
				TS_ASSERT_DELTA( omega0_1, 0.03, 1e-5  );
				TS_ASSERT_DELTA( delta_omega0_1, 1.5 + (i==4 ? 3.141592654 : 0), 1e-5  );
				TS_ASSERT( invert_1 == true );
				TS_ASSERT_DELTA( z1_off, 0, 1e-5 );
				TS_ASSERT_DELTA( z0_off, -0.2, 1e-5 );
			}

		}

		TS_TRACE( "Finished test_params_data_storage unit test." );

		return;
	}

}; //class params_data_storage_Tests
