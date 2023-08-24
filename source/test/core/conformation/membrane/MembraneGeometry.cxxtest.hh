// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/conformation/membrane/MembraneGeometry.cxxtest.hh
///
/// @brief   Unit Test: MembraneGeometry Object
/// @details The membrane geometry is the base class for
///      the different membrane geometries.
///
/// @author  Hope Woods (hope.woods@vanderbilt.edu)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/deriv_funcs.hh>
#include <test/util/pose_funcs.hh>


//Unit Headers
#include <core/energy_methods/FaMPEnvEnergy.hh>
#include <core/energy_methods/EnvSmoothEnergy.hh>

// Package Headers
#include <protocols/membrane/AddMembraneMover.hh>

#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/MembraneGeometry.hh>
#include <core/conformation/membrane/membrane_geometry/Bicelle.hh>
#include <core/conformation/membrane/ImplicitLipidInfo.hh>

#include <core/conformation/Conformation.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <core/types.hh>
#include <numeric/xyzVector.hh>

#include <basic/Tracer.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>

// Utility Headers
#include <utility/vector1.hh>

using namespace core;
using namespace core::conformation;
using namespace core::conformation::membrane;
using namespace protocols::membrane;

static basic::Tracer TR("core.conformation.membrane.MembraneGeometry.cxxtest");

class MembraneGeometryTest : public CxxTest::TestSuite {

public: // test functions

	// Test Setup Functions ///////////////////////////

	/// @brief Setup Test
	void setUp() {
		using namespace core::import_pose;
		using namespace core::pose;


		// Initialize core & options system
		core_init();
		TR << "core_init" << std::endl;

		// Load in pose from pdb
		pose_ = utility::pointer::make_shared< Pose >();
		TR << "initialize pose" << std::endl;
		core::import_pose::pose_from_file( *pose_, "core/conformation/membrane/1PY6_mp_coords.pdb", core::import_pose::PDB_file);
		TR << "pose from file" << std::endl;

		// Initialize Spans from spanfile
		spanfile_ = "core/conformation/membrane/1PY6.span";


		//grab test atoms
		//in membrane
		atom_1_ = pose_->conformation().residue( 206 ).atom( "CA" ).xyz();

		//membrane interface
		atom_2_ = pose_->conformation().residue( 197 ).atom( "CA" ).xyz();

		//extracellular
		atom_3_ = pose_->conformation().residue( 69 ).atom( "CA" ).xyz();

		//intracellular
		atom_4_ = pose_->conformation().residue( 160 ).atom( "CA" ).xyz();

		// Load alpha helical protein test case
		ahelical_pose_ = utility::pointer::make_shared< Pose >();
		pose_from_file( *ahelical_pose_, "core/conformation/membrane/1U19_tr_ignorechain.pdb" , core::import_pose::PDB_file );
		AddMembraneMoverOP add_memb( new AddMembraneMover( "from_structure" ) );
		add_memb->apply( *ahelical_pose_ );

		// Setup some test coordinates
		// outside - water exposed (above the membrane)
		p1_ = ahelical_pose_->conformation().residue( 23 ).atom( "CA" ).xyz();

		// in pore, water exposed
		p2_ = ahelical_pose_->conformation().residue( 125 ).atom( "CA" ).xyz();

		// interface exposed
		p3_ = ahelical_pose_->conformation().residue( 201 ).atom( "CA" ).xyz();

		// lipid exposed
		p4_ = ahelical_pose_->conformation().residue( 209 ).atom( "CA" ).xyz();

		// outside - water exposed (below the membrane)
		p5_ = ahelical_pose_->conformation().residue( 245 ).atom( "CA" ).xyz();


	}

	/// @brief Tear Down Test
	void tearDown() {}

	// Test Methods /////////////////////////////////
	//testing that correct geometry is initialized for different command line options
	void test_default_slab_initialization() {

		AddMembraneMoverOP add_memb( new AddMembraneMover( spanfile_ ));

		add_memb->apply( *pose_ );

		MembraneInfoOP membrane_info = pose_->conformation().membrane_info();
		TR << "geometry string: " << membrane_info->membrane_geometry()->geometry_string() << std::endl;
		TS_ASSERT_EQUALS( membrane_info->membrane_geometry()->geometry_string(), "slab" );
	}

	void test_bicelle_initialization() {
		core_init_with_additional_options( "-mp:restore_lazaridis_imm_behavior -mp:geometry bicelle" );

		AddMembraneMoverOP add_memb( new AddMembraneMover( spanfile_ ));

		add_memb->apply( *pose_ );

		MembraneInfoOP membrane_info = pose_->conformation().membrane_info();
		TR << "geometry string: " << membrane_info->membrane_geometry()->geometry_string() << std::endl;
		TS_ASSERT_EQUALS( membrane_info->membrane_geometry()->geometry_string(), "bicelle" );
	}

	void test_vesicle_initialization() {
		core_init_with_additional_options( "-mp:restore_lazaridis_imm_behavior -mp:geometry vesicle" );

		AddMembraneMoverOP add_memb( new AddMembraneMover( spanfile_ ));

		add_memb->apply( *pose_ );

		MembraneInfoOP membrane_info = pose_->conformation().membrane_info();
		TR << "geometry string: " << membrane_info->membrane_geometry()->geometry_string() << std::endl;
		TS_ASSERT_EQUALS( membrane_info->membrane_geometry()->geometry_string(), "vesicle" );
	}

	void test_doublevesicle_initialization() {
		core_init_with_additional_options( "-mp:restore_lazaridis_imm_behavior -mp:geometry double_vesicle" );

		AddMembraneMoverOP add_memb( new AddMembraneMover( spanfile_ ));

		add_memb->apply( *pose_ );

		MembraneInfoOP membrane_info = pose_->conformation().membrane_info();
		TR << "geometry string: " << membrane_info->membrane_geometry()->geometry_string() << std::endl;
		TS_ASSERT_EQUALS( membrane_info->membrane_geometry()->geometry_string(), "double_vesicle" );
	}

	//testing bicelle radius
	void test_bicelle_radius_cmd_line() {
		core_init_with_additional_options( "-mp:restore_lazaridis_imm_behavior -mp:geometry bicelle -mp:geo:bicelle_radius 40" );

		//initialize Bicelle with steepness, thickness, and protein_diameter
		core::conformation::membrane::membrane_geometry::BicelleOP bicelle( new membrane_geometry::Bicelle( 10, 15, 30 ));
		bicelle->update_radii();

		core::Real radius = 40;
		TR << radius << std::endl;
		TR << bicelle->bicelle_inner_radius() << std::endl;
		TS_ASSERT_EQUALS( bicelle->bicelle_inner_radius(), radius );
	}

	void test_bicelle_radius_from_protein_diameter() {
		core_init_with_additional_options( "-mp:restore_lazaridis_imm_behavior -mp:geometry bicelle" );

		//initialize Bicelle with steepness, thickness, and protein_diameter
		core::conformation::membrane::membrane_geometry::BicelleOP bicelle( new membrane_geometry::Bicelle( 10, 15, 30 ));
		bicelle->set_protein_slice_diameter( 20 );

		TS_ASSERT_EQUALS( bicelle->bicelle_inner_radius(), 30 );
		TS_ASSERT_EQUALS( bicelle->bicelle_outer_radius(), 45 );
	}

	void test_1py6_slab_transition_values() {
		core_init_with_additional_options( "-mp:restore_lazaridis_imm_behavior -has_pore 0" );
		AddMembraneMoverOP add_memb( new AddMembraneMover( spanfile_ ));

		add_memb->apply( *pose_ );

		MembraneInfoOP membrane_info = pose_->conformation().membrane_info();
		core::conformation::membrane::MembraneGeometryCOP mp_geometry( pose_->conformation().membrane_info()->membrane_geometry() );
		core::Real atom_1_transition = mp_geometry->f_transition( pose_->conformation(), 206, 2 );
		core::Real atom_2_transition = mp_geometry->f_transition( pose_->conformation(), 197, 2 );
		core::Real atom_3_transition = mp_geometry->f_transition( pose_->conformation(), 69, 2 );
		core::Real atom_4_transition = mp_geometry->f_transition( pose_->conformation(), 160, 2 );

		TS_ASSERT_DELTA( atom_1_transition, 0.000, 0.001 );
		TS_ASSERT_DELTA( atom_2_transition, 0.022, 0.001 );
		TS_ASSERT_DELTA( atom_3_transition, 0.988, 0.001 );
		TS_ASSERT_DELTA( atom_4_transition, 0.998, 0.001 );

	}

	void test_1py6_bicelle_transition_values() {
		core_init_with_additional_options( "-mp:restore_lazaridis_imm_behavior -mp:geometry bicelle -mp:geo:bicelle_radius 30 -has_pore 0");
		AddMembraneMoverOP add_memb( new AddMembraneMover( spanfile_ ));

		add_memb->apply( *pose_ );

		MembraneInfoOP membrane_info = pose_->conformation().membrane_info();
		core::conformation::membrane::MembraneGeometryCOP mp_geometry( pose_->conformation().membrane_info()->membrane_geometry() );
		core::Real atom_1_transition = mp_geometry->f_transition( pose_->conformation(), 206, 2 );
		core::Real atom_2_transition = mp_geometry->f_transition( pose_->conformation(), 197, 2 );
		core::Real atom_3_transition = mp_geometry->f_transition( pose_->conformation(), 69, 2 );
		core::Real atom_4_transition = mp_geometry->f_transition( pose_->conformation(), 160, 2 );

		TS_ASSERT_DELTA( atom_1_transition, 0.000, 0.001 );
		TS_ASSERT_DELTA( atom_2_transition, 0.022, 0.001 );
		TS_ASSERT_DELTA( atom_3_transition, 0.988, 0.001 );
		TS_ASSERT_DELTA( atom_4_transition, 0.998, 0.001 );

	}

	void test_1py6_vesicle_transition_values() {
		core_init_with_additional_options( "-mp:restore_lazaridis_imm_behavior -mp:geometry vesicle -mp:geo:vesicle_radius 100 -has_pore 0");
		AddMembraneMoverOP add_memb( new AddMembraneMover( spanfile_ ));

		add_memb->apply( *pose_ );

		MembraneInfoOP membrane_info = pose_->conformation().membrane_info();
		core::conformation::membrane::MembraneGeometryCOP mp_geometry( pose_->conformation().membrane_info()->membrane_geometry() );
		core::Real atom_1_transition = mp_geometry->f_transition( pose_->conformation(), 206, 2 );
		core::Real atom_2_transition = mp_geometry->f_transition( pose_->conformation(), 197, 2 );
		core::Real atom_3_transition = mp_geometry->f_transition( pose_->conformation(), 69, 2 );
		core::Real atom_4_transition = mp_geometry->f_transition( pose_->conformation(), 160, 2 );

		TS_ASSERT_DELTA( atom_1_transition, 0.000, 0.001 );
		TS_ASSERT_DELTA( atom_2_transition, 0.046, 0.001 );
		TS_ASSERT_DELTA( atom_3_transition, 0.991, 0.001 );
		TS_ASSERT_DELTA( atom_4_transition, 0.998, 0.001 );

	}

	void test_f_franklin_transition_value() {

		core::Real tau( ahelical_pose_->conformation().membrane_info()->implicit_lipids()-> water_pseudo_thickness());
		core::Real kappa( ahelical_pose_->conformation().membrane_info()->implicit_lipids()-> water_steepness());
		core::conformation::membrane::MembraneGeometryCOP mp_geometry( ahelical_pose_->conformation().membrane_info()->membrane_geometry() );

		TS_ASSERT_DELTA( mp_geometry->f_franklin( p1_.z(), tau, kappa ), 0.963, 0.001 );
		TS_ASSERT_DELTA( mp_geometry->f_franklin( p2_.z(), tau, kappa ), 0.014, 0.001 );
		TS_ASSERT_DELTA( mp_geometry->f_franklin( p3_.z(), tau, kappa ), 0.572, 0.001 );
		TS_ASSERT_DELTA( mp_geometry->f_franklin( p4_.z(), tau, kappa ), 0.018, 0.001 );
		TS_ASSERT_DELTA( mp_geometry->f_franklin( p5_.z(), tau, kappa ), 0.986, 0.001 );

	}

	void test_g_radius() {
		core::conformation::membrane::MembraneGeometryCOP mp_geometry( ahelical_pose_->conformation().membrane_info()->membrane_geometry() );

		TS_ASSERT_DELTA( mp_geometry->g_radius( p1_ ), 0.548, 0.001 );
		TS_ASSERT_DELTA( mp_geometry->g_radius( p2_ ), 1.178, 0.001 );
		TS_ASSERT_DELTA( mp_geometry->g_radius( p3_ ), 8.207, 0.001 );
		TS_ASSERT_DELTA( mp_geometry->g_radius( p4_ ), 8.414, 0.001 );
		TS_ASSERT_DELTA( mp_geometry->g_radius( p5_ ), 1.786, 0.001 );
		//g_radius values changed after adding ellipse_radius_buffer to
		//major and minor radius values in AqueousPoreFinder
		//Needed to add ellipse_radius_buffer to radius values for derivative
		//to be calculated correctly, since ellipse_radius_buffer
		//was added to the min and max radius values.

	}

	void test_f_cavity() {
		core::conformation::membrane::MembraneGeometryCOP mp_geometry( ahelical_pose_->conformation().membrane_info()->membrane_geometry() );

		TS_ASSERT_DELTA( mp_geometry->f_cavity( p1_ ), 0.997, 0.001 );
		TS_ASSERT_DELTA( mp_geometry->f_cavity( p2_ ), 0.162, 0.001 );
		TS_ASSERT_DELTA( mp_geometry->f_cavity( p3_ ), 0.000, 0.001 );
		TS_ASSERT_DELTA( mp_geometry->f_cavity( p4_ ), 0.000, 0.001 );
		TS_ASSERT_DELTA( mp_geometry->f_cavity( p5_ ), 0.003, 0.001 );
	}


	void test_f_hydration() {
		core::Real tau( ahelical_pose_->conformation().membrane_info()->implicit_lipids()-> water_pseudo_thickness());
		core::Real kappa( ahelical_pose_->conformation().membrane_info()->implicit_lipids()-> water_steepness());
		core::conformation::membrane::MembraneGeometryCOP mp_geometry( ahelical_pose_->conformation().membrane_info()->membrane_geometry() );

		core::Real p1_thk(mp_geometry->f_franklin( p1_.z(), tau, kappa ));
		core::Real p2_thk(mp_geometry->f_franklin( p2_.z(), tau, kappa ));
		core::Real p3_thk(mp_geometry->f_franklin( p3_.z(), tau, kappa ));
		core::Real p4_thk(mp_geometry->f_franklin( p4_.z(), tau, kappa ));
		core::Real p5_thk(mp_geometry->f_franklin( p5_.z(), tau, kappa ));

		TS_ASSERT_DELTA( mp_geometry->f_hydration( p1_thk, p1_ ), 0.999, 0.001 );
		TS_ASSERT_DELTA( mp_geometry->f_hydration( p2_thk, p2_ ), 0.174, 0.001 );
		TS_ASSERT_DELTA( mp_geometry->f_hydration( p3_thk, p3_ ), 0.572, 0.001 );
		TS_ASSERT_DELTA( mp_geometry->f_hydration( p4_thk, p4_ ), 0.018, 0.001 );
		TS_ASSERT_DELTA( mp_geometry->f_hydration( p5_thk, p5_ ), 0.986, 0.001 );

	}

	void test_f_franklin_transition_value_gradient() {

		core::Real tau( ahelical_pose_->conformation().membrane_info()->implicit_lipids()-> water_pseudo_thickness());
		core::Real kappa( ahelical_pose_->conformation().membrane_info()->implicit_lipids()-> water_steepness());
		core::conformation::membrane::MembraneGeometryCOP mp_geometry( ahelical_pose_->conformation().membrane_info()->membrane_geometry() );

		TS_ASSERT_DELTA( mp_geometry->f_franklin_gradient( p1_.z(), tau, kappa ), 0.011, 0.001 );
		TS_ASSERT_DELTA( mp_geometry->f_franklin_gradient( p2_.z(), tau, kappa ), 0.004, 0.001 );
		TS_ASSERT_DELTA( mp_geometry->f_franklin_gradient( p3_.z(), tau, kappa ), 0.083, 0.001 );
		TS_ASSERT_DELTA( mp_geometry->f_franklin_gradient( p4_.z(), tau, kappa ), 0.006, 0.001 );
		TS_ASSERT_DELTA( mp_geometry->f_franklin_gradient( p5_.z(), tau, kappa ), 0.004, 0.001 );

	}

	void test_f_cavity_gradient() {
		core::conformation::membrane::MembraneGeometryCOP mp_geometry( ahelical_pose_->conformation().membrane_info()->membrane_geometry() );

		core::Real g_rad1(mp_geometry->g_radius( p1_ ));
		core::Real g_rad2(mp_geometry->g_radius( p2_ ));
		core::Real g_rad3(mp_geometry->g_radius( p3_ ));
		core::Real g_rad4(mp_geometry->g_radius( p4_ ));
		core::Real g_rad5(mp_geometry->g_radius( p5_ ));

		TS_ASSERT_DELTA( mp_geometry->f_cavity_gradient( g_rad1 ), -0.044, 0.001 );
		TS_ASSERT_DELTA( mp_geometry->f_cavity_gradient( g_rad2 ), -1.155, 0.001 );
		TS_ASSERT_DELTA( mp_geometry->f_cavity_gradient( g_rad3 ), 0.000, 0.001 );
		TS_ASSERT_DELTA( mp_geometry->f_cavity_gradient( g_rad4 ), 0.000, 0.001 );
		TS_ASSERT_DELTA( mp_geometry->f_cavity_gradient( g_rad5 ), -0.016, 0.001 );

	}

private:

	core::pose::PoseOP pose_;

	core::pose::PoseOP ahelical_pose_;

	std::string spanfile_;

	//Test atoms
	numeric::xyzVector< core::Real > atom_1_;
	numeric::xyzVector< core::Real > atom_2_;
	numeric::xyzVector< core::Real > atom_3_;
	numeric::xyzVector< core::Real > atom_4_;

	// Test coordinates
	numeric::xyzVector< core::Real > p1_;
	numeric::xyzVector< core::Real > p2_;
	numeric::xyzVector< core::Real > p3_;
	numeric::xyzVector< core::Real > p4_;
	numeric::xyzVector< core::Real > p5_;
	numeric::xyzVector< core::Real > p6_;



}; // MembraneGeometry unit test
