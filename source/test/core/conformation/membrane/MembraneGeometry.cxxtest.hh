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
///    Last Modified: 4/8/2020
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


		core_init_with_additional_options( "-mp:restore_lazaridis_imm_behavior" );

		//grab test atoms
		//in membrane
		atom_1_ = pose_->conformation().residue( 206 ).atom( "CA" ).xyz();

		//membrane interface
		atom_2_ = pose_->conformation().residue( 197 ).atom( "CA" ).xyz();

		//extracellular
		atom_3_ = pose_->conformation().residue( 69 ).atom( "CA" ).xyz();

		//intracellular
		atom_4_ = pose_->conformation().residue( 160 ).atom( "CA" ).xyz();

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
		AddMembraneMoverOP add_memb( new AddMembraneMover( spanfile_ ));

		add_memb->apply( *pose_ );

		MembraneInfoOP membrane_info = pose_->conformation().membrane_info();
		core::conformation::membrane::MembraneGeometryCOP mp_geometry( pose_->conformation().membrane_info()->membrane_geometry() );
		core::Real atom_1_transition = mp_geometry->f_transition( pose_->conformation(), 206, 2 );
		core::Real atom_2_transition = mp_geometry->f_transition( pose_->conformation(), 197, 2 );
		core::Real atom_3_transition = mp_geometry->f_transition( pose_->conformation(), 69, 2 );
		core::Real atom_4_transition = mp_geometry->f_transition( pose_->conformation(), 160, 2 );

		TS_ASSERT_DELTA( atom_1_transition, 0.00, 0.001 );
		TS_ASSERT_DELTA( atom_2_transition, 0.0229, 0.001 );
		TS_ASSERT_DELTA( atom_3_transition, 0.9884, 0.001 );
		TS_ASSERT_DELTA( atom_4_transition,0.9987, 0.001 );

	}

	void test_1py6_bicelle_transition_values() {
		core_init_with_additional_options( "-mp:restore_lazaridis_imm_behavior -mp:geometry bicelle -mp:geo:bicelle_radius 30");
		AddMembraneMoverOP add_memb( new AddMembraneMover( spanfile_ ));

		add_memb->apply( *pose_ );

		MembraneInfoOP membrane_info = pose_->conformation().membrane_info();
		core::conformation::membrane::MembraneGeometryCOP mp_geometry( pose_->conformation().membrane_info()->membrane_geometry() );
		core::Real atom_1_transition = mp_geometry->f_transition( pose_->conformation(), 206, 2 );
		core::Real atom_2_transition = mp_geometry->f_transition( pose_->conformation(), 197, 2 );
		core::Real atom_3_transition = mp_geometry->f_transition( pose_->conformation(), 69, 2 );
		core::Real atom_4_transition = mp_geometry->f_transition( pose_->conformation(), 160, 2 );

		TS_ASSERT_DELTA( atom_1_transition, 0.00, 0.001 );
		TS_ASSERT_DELTA( atom_2_transition, 0.0229, 0.001 );
		TS_ASSERT_DELTA( atom_3_transition, 0.9884, 0.001 );
		TS_ASSERT_DELTA( atom_4_transition, 0.9987, 0.001 );

	}

	void test_fampenvenergy_bicelle_deriv() {
		core_init_with_additional_options( "-mp:restore_lazaridis_imm_behavior -mp:geometry bicelle -mp:geo:bicelle_radius 70" );

		using namespace core::scoring;
		using namespace core::import_pose;
		using namespace core::pose;

		// Load in pose from pdb
		core::pose::PoseOP pose = utility::pointer::make_shared< Pose >();
		TR << "initialize pose" << std::endl;
		core::import_pose::pose_from_file( *pose, "core/conformation/membrane/1AFO_AB.pdb", core::import_pose::PDB_file);
		TR << "pose from file" << std::endl;

		// Initialize Spans from spanfile
		std::string span = "core/conformation/membrane/1AFO_AB.span";

		AddMembraneMoverOP add_memb( new AddMembraneMover( span ));

		add_memb->apply( *pose );
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( FaMPEnv, 0.5 );
		core::kinematics::MoveMap movemap;
		movemap.set_bb( true );
		AtomDerivValidator adv( *pose, sfxn, movemap );
		adv.simple_deriv_check( true, 1e-6);
	}

	void test_fampenvenergy_bicelle_deriv_barrel() {
		core_init_with_additional_options( "-mp:restore_lazaridis_imm_behavior -mp:geometry bicelle -mp:geo:bicelle_radius 40" );

		using namespace core::scoring;
		using namespace core::import_pose;
		using namespace core::pose;

		// Load in pose from pdb
		core::pose::PoseOP pose = utility::pointer::make_shared< Pose >();
		TR << "initialize pose" << std::endl;
		core::import_pose::pose_from_file( *pose, "core/conformation/membrane/3gp6_opm_A.pdb", core::import_pose::PDB_file);
		TR << "pose from file" << std::endl;

		// Initialize Spans from spanfile
		std::string span = "core/conformation/membrane/3gp6_opm_A.span";

		AddMembraneMoverOP add_memb( new AddMembraneMover( span ));

		add_memb->apply( *pose );
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( FaMPEnv, 0.5 );
		core::kinematics::MoveMap movemap;
		movemap.set_bb( true );
		AtomDerivValidator adv( *pose, sfxn, movemap );
		adv.simple_deriv_check( true, 1e-6);
	}

	void test_fampenvenergy_vesicle_deriv() {
		core_init_with_additional_options( "-mp:restore_lazaridis_imm_behavior -mp:geometry vesicle -mp:geo:vesicle_radius 100" );

		using namespace core::scoring;
		using namespace core::import_pose;
		using namespace core::pose;

		// Load in pose from pdb
		core::pose::PoseOP pose = utility::pointer::make_shared< Pose >();
		TR << "initialize pose" << std::endl;
		core::import_pose::pose_from_file( *pose, "core/conformation/membrane/1AFO_AB.pdb", core::import_pose::PDB_file);
		TR << "pose from file" << std::endl;

		// Initialize Spans from spanfile
		std::string span = "core/conformation/membrane/1AFO_AB.span";

		AddMembraneMoverOP add_memb( new AddMembraneMover( span ));

		add_memb->apply( *pose );
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( FaMPEnv, 0.3 );
		core::kinematics::MoveMap movemap;
		movemap.set_bb( true );
		AtomDerivValidator adv( *pose, sfxn, movemap );
		adv.simple_deriv_check( true, 1e-6);

	}

	void test_fampenvenergy_doublevesicle_deriv() {
		core_init_with_additional_options( "-mp:restore_lazaridis_imm_behavior -mp:geometry double_vesicle -mp:geo:vesicle_radius 100" );

		using namespace core::scoring;
		using namespace core::import_pose;
		using namespace core::pose;

		// Load in pose from pdb
		core::pose::PoseOP pose = utility::pointer::make_shared< Pose >();
		TR << "initialize pose" << std::endl;
		core::import_pose::pose_from_file( *pose, "core/conformation/membrane/1AFO_AB.pdb", core::import_pose::PDB_file);
		TR << "pose from file" << std::endl;

		// Initialize Spans from spanfile
		std::string span = "core/conformation/membrane/1AFO_AB.span";

		AddMembraneMoverOP add_memb( new AddMembraneMover( span ));

		add_memb->apply( *pose );
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( FaMPEnv, 0.3 );
		core::kinematics::MoveMap movemap;
		movemap.set_bb( true );
		AtomDerivValidator adv( *pose, sfxn, movemap );
		adv.simple_deriv_check( true, 1e-3);

	}

	void test_fampsolvenergy_bicelle_deriv() {
		core_init_with_additional_options( "-mp:restore_lazaridis_imm_behavior -mp:geometry bicelle -mp:geo:bicelle_radius 70" );

		using namespace core::scoring;
		using namespace core::import_pose;
		using namespace core::pose;

		// Load in pose from pdb
		core::pose::PoseOP pose = utility::pointer::make_shared< Pose >();
		TR << "initialize pose" << std::endl;
		core::import_pose::pose_from_file( *pose, "core/conformation/membrane/1AFO_AB.pdb", core::import_pose::PDB_file);
		TR << "pose from file" << std::endl;

		// Initialize Spans from spanfile
		std::string span = "core/conformation/membrane/1AFO_AB.span";

		AddMembraneMoverOP add_memb( new AddMembraneMover( span ));

		add_memb->apply( *pose );
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( FaMPSolv, 0.3 );
		core::kinematics::MoveMap movemap;
		movemap.set_bb( true );
		AtomDerivValidator adv( *pose, sfxn, movemap );
		adv.simple_deriv_check( true, 1e-6);

	}
private:

	core::pose::PoseOP pose_;

	std::string spanfile_;

	//Test atoms
	numeric::xyzVector< core::Real > atom_1_;
	numeric::xyzVector< core::Real > atom_2_;
	numeric::xyzVector< core::Real > atom_3_;
	numeric::xyzVector< core::Real > atom_4_;


}; // MembraneGeometry unit test
