// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// @brief This file is to test the behavior of "orphan" restypes (that is, ResidueTypes which
// don't have a ResidueTypeSet listed, or residues which aren't in the ResidueTypeSet which the
// do list as a parent.

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ElementSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.hh>

#include <core/chemical/ResidueType.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/residue_io.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/MinimizerOptions.hh> // Needed for MinMover copy constructor

#include <core/types.hh>

#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/minimization_packing/MinMover.hh>

#include <devel/init.hh>

#include <utility/options/FileVectorOption.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#
#include <basic/Tracer.hh>

#include <fstream>
static basic::Tracer TR("protocols.ligand_design.orphan_restypes.cxxtest");


/// @details Note that there really aren't any asserts here - this test is mainly
/// to make sure that things don't crash unnecessarily.

class orphan_restypes_tests : public CxxTest::TestSuite {

	std::string resfile_;
	std::string typeset_;
	std::string offtypeset_;
	core::scoring::ScoreFunctionOP scorefxn_;
	core::chemical::ChemicalManager * chem_mang_;
	core::chemical::AtomTypeSetCAP atom_types_;
	core::chemical::ElementSetCAP elements_;
	core::chemical::MMAtomTypeSetCAP mm_atom_types_;
	core::chemical::orbitals::OrbitalTypeSetCAP orbital_types_;
	protocols::minimization_packing::PackRotamersMoverOP pack_;
	protocols::minimization_packing::MinMoverOP min_;

public:

	void setUp() {
		core_init();

		resfile_ = "protocols/ligand_docking/7cpa.params";
		typeset_ = core::chemical::FA_STANDARD;
		offtypeset_ = core::chemical::CENTROID;

		scorefxn_ = core::scoring::get_score_function();
		scorefxn_->set_weight( core::scoring::cart_bonded, 1.0 );
		scorefxn_->set_weight( core::scoring::pro_close, 0.0 );

		chem_mang_ = core::chemical::ChemicalManager::get_instance();
		atom_types_ = chem_mang_->atom_type_set(typeset_);
		elements_ = chem_mang_->element_set("default");
		mm_atom_types_ = chem_mang_->mm_atom_type_set(typeset_);
		orbital_types_ = chem_mang_->orbital_type_set(typeset_);

		core::pack::task::TaskFactoryOP taskfactory( new core::pack::task::TaskFactory );
		core::pack::task::operation::RestrictToRepackingOP rtrp( new core::pack::task::operation::RestrictToRepacking );
		taskfactory->push_back( rtrp );
		pack_ = protocols::minimization_packing::PackRotamersMoverOP( new protocols::minimization_packing::PackRotamersMover );
		pack_->score_function( scorefxn_ );
		pack_->task_factory( taskfactory );

		core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap );
		min_ = protocols::minimization_packing::MinMoverOP(new protocols::minimization_packing::MinMover);
		min_->score_function( scorefxn_ );
		min_->movemap( movemap );
		min_->tolerance(0.1); // Don't bother to do a good min - just see that we can.
	}

	void tearDown() {}

	void test_load_restype_noset() {
		core::chemical::ResidueTypeOP restype( read_topology_file(
			resfile_, atom_types_, elements_, mm_atom_types_, orbital_types_ ) ); // Null for ResidueTypeSet is deliberate
	}

	void test_score_restype_noset() {
		core::chemical::ResidueTypeOP restype( read_topology_file(
			resfile_, atom_types_, elements_, mm_atom_types_, orbital_types_ ) ); // Null for ResidueTypeSet is deliberate

		core::pose::Pose pose;
		// Make a default residue
		core::conformation::Residue new_rsd( *restype, true );
		pose.append_residue_by_jump(new_rsd, 1);

		TR << "Scoring:" << std::endl;
		(*scorefxn_)( pose );
	}

	void test_pack_restype_noset() {
		core::chemical::ResidueTypeOP restype( read_topology_file(
			resfile_, atom_types_, elements_, mm_atom_types_, orbital_types_ ) ); // Null for ResidueTypeSet is deliberate

		core::pose::Pose pose;
		// Make a default residue
		core::conformation::Residue new_rsd( *restype, true );
		pose.append_residue_by_jump(new_rsd, 1);

		TR << "Packing:" << std::endl;
		pack_->apply( pose );
	}

	void test_min_restype_noset() {
		core::chemical::ResidueTypeOP restype( read_topology_file(
			resfile_, atom_types_, elements_, mm_atom_types_, orbital_types_ ) ); // Null for ResidueTypeSet is deliberate

		core::pose::Pose pose;
		// Make a default residue
		core::conformation::Residue new_rsd( *restype, true );
		pose.append_residue_by_jump(new_rsd, 1);

		TR << "Minimizing" << std::endl;
		min_->apply( pose );
	}

	void test_cartmin_restype_noset() {
		core::chemical::ResidueTypeOP restype( read_topology_file(
			resfile_, atom_types_, elements_, mm_atom_types_, orbital_types_ ) ); // Null for ResidueTypeSet is deliberate

		core::pose::Pose pose;
		// Make a default residue
		core::conformation::Residue new_rsd( *restype, true );
		pose.append_residue_by_jump(new_rsd, 1);

		TR << "Cartesian Minimizing" << std::endl;
		protocols::minimization_packing::MinMoverOP cartmin( new protocols::minimization_packing::MinMover( *min_ ) );
		cartmin->cartesian( true );
		cartmin->min_type( "lbfgs_armijo_nonmonotone" );

		cartmin->apply( pose );
	}

};

