// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// @brief This file is to test the behavior of restypes which are stored in the Conformation,
// and if they can be properly used in many protocols.

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/PoseResidueTypeSet.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ElementSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.hh>

#include <core/chemical/ResidueType.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/residue_io.hh>
#include <core/chemical/MergeBehaviorManager.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/build_pose_as_is.hh>
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
#include <basic/database/open.hh>
#include <basic/Tracer.hh>

#include <fstream>

#ifdef SERIALIZATION
// Cereal headers
#include <cereal/archives/json.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

static basic::Tracer TR("core.conformation.conformation_stored_restypes.cxxtest");

/// @details Note that there really aren't any asserts here - this test is mainly
/// to make sure that things don't crash unnecessarily.

class conformation_stored_restypes_tests : public CxxTest::TestSuite {


	std::string pdbname_;
	std::string resfile_;
	std::string name3_;
	std::string typeset_;
	std::string offtypeset_;
	core::scoring::ScoreFunctionOP scorefxn_;
	protocols::minimization_packing::PackRotamersMoverOP pack_;
	protocols::minimization_packing::MinMoverOP min_;

public:

	void setUp() {
		core_init();

		pdbname_ = "protocols/denovo_design/components/test_enz_remark_input.pdb";
		resfile_ = "devel/denovo_design/D2I.params";
		name3_ = "D2I";

		scorefxn_ = core::scoring::get_score_function();
		scorefxn_->set_weight( core::scoring::cart_bonded, 1.0 );
		scorefxn_->set_weight( core::scoring::pro_close, 0.0 );

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

	void test_load_restype_conf() {
		using namespace core::chemical;
		core::pose::Pose pose;
		TS_ASSERT( ! pose.conformation().residue_type_set_for_conf( FULL_ATOM_t )->has_name(name3_) );
		core::chemical::PoseResidueTypeSetOP mod_restype_set( pose.conformation().modifiable_residue_type_set_for_conf( core::chemical::FULL_ATOM_t ) );
		mod_restype_set->add_base_residue_type( resfile_ );
		pose.conformation().reset_residue_type_set_for_conf( mod_restype_set );
		TS_ASSERT( pose.conformation().residue_type_set_for_conf( FULL_ATOM_t )->has_name(name3_) );
	}

	void test_load_pdb() {
		using namespace core::chemical;
		core::pose::Pose pose;
		core::chemical::PoseResidueTypeSetOP mod_restype_set( pose.conformation().modifiable_residue_type_set_for_conf( core::chemical::FULL_ATOM_t ) );
		mod_restype_set->add_base_residue_type( resfile_ );
		pose.conformation().reset_residue_type_set_for_conf( mod_restype_set );
		TR << "Loading PDB:" << std::endl;
		core::io::pdb::build_pose_from_pdb_as_is( pose, pdbname_ );
		// Just test it doesn't crash.
	}

	void test_load_pdb_with_pose_copy() {
		using namespace core::chemical;
		core::pose::Pose pose;
		core::chemical::PoseResidueTypeSetOP mod_restype_set( pose.conformation().modifiable_residue_type_set_for_conf( core::chemical::FULL_ATOM_t ) );
		mod_restype_set->add_base_residue_type( resfile_ );
		pose.conformation().reset_residue_type_set_for_conf( mod_restype_set );

		// Make sure that we can copy the contained ResidueTypes when we copy the pose
		core::pose::Pose pose_copy( pose );

		TR << "Loading PDB:" << std::endl;
		core::io::pdb::build_pose_from_pdb_as_is( pose_copy, pdbname_ );
		// Just test it doesn't crash.
	}


	void test_score_restype_noset() {
		using namespace core::chemical;
		core::pose::Pose pose;
		core::chemical::PoseResidueTypeSetOP mod_restype_set( pose.conformation().modifiable_residue_type_set_for_conf( core::chemical::FULL_ATOM_t ) );
		mod_restype_set->add_base_residue_type( resfile_ );
		pose.conformation().reset_residue_type_set_for_conf( mod_restype_set );

		core::chemical::ResidueTypeCOP restype( pose.conformation().residue_type_set_for_conf()->name_mapOP(name3_) );
		TS_ASSERT( restype != nullptr );

		// Make a default residue
		core::conformation::Residue new_rsd( *restype, true );
		pose.append_residue_by_jump(new_rsd, 1);

		TR << "Scoring:" << std::endl;
		(*scorefxn_)( pose );
	}

	void test_pack_restype_noset() {
		using namespace core::chemical;
		core::pose::Pose pose;
		core::chemical::PoseResidueTypeSetOP mod_restype_set( pose.conformation().modifiable_residue_type_set_for_conf( core::chemical::FULL_ATOM_t ) );
		mod_restype_set->add_base_residue_type( resfile_ );
		pose.conformation().reset_residue_type_set_for_conf( mod_restype_set );

		core::chemical::ResidueTypeCOP restype( pose.conformation().residue_type_set_for_conf()->name_mapOP(name3_) );
		TS_ASSERT( restype != nullptr );

		// Make a default residue
		core::conformation::Residue new_rsd( *restype, true );
		pose.append_residue_by_jump(new_rsd, 1);

		TR << "Packing:" << std::endl;
		pack_->apply( pose );
	}

	void test_min_restype_noset() {
		using namespace core::chemical;
		core::pose::Pose pose;
		core::chemical::PoseResidueTypeSetOP mod_restype_set( pose.conformation().modifiable_residue_type_set_for_conf( core::chemical::FULL_ATOM_t ) );
		mod_restype_set->add_base_residue_type( resfile_ );
		pose.conformation().reset_residue_type_set_for_conf( mod_restype_set );

		core::chemical::ResidueTypeCOP restype( pose.conformation().residue_type_set_for_conf()->name_mapOP(name3_) );
		TS_ASSERT( restype != nullptr );

		// Make a default residue
		core::conformation::Residue new_rsd( *restype, true );
		pose.append_residue_by_jump(new_rsd, 1);

		TR << "Minimizing" << std::endl;
		min_->apply( pose );
	}

	void test_cartmin_restype_noset() {
		using namespace core::chemical;
		core::pose::Pose pose;
		core::chemical::PoseResidueTypeSetOP mod_restype_set( pose.conformation().modifiable_residue_type_set_for_conf( core::chemical::FULL_ATOM_t ) );
		mod_restype_set->add_base_residue_type( resfile_ );
		pose.conformation().reset_residue_type_set_for_conf( mod_restype_set );

		core::chemical::ResidueTypeCOP restype( pose.conformation().residue_type_set_for_conf()->name_mapOP(name3_) );
		TS_ASSERT( restype != nullptr );

		// Make a default residue
		core::conformation::Residue new_rsd( *restype, true );
		pose.append_residue_by_jump(new_rsd, 1);

		TR << "Cartesian Minimizing" << std::endl;
		protocols::minimization_packing::MinMoverOP cartmin( new protocols::minimization_packing::MinMover( *min_ ) );
		cartmin->cartesian( true );
		cartmin->min_type( "lbfgs_armijo_nonmonotone" );

		cartmin->apply( pose );
	}

	void test_serialization() {
#ifdef SERIALIZATION
		using namespace core::chemical;
		core::pose::Pose pose;
		core::chemical::PoseResidueTypeSetOP mod_restype_set( pose.conformation().modifiable_residue_type_set_for_conf( core::chemical::FULL_ATOM_t ) );
		mod_restype_set->add_base_residue_type( resfile_ );
		mod_restype_set->add_unpatchable_residue_type( "protocols/ligand_docking/7cpa.params" );

		// Add some additional data, to make sure that things get serialized properly.
		utility::vector1< std::string > trunc_patches;
		utility::vector1< std::string > trunc_metapatches;
		trunc_patches.push_back("core/chemical/1pqc_test.patch");
		trunc_metapatches.push_back( basic::database::full_name( "chemical/residue_type_sets/fa_standard/metapatches/connect.txt" ) );
		mod_restype_set->add_patches( trunc_patches, trunc_metapatches );

		MergeBehaviorManagerOP mbm_copy( new MergeBehaviorManager( mod_restype_set->default_rts()->merge_behavior_manager() ) );
		mod_restype_set->set_merge_behavior_manager( mbm_copy );

		pose.conformation().reset_residue_type_set_for_conf( mod_restype_set );

		ResidueTypeSetCOP original( pose.conformation().residue_type_set_for_conf( core::chemical::FULL_ATOM_t ) );
		std::ostringstream oss;
		{
			cereal::JSONOutputArchive arch( oss );
			arch( original );
		}

		TR.Debug << "Serialized Form: \n";
		TR.Debug << oss.str() << std::endl;

		ResidueTypeSetCOP reconstituted;
		std::istringstream iss( oss.str() );
		{
			cereal::JSONInputArchive arch( iss );
			arch( reconstituted );
		}
		TR << std::endl;

		TS_ASSERT( original != reconstituted ); // Shouldn't just be the raw pointer - needs to be a new item.

		ResidueTypeCOPs basetypes( reconstituted->base_residue_types() );
		TS_ASSERT( basetypes.size() >= 1 );
		ResidueTypeCOP reconst_restype( reconstituted->name_mapOP(basetypes[1]->name()) ); // Grabs from the cache

		TS_ASSERT( basetypes[1].get() == reconst_restype.get() ); // It's the same type, so it shouldn't be re-created more than once.

		// Double check that the Atom TypeSet reset works
		TS_ASSERT( reconst_restype->has("O1") );
		TS_ASSERT( reconst_restype->atom("O1").element_type() != nullptr );
		TS_ASSERT_EQUALS( reconst_restype->atom("O1").element_type()->get_chemical_symbol(), "O" );


		TR << "Done with serialization test." << std::endl;

#endif // SERIALIZATION
	}

};

