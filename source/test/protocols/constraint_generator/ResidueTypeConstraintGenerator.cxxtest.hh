// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/constraint_generator/ResidueTypeConstraintGenerator.cxxtest.hh
/// @brief test suite for protocols::constraint_generator::ResidueTypeConstraintGenerator
/// @author Sharon Guffy (guffy@email.unc.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/UTracer.hh>
#include <test/protocols/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/util/pdb1rpb.hh>

// Protocol headers
#include <protocols/constraint_generator/ResidueTypeConstraintGenerator.hh>

// Core headers
#include <core/conformation/Residue.hh>
#include <core/io/pdb/build_pose_as_is.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ResidueTypeConstraint.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/select/residue_selector/OrResidueSelector.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>

// Boost headers
#include <boost/assign.hpp>

// C++ headers

using namespace core::scoring::constraints;
using namespace core::scoring::func;
using namespace protocols::constraint_generator;

static basic::Tracer TR( "protocols.constraint_generator.ResidueTypeConstraintGenerator.cxxtest.hh" );

class ResidueTypeConstraintGeneratorTests : public CxxTest::TestSuite {

private:
	core::pose::PoseOP pose_;
	core::Real const TOLERANCE = 1e-5;
public:
	void setUp()
	{
		protocols_init();
		pose_ = fullatom_poseop_from_string( test_in_pdb_string() );
	}

	void tearDown(){

	}

	void test_constraints_full_pose_native()
	{
		//core::select::residue_selector::ResidueSelectorCOP selector1( new core::select::residue_selector::ResidueIndexSelector( "9" ) );
		//The constraint generator should default to constraining the full pose to its native residue types
		ResidueTypeConstraintGenerator rt_gen;
		rt_gen.set_id( "full_pose_native" );
		//hb_gen.set_residue_selector1( selector );
		TS_ASSERT_EQUALS( rt_gen.class_name(), "ResidueTypeConstraintGenerator" );


		//Use test_in_pdb_string() for pose string
		core::scoring::constraints::ConstraintCOPs const csts = rt_gen.apply( *pose_ );

		//Assert that there is a constraint for every residue
		TS_ASSERT_EQUALS( pose_->total_residue(), csts.size() );
		for ( core::Size i = 1; i <= pose_->total_residue(); ++i ) {
			if ( i > csts.size() ) {
				TR << "Constraints were not added for every residue!" << std::endl;
				continue; //Just so we won't crash (we'll have already failed
			}
			//Assert that each constraint is non-null
			TS_ASSERT( csts.at( i ) );
			//Assert that each constraint is a ResidueTypeConstraint
			TS_ASSERT( utility::pointer::dynamic_pointer_cast< core::scoring::constraints::ResidueTypeConstraint const >( csts.at( i ) ) );
			//Assert that each constraint's rsd_type_name3 is the same as the pose's residue name3
			TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< core::scoring::constraints::ResidueTypeConstraint const >( csts.at( i ) )->get_rsd_type_name3(), pose_->residue_type( i ).name3() );
			//Assert that each constraint's favor_native_bonus is 1
			TS_ASSERT_DELTA( utility::pointer::dynamic_pointer_cast< core::scoring::constraints::ResidueTypeConstraint const >( csts.at( i ) )->get_favor_native_bonus(), 1.0, TOLERANCE );
		}
		ConstraintSet cst_set;
		cst_set.add_constraints( csts );

	}

	void test_constraints_residue_selector()
	{
		core::select::residue_selector::ResidueSelectorCOP select_5( new core::select::residue_selector::ResidueIndexSelector( "5" ) );

		ResidueTypeConstraintGenerator rt_gen;
		rt_gen.set_id( "res_5_native" );
		rt_gen.set_residue_selector( select_5 );
		//hb_gen.set_residue_selector1( selector );


		//Use test_in_pdb_string() for pose string
		core::scoring::constraints::ConstraintCOPs const csts = rt_gen.apply( *pose_ );

		//Assert that there is a constraint only for residue 5
		TS_ASSERT_EQUALS( csts.size(), 1 );
		for ( core::scoring::constraints::ConstraintCOP constraint: csts ) {
			//Assert that the constraint is non-null
			TS_ASSERT( constraint );
			//Assert that the constraint is a ResidueTypeConstraint
			TS_ASSERT( utility::pointer::dynamic_pointer_cast< core::scoring::constraints::ResidueTypeConstraint const >( constraint ) );
			//Assert that the constraint's rsd_type_name3 is the same as the pose's residue name3
			TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< core::scoring::constraints::ResidueTypeConstraint const >( constraint )->get_rsd_type_name3(), pose_->residue_type( 5 ).name3() );
			//Assert that the constraint's favor_native_bonus is 1
			TS_ASSERT_DELTA( utility::pointer::dynamic_pointer_cast< core::scoring::constraints::ResidueTypeConstraint const >( constraint )->get_favor_native_bonus(), 1.0, TOLERANCE );
		}
		ConstraintSet cst_set;
		cst_set.add_constraints( csts );


	}

	void test_constraints_set_type_full_pose(){

		ResidueTypeConstraintGenerator rt_gen;
		rt_gen.set_id( "full_ala" );
		rt_gen.set_rsd_type_name3( "ALA" ); //Constrain all to alanine
		//Use test_in_pdb_string() for pose string
		core::scoring::constraints::ConstraintCOPs const csts = rt_gen.apply( *pose_ );
		for ( core::Size i = 1; i <= pose_->total_residue(); ++i ) {
			if ( i > csts.size() ) {
				TR << "Constraints were not added for every residue!" << std::endl;
				continue; //Just so we won't crash (we'll have already failed
			}
			//Assert that each constraint is non-null
			TS_ASSERT( csts.at( i ) );
			//Assert that each constraint is a ResidueTypeConstraint
			TS_ASSERT( utility::pointer::dynamic_pointer_cast< core::scoring::constraints::ResidueTypeConstraint const >( csts.at( i ) ) );
			//Assert that each constraint's rsd_type_name3 is the same as the pose's residue name3
			TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< core::scoring::constraints::ResidueTypeConstraint const >( csts.at( i ) )->get_rsd_type_name3(), "ALA" );
			//Assert that each constraint's favor_native_bonus is 1
			TS_ASSERT_DELTA( utility::pointer::dynamic_pointer_cast< core::scoring::constraints::ResidueTypeConstraint const >( csts.at( i ) )->get_favor_native_bonus(), 1.0, TOLERANCE );
		}

		ConstraintSet cst_set;
		cst_set.add_constraints( csts );



	}


	void test_constraints_set_type_with_selector()
	{
		core::select::residue_selector::ResidueSelectorCOP select_5( new core::select::residue_selector::ResidueIndexSelector( "5" ) );
		core::select::residue_selector::ResidueSelectorCOP select_3( new core::select::residue_selector::ResidueIndexSelector( "3" ) );
		core::select::residue_selector::ResidueSelectorCOP select_both( new core::select::residue_selector::OrResidueSelector( select_3, select_5 ) );

		ResidueTypeConstraintGenerator rt_gen;
		rt_gen.set_id( "res_3_and_5_ala" );
		rt_gen.set_residue_selector( select_both );
		rt_gen.set_rsd_type_name3( "ALA" );

		//Use test_in_pdb_string() for pose string
		core::scoring::constraints::ConstraintCOPs const csts = rt_gen.apply( *pose_ );

		//Assert that there are constraints for just 2 residues
		TS_ASSERT_EQUALS( csts.size(), 2 );
		for ( core::scoring::constraints::ConstraintCOP constraint: csts ) {
			//Assert that the constraint is non-null
			TS_ASSERT( constraint );
			//Assert that the constraint is a ResidueTypeConstraint
			TS_ASSERT( utility::pointer::dynamic_pointer_cast< core::scoring::constraints::ResidueTypeConstraint const >( constraint ) );
			//Assert that the constraint's rsd_type_name3 is the same as the pose's residue name3
			TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< core::scoring::constraints::ResidueTypeConstraint const >( constraint )->get_rsd_type_name3(), "ALA" );
			//Assert that the constraint's favor_native_bonus is 1
			TS_ASSERT_DELTA( utility::pointer::dynamic_pointer_cast< core::scoring::constraints::ResidueTypeConstraint const >( constraint )->get_favor_native_bonus(), 1.0, TOLERANCE );
		}
		ConstraintSet cst_set;
		cst_set.add_constraints( csts );

	}


	void test_change_favor_native_bonus()
	{
		//core::select::residue_selector::ResidueSelectorCOP selector1( new core::select::residue_selector::ResidueIndexSelector( "9" ) );
		//The constraint generator should default to constraining the full pose to its native residue types
		ResidueTypeConstraintGenerator rt_gen;
		rt_gen.set_id( "full_pose_native" );
		rt_gen.set_favor_native_bonus( 0.5 );
		//hb_gen.set_residue_selector1( selector );
		TS_ASSERT_EQUALS( rt_gen.class_name(), "ResidueTypeConstraintGenerator" );

		//Use test_in_pdb_string() for pose string
		core::scoring::constraints::ConstraintCOPs const csts = rt_gen.apply( *pose_ );

		//Assert that there is a constraint for every residue
		TS_ASSERT_EQUALS( pose_->total_residue(), csts.size() );
		for ( core::Size i = 1; i <= pose_->total_residue(); ++i ) {
			if ( i > csts.size() ) {
				TR << "Constraints were not added for every residue!" << std::endl;
				continue; //Just so we won't crash (we'll have already failed
			}
			//Assert that each constraint is non-null
			TS_ASSERT( csts.at( i ) );
			//Assert that each constraint is a ResidueTypeConstraint
			TS_ASSERT( utility::pointer::dynamic_pointer_cast< core::scoring::constraints::ResidueTypeConstraint const >( csts.at( i ) ) );
			//Assert that each constraint's rsd_type_name3 is the same as the pose's residue name3
			TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< core::scoring::constraints::ResidueTypeConstraint const >( csts.at( i ) )->get_rsd_type_name3(), pose_->residue_type( i ).name3() );
			//Assert that each constraint's favor_native_bonus is 1
			TS_ASSERT_DELTA( utility::pointer::dynamic_pointer_cast< core::scoring::constraints::ResidueTypeConstraint const >( csts.at( i ) )->get_favor_native_bonus(), 0.5, TOLERANCE );
		}
		ConstraintSet cst_set;
		cst_set.add_constraints( csts );

	}



	void test_two_constraints_same_residue()
	{

		//We will favor ALA at positions 3 and 5 and VAL at position 3 only
		core::select::residue_selector::ResidueSelectorCOP select_5( new core::select::residue_selector::ResidueIndexSelector( "5" ) );
		core::select::residue_selector::ResidueSelectorCOP select_3( new core::select::residue_selector::ResidueIndexSelector( "3" ) );
		core::select::residue_selector::ResidueSelectorCOP select_both( new core::select::residue_selector::OrResidueSelector( select_3, select_5 ) );

		ResidueTypeConstraintGenerator rt_gen_1;
		rt_gen_1.set_id( "res_3_and_5_ala" );
		rt_gen_1.set_residue_selector( select_both );
		rt_gen_1.set_favor_native_bonus( 0.5 );
		rt_gen_1.set_rsd_type_name3( "ALA" );
		ResidueTypeConstraintGenerator rt_gen_2;
		rt_gen_2.set_id( "res_3_val" );
		rt_gen_2.set_residue_selector( select_3 );
		rt_gen_2.set_rsd_type_name3( "VAL" );

		core::scoring::constraints::ConstraintCOPs const csts_1 = rt_gen_1.apply( *pose_ );
		core::scoring::constraints::ConstraintCOPs const csts_2 = rt_gen_2.apply( *pose_ );

		TS_ASSERT_EQUALS( csts_1.size(), 2 );
		TS_ASSERT_EQUALS( csts_2.size(), 1 );
		//There should be 3 total constraints
		ConstraintSet cst_set;
		cst_set.add_constraints( csts_1 );
		cst_set.add_constraints( csts_2 );
		core::scoring::constraints::ConstraintCOPs const all_csts = cst_set.get_all_constraints();
		TS_ASSERT_EQUALS( all_csts.size(), 3 );
		core::Size res_3_csts = 0;
		core::Size res_5_csts = 0;
		for ( core::scoring::constraints::ConstraintCOP constraint: all_csts ) {
			//seqpos should be either 3 or 5
			core::Size seqpos = constraint->residues().at( 1 );
			TS_ASSERT( seqpos == 3 || seqpos == 5 );
			if ( seqpos == 5 ) {
				++res_5_csts;
				TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< core::scoring::constraints::ResidueTypeConstraint const >( constraint )->get_rsd_type_name3(), "ALA" );
				TS_ASSERT_DELTA( utility::pointer::dynamic_pointer_cast< core::scoring::constraints::ResidueTypeConstraint const >( constraint )->get_favor_native_bonus(), 0.5, TOLERANCE );
			} else {
				//The constraint on residue 5 should have rsd_type_name3 == "ALA"
				++res_3_csts;
				std::string restype = utility::pointer::dynamic_pointer_cast< core::scoring::constraints::ResidueTypeConstraint const >( constraint )->get_rsd_type_name3();
				//The constraint on residue 3 should have rsd_type_name3 == "ALA" or "VAL"
				TS_ASSERT( restype == "ALA" || restype == "VAL" );
				//The favor native bonus should be 0.5 for the ALA constraints and 1.0 for the VAL constraint
				if ( restype == "ALA" ) {
					TS_ASSERT_DELTA( utility::pointer::dynamic_pointer_cast< core::scoring::constraints::ResidueTypeConstraint const >( constraint )->get_favor_native_bonus(), 0.5, TOLERANCE );
				} else {
					TS_ASSERT_DELTA( utility::pointer::dynamic_pointer_cast< core::scoring::constraints::ResidueTypeConstraint const >( constraint )->get_favor_native_bonus(), 1.0, TOLERANCE );
				}
			}
		}
		//There should be 2 constraints on residue 3 and 1 on residue 5
		TS_ASSERT_EQUALS( res_3_csts, 2 );
		TS_ASSERT_EQUALS( res_5_csts, 1 );
	}

	void test_constraints_full_ref_pose(){

		core::pose::PoseCOP ref_pose = pose_->clone();

		//Mutate a residue in the pose
		ResidueTypeConstraintGenerator rt_gen_1;
		rt_gen_1.set_id( "refpose_mutated_1" );
		rt_gen_1.set_reference_pose( ref_pose );

		//This will effectively mutate residue 1 from ASP to ALA
		core::conformation::Residue replacement( pose_->residue( 2 ) );
		pose_->replace_residue( 5, replacement, true );

		//See if it works
		core::scoring::constraints::ConstraintCOPs const csts_1 = rt_gen_1.apply( *pose_ );
		TS_ASSERT_EQUALS( pose_->total_residue(), csts_1.size() );
		for ( core::Size i = 1; i <= pose_->total_residue(); ++i ) {
			if ( i > csts_1.size() ) {
				TR << "Constraints were not added for every residue!" << std::endl;
				continue; //Just so we won't crash (we'll have already failed
			}
			//Assert that each constraint is non-null
			TS_ASSERT( csts_1.at( i ) );
			//Assert that each constraint is a ResidueTypeConstraint
			TS_ASSERT( utility::pointer::dynamic_pointer_cast< core::scoring::constraints::ResidueTypeConstraint const >( csts_1.at( i ) ) );
			TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< core::scoring::constraints::ResidueTypeConstraint const >( csts_1.at( i ) )->get_rsd_type_name3(), ref_pose->residue_type( i ).name3() );
		}
	}

	void test_constraints_ref_pose_with_selector(){
		core::pose::PoseCOP ref_pose = pose_->clone();
		core::select::residue_selector::ResidueSelectorCOP select_5( new core::select::residue_selector::ResidueIndexSelector( "5" ) );
		core::select::residue_selector::ResidueSelectorCOP select_end( new core::select::residue_selector::ResidueIndexSelector( "110" ) );

		//Replace terminal LEU with non-terminal LEU
		core::conformation::Residue replace_terminus( pose_->residue( 9 ) ); //this is a Leucine too

		pose_->replace_residue( 110, replace_terminus, true );

		//See if it still works
		ResidueTypeConstraintGenerator rt_gen_1;
		rt_gen_1.set_id( "refpose_length_change_pos_control" );
		rt_gen_1.set_reference_pose( ref_pose );
		rt_gen_1.set_residue_selector( select_5 );
		core::scoring::constraints::ConstraintCOPs const csts_1 = rt_gen_1.apply( *pose_ );
		TS_ASSERT_EQUALS( csts_1.size(), 1 );
		std::string old_res = ref_pose->residue_type( 5 ).name3();
		std::string new_res = ref_pose->residue_type( 5 ).name3();
		TR << "Old residue: " << old_res << " New residue: " << new_res << std::endl;
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< core::scoring::constraints::ResidueTypeConstraint const >( csts_1.at( 1 ) )->get_rsd_type_name3(), pose_->residue_type( 5 ).name3() );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< core::scoring::constraints::ResidueTypeConstraint const >( csts_1.at( 1 ) )->get_rsd_type_name3(), ref_pose->residue_type( 5 ).name3() );
		ResidueTypeConstraintGenerator rt_gen_2;
		rt_gen_2.set_id( "refpose_length_change_neg_control" );
		rt_gen_2.set_reference_pose( ref_pose );
		rt_gen_2.set_residue_selector( select_end );
		core::scoring::constraints::ConstraintCOPs const csts_2 = rt_gen_2.apply( *pose_ );
		TS_ASSERT_EQUALS( csts_2.size(), 1 );
		old_res = ref_pose->residue_type( 110 ).name3();
		new_res = ref_pose->residue_type( 110 ).name3();
		TR << "Old residue: " << old_res << " New residue: " << new_res << std::endl;
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< core::scoring::constraints::ResidueTypeConstraint const >( csts_2.at( 1 ) )->get_rsd_type_name3(), ref_pose->residue_type( 110 ).name3() );
		TS_ASSERT( utility::pointer::dynamic_pointer_cast< core::scoring::constraints::ResidueTypeConstraint const >( csts_2.at( 1 ) )->get_rsd_type_name3() != pose_->residue_type( 110 ).name3() );

	}



};

