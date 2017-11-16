// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/constraint_generator/MetalContactsConstraintGenerator.cxxtest.hh
/// @brief test suite for protocols::constraint_generator::MetalContactsConstraintGenerator
/// @author Sharon Guffy (guffy@email.unc.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/UTracer.hh>
#include <test/protocols/init_util.hh>

// Protocol headers
#include <protocols/constraint_generator/MetalContactsConstraintGenerator.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/types.hh>
// Basic headers
#include <basic/datacache/DataMap.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>

// Boost headers


// C++ headers

using namespace core::scoring::constraints;
using namespace core::scoring::func;
using namespace protocols::constraint_generator;

static basic::Tracer TR( "protocols.constraint_generator.MetalContactsConstraintGenerator.cxxtest.hh" );

class MetalContactsConstraintGeneratorTests : public CxxTest::TestSuite {

private:
	core::pose::PoseOP ref_pose_;
	core::Real const TOLERANCE = 1e-5;
public:
	void setUp()
	{
		protocols_init();
		ref_pose_ = core::import_pose::pose_from_file( "core/util/2c9v_stripped.pdb" );
	}

	void tearDown(){

	}





	//Things to test
	//Note: Selectors can only be specified through the datamap

	//DEFAULTS
	void test_default_values(){
		MetalContactsConstraintGenerator cst_gen;
		std::stringstream ss_default;
		utility::tag::TagOP tag( new utility::tag::Tag() );
		basic::datacache::DataMap datamap;
		TS_ASSERT_DELTA( cst_gen.get_dist_cutoff_multiplier(), 1.0, TOLERANCE );
		//Other numeric values have bogus/empty defaults
		TS_ASSERT( !cst_gen.get_score_against_internal_contacts() );
		TS_TRACE( "Default values" );
		ss_default << "<MetalContactsConstraintGenerator name=\"metal\" ligand_resnum=\"310\" ligand_atom_name=\"ZN\" />" << std::endl;
		tag->read( ss_default );
		TS_ASSERT_THROWS_NOTHING( cst_gen.parse_my_tag( tag, datamap ) );
		core::scoring::constraints::ConstraintCOPs constraints = cst_gen.apply( *ref_pose_ );
		//Metal has 4 coordinating residues
		//4 distance constraints
		//4 angle_about_contact constraints
		//4 dihedral_about_contact
		//6 angle about metal
		//6 dihedral about metal
		//6 dihedral3
		//= 30 total
		TS_ASSERT_EQUALS( constraints.size(), 30 );
	}
	//INVALID TAG
	//string and selector for either ligand or contacts
	//Nothing specified for ligand atom name
	//Neither string nor selector for ligand
	void test_invalid_tags(){
		MetalContactsConstraintGenerator cst_gen;
		std::stringstream ss_invalid_both_ligand;
		std::stringstream ss_invalid_both_contact;
		std::stringstream ss_invalid_neither_ligand;
		std::stringstream ss_invalid_no_atom;
		basic::datacache::DataMap datamap;
		core::select::residue_selector::ResidueIndexSelectorOP ligand_310( new core::select::residue_selector::ResidueIndexSelector );
		ligand_310->set_index( "310" );
		datamap.add( "ResidueSelector", "ligand_310", ligand_310 );
		core::select::residue_selector::ResidueIndexSelectorOP contact_218_226( new core::select::residue_selector::ResidueIndexSelector );
		contact_218_226->set_index( "218,226" );
		datamap.add( "ResidueSelector", "contact_218_226", contact_218_226 );

		TS_TRACE( "Invalid both ligand" );
		ss_invalid_both_ligand << "<MetalContactsConstraintGenerator name=\"metal\" ligand_resnum=\"310\" ligand_selector=\"ligand_310\" ligand_atom_name=\"ZN\" />" << std::endl;
		TS_TRACE( "Invalid neither ligand" );
		ss_invalid_neither_ligand << "<MetalContactsConstraintGenerator name=\"metal\" ligand_atom_name=\"ZN\" />" << std::endl;
		TS_TRACE( "Invalid no atom" );
		ss_invalid_no_atom << "<MetalContactsConstraintGenerator name=\"metal\" ligand_resnum=\"310\" />" << std::endl;
		TS_TRACE( "Invalid both contact" );
		ss_invalid_both_contact << "<MetalContactsConstraintGenerator name=\"metal\" ligand_resnum=\"310\" ligand_atom_name=\"ZN\" contact_selector=\"contact_218_226\" contact_resnums=\"218,226\" />" << std::endl;




		TS_TRACE( "Invalid both ligand" );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss_invalid_both_ligand );
		set_throw_on_next_assertion_failure();
		TS_ASSERT_THROWS_ANYTHING( cst_gen.parse_my_tag( tag, datamap ) );

		TS_TRACE( "Invalid both contact" );
		tag = utility::tag::TagOP( new utility::tag::Tag() );
		tag->read( ss_invalid_both_contact );
		set_throw_on_next_assertion_failure();
		TS_ASSERT_THROWS_ANYTHING( cst_gen.parse_my_tag( tag, datamap ) );

		TS_TRACE( "Invalid neither ligand" );
		tag = utility::tag::TagOP( new utility::tag::Tag() );
		tag->read( ss_invalid_neither_ligand );
		set_throw_on_next_assertion_failure();
		TS_ASSERT_THROWS_ANYTHING( cst_gen.parse_my_tag( tag, datamap ) );

		TS_TRACE( "Invalid no atom" );
		tag = utility::tag::TagOP( new utility::tag::Tag() );
		tag->read( ss_invalid_no_atom );
		set_throw_on_next_assertion_failure();
		TS_ASSERT_THROWS_ANYTHING( cst_gen.parse_my_tag( tag, datamap ) );
	}


	//OPTIONS: all in one test function
	//Numeric (use TS_DELTA with TOLERANCE):
	//dist_cutoff_multiplier
	//ideal_angle_about_contact
	//ideal_dihedral_about_contact
	//ideal_angle_about_metal
	//ideal_dihedral_about_metal
	//ideal_dihedral_3
	//ideal_distance
	//Not numeric( use TS_EQUALS ):
	//score_against_internal_contacts (bool)--can test value, not result
	//base_atom_name, base_base_atom_name
	void test_set_options_in_tag(){
		std::stringstream ss_set_multiplier;
		std::stringstream ss_set_ideal_values;
		std::stringstream ss_set_bool;
		std::stringstream ss_set_bases;
		basic::datacache::DataMap datamap;

		TS_TRACE( "Set multiplier" );
		ss_set_multiplier << "<MetalContactsConstraintGenerator name=\"metal\" ligand_resnum=\"310\" ligand_atom_name=\"ZN\" dist_cutoff_multiplier=\"2.0\" />" << std::endl;
		TS_TRACE( "Set ideal values" );
		ss_set_ideal_values << "<MetalContactsConstraintGenerator name=\"metal\" ligand_resnum=\"310\" ligand_atom_name=\"ZN\" ideal_distance=\"2.2\" ideal_angle_about_contact=\"120\" ideal_dihedral_about_contact=\"0,180\" ideal_angle_about_metal=\"109.5\" ideal_dihedral_about_metal=\"30,60,90,120\" ideal_dihedral_3=\"120\"/>" << std::endl;
		TS_TRACE( "Set score against internal contacts" );
		ss_set_bool << "<MetalContactsConstraintGenerator name=\"metal\" ligand_resnum=\"310\" ligand_atom_name=\"ZN\" score_against_internal_contacts=\"true\" />" << std::endl;
		TS_TRACE( "Set bases" );
		ss_set_bases << "<MetalContactsConstraintGenerator name=\"metal\" ligand_resnum=\"310\" ligand_atom_name=\"ZN\" base_atom_name=\"CB\" base_base_atom_name=\"CA\" />" << std::endl;


		MetalContactsConstraintGenerator cst_gen;
		utility::tag::TagOP tag( new utility::tag::Tag() );
		TS_TRACE( "Set multiplier" );
		tag->read( ss_set_multiplier );
		TS_ASSERT_THROWS_NOTHING( cst_gen.parse_my_tag( tag, datamap ) );
		TS_ASSERT_DELTA( cst_gen.get_dist_cutoff_multiplier(), 2.0, TOLERANCE );

		cst_gen = MetalContactsConstraintGenerator();
		tag = utility::tag::TagOP( new utility::tag::Tag() );
		TS_TRACE( "Set ideal values" );
		tag->read( ss_set_ideal_values );
		TS_ASSERT_THROWS_NOTHING( cst_gen.parse_my_tag( tag, datamap ) );
		TS_ASSERT_DELTA( cst_gen.get_ideal_distance(), 2.2, TOLERANCE );
		TS_ASSERT_EQUALS( cst_gen.get_ideal_angle_about_contact().size(), 1 );
		TS_ASSERT_EQUALS( cst_gen.get_ideal_dihedral_about_contact().size(), 2 );
		TS_ASSERT_EQUALS( cst_gen.get_ideal_angle_about_metal().size(), 1 );
		TS_ASSERT_EQUALS( cst_gen.get_ideal_dihedral_about_metal().size(), 4 );
		TS_ASSERT_EQUALS( cst_gen.get_ideal_dihedral_3().size(), 1 );

		cst_gen = MetalContactsConstraintGenerator();
		tag = utility::tag::TagOP( new utility::tag::Tag() );
		tag->read( ss_set_bool );
		TS_TRACE( "Set score against internal contacts" );
		TS_ASSERT_THROWS_NOTHING( cst_gen.parse_my_tag( tag, datamap ) );
		TS_ASSERT_EQUALS( cst_gen.get_score_against_internal_contacts(), true );

		cst_gen = MetalContactsConstraintGenerator();
		tag = utility::tag::TagOP( new utility::tag::Tag() );
		TS_TRACE( "Set atom bases" );
		tag->read( ss_set_bases );
		TS_ASSERT_THROWS_NOTHING( cst_gen.parse_my_tag( tag, datamap ) );
		TS_ASSERT_EQUALS( cst_gen.get_base_atom_name(), "CB" );
		TS_ASSERT_EQUALS( cst_gen.get_base_base_atom_name(), "CA" );

	}

	//SPECIFYING LIGAND (only 1 residue at a time for now)
	//residue selector
	//resnum string
	void test_set_ligand(){
		std::stringstream ss_res_selector;
		std::stringstream ss_string;
		basic::datacache::DataMap datamap;
		core::select::residue_selector::ResidueIndexSelectorOP ligand_310( new core::select::residue_selector::ResidueIndexSelector );
		ligand_310->set_index( "310" );
		datamap.add( "ResidueSelector", "ligand_310", ligand_310 );

		//set up stringstreams
		ss_res_selector << "<MetalContactsConstraintGenerator name=\"metal\" ligand_selector=\"ligand_310\" ligand_atom_name=\"ZN\" />" << std::endl;
		ss_string << "<MetalContactsConstraintGenerator name=\"metal\" ligand_resnum=\"310\" ligand_atom_name=\"ZN\" />" << std::endl;

		//res_selector
		MetalContactsConstraintGenerator cst_gen;
		utility::tag::TagOP tag( new utility::tag::Tag() );
		TS_TRACE( "Ligand specified with selector" );
		tag->read( ss_res_selector );
		TS_ASSERT_THROWS_NOTHING( cst_gen.parse_my_tag( tag, datamap ) );
		TS_ASSERT( cst_gen.get_use_ligand_selector() );
		core::scoring::constraints::ConstraintCOPs constraints = cst_gen.apply( *ref_pose_ );
		TS_ASSERT_EQUALS( constraints.size(), 30 );

		//string
		cst_gen = MetalContactsConstraintGenerator();
		tag = utility::tag::TagOP( new utility::tag::Tag() );
		TS_TRACE( "Ligand specified with string" );
		tag->read( ss_string );
		TS_ASSERT_THROWS_NOTHING( cst_gen.parse_my_tag( tag, datamap ) );
		TS_ASSERT( !cst_gen.get_use_ligand_selector() );
		constraints = cst_gen.apply( *ref_pose_ );
		TS_ASSERT_EQUALS( constraints.size(), 30 );

	}

	//SPECIFYING CONTACTS
	//residue selector
	//resnum string
	void test_set_contacts(){
		std::stringstream ss_res_selector;
		std::stringstream ss_string;
		basic::datacache::DataMap datamap;
		core::select::residue_selector::ResidueIndexSelectorOP contact_218_226( new core::select::residue_selector::ResidueIndexSelector );
		contact_218_226->set_index( "218,226" );
		datamap.add( "ResidueSelector", "contact_218_226", contact_218_226 );

		ss_res_selector << "<MetalContactsConstraintGenerator name=\"metal\" ligand_resnum=\"310\" ligand_atom_name=\"ZN\" contact_selector=\"contact_218_226\" />" << std::endl;
		ss_string << "<MetalContactsConstraintGenerator name=\"metal\" ligand_resnum=\"310\" ligand_atom_name=\"ZN\" contact_resnums=\"218,226\" />" << std::endl;

		MetalContactsConstraintGenerator cst_gen;
		utility::tag::TagOP tag( new utility::tag::Tag() );
		TS_TRACE( "Contacts specified with selector" );
		tag->read( ss_res_selector );
		TS_ASSERT_THROWS_NOTHING( cst_gen.parse_my_tag( tag, datamap ) );
		TS_ASSERT( !cst_gen.get_use_ligand_selector() );
		TS_ASSERT( cst_gen.get_use_contact_selector() );
		core::scoring::constraints::ConstraintCOPs constraints = cst_gen.apply( *ref_pose_ );
		//Only 2 residues selected--should make how many constraints?
		//2 distance/angle/dihedral for each plus 1 angle/2 dihedral between: 9 total
		TS_ASSERT_EQUALS( constraints.size(), 9 );

		cst_gen = MetalContactsConstraintGenerator();
		tag = utility::tag::TagOP( new utility::tag::Tag() );
		TS_TRACE( "Contacts specified with string" );
		tag->read( ss_string );
		TS_ASSERT_THROWS_NOTHING( cst_gen.parse_my_tag( tag, datamap ) );
		TS_ASSERT( !cst_gen.get_use_ligand_selector() );
		TS_ASSERT( !cst_gen.get_use_contact_selector() );
		constraints = cst_gen.apply( *ref_pose_ );
		TS_ASSERT_EQUALS( constraints.size(), 9 );
	}





};

