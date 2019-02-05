// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/select/residue_selector/ResiduePropertySelectorTests.cxxtest.hh
/// @brief  Unit tests for the ResiduePropertySelector.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

#include <test/core/select/residue_selector/utilities_for_testing.hh>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueProperty.hh>
#include <core/select/residue_selector/ResiduePropertySelector.hh>
#include <protocols/rosetta_scripts/XmlObjects.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>

static basic::Tracer TR("ResiduePropertySelectorTests");


class ResiduePropertySelectorTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init_with_additional_options( "-include_sugars" );
		pose_from_file( pose_, "core/chemical/carbohydrates/gp120_2glycans_man5_chainC.pdb" , core::import_pose::PDB_file);
	}

	void test_selector(){
		using namespace core::select::residue_selector;
		using namespace core::chemical;

		utility::vector1< bool > protein_residues(pose_.size(), false);
		utility::vector1< bool > glycan_residues(pose_.size(), false);
		utility::vector1< bool > all_residues(pose_.size(), false);

		for ( core::Size i = 1; i <= pose_.size(); ++i ) {
			if ( pose_.residue_type(i).is_protein() ) {
				protein_residues[i] = true;
				all_residues[i] = true;
			}
			if ( pose_.residue_type(i).is_carbohydrate() ) {
				glycan_residues[i] = true;
				all_residues[i] = true;
			}
		}

		//Test Protein
		ResiduePropertySelector selector = ResiduePropertySelector();
		selector.set_property( PROTEIN );
		utility::vector1< bool > protein_subset = selector.apply(pose_);
		compare_bool_vector( protein_subset, protein_residues);

		//Test Carbohydrates
		selector.set_property( CARBOHYDRATE );
		utility::vector1< bool > glycan_subset = selector.apply(pose_);
		compare_bool_vector( glycan_subset, glycan_residues);

		//Test PROTEIN OR CARBOHDYRATES
		selector.set_selection_logic( core::select::residue_selector::or_logic );
		selector.add_property( PROTEIN );
		utility::vector1< bool > all_subset = selector.apply(pose_);
		compare_bool_vector( all_subset, all_residues);

		//Test CARBOHDYRATE AND VIRTUAL
		core::Size n_to_virt = 0;
		utility::vector1< bool > virtual_glycan_residues( pose_.size(), false);

		for ( core::Size i = 1; i <= pose_.size(); ++i ) {

			//Only virtualize a subset to make sure the AND logic is working correctly.
			if ( pose_.residue_type(i).is_carbohydrate() && n_to_virt < 4 ) {
				pose_.real_to_virtual(i);
				n_to_virt+=1;
				virtual_glycan_residues[ i ] = true;
			}
		}

		utility::vector1< core::chemical::ResidueProperty > combined_properties;
		combined_properties.push_back( CARBOHYDRATE );
		combined_properties.push_back( VIRTUAL_RESIDUE );

		selector.set_selection_logic( core::select::residue_selector::and_logic );
		selector.set_properties( combined_properties );

		utility::vector1< bool > virtual_glycan_subset = selector.apply(pose_);
		compare_bool_vector( virtual_glycan_subset, virtual_glycan_residues );

	}
	void tearDown(){

	}

	void test_xml_interface(){


		using namespace core::select::residue_selector;
		using namespace core::chemical;
		using namespace protocols::rosetta_scripts;



		utility::vector1< bool > protein_residues(pose_.size(), false);
		utility::vector1< bool > glycan_residues(pose_.size(), false);
		utility::vector1< bool > all_residues(pose_.size(), false);

		for ( core::Size i = 1; i <= pose_.size(); ++i ) {
			if ( pose_.residue_type(i).is_protein() ) {
				protein_residues[i] = true;
				all_residues[i] = true;
			}
			if ( pose_.residue_type(i).is_carbohydrate() ) {
				glycan_residues[i] = true;
				all_residues[i] = true;
			}
		}

		//Test Protein
		std::string xml_str = "<ResiduePropertySelector name=\"gen\" properties=\"protein\" />";
		ResidueSelectorOP selector = XmlObjects::static_get_residue_selector(xml_str);
		utility::vector1< bool > protein_subset = selector->apply(pose_);
		compare_bool_vector( protein_subset, protein_residues);

		//Test Carbohydrates
		xml_str = "<ResiduePropertySelector name=\"gen\" properties=\"carbohydrate\" />";
		selector = XmlObjects::static_get_residue_selector(xml_str);
		utility::vector1< bool > glycan_subset = selector->apply(pose_);
		compare_bool_vector( glycan_subset, glycan_residues);

		//Test PROTEIN OR CARBOHDYRATES
		xml_str = "<ResiduePropertySelector name=\"gen\" properties=\"protein,carbohydrate\" logic=\"or_logic\" />";
		selector = XmlObjects::static_get_residue_selector(xml_str);
		utility::vector1< bool > all_subset = selector->apply(pose_);
		compare_bool_vector( all_subset, all_residues);

		//Test CARBOHDYRATE AND VIRTUAL
		core::Size n_to_virt = 0;
		utility::vector1< bool > virtual_glycan_residues( pose_.size(), false);

		for ( core::Size i = 1; i <= pose_.size(); ++i ) {

			//Only virtualize a subset to make sure the AND logic is working correctly.
			if ( pose_.residue_type(i).is_carbohydrate() && n_to_virt < 4 ) {
				pose_.real_to_virtual(i);
				n_to_virt+=1;
				virtual_glycan_residues[ i ] = true;
			}
		}

		xml_str = "<ResiduePropertySelector name=\"gen\" properties=\"virtual_residue,carbohydrate\" logic=\"and_logic\" />";
		selector = XmlObjects::static_get_residue_selector(xml_str);

		utility::vector1< bool > virtual_glycan_subset = selector->apply(pose_);
		compare_bool_vector( virtual_glycan_subset, virtual_glycan_residues );
	}


private:

	core::pose::Pose pose_;




};
