// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   qsar/scoring_grid/ScoreNormalization.cxxtest.hh
/// @brief  test suite for Score Normalization
/// @author Sam DeLuca (samuel.l.deluca@vanderbilt.edu)

#include <cxxtest/TestSuite.h>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <protocols/qsar/scoring_grid/ScoreNormalization.hh>
#include <utility/excn/Exceptions.hh>
#include <test/core/init_util.hh>

class ScoreNormalizationTest : public CxxTest::TestSuite {
public:

	ScoreNormalizationTest()
	{

	}

	void setUp() {
		core_init();

	}

	void tearDown() {}

	void test_get_score_normalization_function() {
		using namespace protocols::qsar::scoring_grid;

		ScoreNormalizationOP normalization_class =
			protocols::qsar::scoring_grid::get_score_normalization_function("HeavyAtomNormalization");
		TS_ASSERT(normalization_class->get_name() == "HeavyAtomNormalization");

		normalization_class =
			protocols::qsar::scoring_grid::get_score_normalization_function("AllAtomNormalization");
		TS_ASSERT(normalization_class->get_name() == "AllAtomNormalization");

		normalization_class =
			protocols::qsar::scoring_grid::get_score_normalization_function("ChiAngleNormalization");
		TS_ASSERT(normalization_class->get_name() == "ChiAngleNormalization");

		normalization_class =
			protocols::qsar::scoring_grid::get_score_normalization_function("MolecularWeightNormalization");
		TS_ASSERT(normalization_class->get_name() == "MolecularWeightNormalization");


		try {
			normalization_class =
				protocols::qsar::scoring_grid::get_score_normalization_function("NonExistantNormalization");
			TS_ASSERT(false); // last line ought to throw an exception
		}catch(utility::excn::EXCN_RosettaScriptsOption const &)
		{
			TS_ASSERT(true);
		}
	}

	void test_heavy_atom_normalization(){
		using namespace protocols::qsar::scoring_grid;

		core::conformation::ResidueCOPs residue_list;
		core::chemical::ResidueTypeSetCOP fa_std_set = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");
		core::chemical::ResidueType ser_type = fa_std_set->name_map("SER");
		core::chemical::ResidueType gly_type = fa_std_set->name_map("GLY");

		residue_list.push_back(new core::conformation::Residue(ser_type,false));
		residue_list.push_back(new core::conformation::Residue(gly_type,false));

		ScoreNormalizationOP normalization_class =
			protocols::qsar::scoring_grid::get_score_normalization_function("HeavyAtomNormalization");

		core::Real input_score = 1.0;
		core::Real ser_score = (*normalization_class)(input_score,*residue_list[1]);
		core::Real ser_gly_score = (*normalization_class)(input_score,residue_list);

		TS_ASSERT_DELTA(ser_score,0.166,0.001); // SER has 6 heavy atoms

		TS_ASSERT_DELTA(ser_gly_score,0.1,0.001); // SER-GLY has 10 heavy atoms
	}

	void test_all_atom_normalization(){
		using namespace protocols::qsar::scoring_grid;

		core::conformation::ResidueCOPs residue_list;
		core::chemical::ResidueTypeSetCOP fa_std_set = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");
		core::chemical::ResidueType ser_type = fa_std_set->name_map("SER");
		core::chemical::ResidueType gly_type = fa_std_set->name_map("GLY");

		residue_list.push_back(new core::conformation::Residue(ser_type,false));
		residue_list.push_back(new core::conformation::Residue(gly_type,false));

		ScoreNormalizationOP normalization_class =
			protocols::qsar::scoring_grid::get_score_normalization_function("AllAtomNormalization");

		core::Real input_score = 1.0;
		core::Real ser_score = (*normalization_class)(input_score,*residue_list[1]);
		core::Real ser_gly_score = (*normalization_class)(input_score,residue_list);

		TS_ASSERT_DELTA(ser_score,0.090,0.001); // SER has 11 atoms

		TS_ASSERT_DELTA(ser_gly_score,0.055,0.001); // SER-GLY has 18 atoms
	}

	void test_chi_angle_normalization(){
		using namespace protocols::qsar::scoring_grid;

		core::conformation::ResidueCOPs residue_list;
		core::chemical::ResidueTypeSetCOP fa_std_set = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");
		core::chemical::ResidueType ser_type = fa_std_set->name_map("SER");
		core::chemical::ResidueType gly_type = fa_std_set->name_map("GLY");

		residue_list.push_back(new core::conformation::Residue(ser_type,false));
		residue_list.push_back(new core::conformation::Residue(gly_type,false));

		ScoreNormalizationOP normalization_class =
			protocols::qsar::scoring_grid::get_score_normalization_function("ChiAngleNormalization");

		core::Real input_score = 1.0;
		core::Real ser_score = (*normalization_class)(input_score,*residue_list[1]);
		core::Real ser_gly_score = (*normalization_class)(input_score,residue_list);

		TS_ASSERT_DELTA(ser_score,0.5,4); // SER has 2 chi angles

		TS_ASSERT_DELTA(ser_gly_score,0.5,4); // SER-GLY has 2
	}

	void test_mol_weight_normalization(){
		using namespace protocols::qsar::scoring_grid;

		core::conformation::ResidueCOPs residue_list;
		core::chemical::ResidueTypeSetCOP fa_std_set = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");
		core::chemical::ResidueType ser_type = fa_std_set->name_map("SER");
		core::chemical::ResidueType gly_type = fa_std_set->name_map("GLY");

		residue_list.push_back(new core::conformation::Residue(ser_type,false));
		residue_list.push_back(new core::conformation::Residue(gly_type,false));

		ScoreNormalizationOP normalization_class =
			protocols::qsar::scoring_grid::get_score_normalization_function("MolecularWeightNormalization");

		core::Real input_score = 1.0;
		core::Real ser_score = (*normalization_class)(input_score,*residue_list[1]);
		core::Real ser_gly_score = (*normalization_class)(input_score,residue_list);

		TS_ASSERT_DELTA(ser_score,0.0114,0.001); // SER has MW 87.08

		TS_ASSERT_DELTA(ser_gly_score,0.0069,0.0001); // SER-GLY has MW 144.13
	}

};
