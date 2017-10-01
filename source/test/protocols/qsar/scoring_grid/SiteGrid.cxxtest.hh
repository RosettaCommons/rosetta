// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   qsar/scoring_grid/SiteGrid.cxxtest.hh
/// @brief  test suite for SiteGrid setup and scoring
/// @author Darwin Fu (darwinyfu@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/UTracer.hh>

#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/UltraLightResidue.hh>
#include <core/chemical/ResidueType.hh>
#include <protocols/qsar/scoring_grid/SiteGridCreator.hh>
#include <protocols/qsar/scoring_grid/GridFactory.hh>
#include <protocols/qsar/scoring_grid/GridBase.hh>
#include <protocols/qsar/qsarMap.hh>
#include <utility/excn/Exceptions.hh>
#include <test/core/init_util.hh>
#include <utility/tag/Tag.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>

static basic::Tracer TR("protocols.qsar.scoring_grid.SiteGrid.cxxtest");

class SiteGridTests : public CxxTest::TestSuite {
public:

	void setUp() {
		core_init_with_additional_options("-extra_res_fa protocols/qsar/scoring_grid/LIG.params");
	}

	void tearDown() {}

	void test_score_function() {
		using namespace protocols::qsar::scoring_grid;

		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, "protocols/qsar/scoring_grid/2JJC.pdb" , core::import_pose::PDB_file);

		//148 Protein residues, ligand number 149
		core::conformation::UltraLightResidue ligand(pose.residue(149).get_self_ptr());
		core::Vector center(ligand.center());

		std::string xmlfile = " <SiteGrid grid_name=site residues=11,52/>\n";
		std::istringstream resource_options_stream( xmlfile );
		utility::tag::TagCOP grid_tag = utility::tag::Tag::create( resource_options_stream );

		GridBaseOP grid(GridFactory::get_instance()->new_grid(grid_tag));
		grid->initialize(center,30,0.25);
		grid->set_chain('X');
		grid->refresh(pose,center,'X');

		protocols::qsar::qsarMapOP map;
		TS_ASSERT_EQUALS(grid->score(ligand, 100, map), -4);

		utility::vector1<core::Real> correct_values;

		correct_values.push_back(-1);
		correct_values.push_back(-1);
		correct_values.push_back(0);
		correct_values.push_back(0);
		correct_values.push_back(0);
		correct_values.push_back(0);
		correct_values.push_back(-1);
		correct_values.push_back(-1);
		correct_values.push_back(0);
		correct_values.push_back(0);
		correct_values.push_back(0);
		correct_values.push_back(0);

		for ( core::Size i = 1; i<=ligand.natoms(); i++ ) {
			TR << "atom score for " << i << " is " << grid->atom_score(ligand,i,0) << std::endl;
			TS_ASSERT_EQUALS(grid->atom_score(ligand,i,map), correct_values[i]);
		}
	}

};
