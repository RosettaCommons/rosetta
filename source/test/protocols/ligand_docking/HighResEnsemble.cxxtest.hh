// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/UTracer.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/UltraLightResidue.hh>
#include <protocols/qsar/scoring_grid/GridManager.hh>
#include <core/pose/util.hh>
#include <protocols/rigid/RB_geometry.hh>

#include <protocols/ligand_docking/HighResEnsemble.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/random/random.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <utility/vector1.hh>
//#include <test/core/init_util.hh>


static basic::Tracer TR("protocols.ligand_docking.HighResEnsemble.cxxtest");


class HighResEnsembleTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init_with_additional_options("-extra_res_fa protocols/ligand_docking/ZNx.params protocols/ligand_docking/7cpa.params");

	}

	void tearDown() {}

	//inline void core_init_with_additional_options( std::string const & commandline_in );

	void test_spearman() {

		utility::vector1<std::pair<core::Size, core::Real> > vector_1;
		utility::vector1<std::pair<core::Size, core::Real> > vector_2;

		vector_1.push_back(std::make_pair(1,1));
		vector_1.push_back(std::make_pair(2,2));
		vector_1.push_back(std::make_pair(3,3));
		vector_1.push_back(std::make_pair(4,4));
		vector_1.push_back(std::make_pair(5,5));

		vector_2.push_back(std::make_pair(1,1));
		vector_2.push_back(std::make_pair(2,4));
		vector_2.push_back(std::make_pair(3,3));
		vector_2.push_back(std::make_pair(4,2));
		vector_2.push_back(std::make_pair(5,5));

		TS_ASSERT_EQUALS(protocols::ligand_docking::spearman(vector_1, vector_2), 0.6);
	}

	void test_vector_to_rank() {

		utility::vector1<std::pair<core::Size, core::Real> > data_vector;
		utility::vector1<std::pair<core::Size, core::Real> > correct_ranks;

		data_vector.push_back(std::make_pair(1,1.0));
		data_vector.push_back(std::make_pair(2,2.0));
		data_vector.push_back(std::make_pair(3,2.0));
		data_vector.push_back(std::make_pair(4,2.0));
		data_vector.push_back(std::make_pair(5,3.5));
		data_vector.push_back(std::make_pair(6,4.0));
		data_vector.push_back(std::make_pair(7,4.5));
		data_vector.push_back(std::make_pair(8,5.0));
		data_vector.push_back(std::make_pair(9,5.0));

		correct_ranks.push_back(std::make_pair(1,1));
		correct_ranks.push_back(std::make_pair(2,3));
		correct_ranks.push_back(std::make_pair(3,3));
		correct_ranks.push_back(std::make_pair(4,3));
		correct_ranks.push_back(std::make_pair(5,5));
		correct_ranks.push_back(std::make_pair(6,6));
		correct_ranks.push_back(std::make_pair(7,7));
		correct_ranks.push_back(std::make_pair(8,8.5));
		correct_ranks.push_back(std::make_pair(9,8.5));

		protocols::ligand_docking::vector_to_rank(data_vector);

		for ( core::Size i=1; i <= correct_ranks.size(); i++ ) {
			TS_ASSERT_EQUALS(data_vector[i].second, correct_ranks[i].second);
		}
	}

};

