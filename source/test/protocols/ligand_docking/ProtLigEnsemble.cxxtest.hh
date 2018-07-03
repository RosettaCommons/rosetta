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

#include <protocols/ligand_docking/ProtLigEnsemble.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/random/random.hh>
#include <ObjexxFCL/FArray2D.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <utility/vector1.hh>
//#include <test/core/init_util.hh>


static basic::Tracer TR("protocols.ligand_docking.ProtLigEnsemble.cxxtest");

class ProtLigEnsembleTests : public CxxTest::TestSuite {

public:

	void setUp() {
		// core_init_with_additional_options("-extra_res_fa protocols/ligand_docking/ZNx.params protocols/ligand_docking/7cpa.params");

	}

	void tearDown() {}

	//inline void core_init_with_additional_options( std::string const & commandline_in );

	void test_process_line(){

		protocols::ligand_docking::ProtLigPair_info pair_info;

		std::string wt_line = "WT X 2.0";

		pair_info = protocols::ligand_docking::process_line(wt_line);
		TS_ASSERT_EQUALS(pair_info.wild_type, true);
		TS_ASSERT_EQUALS(pair_info.has_bind, true);
		TS_ASSERT_EQUALS(pair_info.lig_chain, 'X');
		TS_ASSERT_EQUALS(pair_info.bind_data, 2.0);

		std::string mut_line = "100 A X";

		pair_info = protocols::ligand_docking::process_line(mut_line);
		TS_ASSERT_EQUALS(pair_info.wild_type, false);
		TS_ASSERT_EQUALS(pair_info.has_bind, false);
		TS_ASSERT_EQUALS(pair_info.mut_resid, 100);
		TS_ASSERT_EQUALS(pair_info.mut_target, 'A');
		TS_ASSERT_EQUALS(pair_info.lig_chain, 'X');
		TS_ASSERT_EQUALS(pair_info.bind_data, 0.0);
	}

	void test_sort_by_binding() {

		utility::vector1<protocols::ligand_docking::ProtLigPair_info> infos_to_sort;
		protocols::ligand_docking::ProtLigPair_info pair_info;

		pair_info.has_bind = true;
		pair_info.bind_data = -3.0;
		infos_to_sort.push_back(pair_info);

		pair_info.has_bind = false;
		pair_info.bind_data = 0.0;
		infos_to_sort.push_back(pair_info);

		pair_info.has_bind = true;
		pair_info.bind_data = -1.0;
		infos_to_sort.push_back(pair_info);

		pair_info.has_bind = true;
		pair_info.bind_data = -2.0;
		infos_to_sort.push_back(pair_info);

		pair_info.has_bind = false;
		pair_info.bind_data = 0.0;
		infos_to_sort.push_back(pair_info);

		pair_info.has_bind = true;
		pair_info.bind_data = -1.0;
		infos_to_sort.push_back(pair_info);

		std::sort(infos_to_sort.begin(), infos_to_sort.end(), protocols::ligand_docking::sort_by_binding);

		utility::vector1<std::pair<core::Size, core::Real> > correct_vector;

		correct_vector.push_back(std::make_pair(true,-3));
		correct_vector.push_back(std::make_pair(true,-2));
		correct_vector.push_back(std::make_pair(true,-1));
		correct_vector.push_back(std::make_pair(true,-1));
		correct_vector.push_back(std::make_pair(false,0));
		correct_vector.push_back(std::make_pair(false,0));

		for ( core::Size i=1; i<= infos_to_sort.size(); i++ ) {
			TS_ASSERT_EQUALS(infos_to_sort[i].has_bind, correct_vector[i].first);
			TS_ASSERT_EQUALS(infos_to_sort[i].bind_data, correct_vector[i].second);
		}

	}


};

