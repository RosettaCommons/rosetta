// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  protocols/antibody/residue_selector/AntibodyResidueSelectorsTest.cxxtest.hh
/// @brief  Test residue_selectors for antibodies.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/core/select/residue_selector/utilities_for_testing.hh>

// Project Headers
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/util.hh>
#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/antibody/clusters/CDRClusterEnum.hh>

#include <protocols/antibody/task_operations/DisableCDRsOperation.hh>
#include <protocols/antibody/task_operations/DisableAntibodyRegionOperation.hh>
#include <protocols/antibody/task_operations/RestrictToCDRsAndNeighbors.hh>
#include <protocols/antibody/task_operations/AddCDRProfilesOperation.hh>
#include <protocols/antibody/task_operations/AddCDRProfileSetsOperation.hh>

#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/residue_selector/CDRResidueSelector.hh>
#include <protocols/antibody/residue_selector/AntibodyRegionSelector.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// Protocol Headers
#include <basic/Tracer.hh>
#include <boost/foreach.hpp>

#include <utility/vector1.hh>
#include <utility/string_util.hh>

// Project Headers

static THREAD_LOCAL basic::Tracer TR("protocols.antibody.residue_selector.AntibodyResidueSelectorsTest");

using namespace protocols::antibody;
using namespace protocols::antibody::clusters;
using namespace protocols::antibody::task_operations;
using namespace protocols::antibody::design;
using namespace protocols::antibody::residue_selector;
using namespace core::pack::task::operation;

using namespace core::pack::task;
using utility::vector1;
using utility::to_string;

class AntibodyResidueSelectorsTest : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init();
		core::import_pose::pose_from_file(pose_, "protocols/antibody/aho_with_antigen.pdb", core::import_pose::PDB_file); //AHO renumbered pose
		core::import_pose::pose_from_file(pose_chothia_, "protocols/antibody/1bln_AB_chothia.pdb", core::import_pose::PDB_file);

		//TenA NeighborGraph setup
		core::scoring::ScoreFunctionOP score = core::scoring::get_score_function();
		score->score(pose_);
		score->score(pose_chothia_);

		ab_info_ = AntibodyInfoOP( new AntibodyInfo(pose_, AHO_Scheme, North) );
		ab_info_chothia_ = AntibodyInfoOP( new AntibodyInfo(pose_chothia_));
	}

	void tearDown(){

	}

	void test_selectors() {
		utility::vector1< bool > cdrs(8, false);
		utility::vector1< bool > all_cdrs(8, true);

		cdrs[ l1 ] = true;
		cdrs[ l2 ] = true;

		all_cdrs[ l4 ] = false;
		all_cdrs[ h4 ] = false;


		//Setup what we are testing on. These functions themselves are tested elsewhere.
		utility::vector1< bool > correct_cdr_residues = get_cdr_residues( cdrs );
		utility::vector1< bool > correct_all_cdr_residues = get_cdr_residues( all_cdrs );


		//Test CDRResidueSelector
		CDRResidueSelector selector = CDRResidueSelector(ab_info_, all_cdrs);
		utility::vector1< bool > cdr_residues = selector.apply(pose_);
		//TR << "cor" << to_string( correct_cdr_residues) << std::endl;
		//TR << "new" << to_string( cdr_residues) << std::endl;
		//TR.flush();
		TS_ASSERT( cdr_residues == correct_all_cdr_residues);

		// ->Use task to test.

		selector.set_cdrs( cdrs );
		cdr_residues = selector.apply(pose_);
		//TR << "cor" << to_string( correct_all_cdr_residues) << std::endl;
		//TR << "new" << to_string( cdr_residues) << std::endl;
		//TR.flush();
		TS_ASSERT( cdr_residues == correct_cdr_residues);

		//Test AntibodyRegionSelector
		AntibodyRegionSelector region_selector = AntibodyRegionSelector(ab_info_);
		utility::vector1< bool > region_residues;
		utility::vector1< bool > correct_region_residues;

		for (core::Size i = 1; i <=3; ++i){

			AntibodyRegionEnum region = static_cast<AntibodyRegionEnum>( i );
			correct_region_residues = get_region_residues( region );
			region_selector.set_region(region);
			region_residues = region_selector.apply( pose_ );
			TS_ASSERT( region_residues ==  correct_region_residues);
		}


	}

	utility::vector1<bool> get_cdr_residues(utility::vector1< bool > const & cdrs){

		utility::vector1<bool> residues(pose_.total_residue(), false);
		for (core::Size i = 1; i <= 8; ++i){
			CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
			if (cdrs[ i ]){
				for (core::Size resnum = ab_info_->get_CDR_start(cdr, pose_); resnum <= ab_info_->get_CDR_end(cdr, pose_); ++resnum){
					//TR << resnum << std::endl;
					//TR.flush();
					residues[ resnum ] = true;
				}

			}

		}
		return residues;

	}

	utility::vector1<bool> get_region_residues(AntibodyRegionEnum region){

		utility::vector1< bool > region_residues(pose_.total_residue(), false);
		for (core::Size i = 1; i <= pose_.total_residue(); ++i){
			if (region == ab_info_->get_region_of_residue(pose_, i)){
				region_residues[ i ] = true;
			}
		}
		return region_residues;
	}
private:

	core::pose::Pose pose_;
	core::pose::Pose pose_chothia_;

	AntibodyInfoOP ab_info_;
	AntibodyInfoOP ab_info_chothia_;





};



