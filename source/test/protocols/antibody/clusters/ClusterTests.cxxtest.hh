// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/antibody/clusters/ClusterTests.cxxtest.hh
/// @brief  tests for antibody cluster functions
/// @author Jared Adolf-Bryfogle

// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Project Headers
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/util.hh>
#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/antibody/CDRClusterEnum.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/rms_util.hh>

// Protocol Headers
#include <basic/Tracer.hh>

using namespace protocols::antibody;

static basic::Tracer TR("protocols.antibody.clusters.ClusterTests");
class ClusterTests : public CxxTest::TestSuite {
	core::pose::Pose ab_pose; //Full PDB

	AntibodyInfoOP ab_info;

	
public:
	
	void setUp(){
		
		core_init();
		core::import_pose::pose_from_pdb(ab_pose, "protocols/antibody/2j88.pdb");

		ab_info = new AntibodyInfo(ab_pose, Modified_AHO);
		ab_info->setup_CDR_clusters(ab_pose);
		TR <<"Setup"<<std::endl;
	}
	
	void tearDown(){
		ab_pose.clear();
	}
    
	void test_cluster_identification(){
		
		TS_ASSERT(check_if_pose_renumbered_for_clusters(ab_pose)); 
		
		TS_ASSERT_EQUALS("L1-11-1", ab_info->get_cluster_name(ab_info->get_CDR_cluster(l1).first));
		TS_ASSERT_EQUALS("L2-8-1", ab_info->get_cluster_name(ab_info->get_CDR_cluster(l2).first));
		TS_ASSERT_EQUALS("L3-8-1", ab_info->get_cluster_name(ab_info->get_CDR_cluster(l3).first));
	}
	
	void test_basic_cluster_functions(){
		TS_ASSERT_EQUALS(l1, ab_info->get_cdr_enum_for_cluster(L1_11_1));
		TS_ASSERT_EQUALS(11, ab_info->get_cluster_length(L1_11_1));
		TS_ASSERT_EQUALS(L1_11_1, ab_info->get_cluster_enum("L1-11-1"));
	}
	
	void test_cluster_constraints(){
		
		std::map< CDRNameEnum, bool > result = add_harmonic_cluster_constraints(ab_info, ab_pose);
		for (core::Size i = 1; i <= ab_info->get_total_num_CDRs(); ++i){
			CDRNameEnum cdr_name = static_cast<CDRNameEnum>(i);
			TS_ASSERT(result[cdr_name]);
		}
		
		ab_pose.remove_constraints();
		
		utility::vector1< core::scoring::constraints::ConstraintCOP > constraints;
		result = protocols::antibody::add_harmonic_cluster_constraints(ab_info, ab_pose, constraints);
		for (core::Size i = 1; i <= ab_info->get_total_num_CDRs(); ++i){
			CDRNameEnum cdr_name = static_cast<CDRNameEnum>(i);
			TS_ASSERT(result[cdr_name]);
		}
		
		ab_pose.remove_constraints();
		
		bool one_constraint_result = add_harmonic_cluster_constraint(ab_info, ab_pose, L1_11_1);
		TS_ASSERT(one_constraint_result);
		
	}
    
};

