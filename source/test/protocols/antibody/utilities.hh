// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/antibody/constraints/AntibodyConstraintTests.cxxtest.hh
/// @brief  Utility functions for Antibody Unit tests
/// @author Jared Adolf-Bryfogle

#ifndef INCLUDED_protocols_antibody_utilities_HH
#define INCLUDED_protocols_antibody_utilities_HH

// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

#include <protocols/antibody/AntibodyInfo.hh>

// Core headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>

// C++ headers
#include <string>
#include <iostream>



inline
void
output_or_test(core::pack::task::TaskFactoryCOP tf,core::pose::Pose const & pose, bool first_run, std::string name, std::string inpath, std::string first_run_outpath) {
	

	std::string inname = inpath+"/"+name+".u";
	std::string outname = first_run_outpath+"/"+name+".u";
	if (first_run){
		std::cout <<"////"<<  std::endl << inname << std::endl << std::endl;
		core::pack::task::PackerTaskOP task = tf->create_task_and_apply_taskoperations(pose);
		//task->show(std::cout);
			
		//Write to a file, which we can manually check, diff, and rename if needed.
		std::ofstream OUT;
			
		std::cout << outname << std::endl;
		OUT.open(outname.c_str());
		task->show(OUT);
		OUT << std::endl << std::endl;
		OUT.close();
			
	}
	else{
		test::UTracer UT(inname, true);
		std::cout <<"////"<<  std::endl << inname << std::endl << std::endl;
		//tf->create_task_and_apply_taskoperations(pose)->show(std::cout);
		tf->create_task_and_apply_taskoperations(pose)->show(UT);
	}
}

inline
void
assert_region_design_is_disabled(
	core::pose::Pose const & pose,
	core::pack::task::PackerTaskOP task,
	protocols::antibody::AntibodyInfoCOP ab_info,
	protocols::antibody::AntibodyRegionEnum region) {
	
	
	std::string r;
	if (region == protocols::antibody::antigen_region) {r = "antigen";};
	if (region == protocols::antibody::framework_region) {r = "framework";}
	if (region == protocols::antibody::cdr_region) {r = "cdr";}
		
	std::cout <<"Checking region: " << r << std::endl;
	for (core::Size i = 1; i <= pose.total_residue(); ++i ){
		if (ab_info->get_region_of_residue(pose, i) == region){
			TS_ASSERT_EQUALS(task->design_residue( i ), false );
		}
	}
}

inline
void
assert_region_packing_is_disabled(
	core::pose::Pose const & pose,
	core::pack::task::PackerTaskOP task,
	protocols::antibody::AntibodyInfoCOP ab_info,
	protocols::antibody::AntibodyRegionEnum region) {
	
	
	std::string r;
	if (region == protocols::antibody::antigen_region) {r = "antigen";};
	if (region == protocols::antibody::framework_region) {r = "framework";}
	if (region == protocols::antibody::cdr_region) {r = "cdr";}
		
	std::cout <<"Checking region: " << r << std::endl;
	for (core::Size i = 1; i <= pose.total_residue(); ++i ){
		if (ab_info->get_region_of_residue(pose, i) == region){
			TS_ASSERT_EQUALS(task->pack_residue( i ), false );
		}
	}
}

inline
void
assert_region_design_is_disabled_rr(
	core::pose::Pose const & pose,
	core::pack::task::operation::RestrictResidueToRepackingOP disable,
	protocols::antibody::AntibodyInfoCOP ab_info,
	protocols::antibody::AntibodyRegionEnum region){

	
	core::pack::task::TaskFactoryOP tf( new core::pack::task::TaskFactory() );
	tf->push_back(disable);
	assert_region_design_is_disabled(pose, tf->create_task_and_apply_taskoperations(pose), ab_info, region);
}

inline
void
assert_cdr_design_is_enabled_or_disabled(
	core::pose::Pose const & pose,
	core::pack::task::PackerTaskOP task,
	protocols::antibody::AntibodyInfoCOP ab_info,
	utility::vector1<bool> cdrs_to_check_disabled )
{
	
	//Checks to make sure that the cdrs that are set to be disabled are disabled and that the other cdr residues are ALL enabled. 
	std::cout << "Checking CDR Design" << std::endl;
	assert(cdrs_to_check_disabled.size() == 6);

	for (core::Size i = 1; i <= 6; ++i ){
		protocols::antibody::CDRNameEnum cdr = static_cast<protocols::antibody::CDRNameEnum>( i );
		core::Size start = ab_info->get_CDR_start(cdr, pose);
		core::Size end = ab_info->get_CDR_end(cdr, pose);
		std::cout <<"CDR: " << ab_info->get_CDR_name(cdr) << std::endl;
			
		for (core::Size res = start; res <= end; ++res ){
			if (cdrs_to_check_disabled [ i ]) {
				TS_ASSERT_EQUALS(task->design_residue( res ) , false);
			}
			else {
				TS_ASSERT_EQUALS(task->design_residue( res ), true);
			}	
		}
	}
		
}

inline
void
assert_cdr_design_disabled(
	core::pose::Pose const & pose,
	core::pack::task::PackerTaskOP task,
	protocols::antibody::AntibodyInfoCOP ab_info,
	utility::vector1<bool> cdrs_to_check_disabled)
{
	
	assert(cdrs_to_check_disabled.size() == 6);
	std::cout << "Checking CDR Design" << std::endl;
	for (core::Size i = 1; i <= 6; ++i ){
		protocols::antibody::CDRNameEnum cdr = static_cast<protocols::antibody::CDRNameEnum>( i );
		core::Size start = ab_info->get_CDR_start(cdr, pose);
		core::Size end = ab_info->get_CDR_end(cdr, pose);
		for (core::Size res = start; res <= end; ++res ){
			if (cdrs_to_check_disabled [ i ]) {
				TS_ASSERT_EQUALS(task->design_residue( res ) , false);
			}
		}
	}
}
	
inline	
void
assert_cdr_packing_is_enabled_or_disabled(
	core::pose::Pose const & pose,
	core::pack::task::PackerTaskOP task,
	protocols::antibody::AntibodyInfoCOP ab_info,
	utility::vector1<bool> cdrs_to_check_disabled)
{
	
	assert(cdrs_to_check_disabled.size() == 6);
	std::cout << "Checking packing " << std::endl;
	for (core::Size i = 1; i <= 6; ++i ){
		protocols::antibody::CDRNameEnum cdr = static_cast<protocols::antibody::CDRNameEnum>( i );
		core::Size start = ab_info->get_CDR_start(cdr, pose);
		core::Size end = ab_info->get_CDR_end(cdr, pose);
		for (core::Size res = start; res <= end; ++res ){
			if (cdrs_to_check_disabled [ i ]) {
				TS_ASSERT_EQUALS(task->pack_residue( res ) , false);
			}
			else {
				TS_ASSERT_EQUALS(task->pack_residue( res ), true);
			}	
		}
	}
		
}

#endif

