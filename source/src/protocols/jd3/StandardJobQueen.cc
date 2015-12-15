// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/StandardJobQueen.cc
/// @brief  StandardJobQueen class's method definitions
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

//unit headers
#include <protocols/jd3/StandardJobQueen.hh>

// package headers
#include <protocols/jd3/InnerLarvalJob.hh>
#include <protocols/jd3/LarvalJob.hh>
#include <protocols/jd3/MoverAndPoseJob.hh>
#include <protocols/jd3/pose_inputters/PDBPoseInputter.hh>
#include <protocols/jd3/pose_outputters/PDBPoseOutputter.hh>

//project headers
#include <core/pose/Pose.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

//basic headers
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

namespace protocols {
namespace jd3 {


StandardJobQueen::StandardJobQueen() :
	pose_inputter_( new pose_inputters::PDBPoseInputter ), // TEMP!  Ask for PoseInputter from the PoseInputterFactory.
	pose_outputter_( new pose_outputters::PDBPoseOutputter ) // TEMP!  Ask for PoseOutputter from the PoseOutputterFactory.
{}

StandardJobQueen::~StandardJobQueen() {}

std::string StandardJobQueen::job_definition_xsd() const {
	return "";
}

InnerLarvalJobs
StandardJobQueen::prepare_preliminary_job_list()
{
	// temp
	PoseInputSources input_sources = pose_inputter_->initialize_pose_input_sources();
	InnerLarvalJobs inner_jobs;
	inner_jobs.reserve( input_sources.size() );
	for ( core::Size ii = 1; ii <= input_sources.size(); ++ii ) {
		InnerLarvalJobOP inner_job( new InnerLarvalJob );
		inner_job->input_source( *input_sources[ii] );
		inner_jobs.push_back( inner_job );
	}
	return inner_jobs;
}

LarvalJobs
StandardJobQueen::expand_job_list( InnerLarvalJobs const & inner_jobs ) const {
	// temp
	core::Size nstruct = basic::options::option[ basic::options::OptionKeys::out::nstruct ];
	LarvalJobs jobs;
	jobs.reserve( inner_jobs.size() * nstruct );
	for ( core::Size ii = 1; ii <= inner_jobs.size(); ++ii ) {
		for ( core::Size jj = 1; jj <= nstruct; ++jj ) {
			LarvalJobOP job = create_job( inner_jobs[ii], jj );
			jobs.push_back( job );
		}
	}
	return jobs;
}

LarvalJobOP
StandardJobQueen::create_job( InnerLarvalJobOP job, core::Size nstruct_index ) const
{
	return LarvalJobOP( new LarvalJob( job, nstruct_index ));
}

void StandardJobQueen::add_option( utility::options::BooleanOptionKey const & /*key*/ ) {}
void StandardJobQueen::add_option( utility::options::BooleanVectorOptionKey const & /*key*/ ) {}
void StandardJobQueen::add_option( utility::options::FileOptionKey const & /*key*/ ) {}
void StandardJobQueen::add_option( utility::options::FileVectorOptionKey const & /*key*/ ) {}
void StandardJobQueen::add_option( utility::options::IntegerOptionKey const & /*key*/ ) {}
void StandardJobQueen::add_option( utility::options::IntegerVectorOptionKey const & /*key*/ ) {}
void StandardJobQueen::add_option( utility::options::PathOptionKey const & /*key*/ ) {}
void StandardJobQueen::add_option( utility::options::PathVectorOptionKey const & /*key*/ ) {}
void StandardJobQueen::add_option( utility::options::RealOptionKey const & /*key*/ ) {}
void StandardJobQueen::add_option( utility::options::RealVectorOptionKey const & /*key*/ ) {}
void StandardJobQueen::add_option( utility::options::StringOptionKey const & /*key*/ ) {}
void StandardJobQueen::add_option( utility::options::StringVectorOptionKey const & /*key*/ ) {}
void StandardJobQueen::remove_default_input_element() {}

//void StandardJobQueen::add_all_sub_element(
// std::string const & nesting_level,
// std::string const & element_name
//)
//{}
//
//void StandardJobQueen::add_choice_sub_element(
// std::string const & nesting_level,
// core::Size choice_block_index,
// std::string const & element_name
//)
//{}
//
//void StandardJobQueen::add_attribute(
// std::string const & nested_element_name,
// std::string const & attribute_nanme,
// XMLAttributeType attribute_type
//)
//{}
//
//void StandardJobQueen::add_attribute(
// std::string const & nested_element_name,
// std::string const & attribute_nanme,
// XMLAttributeType attribute_type,
// std::string const & default_value
//)
//{}
//
//void StandardJobQueen::add_managed_resource_attribute(
// std::string const & nested_element_name,
// std::string const & attribute_name
//)
//{}

core::pose::PoseOP
StandardJobQueen::pose_for_job( LarvalJobCOP job )
{
	// either read the Pose in using the pose_inputter (and then keep a copy
	// in the resource manager), or retrieve the Pose from the resource manager
	// initial version: just read the pose in for each job.
	return pose_inputter_->pose_from_input_source( job->inner_job()->input_source() );
}

//ResourceManagerOP StandardJobQueen::resource_manager()
//{}

/// @brief Access the pose inputter
PoseInputter &
StandardJobQueen::pose_inputter()
{
	return *pose_inputter_;
}

/// @brief Access the pose outputter
PoseOutputter &
StandardJobQueen::pose_outputter()
{
	return *pose_outputter_;
}


} // namespace jd3
} // namespace protocols
