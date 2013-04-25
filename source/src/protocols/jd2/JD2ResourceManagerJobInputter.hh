// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/resource_manager/planner/JD2ResourceManagerJobInputter.hh
/// @brief
/// @author

#ifndef INCLUDED_protocols_jd2_JD2ResourceManagerJobInputter_hh
#define INCLUDED_protocols_jd2_JD2ResourceManagerJobInputter_hh

// Package headers
#include <protocols/jd2/JobInputter.hh>

// Basic headers
#include <basic/options/keys/OptionKeys.hh>
#include <basic/resource_manager/JobOptions.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>

//C++ headers
#include <istream>
#include <string>

namespace protocols {
namespace jd2 {


class JD2ResourceManagerJobInputter : public JobInputter
{
public:
	JD2ResourceManagerJobInputter();
	virtual ~JD2ResourceManagerJobInputter();

	///@brief this function is responsible for filling the pose reference with the pose indicated by the job.  The Job object (within its InnerJob) contains a PoseCOP.  This function needs to either fill the pose reference from the InnerJob or, on first demand of a pose from that InnerJob, instantiate the pose, hand off a COP to the InnerJob, and fill the reference.
 	virtual void pose_from_job( core::pose::Pose & pose, JobOP job );

	///@brief this function determines what jobs exist.  This function neither knows nor cares what jobs are already complete on disk/memory - it just figures out what ones should exist given the input.  NOTE: your JobInputter should order Job objects in the Jobs vector to have as few "transitions" between inputs as possible (group all Jobs of the same input next to each other).  This improves efficiency of the "FAIL_BAD_INPUT" functionality.  Note I said "should", not "must".
	virtual void fill_jobs( Jobs & jobs );

	/// @brief return the type of input source that the JobInputter is currently
	///  using
	virtual JobInputterInputSource::Enum input_source() const;

	void
	fill_jobs_from_stream(
		std::istream & instream,
		Jobs & jobs
	);

private:

	/// @brief Delete the Resources associated with a job
	virtual
	void
	cleanup_after_job_completion(
		std::string const & job_tag);

	void
	parse_jobs_tags(
		utility::tag::TagPtr tag,
		Jobs & jobs
	);

	void
	parse_job_tag(
		utility::tag::TagPtr jobs_tags,
		std::map< std::string, std::string > const & generic_resources_for_job,
		basic::resource_manager::JobOptions const & generic_job_options,
		Jobs & jobs
	);

	void
	parse_jobs_table_tag(
		utility::tag::TagPtr tag,
		std::map< std::string, std::string > const & generic_resources_for_job,
		basic::resource_manager::JobOptions const & generic_job_options,
		Jobs & jobs
	);

	void
	record_job(
		std::string const & job_name,
		std::map< std::string, std::string > const & resources_for_job,
		basic::resource_manager::JobOptionsOP job_options,
		Jobs & jobs
	);

	void
	read_Option_subtag_for_job(
		utility::tag::TagPtr options_tag,
		basic::resource_manager::JobOptionsOP job_options
	);

	void
	parse_options_name_and_value(
		std::string const & optname,
		std::string const & value,
		basic::resource_manager::JobOptionsOP job_options
	);

	void
	read_BooleanOption_subtag_for_job(
		basic::options::BooleanOptionKey const & boolopt,
		std::string const & optname,
		std::string const & val,
		basic::resource_manager::JobOptionsOP job_options
	);

	void
	read_FileOption_subtag_for_job(
		basic::options::FileOptionKey const & fileopt,
		std::string const & optname,
		std::string const & val,
		basic::resource_manager::JobOptionsOP job_options
	);

	void
	read_IntegerOption_subtag_for_job(
		basic::options::IntegerOptionKey const & intopt,
		std::string const & optname,
		std::string const & val,
		basic::resource_manager::JobOptionsOP job_options
	);

	void
	read_PathOption_subtag_for_job(
		basic::options::PathOptionKey const & pathopt,
		std::string const & optname,
		std::string const & val,
		basic::resource_manager::JobOptionsOP job_options
	);

	void
	read_RealOption_subtag_for_job(
		basic::options::RealOptionKey const & realopt,
		std::string const & optname,
		std::string const & val,
		basic::resource_manager::JobOptionsOP job_options
	);

	void
	read_StringOption_subtag_for_job(
		basic::options::StringOptionKey const & stringopt,
		std::string const & optname,
		std::string const & val,
		basic::resource_manager::JobOptionsOP job_options
	);

	void
	read_BooleanVectorOption_subtag_for_job(
		basic::options::BooleanVectorOptionKey const & boolvectopt,
		std::string const & optname,
		std::string const & val,
		utility::vector1< std::string > const & vals,
		basic::resource_manager::JobOptionsOP job_options
	);

	void
	read_FileVectorOption_subtag_for_job(
		basic::options::FileVectorOptionKey const & filevectopt,
		std::string const & optname,
		std::string const & val,
		utility::vector1< std::string > const & vals,
		basic::resource_manager::JobOptionsOP job_options
	);

	void
	read_IntegerVectorOption_subtag_for_job(
		basic::options::IntegerVectorOptionKey const & intvectopt,
		std::string const & optname,
		std::string const & val,
		utility::vector1< std::string > const & vals,
		basic::resource_manager::JobOptionsOP job_options
	);

	void
	read_PathVectorOption_subtag_for_job(
		basic::options::PathVectorOptionKey const & pathvectopt,
		std::string const & optname,
		std::string const & val,
		utility::vector1< std::string > const & vals,
		basic::resource_manager::JobOptionsOP job_options
	);

	void
	read_RealVectorOption_subtag_for_job(
		basic::options::RealVectorOptionKey const & realvectopt,
		std::string const & optname,
		std::string const & val,
		utility::vector1< std::string > const & vals,
		basic::resource_manager::JobOptionsOP job_options
	);

	void
	read_StringVectorOption_subtag_for_job(
		basic::options::StringVectorOptionKey const & strinvectopt,
		std::string const & optname,
		std::string const & val,
		utility::vector1< std::string > const & vals,
		basic::resource_manager::JobOptionsOP job_options
	);


	void
	read_Data_for_subtag(
		utility::tag::TagPtr options_tag,
		std::string const & jobname,
		std::string & input_tag,
		std::map< std::string, std::string > & resources_for_job
	);

	void
	read_ResidueType_for_subtag(
		utility::tag::TagPtr options_tag,
		std::map< std::string, std::string > & resources_for_job
	);

	void
	check_each_job_has_startstruct(
		Jobs const & jobs
	) const;

private:
	///@brief save the last input tag so the resources loaded for it can
	///be unloaded when we see a new input tag
	std::string last_input_tag_;

};

} // namespace jd2
} // namespace protocols

#endif
