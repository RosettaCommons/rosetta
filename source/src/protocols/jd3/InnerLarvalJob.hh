// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/InnerLarvalJob.hh
/// @brief  Class definition for the InnerLarvalJob, which describes the main details of a collection
///         of Jobs that differ only in their random number seed.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_jd3_InnerLarvalJob_hh
#define INCLUDED_protocols_jd3_InnerLarvalJob_hh

// Unit headers
#include <protocols/jd3/InnerLarvalJob.fwd.hh>

// Package headers
#include <protocols/jd3/PoseInputSource.fwd.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>

// Basic headers
#include <basic/datacache/ConstDataMap.fwd.hh>
#include <basic/resource_manager/JobOptions.fwd.hh>

//C++ headers
#include <string>

#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace jd3 {

/// @details The %InnerLarvalJob class is responsible for holding the input requirements
/// for a given job - how many nstruct, and what the input is.  %InnerLarvalJobs are
/// relatively heavy; there is no need to duplicate a series of %InnerLarvalJobs for
/// each index into nstruct.  The companion Job class handles the nstruct index
/// and has a pointer to an InnerLarvalJob (which is shared across many Jobs).
/// The JobQueen has considerable leeway in how she shares data held in the %InnerLarvalJobs
/// she creates: if all of the options for a set of %InnerLarvalJob are identical, and
/// they differ only in their input structure, then she can, e.g. share a single
/// JobOptions object between all of them.
/// InnerLarvalJobs are serialized by the JobDistributor and shipped between nodes.  For this
/// reason, they should not be loaded with ResourceManager-managed data (that data should
/// not get shipped between nodes, though it should be instantiatable on any node), but they
/// may need to store Poses, as would be necessary in any multi-round protocol, where the
/// JobResults from round i are the starting points for round i+1.
class InnerLarvalJob : public utility::pointer::ReferenceCount {
public:

	InnerLarvalJob();
	InnerLarvalJob( core::Size nstruct );
	InnerLarvalJob( InnerLarvalJob const & src );

	~InnerLarvalJob() override;

	/// @brief Determine if two jobs are defined with the same set of PoseInputSources
	bool sources_same( InnerLarvalJob const & other ) const;

	/// @brief Mutual comparison of this inner job to the other inner job
	/// so that if either one thinks it's not the same as the other, then
	/// it returns false.  Invokes the same() function on both this and other
	///
	/// @details Note: only compare if the pointers to the poses are to the
	/// same location
	bool
	operator == ( InnerLarvalJob const & other ) const;

	bool
	operator != (InnerLarvalJob const & other ) const;

	/// @brief returns true if this is the same type as other;
	/// does not call other.same_type()
	virtual
	bool
	same_type( InnerLarvalJob const & other ) const;

	virtual
	void
	show( std::ostream & out ) const;

	friend
	std::ostream &
	operator<< ( std::ostream & out, const InnerLarvalJob & inner_job );

	/// @brief Set the single PoseInputSource for this job
	void input_source( PoseInputSourceCOP setting );

	/// @brief Clear the vector of PoseInputSources
	void clear_input_sources();

	/// @brief Append a PoseInputSource to the vector of them
	void append_input_source( PoseInputSourceCOP setting );

	/// @brief Return the input tag -- a string that describes the input structure, but that
	/// is in no way a complete description of the job.  The job_tag instead should be looked
	/// to as the name to use to identify information about this job: a single input structure
	/// may be used for multiple jobs.
	std::string input_tag() const;

	/// @brief Set the input tag -- if it is not set, then the PoseInputSources will be used,
	/// the file names / tags being concatenated together, separated by underscores.
	void input_tag( std::string const & setting );

	/// @brief Return the job tag -- a string that describes what task is being performed and
	/// and (usually) what input structure that job is being performed upon.  This is the tag
	/// can be defined by the user in the job-definition XML file.
	std::string job_tag() const;

	/// @brief Set the job tag.  If not set, then the input tag will be returned.
	void job_tag( std::string const & setting );

	/// @brief Return the outputter class; i.e. the key used by the PoseOutputterFactory
	std::string outputter() const;

	/// @brief Set the outputter class for this job; i.e. the key used by the PoseOutputterFactory
	void outputter( std::string const & setting );

	/// @brief The number of job replicates to be performed for this %InnerLarvalJob.
	core::Size nstruct_max() const;

	/// @brief Read access to the ConstDataMap that the %InnerLarvalJob holds.
	basic::datacache::ConstDataMap const & const_data_map() const;

	utility::tag::TagCOP jobdef_tag() const;

	void nstruct_max( core::Size setting );

	void jobdef_tag( utility::tag::TagCOP setting );

	/// @brief Return the number of pose input sources for this job.
	core::Size n_input_sources() const;

	/// @brief Read access to the PoseInputSource object that describes how the Pose for
	/// this job should be constructed.  Requires that there is only a single PoseInputSource
	/// for this job.
	PoseInputSource const & input_source() const;

	/// @brief Retrieve a particular PoseInputSource
	PoseInputSource const & input_source( Size index ) const;

	/// @brief The list of the JobResults required to mature this %LarvalJob, by global index of the
	/// already-executed (Lavral)Jobs
	utility::vector1< core::Size > const &
	input_job_result_indices() const;

	void add_input_job_result_index( core::Size job_index );


	/// @brief Store the fact that this job is bad, for some reason or another, e.g. it has
	/// malformed inputs.
	void set_bad( bool setting = true );

	/// @brief Has this job been labeled as bad?
	bool bad() const;

	/// @brief Write access to the ConstDataMap that the %InnerLarvalJob holds.  The objects pointed
	/// to by the const data map can be shared between two %InnerLarvalJobs, but the data maps
	/// themselves should not be shared.
	basic::datacache::ConstDataMap & const_data_map();

private:

	/// @brief construct a string from the input_sources_
	std::string concatenated_input_source_names() const;

private:

	std::string input_tag_;
	std::string job_tag_;
	std::string outputter_;

	utility::vector1< PoseInputSourceCOP > input_sources_;

	core::Size nstruct_;

	/// @brief What set of JobResults are required as inputs for the JobQueen to mature larval jobs that point at
	/// this %InnerLarvalJob into a full job?
	utility::vector1< core::Size > input_job_result_indices_;

	bool bad_;

	/// @brief The %InnerLarvalJob holds a const data cache to ensure that if two %InnerLarvalJobs point to the same
	/// piece of data, modification to one of the %InnerLarvalJobs objects cannot disrupt the other %InnerLarvalJob.
	/// Arbitrary data can be cached here by the JobQueen as she initializes the %InnerLarvalJobs.
	basic::datacache::ConstDataMapOP const_data_cache_;
	utility::tag::TagCOP jobdef_tag_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // InnerLarvalJob

} // namespace jd3
} // namespace protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_jd3_InnerLarvalJob )
#endif // SERIALIZATION


#endif //INCLUDED_protocols_jd3_InnerLarvalJob_HH
