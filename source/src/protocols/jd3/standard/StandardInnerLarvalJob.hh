// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/standard/StandardInnerLarvalJob.hh
/// @brief  Class definition for the StandardInnerLarvalJob, which holds the
///         index of the preliminary job node from which it was derived,
///         in addition to the other data held by the InnerLarvalJob base
///         class.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_jd3_standard_StandardInnerLarvalJob_hh
#define INCLUDED_protocols_jd3_standard_StandardInnerLarvalJob_hh

// Unit headers
#include <protocols/jd3/standard/StandardInnerLarvalJob.fwd.hh>

// Package headers
#include <protocols/jd3/InnerLarvalJob.hh>
#include <protocols/jd3/strong_types.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace jd3 {
namespace standard {

/// @details The %StandardInnerLarvalJob class is responsible for holding the input requirements
/// for a given job - how many nstruct, and what the input is.  %StandardInnerLarvalJobs are
/// relatively heavy; there is no need to duplicate a series of %StandardInnerLarvalJobs for
/// each index into nstruct.  The companion Job class handles the nstruct index
/// and has a pointer to an StandardInnerLarvalJob (which is shared across many Jobs).
/// The JobQueen has considerable leeway in how she shares data held in the %StandardInnerLarvalJobs
/// she creates: if all of the options for a set of %StandardInnerLarvalJob are identical, and
/// they differ only in their input structure, then she can, e.g. share a single
/// JobOptions object between all of them.
/// StandardInnerLarvalJobs are serialized by the JobDistributor and shipped between nodes.  For this
/// reason, they should not be loaded with ResourceManager-managed data (that data should
/// not get shipped between nodes, though it should be instantiatable on any node), but they
/// may need to store Poses, as would be necessary in any multi-round protocol, where the
/// JobResults from round i are the starting points for round i+1.
class StandardInnerLarvalJob : public jd3::InnerLarvalJob {
public:
	typedef jd3::InnerLarvalJob parent;

public:

	StandardInnerLarvalJob();
	StandardInnerLarvalJob( core::Size nstruct, JobDAGNodeID job_node, PrelimJobNodeID preliminary_job_node );
	StandardInnerLarvalJob( StandardInnerLarvalJob const & src );

	~StandardInnerLarvalJob() override;

public:


	/// @brief Mutual comparison of this inner job to the other inner job
	/// so that if either one thinks it's not the same as the other, then
	/// it returns false.  Invokes the same_type() function on both this and other
	bool
	operator == ( InnerLarvalJob const & other ) const override;

	/// @brief returns true if this is the same type as other;
	/// does not call other.same_type()
	bool
	same_type( InnerLarvalJob const & other ) const override;

	void
	show( std::ostream & out ) const override;

	friend
	std::ostream &
	operator<< ( std::ostream & out, const StandardInnerLarvalJob & inner_job );

public:

	///@brief Return the job node this StandardInnerLarvalJob belongs to
	/// If the Job Node has not been set, we return 0.
	PrelimJobNodeID
	preliminary_job_node() const;

private:

	PrelimJobNodeID preliminary_job_node_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // StandardInnerLarvalJob

} // namespace standard
} // namespace jd3
} // namespace protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_jd3_standard_StandardInnerLarvalJob )
#endif // SERIALIZATION


#endif //INCLUDED_protocols_jd3_standard_StandardInnerLarvalJob_HH
