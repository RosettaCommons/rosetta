// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/jd3/standard/PreliminaryLarvalJobTracker.hh
/// @brief A class that tracks PreliminaryLarvalJobs in the SJQ.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_jd3_standard_PreliminaryLarvalJobTracker_hh
#define INCLUDED_protocols_jd3_standard_PreliminaryLarvalJobTracker_hh

#include <protocols/jd3/standard/PreliminaryLarvalJobTracker.fwd.hh>
#include <protocols/jd3/standard/PreliminaryLarvalJob.fwd.hh>

#include <protocols/jd3/LarvalJob.fwd.hh>

#include <core/types.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <map>
#include <set>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace jd3 {
namespace standard {

///@brief A class for tracking the progress of all PreliminaryLarvalJobs.
/// This tracking is partially responsible for deallocation of the input poses.
class PreliminaryLarvalJobTracker : public utility::pointer::ReferenceCount {

public:

	PreliminaryLarvalJobTracker();

public:

	///@brief Initialize the PreliminaryLarvalJobTracker's counters.
	void
	initialize_tracker( utility::vector1< PreliminaryLarvalJob > const & prelim_jobs );

	///@brief Track a completed LarvalJob.
	/// If the LarvalJob is NOT a preliminary job, we skip it.
	void
	track_job_completed( LarvalJobCOP job );

	///@brief Called when we have given out ALL Larval Jobs for a particular PJN.
	void
	track_job_node_assigned( core::Size job_dag_node_index, core::Size last_global_nstruct_index_for_node);

	///@brief Called during assignment of the PJN.
	void
	track_job_node_being_assigned( core::Size job_dag_node_index, core::Size global_index_start, core::Size global_index_end);

public:

	///@brief Are ALL jobs for this job_node_index complete?
	bool
	get_job_node_complete( core::Size job_dag_node_index ) const;

	///@brief Are jobs for this job_node_index assigned (IE - LarvalJobs created)?
	bool
	get_job_node_assigned( core::Size job_dag_node_index ) const;

	///@brief Get the job dag node for the global job index if tracked through
	///  track_job_node_being_assigned during assignment of the PJN
	///
	/// Returns 0 if the job_id is NOT a preliminary job.
	core::Size
	get_job_node_for_job_index( core::Size job_id ) const;



public:

	///@brief Get the starting global job index that starts for a particular preliminary job node
	core::Size
	get_job_index_starting_job_node(core::Size job_node);

	///@brief Get the ending global job index that starts for a particular preliminary job node
	core::Size
	get_job_index_ending_job_node(core::Size job_node);

public:

	///@brief Get a list of input poses that the JQ has not deallocated and that are no longer needed by PJNs.
	utility::vector1< core::Size >
	get_input_poses_to_deallocate() const;

	//Deallocation functions
	///@details Is the input pose still needed by any Preliminary Jobs or Preliminary Job Nodes?
	bool
	input_pose_no_longer_needed( core::Size pose_index ) const;

	///@brief Track a deallocated input pose. (indicate that this pose id has been deallocated)
	void
	track_deallocated_input_pose( core::Size pose_id );

public:
	//(used for unit testing)

	utility::vector1< core::Size >
	get_preliminary_job_node_indices() const;

private:

	///@bref All Preliminary Job Node indices.
	utility::vector1< core::Size > preliminary_job_node_inds_;

	// If the DJQ uses the SJQ's version of next_batch_of_larval_jobs_from_prelim in its
	// determine_job_list method when handling those jobs from prelimary job nodes (e.g.
	// as would happen automatically if the DJQ does not override determine_job_list, but
	// rather only overrides the next_batch_of_larval_jobs_for_job_node method -- this is
	// recommended!), then the SJQ will keep track of the starting and ending job indices
	// for each preliminary job node, and will thus be able to decide when to deallocate
	// the input Poses that will no longer be needed.

	///@brief The vector of global job indexes that start for a particular preliminary job node
	utility::vector1< core::Size > pjn_job_ind_begin_;

	///@brief The vector of global job indexes that end for a particular preliminary job node.
	utility::vector1< core::Size > pjn_job_ind_end_;

	///@brief Tracker of which preliminary job nodes have been assigned.
	utility::vector1< bool > preliminary_job_nodes_assigned_;

	///@brief The nstruct left for each preliminary job node.
	utility::vector1< core::Size > outstanding_job_count_;

	///@brief The max nstruct for each preliminary job node.
	utility::vector1< core::Size > pjn_max_nstruct_;

	///@brief A map of pjns that use the same input pose that are not yet complete.
	/// Input poses are indexed by size.  This is updated on completed_job.
	std::map< core::Size, std::set< core::Size > > incomplete_pjns_using_input_pose_;

	///A list of poses that are not deallocated.  These may be complete and not needed, but they have not been deallocated
	utility::vector1< core::Size > non_deallocated_input_poses_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} //standard
} //jd3
} //protocols



#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_jd3_standard_PreliminaryLarvalJobTracker )
#endif // SERIALIZATION


#endif //INCLUDED_protocols_jd3_standard_PreliminaryLarvalJobTracker_hh





