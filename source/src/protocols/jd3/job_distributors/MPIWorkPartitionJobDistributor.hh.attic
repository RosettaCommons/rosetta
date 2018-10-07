// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/job_distributors/MPIWorkPartitionJobDistributor.hh
/// @brief  jd3 header for MPIWorkPartitionJobDistributor - intended for MPI jobs on small numbers of nodes where the load can be balanced equally by the user
/// @author P. Douglas Renfrew (renfrew@nyu.edu)
/// @author Andy Watkins (amw579@nyu.edu)

#ifndef INCLUDED_protocols_jd3_job_distributors_MPIWorkPartitionJobDistributor_hh
#define INCLUDED_protocols_jd3_job_distributors_MPIWorkPartitionJobDistributor_hh

// Unit headers
#include <protocols/jd3/job_distributors/MPIWorkPartitionJobDistributor.fwd.hh>

// Package headers
#include <protocols/jd3/JobDistributor.hh>
#include <protocols/jd3/Job.fwd.hh>

#include <protocols/moves/Mover.fwd.hh>

// Utility headers
#include <core/types.hh>

// C++ headers
#include <string>

#include <platform/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/jd3/InnerLarvalJob.fwd.hh>
#include <protocols/jd3/LarvalJob.fwd.hh>
#include <protocols/jd3/Job.hh>
#include <protocols/jd3/JobDistributor.fwd.hh>
#include <protocols/jd3/JobQueen.fwd.hh>
#include <utility/down_cast.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <numeric/xyzVector.fwd.hh>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iosfwd>
#include <limits>
#include <list>
#include <map>
#include <vector>


namespace protocols {
namespace jd3 {
namespace job_distributors {

/// @details This job distributor is meant for running jobs where the number of jobs is equal to the number of processors
///(or, similarly, the jobs % processors calculation is very close to the number of processors and NOT a small number).
///It will blindly divide up jobs across processors and then start running them; it will NOT attempt load-balancing by
///giving more jobs to the processors that finished their original jobs.  This is intended for use on smaller numbers of
///processors, and/or where the jobs are known to be equal in runtime.  (The WorkPool implementation is meant for when
///runtimes are uncertain, or you have many many processors). It does not "waste" a processor as a master node, instead
///all processors run jobs.
class MPIWorkPartitionJobDistributor : public JobDistributor
{
protected:
	virtual void handle_interrupt() {}

public:
	///@brief ctor is protected; singleton pattern
	MPIWorkPartitionJobDistributor();

	virtual ~MPIWorkPartitionJobDistributor();

	virtual
	void
	go( JobQueenOP queen );

	virtual
	core::Size
	get_new_job_id();

	//virtual
	//void
	//mark_current_job_id_for_repetition();

	//virtual
	//void
	//remove_bad_inputs_from_job_list();

	friend class JobDistributorFactory;  //singleton management

private:
	///@brief ctor helper function splits up job list
	void
	determine_job_ids_to_run();

	bool
	more_jobs_in_current_round();

	LarvalJobOP
	select_next_job();

	///@brief total number of processing elements
	core::Size npes_;

	///@brief rank of the "local" instance
	core::Size rank_;

	//@brief start of Jobs vector slice
	core::Size job_id_start_;

	//@brief end of Jobs vector slice
	core::Size job_id_end_;

	core::Size next_job_to_try_assigning_;
};

}//job_distributors
}//jd3
}//protocols

#endif //INCLUDED_protocols_jd2_MPIWorkPartitionJobDistributor_HH
