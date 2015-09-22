// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/mmt_msd/MMTPackingJob.cc
/// @brief  Implementation for class MMTPackingJob
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <devel/mmt_msd/MMTPackingJob.hh>

// core headers
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/PackerTask.hh>

// Basic headers
#include <basic/Tracer.hh>

#ifdef MULTI_THREADED
#ifdef CXX11
// C++11 headers
#include <chrono>
#endif
#endif

namespace devel {
namespace mmt_msd {

static THREAD_LOCAL basic::Tracer TR( "devel.mmt_msd.MMTPackingJob" );

MMTPackingJob::MMTPackingJob() : running_time_( -1.0 ) {}
MMTPackingJob::~MMTPackingJob() {}

void MMTPackingJob::set_pose( core::pose::Pose const & pose )
{
	pose_ = pose.clone();
}

void MMTPackingJob::set_sfxn( core::scoring::ScoreFunction const & sfxn )
{
	sfxn_ = sfxn.clone();
}

void MMTPackingJob::set_packer_task( core::pack::task::PackerTask const & task )
{
	task_ = task.clone();
}

void MMTPackingJob::set_npd_properties( required_npds const & npd_to_do_list )
{
	npd_to_do_list_ = npd_to_do_list;
}


core::pose::Pose const &
MMTPackingJob::get_pose() const { return *pose_; }

core::scoring::ScoreFunction const &
MMTPackingJob::get_sfxn() const { return *sfxn_; }

core::pack::task::PackerTaskCOP
MMTPackingJob::get_task() const { return task_; }

bool MMTPackingJob::has_pose() const { return pose_ != 0; }
bool MMTPackingJob::has_sfxn() const { return sfxn_ != 0; }
bool MMTPackingJob::has_task() const { return task_ != 0; }

void
MMTPackingJob::go()
{
	TR << "Starting packing job" << std::endl;

	core::Real running_time; // local copy
#ifdef MULTI_THREADED
#ifdef CXX11
	auto starttime = std::chrono::system_clock::now();
#endif
#else
	clock_t starttime = clock();
#endif

	setup();
	optimize();
	compute_npd_properties();


#ifdef MULTI_THREADED
#ifdef CXX11
	auto stoptime = std::chrono::system_clock::now();
	running_time = std::chrono::duration_cast< std::chrono::seconds >( stoptime - starttime ).count();
#endif
#else
	clock_t stoptime = clock();
	running_time = ((double) stoptime - starttime ) / CLOCKS_PER_SEC;
#endif

	clean_up();

	// this must be the last statement in this function to avoid race conditions
	// since it is the signal that the MMTReceiver will use to determine if the
	// job has finished.
	running_time_ = running_time;

}

core::Real
MMTPackingJob::running_time() const {
	return running_time_;
}

bool MMTPackingJob::optimization_complete() const { return running_time_ >= 0; }

core::Size MMTPackingJob::n_npd_properties() const
{
	return npd_to_do_list_.size();
}

MMTPackingJob::npd_properties::const_iterator
MMTPackingJob::npd_properties_begin() const
{
	return computed_npd_properties_.begin();
}

MMTPackingJob::npd_properties::const_iterator
MMTPackingJob::npd_properties_end() const
{
	return computed_npd_properties_.end();
}


core::pose::Pose & MMTPackingJob::pose() { return *pose_; }
core::scoring::ScoreFunction & MMTPackingJob::sfxn() { return *sfxn_; }
core::pack::task::PackerTaskOP MMTPackingJob::task() { return task_; }

void MMTPackingJob::compute_npd_properties()
{
	if ( npd_to_do_list_.empty() ) return;

	update_pose( *pose_ );
	// ok -- code will be inserted here in the near future
	// to create the npd evaluators and to invoke them

	// junk code temporarily needed to avoid any MPI_Receive calls from hanging!
	for ( required_npds::const_iterator iter = npd_to_do_list_.begin(), iter_end = npd_to_do_list_.end();
			iter != iter_end; ++iter ) {
		computed_npd_properties_.push_back( std::make_pair( iter->first, -1234.5 ) );
	}
}


void MMTPackingJob::clean_up() {
	pose_.reset();
	sfxn_.reset();
	task_.reset();
}

}
}
