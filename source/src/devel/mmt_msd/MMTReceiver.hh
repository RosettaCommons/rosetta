// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   devel/mmt_msd/MMTReceiver.hh
/// @brief  declaration for class MMTReceiver
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_devel_mmt_msd_MMTReceiver_HH
#define INCLUDED_devel_mmt_msd_MMTReceiver_HH

// Unit headers
#include <devel/mmt_msd/MMTReceiver.fwd.hh>

#include <devel/mmt_msd/MMTPackingJob.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/ResfileReader.fwd.hh>

#include <protocols/pack_daemon/EntityCorrespondence.fwd.hh>

// C++ headers
#include <algorithm>
#include <list>
#include <map>
#include <string>

namespace devel {
namespace mmt_msd {

enum mmt_message
{
	handshake_begin = 1,
	handshake_acknowledged,
	handshake_failure,
	new_generation,
	requesting_new_job,
	test_state_inputs,
	new_job_ready,
	no_ready_job,
	job_complete,
	generation_complete,
	result_from_last_generation_needs_saving,
	old_result_can_be_discarded,
	recover_pose_for_result,
	recovery_successful,
	error,
	spin_down
};

struct StateInputData {
	typedef std::list< std::pair< core::Size, std::string > > required_npds;

	core::Size  state_index;
	std::string pdb_name;
	std::string pdb_string;
	std::string corr_fname;
	std::string corr_file;
	std::string sec_resfname;
	std::string sec_resfile;
	required_npds npd_to_do_list;
};

struct StateData {
	StateData();
	~StateData();

	core::pose::PoseOP pose;
	protocols::pack_daemon::EntityCorrespondenceOP ent_corr;
	core::pack::task::ResfileContentsOP secondary_resfile_contents;
	core::pack::task::PackerTaskOP task;
};

class MMTReceiver : public utility::pointer::ReferenceCount
{
public:

	typedef std::pair< core::Size, std::string > StateAndSequencePair;
	typedef std::map< StateAndSequencePair, MMTPackingJobOP > JobMap;
	typedef std::pair< StateAndSequencePair, MMTPackingJobOP > RunningJob;
	typedef std::list< RunningJob > RunningJobsList;

public:

	MMTReceiver();
	~MMTReceiver() override;

	void set_max_capacity( core::Size nthreads_max );

	bool initial_handshake();
	void main_optimization_loop();

private:

	void initialize_entity_resfile();

	bool start_new_generation();
	mmt_message request_new_job();
	void send_job_result_to_node0( RunningJob const & job );

	bool save_result_from_last_generation();
	bool discard_old_result();
	bool recreate_previously_generated_result_pose();

	/// @brief put the receiver thread to sleep for a little while
	/// and for an increasing amount of time if it repeatedly
	/// gets through the start_new_generation loop without anything
	/// having changed
	void sleep_a_bit();

	/// @brief Reset the amount of time the master thread is supposed to sleep
	/// back to it's minimum
	void reset_sleep_progression();

	StateInputData receive_state_input_data() const;
	StateData initialize_state_data( StateInputData const & sid );

	void receive_new_job();

	protocols::pack_daemon::EntityCorrespondenceOP
	entity_correspondence_from_string(
		std::string const & corr_file,
		core::pose::PoseOP pose) const;

	core::pack::task::PackerTaskOP
	initialize_packer_task(
		core::pose::PoseOP pose,
		protocols::pack_daemon::EntityCorrespondenceOP ent_corr,
		core::pack::task::ResfileContents const & secondary_resfile_contents
	) const;

	void restrict_task_for_job(
		core::pose::PoseOP pose,
		protocols::pack_daemon::EntityCorrespondenceOP ent_corr,
		std::string const & seqstring,
		core::pack::task::PackerTaskOP task
	) const;

	void
	delete_running_job( RunningJobsList::iterator & iter );


private:

	RunningJobsList running_jobs_;
	JobMap curr_gen_jobs_;
	JobMap best_jobs_;

	core::Size my_mpi_rank_;

	core::Size curr_njobs_running_;
	core::Size max_capacity_;

	core::scoring::ScoreFunctionOP sfxn_;

	core::pack::task::PackerTaskOP entity_task_;
	core::pack::task::ResfileContentsOP entity_resfile_contents_;
	core::Size num_entities_;

	core::Size const sleep_min_;
	core::Size const sleep_max_;
	core::Size sleep_nsecs_;
};

}
}

#endif
