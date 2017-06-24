// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/MPIWorkPoolJobDistributor.hh
/// @brief  header for MPIWorkPoolJobDistributor - intended for continuous resamplig jobs  that spawn new jobs based on a pool/archive of
///         structures
/// @author Oliver Lange olange@u.washington.edu

#ifndef INCLUDED_protocols_jd2_archive_MPIArchiveJobDistributor_hh
#define INCLUDED_protocols_jd2_archive_MPIArchiveJobDistributor_hh

// Unit headers
#include <protocols/jd2/archive/MPIArchiveJobDistributor.fwd.hh>
#include <protocols/jd2/archive/ArchiveBase.hh>
#include <protocols/jd2/MPIFileBufJobDistributor.hh>
// Package headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.fwd.hh>
#include <protocols/jd2/JobDistributorFactory.fwd.hh>

#include <protocols/moves/Mover.fwd.hh>

// Utility headers
#include <core/types.hh>
#include <utility/exit.hh>

// C++ headers
#include <string>
#include <deque>

#include <platform/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/jd2/InnerJob.fwd.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.fwd.hh>
#include <protocols/jd2/JobInputter.fwd.hh>
#include <protocols/jd2/JobOutputter.fwd.hh>
#include <protocols/jd2/MPIFileBufJobDistributor.fwd.hh>
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
#include <utility/assert.hh>
#include <cstddef>
#include <iosfwd>
#include <limits>
#include <list>
#include <map>
#include <vector>


namespace protocols {
namespace jd2 {
namespace archive {

//Archive has numbers 100+
core::Size const BATCH_SYNC = 101;
core::Size const QUEUE_EMPTY = 102;
core::Size const ADD_BATCH = 103;
core::Size const CANCEL_BATCH = 104;
core::Size const JOB_COMPLETION = 105;

core::Size const MPI_ARCHIVE_TAG = 12310925; //keep unique TAG to communicate between ArchiveManager and ArchiveJobDistributor
///Tags used to tag messeges sent by MPI functions used to decide whether a slave is requesting a new job id or
///flagging as job as being a bad input

/// @brief JobDistributor for the iterative ArchiveManager/Archive Framework
/// @details This job distributor is meant for running iterative jobs with the ArchiveManager/Archive Framework.
///could vary greatly. In this configuration the three first nodes are dedicated processes (JobDistributor, FileBuffer, and ArchiveManger )
///and the remaining CPUs form slave or worker nodes. This JD will not work at all
///without MPI and the implementations of all but the interface functions have been put inside of ifdef directives.
///Generally each function has a master and slave version, and the interface functions call one or the other depending
///on processor rank.

class MPIArchiveJobDistributor : public MPIFileBufJobDistributor
{
public:

	/// @brief CompletionMessage(s) are send to the ArchiveManager whenever more than nr_notify decoys have been finished
	//// or when the full batch is finished.
	struct CompletionMessage {
	public:
		CompletionMessage() : batch_id( 0 ), final( false ), bad( 0 ), good( 0 ), njobs( 0 ), msg_tag( JOB_COMPLETION ) {};
		CompletionMessage( core::Size id, bool fi, core::Size bad_in, core::Size good_in, core::Size total_in )
		: batch_id( id ), final( fi), bad( bad_in ),good( good_in ), njobs( total_in ), msg_tag( JOB_COMPLETION ) {};
		CompletionMessage( core::Size batch_id, core::Size tag )
		: batch_id( batch_id ), final( false ), bad( 0 ), good( 0 ), njobs( 0 ), msg_tag( QUEUE_EMPTY )
		{ runtime_assert( tag == QUEUE_EMPTY ); };
		core::Size batch_id;
		bool final;
		core::Size bad;
		core::Size good;
		core::Size njobs;
		core::Size msg_tag;
	};

protected:
	typedef MPIFileBufJobDistributor Parent;

	/// @brief ctor is protected; singleton pattern
	MPIArchiveJobDistributor();
	friend class protocols::jd2::JobDistributorFactory; //ctor access

	virtual void handle_interrupt() {}

public:

	/// @brief overloaded to also start the ArchiveManager process
	virtual
	void
	go( protocols::moves::MoverOP mover );

	void
	set_archive( archive::ArchiveBaseOP );

	bool is_archive_rank() const {
		return archive_rank() == rank();
	}

protected:
	/// @brief triggered in slave if new batch_ID comes in.
	virtual void batch_underflow();

	/// @brief act on a message, return true if message was understood
	virtual bool process_message(
		core::Size msg_tag,
		core::Size slave_rank,
		core::Size slave_job_id,
		core::Size slave_batch_id,
		core::Real run_time
	);

	/// @brief overloaded to allow statistics and sending of CompletionMessages
	virtual void mark_job_as_completed( core::Size job_id, core::Size batch_id, core::Real run_time );

	/// @brief overloaded to allow statistics and sending of CompletionMessages
	virtual void mark_job_as_bad( core::Size job_id, core::Size batch_id );

	/// @brief overloaded to start new entries in nr_new_completed_, nr_completed_, nstruct_ and nr_bad_ ...
	virtual void load_new_batch();

	/// @brief rank of ArchiveManger process
	core::Size archive_rank() const {
		return archive_rank_;
	}

private:


	//actually transmit a notify msg -- this should only be called out of the process_message method
	void _notify_archive();

	/// @brief receive a new Batch from ArchiveManager
	bool receive_batch( core::Size source_rank );

	/// @brief sync batch queue with slave node
	void sync_batches( core::Size slave_rank );

	/// @brief send message to ArchiveManager
	void master_to_archive( core::Size tag );

	//some statistics about completion for ArchiveManager to query
	utility::vector1< core::Size > nr_jobs_; //nr_jobs per batch
	utility::vector1< core::Size > nr_completed_; // completed per batch
	utility::vector1< core::Size > nr_new_completed_; // completed per batch
	utility::vector1< core::Size > nr_bad_; // bad per batch
	utility::vector1< core::Size > nstruct_;

	//after how many completed decoys should we tell the Archive ?
	core::Size nr_notify_;
	core::Size archive_rank_;

	/// @brief add a notifcation (CompletionMessage) to the msg queue ...
	// these are send out by _notify_archive() at beginning of process_message()
	void notify_archive( CompletionMessage const& );

	/// @brief work out if a notifcation should be send (using above method)
	void notify_archive( core::Size batch_id );

	/// @brief unsent notifications
	std::deque< CompletionMessage > pending_notifications_;


	ArchiveBaseOP theArchive_;
};

} //archive
} //jd2
} //protocols

#endif //INCLUDED_protocols_jd2_MPIArchiveJobDistributor_HH
