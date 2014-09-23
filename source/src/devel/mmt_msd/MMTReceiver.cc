// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/mmt_msd/MMTReceiver.cc
/// @brief  Implementation for class MMTReceiver
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <devel/mmt_msd/MMTReceiver.hh>

// Package headers
#include <devel/mmt_msd/MMTPackingJob.hh>
#include <devel/mmt_msd/MMTMinPackingJob.hh>
#include <devel/mmt_msd/MMTOffRotamerPackingJob.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/import_pose/import_pose.hh>
#include <core/import_pose/import_pose_options.hh>

// protocols headers
#include <protocols/pack_daemon/EntityCorrespondence.hh>
#include <protocols/pack_daemon/util.hh>

// Utility headers
#include <utility/mpi_util.hh>
#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>

// Basic headers
#include <basic/Tracer.hh>

// C++ headers
#include <sstream>

#ifdef MULTI_THREADED
#ifdef CXX11

// C++11 headers
#include <chrono>
#include <thread>

#endif
#endif


namespace devel {
namespace mmt_msd {

static thread_local basic::Tracer TR( "devel.mmt_msd.MMTReceiver" );

StateData::StateData() {}
StateData::~StateData() {}

MMTReceiver::MMTReceiver() :
	my_mpi_rank_( 0 ),
	curr_njobs_running_( 0 ),
	max_capacity_( 0 ),
	sleep_min_( 10 ),
	sleep_max_( 2000 ),
	sleep_nsecs_( sleep_min_ )
{
	sfxn_ = core::scoring::get_score_function();
	my_mpi_rank_ = utility::mpi_rank();
}

MMTReceiver::~MMTReceiver() {}

void MMTReceiver::set_max_capacity( core::Size nthreads_max )
{
	max_capacity_ = nthreads_max;
}

bool
MMTReceiver::initial_handshake()
{
	core::Size handshake = utility::receive_integer_from_node( 0 );
	if ( handshake == handshake_begin ) {
		utility::send_integer_to_node( 0, handshake_acknowledged );
		utility::send_integer_to_node( 0, max_capacity_ );
		initialize_entity_resfile();
		return true;
	} else {
		utility::send_integer_to_node( 0, handshake_failure );
		return false;
	}
	TR << "Handshake successful" << std::endl;
}

void
MMTReceiver::main_optimization_loop()
{
	bool keep_going = true;
	while( keep_going ) {
		core::Size message = utility::receive_integer_from_node( 0 );
		if ( message == new_generation ) {
			keep_going = start_new_generation();
		} else if ( message == result_from_last_generation_needs_saving ) {
			keep_going = save_result_from_last_generation();
		} else if ( message == old_result_can_be_discarded ) {
			keep_going = discard_old_result();
		} else if ( message == recover_pose_for_result ) {
			keep_going = recreate_previously_generated_result_pose();
		} else if ( message == spin_down ) {
			keep_going = false;
		}
	}
}

void
MMTReceiver::initialize_entity_resfile()
{
	std::string entity_resfile_name = utility::receive_string_from_node( 0 );
	std::string entity_resfile = utility::receive_string_from_node( 0 );
	std::istringstream entity_resfile_stream( entity_resfile );

	protocols::pack_daemon::create_entity_resfile_contents(
		entity_resfile_stream, entity_resfile_name,
		entity_resfile_contents_, entity_task_, num_entities_ );

}

bool
MMTReceiver::start_new_generation()
{
	//TR << "Start new generation" << std::endl;

	curr_gen_jobs_.clear();

	bool generation_is_complete = false;
	while( true ) {
		bool something_happened = false;
		if ( ! generation_is_complete && curr_njobs_running_ < max_capacity_ ) {
			//TR << "Requesting new job" << std::endl;
			core::Size message = request_new_job();
			if ( message == new_job_ready ) {
				receive_new_job();
				//TR << "Received new job" << std::endl;
				something_happened = true;
			} else if ( message == generation_complete ) {
				// just wait until the jobs on this node have finished running
				generation_is_complete = true;
				//TR << "Received generation complete signal" << std::endl;
			}
		}

		if ( generation_is_complete && running_jobs_.empty() ) {
			return true;
		}

		for ( RunningJobsList::iterator iter = running_jobs_.begin(); iter != running_jobs_.end(); /* no increment */ ) {
			if ( iter->second->optimization_complete() ) {
				//TR << "Sending completed job to node 0" << std::endl;
				send_job_result_to_node0( *iter );
				delete_running_job( iter );
				something_happened = true;
			} else {
				++iter;
			}
		}

		if ( ! something_happened ) { sleep_a_bit(); }
		else { reset_sleep_progression(); }
	}
}

mmt_message
MMTReceiver::request_new_job()
{
	utility::send_integer_to_node( 0, my_mpi_rank_ );
	utility::send_integer_to_node( 0, requesting_new_job );
	return mmt_message( utility::receive_integer_from_node( 0 ) );
}

void MMTReceiver::send_job_result_to_node0( RunningJob const & job )
{

	utility::send_integer_to_node( 0, my_mpi_rank_ );
	utility::send_integer_to_node( 0, job_complete );
	utility::send_integer_to_node( 0, job.first.first );
	utility::send_string_to_node( 0, job.first.second );
	utility::send_double_to_node( 0, job.second->final_energy() );
	utility::send_double_to_node( 0, job.second->running_time() );
	utility::send_integer_to_node( 0, job.second->n_npd_properties() );
	for ( MMTPackingJob::npd_properties::const_iterator
			iter = job.second->npd_properties_begin(),
			iter_end = job.second->npd_properties_end();
			iter != iter_end; ++iter ) {
		utility::send_integer_to_node( 0, iter->first );
		utility::send_double_to_node(  0, iter->second );
	}
}

bool MMTReceiver::save_result_from_last_generation()
{
	core::Size  state_index  = utility::receive_integer_from_node( 0 );
	std::string seqstring    = utility::receive_string_from_node( 0 );
	StateAndSequencePair job_ssp = std::make_pair( state_index, seqstring );
	MMTPackingJobOP job = curr_gen_jobs_[ job_ssp ];
	best_jobs_[ job_ssp ] = job;
	return true;
}

bool MMTReceiver::discard_old_result()
{
	core::Size  state_index  = utility::receive_integer_from_node( 0 );
	std::string seqstring    = utility::receive_string_from_node( 0 );
	StateAndSequencePair job_ssp = std::make_pair( state_index, seqstring );
	best_jobs_[ job_ssp ] = 0; // set the job pointer to 0; leave the element in the map
	return true;
}

bool MMTReceiver::recreate_previously_generated_result_pose()
{
	StateInputData sid = receive_state_input_data();
	std::string seqstring    = utility::receive_string_from_node( 0 );

	StateData sd = initialize_state_data( sid );
	restrict_task_for_job( sd.pose, sd.ent_corr, seqstring, sd.task );

	StateAndSequencePair job_ssp = std::make_pair( sid.state_index, seqstring );

	MMTPackingJobOP job = best_jobs_[ job_ssp ];

	if ( job == 0 ) {
		utility::send_integer_to_node( 0, error );
		std::string emessage = "Could not find result for previously generated pose. Node " + utility::to_string( utility::mpi_rank() ) + "\n";
		utility::send_string_to_node( 0, emessage );
		return false;
	}

	// reinitialize the job data
	job->set_pose( *sd.pose );
	job->set_sfxn( *sfxn_ );
	job->set_packer_task( *sd.task );
	job->setup();

	try {
		job->update_pose( *sd.pose );
	} catch ( utility::excn::EXCN_Msg_Exception e ) {
		utility::send_integer_to_node( 0, error );
		std::string emessage = "Could not recover previously generated pose. Node " + utility::to_string( utility::mpi_rank() ) + "\n";
		utility::send_string_to_node( 0, emessage + e.msg() );
		return false;
	}

	std::ostringstream pdb_stream;
	sd.pose->dump_pdb( pdb_stream );
	std::string pdb_string = pdb_stream.str();
	utility::send_integer_to_node( 0, recovery_successful );
	utility::send_string_to_node( 0, pdb_string );

	return true;
}

void MMTReceiver::sleep_a_bit()
{
#ifdef MULTI_THREADED
#ifdef CXX11

	if ( sleep_nsecs_ < sleep_max_ ) {
		sleep_nsecs_ *= 2;
	}

	std::chrono::nanoseconds dur( sleep_nsecs_ );
	std::this_thread::sleep_for( dur );
#endif
#else
	/// KAB - casting variable to void to avoid unused variable error
	(void) sleep_max_;
#endif
}

void MMTReceiver::reset_sleep_progression()
{
	sleep_nsecs_ = sleep_min_;
}

/// @brief Retrieve data from node0 describing a "state" for
/// multistate design.
StateInputData
MMTReceiver::receive_state_input_data() const
{
	StateInputData input_data;
	input_data.state_index  = utility::receive_integer_from_node( 0 );
	input_data.pdb_name     = utility::receive_string_from_node( 0 );
	input_data.pdb_string   = utility::receive_string_from_node( 0 );
	input_data.corr_fname   = utility::receive_string_from_node( 0 );
	input_data.corr_file    = utility::receive_string_from_node( 0 );
	input_data.sec_resfname = utility::receive_string_from_node( 0 );
	input_data.sec_resfile  = utility::receive_string_from_node( 0 );

	core::Size n_npd_props  = utility::receive_integer_from_node( 0 );
	for ( core::Size ii = 1; ii <= n_npd_props; ++ii ) {
		core::Size npd_ind   = utility::receive_integer_from_node( 0 );
		std::string npd_name = utility::receive_string_from_node( 0 );
		input_data.npd_to_do_list.push_back( std::make_pair( npd_ind, npd_name ) );
	}

	return input_data;
}

/// @details Convert the data that can be used to construct the objects needed
/// in a repacking job to those objects.
StateData
MMTReceiver::initialize_state_data( StateInputData const & sid )
{
	StateData sd;

	// initialize the pose
	sd.pose = core::pose::PoseOP( new core::pose::Pose );
	core::import_pose::ImportPoseOptions import_opts;
	core::import_pose::pose_from_pdbstring( *sd.pose, sid.pdb_string, import_opts, sid.pdb_name );

	sd.ent_corr = entity_correspondence_from_string( sid.corr_file, sd.pose );

	// initialize the secondary resfile contents
	std::istringstream sec_resfile_stream( sid.sec_resfile );
	sd.secondary_resfile_contents = core::pack::task::ResfileContentsOP( new core::pack::task::ResfileContents( *sd.pose, sec_resfile_stream ) );

	sd.task =	initialize_packer_task( sd.pose, sd.ent_corr, *sd.secondary_resfile_contents );
	return sd;
}

void
MMTReceiver::receive_new_job()
{
	StateInputData sid    = receive_state_input_data();
	std::string seqstring = utility::receive_string_from_node( 0 );

	StateData sd = initialize_state_data( sid );
	restrict_task_for_job( sd.pose, sd.ent_corr, seqstring, sd.task );

	// replace the following line with a call to some factory to produce a generic MMTPackingJob
	MMTOffRotamerPackingJobOP job( new MMTOffRotamerPackingJob );

	job->set_pose( *sd.pose );
	job->set_sfxn( *sfxn_ );
	job->set_packer_task( *sd.task );
	job->set_npd_properties( sid.npd_to_do_list );

#ifdef MULTI_THREADED
#ifdef CXX11

	//TR << "Launching job in separate thread" << std::endl;

	// launch the job in a new thread and let that thread run in the background
	std::thread jobthread( &MMTPackingJob::go, job() );
	jobthread.detach();

	//TR << "Launched job in separate thread" << std::endl;

#endif
#else
	//TR << "Running job in main thread" << std::endl;
	job->go(); // run the packing job in the (single) MPI-listener thread
#endif

	StateAndSequencePair job_ssp = std::make_pair( sid.state_index, seqstring );
	curr_gen_jobs_[ job_ssp ] = job;
	running_jobs_.push_back( std::make_pair( job_ssp, job ) );
	++curr_njobs_running_;

}


protocols::pack_daemon::EntityCorrespondenceOP
MMTReceiver::entity_correspondence_from_string(
	std::string const & corr_file,
	core::pose::PoseOP pose
) const
{

	// initialize the entity correspondence
	protocols::pack_daemon::EntityCorrespondenceOP ent_corr( new protocols::pack_daemon::EntityCorrespondence );
	ent_corr->set_pose( pose );
	ent_corr->set_num_entities( num_entities_ );
	std::istringstream corrfile_stream( corr_file );
	ent_corr->initialize_from_correspondence_file( corrfile_stream );

	return ent_corr;
}

core::pack::task::PackerTaskOP
MMTReceiver::initialize_packer_task(
	core::pose::PoseOP pose,
	protocols::pack_daemon::EntityCorrespondenceOP ent_corr,
	core::pack::task::ResfileContents const & secondary_resfile_contents
) const
{
	// initialize the packer task
	core::pack::task::PackerTaskOP task = core::pack::task::TaskFactory::create_packer_task( *pose );
	protocols::pack_daemon::initialize_task_from_entity_resfile_and_secondary_resfile(
		*pose, ent_corr, *entity_resfile_contents_, secondary_resfile_contents, task );

	return task;
}

/// @throws Throws an EXCN_Msg_Exception if one of the residues in the entity correspondence file becomes
/// completely disabled or if the sequence that's broadcast to contains a non-canonical amino acid
void
MMTReceiver::restrict_task_for_job(
	core::pose::PoseOP pose,
	protocols::pack_daemon::EntityCorrespondenceOP ent_corr,
	std::string const & seqstring,
	core::pack::task::PackerTaskOP task
) const
{
	//restrict task for the given sequence string
	for ( core::Size ii = 1; ii <= pose->total_residue(); ++ii ) {
		core::Size ii_entity = ent_corr->entity_for_residue( ii );
		if ( ii_entity == 0 ) continue;
		utility::vector1< bool > aas_present( core::chemical::num_canonical_aas, false );
		char iiaachar = seqstring[ ii_entity-1 ];
		if ( ! core::chemical::oneletter_code_specifies_aa( iiaachar ) ) {
			throw utility::excn::EXCN_Msg_Exception( "Bad aa received in seqstring " + utility::to_string( ii_entity ) + " " + seqstring );
		}
		core::chemical::AA iiaa = core::chemical::aa_from_oneletter_code( iiaachar );
		aas_present[ iiaa ] = true;
		task->nonconst_residue_task( ii ).restrict_absent_canonical_aas( aas_present );
		if ( ! task->residue_task(ii).being_packed() ) {
			// sanity check; residue ii should be repacked if it's listed in the correspondence file
			throw utility::excn::EXCN_Msg_Exception( "Completely disabled a residue listed in the correspondence file.  Something has gone very wrong" );
		}
	}

}

void
MMTReceiver::delete_running_job( RunningJobsList::iterator & iter )
{
	--curr_njobs_running_;
	RunningJobsList::iterator iternext = iter;
	++iternext;
	running_jobs_.erase( iter );
	iter = iternext;
}


}
}
