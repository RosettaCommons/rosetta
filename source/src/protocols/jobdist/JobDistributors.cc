// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jobdist/JobDistributors.cc
/// @brief
/// @author Sergey Lyskov

#include <utility/io/izstream.hh>
#include <protocols/jobdist/JobDistributors.hh>
#include <basic/Tracer.hh>

#include <core/svn_version.hh>
#include <core/types.hh>
#include <core/import_pose/atom_tree_diffs/atom_tree_diff.hh>
#include <core/conformation/Residue.hh>
#include <basic/options/option.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>

#include <core/io/silent/util.hh>
#include <core/io/raw_data/DecoyFileData.hh>
#include <core/io/raw_data/ScoreFileData.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/pdb/pdb_writer.hh>

#include <numeric/random/random.hh>
#include <numeric/numeric.functions.hh>

#include <utility/exit.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/file/FileName.hh>
#include <utility/io/ozstream.hh>
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>
#ifdef USEMPI
#include <utility/string_util.hh>
#endif

// ObjexxFCL headers
#include <ObjexxFCL/string.functions.hh>

#include <map>
#include <set>
#include <sstream>
#include <string>


// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

#include <protocols/checkpoint/Checkpoint.hh>
#include <numeric/random/random.fwd.hh>


namespace protocols {
namespace jobdist {

static THREAD_LOCAL basic::Tracer JobDistributorTracer( "protocols.jobdist.JobDistributors" );

BaseJobDistributor::BaseJobDistributor(JobVector jobs):
	overwrite_( basic::options::option[ basic::options::OptionKeys::out::overwrite ] ),
	jobs_(jobs),
	current_job_(1),
	current_nstruct_(0),
	is_started_(false),
	nproc_( 0 ),
	proc_id_( 0 ),
	curr_jobid_( 0 ),
#ifdef USEMPI
	mpi_rank_( 0 ),
	mpi_nprocs_( 0 ),
	tag_( 1 ),
#endif
	ignorefinished_( false ),
	nooutput_( false ),
	inprogress_( false )
{
	start_time_ = time(nullptr);
}

BaseJobDistributor::BaseJobDistributor( BaseJobDistributor const & src ) :
	ReferenceCount(),
	overwrite_( src.overwrite_ ),
	jobs_( src.jobs_ ),
	current_job_( src.current_job_ ),
	current_nstruct_( src.current_nstruct_ ),
	is_started_( src.is_started_ ),
	nproc_( src.nproc_ ),
	proc_id_( src.proc_id_ ),
	curr_jobid_( src.curr_jobid_ ),
#ifdef USEMPI
	mpi_rank_( src.mpi_rank_ ),
	mpi_nprocs_( src.mpi_nprocs_ ),
	tag_( src.tag_ ),
#endif
	ignorefinished_( src.ignorefinished_ ),
	nooutput_( src.nooutput_ ),
	inprogress_( src.inprogress_ ),
	start_time_( src.start_time_ ),
	random_counter_( src.random_counter_ ),
	random_store_( src.random_store_ )
{

}


BaseJobDistributor::~BaseJobDistributor()
{
	if ( is_started_ ) {
		basic::Error() << "Must call shutdown() when finished using job distributor!" << std::endl;
	}
}


/// @details
/// Deliberately not virtual:  should not be overriden.
/// This is where to insert ifdefs and code for different cluster architectures!
bool BaseJobDistributor::next_job(BasicJobOP & job, int & struct_n)
{
#ifdef USEMPI
	JobDistributorTracer << "Node: " << mpi_rank_ << " next_job()" << std::endl;
#endif

	int elapsedtime = time(nullptr) - start_time_;
	if ( ( basic::options::option[ basic::options::OptionKeys::run::maxruntime ].user() ) &&
			( basic::options::option[ basic::options::OptionKeys::run::maxruntime ]() > 0 ) &&
			( basic::options::option[ basic::options::OptionKeys::run::maxruntime ]() < elapsedtime ) ) {
		std::cerr << "JobTerminated because runtime of " << elapsedtime << " s exceeded maxruntime of " << basic::options::option[ basic::options::OptionKeys::run::maxruntime ]() << " s " << std::endl;
		return false;
	}


	if ( !is_started_ ) {
		basic::Error() << "Must call startup() before using job distributor!" << std::endl;
		return false;
	}

#ifdef USEMPI
	if ( mpi_rank() == 0 ) {
		master_node_distribute_jobs();
		return false;
	} else {
		bool const job_recieved = request_job_from_master_node();
		if ( job_recieved ) {
			job = jobs_[ current_job_ ];
			struct_n = current_nstruct_;
		}
		job->set_output_file_name( get_output_filename() );
		return job_recieved;
	}
#else // one machine, Condor cluster, or BOINC
	bool job_found = find_available_job();
	if ( ! job_found ) return job_found;

	job = jobs_[ current_job_ ];
	struct_n = current_nstruct_;

	checkpoint_write();
	//current_nstruct_ += 1; this should -- and is now -- handled by find_available_job
	job->set_output_file_name( get_output_filename() );
	return true;
#endif // BOINC vs MPI vs etc

}

/// @details
/// Deliberately not virtual:  should not be overriden.
/// Iterate over the jobs that exist, check that nothing has been started
/// for each, and point the private data current_nstruct_ and current_job_ at
/// whichever available job it finds first.  Returns false if no job can be found.
bool BaseJobDistributor::find_available_job()
{
	bool shuffle_mode = basic::options::option[ basic::options::OptionKeys::run::shuffle ].user() ;

	while ( current_job_ <= jobs_.size() ) {

		// if shuffle mode then choose a random job and reset current_nstruct_ to 0. It should then work out which ones it's already done.
		if ( shuffle_mode ) {

			do{
				current_job_ = get_next_random_range( 1, jobs_.size() );
			}while( (current_job_ < 1) || (current_job_ > jobs_.size()) );
			current_nstruct_ = 0;

			// in shuffle mode nstruct means how many structures **in total** not per job.
			if ( (int)curr_jobid_ >= (int)basic::options::option[ basic::options::OptionKeys::out::nstruct]() ) return false;
		}
		while ( current_nstruct_ < jobs_[ current_job_ ]->nstruct() ) {
			++current_nstruct_; //running number within  current job
			++curr_jobid_; //running number across all jobs
			if ( shuffle_mode ) current_nstruct_ = curr_jobid_;

			JobDistributorTracer << "Looking for an available job: "
				<< current_nstruct_ << " "
				<< current_job_ << " "
				<< jobs_[ current_job_ ]->input_tag() << " "
				<< curr_jobid_
				<< std::endl;

			bool processed( !overwrite_ && is_finished( jobs_[ current_job_ ], current_nstruct_ ) );
			bool skipped(  nproc_ && numeric::mod( curr_jobid_, nproc_ ) != ( proc_id_ - 1 ) );
			if ( shuffle_mode && processed ) { break; }
			if ( !processed && !skipped ) {
#ifdef BOINC
				if( shuffle_mode ) {
				 	if (protocols::boinc::Boinc::worker_is_finished( curr_jobid_-1 )) return false;
				}else{
					if (protocols::boinc::Boinc::worker_is_finished( current_nstruct_-1 )) return false;
				}
#endif
				return true;
			}
		}
		++current_job_;
		current_nstruct_ = 0;
	}
	return false;
}

/// @details
/// If overriden by a subclass, it MUST call the superclass implementation.
void BaseJobDistributor::startup()
{
	if ( is_started_ ) {
		basic::Error() << "Distributor already started, don't call startup() again!" << std::endl;
	} else {
		checkpoint_read();

#ifdef USEMPI
		/// assume that a call to core_init(argv,argc) has already happened, and that MPI_Init has happened;
		runtime_assert( MPI_has_been_initialized() );
		MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank_);/* get current process id */
		MPI_Comm_size (MPI_COMM_WORLD, &mpi_nprocs_);/* get number of processes */
#endif

#ifdef BOINC
		protocols::boinc::Boinc::worker_startup();
#endif

		is_started_ = true;
	}
}

/// @details
/// If overriden by a subclass, it MUST call the superclass implementation.
void BaseJobDistributor::shutdown()
{
	if ( !is_started_ ) {
		basic::Error() << "Distributor not started or already stopped, don't call shutdown() again!" << std::endl;
	} else {
		checkpoint_clear();
		is_started_ = false;
	}

#ifdef USEMPI
	JobDistributorTracer << "Node " << mpi_rank_ << " -- ready to call mpi finalize" << std::endl;
	/// No other MPI calls may be made after this call
	MPI_Barrier( MPI_COMM_WORLD ); // make all nodes reach this point together.
	MPI_Finalize();
#endif

#ifdef BOINC
	bool shuffle_mode = basic::options::option[ basic::options::OptionKeys::run::shuffle ].user() ;
	if( shuffle_mode ){
		protocols::boinc::Boinc::worker_finish_summary( curr_jobid_-1, curr_jobid_-1, jobs_.size() );
	}else{
		protocols::boinc::Boinc::worker_finish_summary( current_nstruct_-1, current_nstruct_-1, jobs_.size() );
	}
	protocols::boinc::Boinc::worker_shutdown(); // Does not return.
#endif

}


/// @details
/// If overriden by a subclass, it MUST call the superclass implementation.
void BaseJobDistributor::begin_critical_section()
{
#ifdef BOINC
	boinc_begin_critical_section();
#endif // BOINC
}


/// @details
/// If overriden by a subclass, it MUST call the superclass implementation.
void BaseJobDistributor::end_critical_section()
{
#ifdef BOINC
	boinc_end_critical_section();
#endif // BOINC
}


/// @details
/// Needed to output .in_progress files for multiple processor jobs
void BaseJobDistributor::temp_file( std::string const & )
{
}

/// @details
/// Needed if a score_map is to be output by the derived class
void BaseJobDistributor::dump_pose_and_map( std::string const &, core::pose::Pose & )
{
}


/// @details
/// Needed if a score_map is to be output by the derived class
void BaseJobDistributor::score_map( std::map < std::string, core::Real > & )
{
}

std::string BaseJobDistributor::get_current_output_tag(){
	return current_job()->output_tag( current_nstruct() );
}

std::string BaseJobDistributor::get_output_filename() {
	return "STUBBED_IN_BASEJOBDISTRIBUTOR";
}


void BaseJobDistributor::checkpoint_read()
{
	begin_critical_section();
#ifdef BOINC
	if( utility::file::file_exists("rng.state.gz") ) {
		utility::io::izstream izs("rng.state.gz");
		numeric::random::rg().restoreState(izs);
		izs.close();
	}
#endif // BOINC
	end_critical_section();
}


/// @details Calls to this function are only a suggestion.
/// If checkpointing is expensive, this function must track time between calls
/// to avoid excessive disk activity, etc.
void BaseJobDistributor::checkpoint_write()
{
	begin_critical_section();
	static time_t last_chkpt_time = time(nullptr);
	time_t time_now = time(nullptr);
	// Refuse to checkpoint more than once a minute, no matter what BOINC wants.
	// Random number checkpoint files can be large (100k or more uncompressed).
	if ( time_now - last_chkpt_time > 60 ) {
#ifdef BOINC
		// BOINC automatically handles begin/end_critical_section() calls.
		utility::io::ozstream ozs("rng.state.gz");
		numeric::random::rg().saveState(ozs);
		ozs.close();
#endif // BOINC
		protocols::checkpoint::Timer::reset();
		last_chkpt_time = time_now;
	}
	end_critical_section();
}


void BaseJobDistributor::checkpoint_clear()
{
#ifdef BOINC
	if( utility::file::file_exists("rng.state.gz") ) {
		utility::file::file_delete("rng.state.gz");
	}
#endif // BOINC
}

BasicJobOP BaseJobDistributor::current_job() { return jobs_[current_job_]; }


#ifdef USEMPI

/// @brief Check that a call to MPI_Init has occurred already -- for use in runtime_assert statements
bool BaseJobDistributor::MPI_has_been_initialized() const {
	int already_initialized = 0;
	MPI_Initialized( & already_initialized );
	return already_initialized != 0;
}

/// @brief read access to derived classes for MPI related (const) data
/// @details must not be called until startup has been called
int BaseJobDistributor::mpi_rank() const {
	runtime_assert( is_started_ );
	return mpi_rank_;
}

/// @brief read access to derived classes for MPI related (const) data
/// @details must not be called until startup has been called.
int BaseJobDistributor::mpi_nprocs() const {
	runtime_assert( is_started_ );
	return mpi_nprocs_;
}

void BaseJobDistributor::master_node_distribute_jobs()
{
	// 1. Accept job requests from other processors, and distribute those jobs
	// in the form of unsigned integers: job_id and struct_n util all jobs have run out

	// 2. Accept job requests from other processors, but return job_id "0" as a signal
	// that no jobs remain, until all mpi_nprocs_ - 1 processors have been informed that no
	// jobs remain. ie. All of them have finished their jobs.

	while ( true ) {
		int node_requesting_job( 0 );

		JobDistributorTracer << "Master Node -- Waiting for job request; tag_ = " << tag_ << std::endl;
		MPI_Recv( & node_requesting_job, 1, MPI_INT, MPI_ANY_SOURCE, tag_, MPI_COMM_WORLD, & stat_ );
		bool const available_job_found = find_available_job();

		JobDistributorTracer << "Master Node --available job? " << available_job_found << std::endl;

		Size job_index = ( available_job_found ? current_job_ : 0 );
		int struct_n  = ( available_job_found ? current_nstruct_ : 0 );
		if ( ! available_job_found ) {
			JobDistributorTracer << "Master Node -- Spinning down node " << node_requesting_job << std::endl;
			MPI_Send( & job_index, 1, MPI_UNSIGNED_LONG, node_requesting_job, tag_, MPI_COMM_WORLD );
			break;
		} else {
			JobDistributorTracer << "Master Node -- Assigning job " << job_index << " " << struct_n << " to node " << node_requesting_job << std::endl;
			MPI_Send( & job_index, 1, MPI_UNSIGNED_LONG, node_requesting_job, tag_, MPI_COMM_WORLD );
			MPI_Send( & struct_n,  1, MPI_INT, node_requesting_job, tag_, MPI_COMM_WORLD );
			//		++current_nstruct_; handled now by find_available_job
		}
	}

	// we've just told one node to spin down, and
	// we don't have to spin ourselves down.
	Size nodes_left_to_spin_down( mpi_nprocs() - 1 - 1);

	while ( nodes_left_to_spin_down > 0 ) {
		int node_requesting_job( 0 );
		int recieve_from_any( MPI_ANY_SOURCE );
		MPI_Recv( & node_requesting_job, 1, MPI_INT, recieve_from_any, tag_, MPI_COMM_WORLD, & stat_ );
		Size job_index( 0 ); // No job left.
		MPI_Send( & job_index, 1, MPI_UNSIGNED_LONG, node_requesting_job, tag_, MPI_COMM_WORLD );
		JobDistributorTracer << "Master Node -- Spinning down node " << node_requesting_job << " with " << nodes_left_to_spin_down << " remaining nodes."  << std::endl;
		--nodes_left_to_spin_down;
	}

}

bool BaseJobDistributor::request_job_from_master_node()
{
	JobDistributorTracer << "Slave Node " << mpi_rank_ << " -- requesting job from master node; tag_ " << tag_ << std::endl;

	/// Tell the master node that I'm ready to recieve a new job.
	MPI_Send( & mpi_rank_, 1, MPI_INT, 0, tag_, MPI_COMM_WORLD );

	/// Get the index of the job that I'm to perform,
	MPI_Recv( & current_job_, 1, MPI_UNSIGNED_LONG, 0, tag_, MPI_COMM_WORLD, & stat_ );

	/// ... and if that job has a valid index, also get the index of the struct_n I'm to work on.
	/// otherwise, spin down.
	if ( current_job_ == 0 ) return false;

	MPI_Recv( & current_nstruct_, 1, MPI_INT, 0, tag_, MPI_COMM_WORLD, & stat_ );
	runtime_assert( current_nstruct_ != 0 ); // Master node should only return a valid struct_n.
	return true;
}
#endif

int BaseJobDistributor::get_next_random_range(int low, int high)
{
	// save yourself some random numbers.
	if ( random_store_.size() == 0 ) {
		for ( int k = 0; k < 1000; k ++ ) {
			random_store_.push_back( numeric::random::uniform() );
		}
		random_counter_ = 1;
	}

	if ( random_counter_ > 1000 ) random_counter_ = 1;

	if ( low > high ) {
		int temp;
		temp = low;
		low = high;
		high = temp;
	}

	int const range( high - low + 1 );

	int range_result = static_cast< int >( range * random_store_[random_counter_]) + low;
	random_counter_ ++;
	return range_result;
}

AtomTreeDiffJobDistributor::AtomTreeDiffJobDistributor(JobVector jobs, std::string outfile_name):
	BaseJobDistributor(jobs),
	out_(),
	used_tags_(),
	last_ref_pose_(nullptr),
	bb_precision_(6),
	sc_precision_(4),
	bondlen_precision_(2)
{
	// Add directory path and prefix/suffix (if specified) to plain file name:
	{
		utility::file::FileName outfile(outfile_name);
		std::ostringstream oss;
		oss << basic::options::option[ basic::options::OptionKeys::out::prefix ]() << outfile.base()
			<< basic::options::option[ basic::options::OptionKeys::out::suffix ]();
		outfile.base( oss.str() );
		outfile.path( basic::options::option[ basic::options::OptionKeys::out::path::pdb ]().path() );
		outfile.vol( basic::options::option[ basic::options::OptionKeys::out::path::pdb ]().vol() );
		if ( basic::options::option[ basic::options::OptionKeys::out::pdb_gz ] && outfile.ext() != "gz" ) {
			outfile.ext( ".gz" ); // else use input extension
		}
		outfile_name = outfile.name();
	}

	if ( utility::file::file_exists(outfile_name) ) {
		// load list of tags
		utility::io::izstream in( outfile_name.c_str() );
		if ( !in.good() ) {
			utility_exit_with_message( "Unable to open file: " + outfile_name + "\n" );
		}
		while ( !in.eof() ) {
			std::string tmp_tag;
			std::map< std::string, core::Real > tmp_scores;
			if ( ! core::import_pose::atom_tree_diffs::header_from_atom_tree_diff(in, tmp_tag, tmp_scores) ) break;
			used_tags_.insert( tmp_tag );
		}
		in.close();
		// re-open for appending (will fail & exit for gzipped files)
		out_.open_append( outfile_name.c_str() );
	} else {
		out_.open( outfile_name.c_str() );
		if ( basic::options::option[ basic::options::OptionKeys::run::version ] ) {
			out_ << "# Mini-Rosetta version " << core::minirosetta_svn_version() << " from " << core::minirosetta_svn_url() << "\n";
		}
	}
	if ( !out_.good() ) {
		utility_exit_with_message( "Unable to open file: " + outfile_name + "\n" );
	}
}

AtomTreeDiffJobDistributor::~AtomTreeDiffJobDistributor() = default;

void AtomTreeDiffJobDistributor::dump_pose(
	std::string const & tag,
	std::map< std::string, core::Real > const & scores,
	core::pose::Pose const & ref_pose,
	core::pose::Pose const & pose
)
{
	this->begin_critical_section();
	if ( last_ref_pose_ != &ref_pose ) {
		std::map< std::string, core::Real > empty_scores;
		core::import_pose::atom_tree_diffs::dump_reference_pose(out_, "%REF%_"+tag, empty_scores, ref_pose);
		last_ref_pose_ = &ref_pose;
	}
	if ( used_tags_.find(tag) != used_tags_.end() ) {
		basic::Error() << "Tag " << tag << " already exists in silent file; writing structure anyway..." << std::endl;
	}
	core::import_pose::atom_tree_diffs::dump_atom_tree_diff(out_, tag, scores, ref_pose, pose, bb_precision_, sc_precision_, bondlen_precision_);
	used_tags_.insert(tag);
	// Can't flush compressed streams -- results in file truncation
	if ( out_.uncompressed() ) out_.flush();
	this->end_critical_section();
}

/// @brief Sets number of digits used in writing atomtree diff.
void AtomTreeDiffJobDistributor::set_precision(
	int bb_precision,
	int sc_precision,
	int bondlen_precision
)
{
	bb_precision_ = bb_precision;
	sc_precision_ = sc_precision;
	bondlen_precision_ = bondlen_precision;
}

void AtomTreeDiffJobDistributor::shutdown()
{
	out_.close();
	BaseJobDistributor::shutdown();
}

bool AtomTreeDiffJobDistributor::is_finished(BasicJobOP const & job, int struct_n)
{
	return ( used_tags_.find(job->output_tag(struct_n)) != used_tags_.end() );
}


PlainPdbJobDistributor::PlainPdbJobDistributor(JobVector jobs, std::string outfile_name):
	BaseJobDistributor(jobs),
	scorefile_name_()
{
	if ( outfile_name != "none" ) {
		scorefile_ = true;
		// set up all the information for the scorefile
		utility::file::FileName outfile("");
		std::ostringstream oss;
		oss << basic::options::option[ basic::options::OptionKeys::out::prefix ]() << outfile.base()
			<< basic::options::option[ basic::options::OptionKeys::out::suffix ]();
		outfile.base( oss.str() );
		outfile.path( basic::options::option[ basic::options::OptionKeys::out::path::pdb ]().path() );
		outfile.vol( basic::options::option[ basic::options::OptionKeys::out::path::pdb ]().vol() );
		// determine the extension based on fullatom or centroid
		if ( basic::options::option[ basic::options::OptionKeys::out::file::fullatom ] && outfile.ext() != "fasc" ) {
			outfile.ext( ".fasc" ); // else use input extension
		} else {
			outfile.ext( ".sc" );
		}
		outfile_name = outfile.name();
		scorefile_name_ = outfile.name();
	} else {
		scorefile_ = false;
	}
}

PlainPdbJobDistributor::~PlainPdbJobDistributor() = default;


/// @details
/// Over riding baseclass to enable inprogress file.
void PlainPdbJobDistributor::startup()
{
	/// this could probably be done with a commandline flag to allow disabling
	this->enable_inprogress();

	// calling base class startup
	BaseJobDistributor::startup();
}

std::string PlainPdbJobDistributor::get_output_filename(std::string const & tag)
{
#ifndef USEMPI
	utility::file::FileName output_pdb_name;
	output_pdb_name.path( basic::options::option[ basic::options::OptionKeys::out::path::pdb ]().path() );
	output_pdb_name.vol( basic::options::option[ basic::options::OptionKeys::out::path::pdb ]().vol() );
	output_pdb_name.base( tag ); // could contain embedded dots that look ~ like extensions
	if ( basic::options::option[ basic::options::OptionKeys::out::pdb_gz ] ) {
		output_pdb_name.ext( ".pdb.gz" );
	} else {
		output_pdb_name.ext( ".pdb" );
	}
	return output_pdb_name.name();

#endif
#ifdef USEMPI
  /// Requires that outdir_{0..(nprocs-1)} exist
  /// burden is on MPI user to create those directories before launching.
  /// note: node0 does no work, so that dir never gets written to.
	utility::file::FileName output_pdb_name(tag);
  output_pdb_name.path(
											 basic::options::option[ basic::options::OptionKeys::out::path::pdb ]().path() +
											 "/outdir_" + utility::to_string( parent::mpi_rank() ) + "/" );
  output_pdb_name.vol( basic::options::option[ basic::options::OptionKeys::out::path::pdb ]().vol() );
	if( basic::options::option[ basic::options::OptionKeys::out::pdb_gz ] ) {
		output_pdb_name.ext( ".pdb.gz" );
	} else {
		output_pdb_name.ext( ".pdb" );
	}
	std::cout << "output file name: " << output_pdb_name.name() << std::endl;
  return output_pdb_name.name();
#endif
}


void PlainPdbJobDistributor::temp_file(std::string const & tag)
{
	utility::io::ozstream tempfile;
	std::string output_tag = get_output_filename( tag );
	tempfile.open( output_tag + ".in_progress" );
}

void PlainPdbJobDistributor::dump_pose_and_map(
	std::string const & tag,
	core::pose::Pose & pose
)
{
	if ( BaseJobDistributor::nooutput() ) return;

	this->begin_critical_section();
	std::string outfile_name = this->get_output_filename(tag);
	utility::io::ozstream out( outfile_name.c_str() );
	if ( !out.good() ) {
		utility_exit_with_message( "Unable to open file: " + outfile_name + "\n" );
	}
	core::io::pdb::dump_pdb( pose, out );
	dump_scores(out, tag, pose);
	out.close();

	if ( utility::file::file_exists( outfile_name + ".in_progress" ) ) {
		utility::file::file_delete( outfile_name + ".in_progress" );
	}

	if ( scorefile_ ) {
		std::string name ("");
		if ( basic::options::option[ basic::options::OptionKeys::out::file::o ].user() ) {
			name = scorefile_name_.base()+basic::options::option[ basic::options::OptionKeys::out::file::o ]()+"."+scorefile_name_.ext();
		} else {
			name = scorefile_name_.base()+"score."+scorefile_name_.ext();
		}
		core::io::raw_data::ScoreFileData sfd(name);
		sfd.write_pose( pose, score_map_, outfile_name );
	}
	this->end_critical_section();
}

bool PlainPdbJobDistributor::is_finished( BasicJobOP const & job, int struct_n )
{
	if ( BaseJobDistributor::ignorefinished() ) return false;
	bool file_exists (false);
	std::string filename ( get_output_filename(job->output_tag(struct_n))+".in_progress" );
	if ( utility::file::file_exists(filename) ) {
		file_exists = true;
	} else if ( utility::file::file_exists( get_output_filename(job->output_tag(struct_n)) ) ) {
		file_exists = true;
	}
	return file_exists;
}


/// @details In order for this to work as expected, the Pose's cached energies
/// must match up with the (current) conformation.
/// A good time to do this is at the end of your protocol's apply() method:
///   scorefxn( pose );
///   scorefxn.accumulate_residue_total_energies( pose );
void PlainPdbJobDistributor::dump_scores(
	utility::io::ozstream & out,
	std::string const & tag,
	core::pose::Pose & pose
)
{
	// Which score terms to use
	core::scoring::EnergyMap weights = pose.energies().weights();
	typedef utility::vector1<core::scoring::ScoreType> ScoreTypeVec;
	ScoreTypeVec score_types;
	for ( int i = 1; i <= core::scoring::n_score_types; ++i ) {
		core::scoring::ScoreType ii = core::scoring::ScoreType(i);
		if ( weights[ii] != 0 ) score_types.push_back(ii);
	}
	// This version is formatted for easy parsing by R, Excel, etc.
	out << "# All scores below are weighted scores, not raw scores.\n";
	out << "#BEGIN_POSE_ENERGIES_TABLE " << tag << "\n";
	out << "label";
	for ( auto & score_type : score_types ) {
		out << " " << name_from_score_type(score_type);
	}
	out << " total\n";
	out << "weights";
	for ( auto & score_type : score_types ) {
		out << " " << weights[score_type];
	}
	out << " NA\n";
	out << "pose";
	core::Real pose_total = 0.0;
	for ( auto & score_type : score_types ) {
		core::Real score = (weights[score_type] * pose.energies().total_energies()[ score_type ]);
		out << " " << score;
		pose_total += score;
	}
	out << " " << pose_total << "\n";
	for ( core::Size j = 1, end_j = pose.size(); j <= end_j; ++j ) {
		core::Real rsd_total = 0.0;
		out << pose.residue(j).name() << "_" << j;
		for ( auto & score_type : score_types ) {
			core::Real score = (weights[score_type] * pose.energies().residue_total_energies(j)[ score_type ]);
			out << " " << score;
			rsd_total += score;
		}
		out << " " << rsd_total << "\n";
	}
	out << "#END_POSE_ENERGIES_TABLE " << tag << "\n";
	// This version uses YAML instead -- hard to read by eye, requires scripts to process.
	//basic::EmitterOP yaml = new basic::YamlEmitter(out);
	//yaml->start_doc(); // This breaks the output into chunks that fit in memory, yay!
	//yaml->start_map(tag);
	//yaml->write("tag", tag); // redundant, but Charlie prefers it here too
	//yaml->write("total_score", total_score);
	//yaml->start_list("score_names", false);
	//for(ScoreTypeVec::iterator ii = score_types.begin(), end_ii = score_types.end(); ii != end_ii; ++ii)
	// yaml->write(name_from_score_type(*ii));
	//yaml->end_list();
	//yaml->start_list("score_weights", false);
	//for(ScoreTypeVec::iterator ii = score_types.begin(), end_ii = score_types.end(); ii != end_ii; ++ii)
	// yaml->write(scorefxn.get_weight(*ii));
	//yaml->end_list();
	//yaml->start_list("scores_raw", false);
	//for(ScoreTypeVec::iterator ii = score_types.begin(), end_ii = score_types.end(); ii != end_ii; ++ii)
	// yaml->write(pose.energies().total_energies()[ *ii ]);
	//yaml->end_list();
	//yaml->start_list("scores_weighted", false);
	//for(ScoreTypeVec::iterator ii = score_types.begin(), end_ii = score_types.end(); ii != end_ii; ++ii)
	// yaml->write(scorefxn.get_weight(*ii) * pose.energies().total_energies()[ *ii ]);
	//yaml->end_list();
	//yaml->start_map("per_res_weighted");
	//for(core::Size j = 1, end_j = pose.size(); j <= end_j; ++j) {
	// std::ostringstream resname; resname << pose.residue(j).name() << " " << j;
	// yaml->start_list(resname.str(), false);
	// for(ScoreTypeVec::iterator ii = score_types.begin(), end_ii = score_types.end(); ii != end_ii; ++ii)
	//  yaml->write(scorefxn.get_weight(*ii) * pose.energies().residue_total_energies(j)[ *ii ]);
	// yaml->end_list();
	//}
	//yaml->end_map();
	//yaml->end_map();
	////yaml->end();
}

PlainRawJobDistributor::PlainRawJobDistributor(JobVector jobs, std::string outfile_name):
	BaseJobDistributor(jobs),
	rawfile_name_(), used_tags_()
{
	utility::file::FileName outfile(outfile_name);
	std::ostringstream oss;
	oss << basic::options::option[ basic::options::OptionKeys::out::prefix ]() << outfile.base()
		<< basic::options::option[ basic::options::OptionKeys::out::suffix ]();
	outfile.base( oss.str() );
	outfile.path( basic::options::option[ basic::options::OptionKeys::out::path::pdb ]().path() );
	outfile.vol( basic::options::option[ basic::options::OptionKeys::out::path::pdb ]().vol() );
	outfile_name = outfile.name();
	rawfile_name_ = outfile.name();

	if ( !utility::file::file_exists( rawfile_name_ ) ) return;
	core::io::raw_data::DecoyFileData dfd(rawfile_name_);
	core::io::raw_data::StructureMap::const_iterator iter;
	used_tags_ = dfd.read_tags_fast( rawfile_name_ );
	utility::vector1< std::string >::const_iterator i;
	for ( i = used_tags_.begin(); i != used_tags_.end(); ++i ) {
		JobDistributorTracer << *i << std::endl;
	}
}

PlainRawJobDistributor::~PlainRawJobDistributor() = default;

void PlainRawJobDistributor::dump_pose_and_map(
	std::string const & tag,
	core::pose::Pose & pose
)
{
	this->begin_critical_section();
	std::string output_tag = get_output_filename( tag );
	bool fa = basic::options::option[ basic::options::OptionKeys::out::file::fullatom ];

	core::io::raw_data::DecoyFileData dfd(rawfile_name_);
	dfd.write_pose( pose, score_map_, output_tag, fa );
	this->end_critical_section();
}

std::string PlainRawJobDistributor::get_output_filename(std::string const & tag)
{

#ifndef USEMPI
	utility::file::FileName output_pdb_name(tag);
	output_pdb_name.path( basic::options::option[ basic::options::OptionKeys::out::path::pdb ]().path() );
	output_pdb_name.vol( basic::options::option[ basic::options::OptionKeys::out::path::pdb ]().vol() );
	return output_pdb_name.name();
#endif
#ifdef USEMPI
	/// Requires that outdir_{0..(nprocs-1)} exist
	/// burden is on MPI user to create those directories before launching.
	/// note: node0 does no work, so that dir never gets written to.
	utility::file::FileName output_pdb_name(tag);
	output_pdb_name.path(
		basic::options::option[ basic::options::OptionKeys::out::path::pdb ]().path() +
		"/outdir_" + utility::to_string( parent::mpi_rank() ) + "/" );
	output_pdb_name.vol( basic::options::option[ basic::options::OptionKeys::out::path::pdb ]().vol() );
	std::cout << "output file name: " << output_pdb_name.name() << std::endl;
	return output_pdb_name.name();
#endif
}

std::string PlainRawJobDistributor::get_output_tag( int const & struct_n )
{ return "D_" + ObjexxFCL::right_string_of(struct_n, TAG_NUM_FORMAT_LENGTH, '0'); }

bool PlainRawJobDistributor::is_finished(BasicJobOP const & job, int struct_n )
{
	std::string output_tag = get_output_filename(job->output_tag(struct_n));
	bool already_processed = false;
	utility::vector1< std::string >::const_iterator i;
	for ( i = used_tags_.begin(); i != used_tags_.end(); ++i ) {
		if ( (*i) == output_tag ) {
			already_processed = true;
			JobDistributorTracer << "Tag: " << output_tag << " " << job->output_tag(struct_n) <<
				" - already processed"  << std::endl;
			break;
		}
	}
	return already_processed;
}


PlainSilentFileJobDistributor::PlainSilentFileJobDistributor(JobVector jobs):
	BaseJobDistributor(jobs),
	used_tags_()
{
	///... relocated reading of silent files to startup()
}

PlainSilentFileJobDistributor::~PlainSilentFileJobDistributor() = default;

void PlainSilentFileJobDistributor::dump_pose(
	BasicJobOP const & job,
	int const & struct_n,
	bool const & /*fullatom*/,
	core::pose::Pose & pose
)
{
	if ( BaseJobDistributor::nooutput() ) return;

	this->begin_critical_section();
	std::string silent_file = get_output_filename();
	std::string output_tag = get_output_tag( job, struct_n );

	using namespace core::io::silent;
	SilentFileOptions opts;
	SilentFileDataOP sfd( new core::io::silent::SilentFileData( opts ) );
	SilentStructOP ss = SilentStructFactory::get_instance()->get_instance()->get_silent_struct_out( opts );
	ss->fill_struct( pose, output_tag );
	sfd->write_silent_struct( *ss, silent_file );
	this->end_critical_section();
}

std::string PlainSilentFileJobDistributor::get_output_filename() {
	std::string silent_file = basic::options::option[ basic::options::OptionKeys::out::file::silent ]();

#ifdef USEMPI
	// attach mpi rank to out files
	size_t lastslash = silent_file.find_last_of("/\\");
	size_t lastdot   = silent_file.find_last_of('.');

	if ( !basic::options::option[ basic::options::OptionKeys::out::path::mpi_rank_dir ]() ) {
		if ( lastdot == silent_file.npos || (lastslash > lastdot && lastslash != silent_file.npos) ) {
			silent_file = silent_file+"_"+utility::to_string( parent::mpi_rank() );
		} else {
			silent_file = silent_file.substr( 0,lastdot )+"_"+utility::to_string( parent::mpi_rank() )+silent_file.substr( lastdot );
		}
	} else {
		if ( lastslash == silent_file.npos)  {
			silent_file = utility::to_string( parent::mpi_rank() )+"/"+silent_file;
		} else {
			silent_file = silent_file.substr( 0,lastslash )+"/"+utility::to_string( parent::mpi_rank() )+silent_file.substr( lastslash );
		}
	}
#endif

	return silent_file;
}


void PlainSilentFileJobDistributor::startup() {

	// calling base class startup
	BaseJobDistributor::startup();

	//this call requires that startup() has been called.... moved from constructor
	std::string const silent_filename = get_output_filename();

	if ( !utility::file::file_exists( silent_filename ) ) return;
	core::io::silent::SilentFileOptions opts;
	core::io::silent::SilentFileData sfd( opts );
	core::io::silent::Structure_Map::const_iterator iter;
	used_tags_ = sfd.read_tags_fast( silent_filename );
	utility::vector1< std::string >::const_iterator i;
	for ( i = used_tags_.begin(); i != used_tags_.end(); ++i ) {
		JobDistributorTracer << *i << std::endl;
	}
}

void PlainSilentFileJobDistributor::shutdown()
{
	core::io::silent::gzip();
	BaseJobDistributor::shutdown();
}

bool PlainSilentFileJobDistributor::is_finished(BasicJobOP const & job, int struct_n )
{
	std::string querytag = job->output_tag(struct_n);
	bool found = ( std::find( used_tags_.begin(), used_tags_.end(), querytag ) != used_tags_.end() );
	if ( found ) return found;
	// compare string without the first two chars as they might be changed from "S_" to "F_"
	querytag = get_output_tag( job, struct_n );
	// this might cause problems ... if we don't have tags longer than 2 characters...
	querytag = querytag.substr(2);
	for ( core::Size i = 1; i <= used_tags_.size(); ++i ) {
		if ( used_tags_[i].substr(2) == querytag ) {
			found = true;
			return found;
		}
	}
	return found;
}


void PlainSilentFileJobDistributor::dump_silent(
	int const & /*struct_n*/,
	core::io::silent::SilentStruct & silent_struct
)
{
	this->begin_critical_section();
	std::string silent_file = get_output_filename();

	// std::string output_tag = get_output_tag( struct_n );
	// silent_struct.decoy_tag( output_tag );
	core::io::silent::SilentFileOptions opts;
	core::io::silent::SilentFileDataOP sfd( new core::io::silent::SilentFileData( opts ) );
	sfd->write_silent_struct( silent_struct, silent_file );
	this->end_critical_section();
}

void PlainSilentFileJobDistributor::dump_silent(
	core::io::silent::SilentFileData const& sfd
)
{
	this->begin_critical_section();
	std::string silent_file = get_output_filename();

	// std::string output_tag = get_output_tag( struct_n );
	// silent_struct.decoy_tag( output_tag );
	sfd.write_all( silent_file );
	this->end_critical_section();
}

/// @details override base class function to
std::string PlainSilentFileJobDistributor::get_current_output_tag( ){
	return get_output_tag( parent::current_job(), parent::current_nstruct() );
}

/// @details returns an output tag generated as follows
/// "S_" + the current jobs tag + "option[ user_tag ]" + nstruct
std::string PlainSilentFileJobDistributor::get_output_tag( BasicJobOP const & job, int const & struct_n ) const
{
	std::string const prefix( "S" );

	std::string number("");
	//why not:
	// if ( job->nstruct() > 1 ) { causes problems ... probably in query_tag
	number = "_"+ ObjexxFCL::right_string_of(struct_n, TAG_NUM_FORMAT_LENGTH, '0');
	// }

	std::string itag( job->input_tag() );
	if ( itag.size() ) itag = "_" + itag;

	std::string user_tag("");
	if ( basic::options::option[ basic::options::OptionKeys::out::user_tag ].user() ) {
		user_tag = "_" + basic::options::option[ basic::options::OptionKeys::out::user_tag ];
	}

	return prefix + itag + user_tag + number;

}


}
}
