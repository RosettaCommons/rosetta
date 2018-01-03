// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide_predict/SimpleCycpepPredictApplication_MPI.hh
/// @brief Application-level code for the simple_cycpep_predict app, MPI version.
/// @details  This application predicts structures of simple backbone-cyclized peptides made of alpha-, beta-, or gamma-amino acids (of any chirality)
/// using generalized kinematic closure (GenKIC) for cyclization, and enforcing user-defined requiresments for numbers of mainchain hydrogen bonds.  This
/// version uses MPI for cross-communication with parallel processes.
/// On 4 Aug 2017, work commenced to make this appliction compatible with a job distribution scheme in which the last level of the hierarchy splits
/// jobs over many threads (hierarchical process-level parallelism plus thread-level parallelism).
/// @author Vikram K. Mulligan, Baker laboratory (vmullig@uw.edu)

//#define USEMPI //DELETE ME -- this is temporarily here for my IDE
//#define MULTI_THREADED //DELETE ME -- this is temporarily here for my IDE

#ifdef USEMPI

#ifndef INCLUDED_protocols_cyclic_peptide_predict_SimpleCycpepPredictApplication_MPI_hh
#define INCLUDED_protocols_cyclic_peptide_predict_SimpleCycpepPredictApplication_MPI_hh

// Unit Headers
#include <protocols/cyclic_peptide_predict/SimpleCycpepPredictApplication_MPI.fwd.hh>
#include <protocols/cyclic_peptide_predict/SimpleCycpepPredictApplication_MPI_JobResultsSummary.fwd.hh>
#include <protocols/cyclic_peptide_predict/util.hh>

// Package Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Project Headers
#include <protocols/cyclic_peptide/DeclareBond.fwd.hh>
#include <protocols/filters/BasicFilters.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/Ramachandran.hh>
#include <core/io/silent/SilentStruct.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C++ headers
#include <stdio.h>
#include <time.h>

#ifdef MULTI_THREADED
#include <mutex>
#endif

namespace protocols {
namespace cyclic_peptide_predict {

enum SIMPLE_CYCPEP_PREDICT_MPI_COMMUNICATION_TYPE {
	NULL_MESSAGE = 1,
	REQUEST_NEW_JOBS_BATCH_UPWARD,
	OFFER_NEW_JOBRESULTSSUMMARY_BATCH_UPWARD,
	OFFER_NEW_JOBS_ATTEMPTED_COUNT_UPWARD,
	OFFER_NEW_POSE_BATCH_UPWARD,
	SILENT_STRUCT_TRANSMISSION,
	GIVE_COMPLETION_SIGNAL_UPWARD,
	OFFER_NEW_JOBS_BATCH_DOWNWARD,
	REQUEST_NEW_POSE_BATCH_DOWNWARD,
	HALT_SIGNAL
};

enum SIMPLE_CYCPEP_PREDICT_MPI_TAG_TYPE {
	GENERAL_REQUEST=1,
	NEW_JOBS_DOWNWARD,
	RESULTS_SUMMARY_UPWARD,
	JOBS_ATTEMPTED_COUNT_UPWARD
};

/// @brief Application-level code for simple_cycpep_predict application, MPI version.
/// @author Vikram K. Mulligan, Baker laboratory (vmullig@uw.edu)
class SimpleCycpepPredictApplication_MPI : public utility::pointer::ReferenceCount
{
public:
	/// @brief Constructor
	///
	SimpleCycpepPredictApplication_MPI();

	/// @brief Constructor with options
	///
	SimpleCycpepPredictApplication_MPI(
		int const MPI_rank,
		int const MPI_n_procs,
		core::scoring::ScoreFunctionCOP sfxn_in,
		core::Size const total_hierarchy_levels,
		utility::vector1 < core::Size > const & procs_per_hierarchy_level,
		utility::vector1 < core::Size > const &batchsize_per_level,
		std::string const &sort_type,
		bool const select_highest,
		core::Real const &output_fraction,
		std::string const &output_filename,
		core::Real const &lambda,
		core::Real const &kbt,
		core::Size const threads_per_slave_proc //Only used in multi-threaded build.
	);

	/// @brief Explicit virtual destructor.
	///
	virtual ~SimpleCycpepPredictApplication_MPI();

	/// @brief Explicit copy constructor.
	///
	SimpleCycpepPredictApplication_MPI( SimpleCycpepPredictApplication_MPI const &src );

#ifdef MULTI_THREADED
	/// @brief Initialize private member variables related to thread random seeds from the options
	/// system.  Does nothing if this isn't a multi-threaded compilation.
	void init_random_info_from_options();
#endif

	/// @brief Actually run the application.
	/// @details On slave nodes, this creates an instance of SimpleCycpepPredictApplication and runs that.  Nodes higher
	/// in the communications hierarchy are just involved in sending out jobs and pulling in results.
	/// @note The run() function is nonconst due to some setup steps that it performs.  It then calls const run functions
	/// for emperor, master, and slave nodes.
	void run();

public:
	/// ------------- Public Methods ---------------------

	/// @brief Set the sort type, by string.
	///	
	void set_sort_type( std::string const &sort_type );

	/// @brief Set the sort type, by enum.
	///	
	void set_sort_type( SIMPLE_CYCPEP_PREDICT_MPI_SORT_TYPE const sort_type );

	/// @brief Set the ouput fraction.
	/// @details Checks that this is between 0 and 1.
	void set_output_fraction( core::Real const &val );


private:
	/// ------------- Methods ----------------------------

	/// @brief The bredth of the Gaussian used to determine whether something is
	/// native or not, for calculating PNear (the funnel quality metric).
	inline core::Real const & lambda() const { return lambda_; }

	/// @brief The Boltzmann temperature used for calculating PNear (the
	/// funnel quality metric).
	inline core::Real const & kbt() const { return kbt_; }

	/// @brief Check the current time and determine whether it's past the timeout time.
	///
	bool halting_condition( clock_t const start_time, core::Size const timeout ) const;

	/// @brief Set the number of processes at each level of the communications hierarchy.
	/// @details The total_hierarchy_levels, MPI_rank, and MPI_n_procs variables must be set first.  This does
	/// some checks to make sure that data_in[1] is 1, that data_in[n] >= data_in[n-1], and that the sum of entries
	/// equals the total number of processes.  The total number of hierarchy levels must also match total_hierarchy_levels.
	void set_procs_per_hierarchy_level( utility::vector1 < core::Size > const &data_in );

	/// @brief Sets the batch size for this proc.
	/// @details The total_hierarchy_levels, MPI_rank, MPI_n_procs, my_children_, my_parent_, and hierarchy_level_ variables
	/// must be set first.  This does some checks to ensure that the batch sizes of subsequent levels are smaller and previous
	/// levels are larger.
	void set_batchsize( utility::vector1< core::Size > const &sizes_in );

	/// @brief Figure out which processes are my children, and which is my parent.  Also, figure out the level in the
	/// communications hierarchy that I'm in.
	/// @details This must be done AFTER the MPI_rank_, MPI_n_procs_, total_hierarchy_levels_, and procs_per_hierarchy_level_,
	/// variables are set.
	void assign_level_children_and_parent();

	/// @brief Get the amino acid sequence of the peptide we're going to predict; set the sequence_ private member variable.
	/// @details The emperor reads this from disk and broadcasts it to all other nodes.  This function should be called from all nodes;
	/// it figures out which behaviour it should be performing.
	/// @note This function necessarily uses new and delete[].
	void get_sequence();

	/// @brief Get the native structure of the peptide we're going to predict (if the user has specified one with the -in:file:native flag).
	/// @details The emperor reads this from disk and broadcasts it to all other nodes.  This function should be called from all nodes;
	/// it figures out which behaviour it should be performing.
	void get_native();

	/// @brief Get the allowed amino acids at each position, if we're designing.
	/// @details The emperor reads these from disk and broadcasts them to all other nodes.  This function should be called from all nodes;
	/// it figures out which behaviour it should be performing.
	void get_design_settings();

	/// @brief Given a map of indicies to lists of residue names, broadcast it to all MPI ranks.
	/// @note This necessarily uses new and delete for the data sent via MPI_Bcast.
	void broadcast_res_list( std::map< core::Size, utility::vector1 < std::string > > &res_list) const;

	/// @brief The emperor sends a message to everyone letting them know it's time to start.
	/// @details Following the go signal, slaves send requests for jobs upward.
	void go_signal() const;

	/// @brief The emperor sends a message to everyone letting them know it's time to stop.
	/// @details Following the stop signal, SimpleCycpepPredictApplication_MPI::run() terminates.
	void stop_signal() const;

	/// @brief Any non-slave node can wait for a node, above or below in the hierarchy, to send it some sort of request.
	/// @details Only messags with tag GENERAL_REQUEST.
	/// @param[out] requesting_node The node from which the request came.
	/// @param[out] message The type of request received.
	void wait_for_request( int &requesting_node, SIMPLE_CYCPEP_PREDICT_MPI_COMMUNICATION_TYPE &message ) const;

	/// @brief Send a signal to stop job distribution to a set of intermediate masters.
	///
	void send_halt_signal( utility::vector1 < int > const &ranks_to_target ) const;

	/// @brief Any node can send some sort of request to a specific node in the hierarchy (parent or child).
	/// @details Sends messags with tag GENERAL_REQUEST.
	/// @param[in] target_node The node to which we're sending the request.
	/// @param[in] message The type of request to send.
	void send_request( int const target_node, SIMPLE_CYCPEP_PREDICT_MPI_COMMUNICATION_TYPE const message ) const;

	/// @brief Send an integer to a child indicating that the child is now holding this number of jobs.
	/// @details Sends messags with tag NEW_JOBS_DOWNWARD.  Sends zero if no jobs remain to send.
	/// @param[in] target_node The rank of the node to which we're sending the request.
	/// @param[in,out] structures_remaining_to_send The number of structures remaining in the send buffer.  Decremented by this operation to a minimum of zero.
	/// @param[in] number_to_send How many jobs should I send down?
	void send_new_jobs_downward( int const target_node, core::Size &structures_remaining_to_send, core::Size const number_to_send ) const;

	/// @brief Receive an integer from my parent indicating that I should add N jobs to my set to carry out or to pass to my children.
	/// @details Recieves message with tag NEW_JOBS_DOWNARD.  Receives zero if no jobs remain in parent to send.
	/// @param[in,out] njobs The number of jobs held on this node that are to be done.  Incremented by this function with however many are recieved from above.
	void receive_njobs_from_above( core::Size &njobs ) const;

	/// @brief Non-emperor nodes must call this when the emperor calls emperor_broadcast_silent_struct.
	/// @details This will build a pose and return an owning pointer to it.
	/// @note This function necessarily uses new and delete[].
	core::pose::PoseCOP receive_broadcast_silent_struct_and_build_pose() const;

	/// @brief Convert a vector of silent structs into a character string and send it to a node.
	/// @details Intended to be used with receive_pose_batch_as_string() to allow another node to receive the transmission.
	/// Message is tagged with SILENT_STRUCT_TRANSMISSION
	/// @note This function necessarily uses new and delete[].
	void send_silent_structs( utility::vector1 < core::io::silent::SilentStructOP > const &ss_vect, int const target_node) const;

	/// @brief Receive a transmitted set of poses (as silent strings).
	/// @details Appends received information to results_string.
	void receive_pose_batch_as_string( int const originating_node, std::string &results_string ) const;

	/// @brief Receive the number of jobs attempted by all of the nodes below a child node, and add this to the total.
	///
	void receive_jobs_attempted_count( core::Size &total_jobs_attempted,	int const originating_node ) const;

	/// @brief Send the number of jobs attempted by this nodes or all of the nodes below this node.
	///
	void send_jobs_attempted_count( core::Size const total_jobs_attempted,	int const target_node) const;	

	/// @brief Recieve a sorted list of job summaries, and merge them with an existing sorted list to make a combined sorted list.
	/// @details To be used in conjunction with send_job_summaries().  Sending and receiving procs must send messages to synchronize, first.
	/// @param[in,out] original_summary_list The current list of job summaries on this node, sorted.  A new list will be appended and sorted into this one.  The result is a sorted list.
	/// @param[in] originating_node The MPI rank from which we're receiving summaries.
	/// @param[in] append_to_handler_list If true, the rank of the originating node is appended to the list of nodes that handled each JobSummary.  This facilitates sending a message back down
	/// the hierarchy to the original node that carried out the job.
	void receive_and_sort_job_summaries( utility::vector1< SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP > &original_summary_list, int const originating_node, bool const append_to_handler_list ) const;

	/// @brief Recieve a list of job summaries.
	/// @details To be used in conjunction with receive_and_sort_job_summaries().  Sending and receiving procs must send messages to synchronize, first.
	void send_job_summaries( utility::vector1< SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP > const &summary_list, int const target_node ) const;

	/// @brief Given a short list of job summaries, split the list by the index of the node that I'd have to send the request to, and send requests for full poses down the hierarchy.
	/// @details Throws an error if any of the jobs in the list cannot be reached by sending a request by way of my_children_.  To be used with receive_pose_requests_from_above().
	/// @param[in] summary_shortlist The list of job summaries to split up and send downard.
	/// @param[out] children_receiving_nonzero_requests The number of children receiving a request for at least one pose.
	void request_poses_from_below( utility::vector1< SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP > const &summary_shortlist, core::Size &children_receiving_nonzero_requests ) const;

	/// @brief Recieve a short list of job summaries from above.
	/// @param[out] summary_shortlist The job summaries list, cleared and populated by this function.
	void receive_pose_requests_from_above( utility::vector1< SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP > &summary_shortlist ) const;

	/// @brief Given a string on the emperor node, send it to all nodes.
	/// @details The "mystring" string is the input on the emperor node, and the output on all other nodes.
	void broadcast_string_from_emperor( std::string &mystring ) const;

	/// ------------- Emperor Methods --------------------

	/// @brief Is this an emperor (root) node?
	/// @details The emperor is responsible for sending out and retrieving all jobs, and for all file output.
	bool i_am_emperor() const;

	/// @brief The jobs done by the emperor during parallel execution.
	/// @details The emperor is responsible for sending out and retrieving all jobs, and for all file output.
	void run_emperor() const;

	/// @brief Convert a silent struct into a character string and broadcast it to all nodes.
	/// @details Intended to be used with receive_broadcast_silent_struct_and_build_pose() to allow all other nodes to receive the broadcast.
	/// @note This function necessarily uses new and delete[].
	void emperor_broadcast_silent_struct( core::io::silent::SilentStructOP ss ) const;

	/// @brief Write out a summary of the jobs completed (node, job index on node, total energy, rmsd, handler path) to the summary tracer.
	///
	void emperor_write_summaries_to_tracer( utility::vector1< SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP > const &summary_list ) const;

	/// @brief Based on the sorted list of summaries, populate a short list of jobs, the results of which will be collected from below for output to disk.
	/// @param[out] summary_shortlist The short list of job summaries populated by this function.
	/// @param[in] summary_full_sorted_list The full list of job summaries collected from below, sorted.
	/// @param[in] output_fraction The fraction of total jobs to collect.
	/// @param[in] select_highest Should we select from the top of the summary list (lowest values) or from the bottom (highest)?
	void emperor_select_best_summaries(
		utility::vector1< SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP > &summary_shortlist,
		utility::vector1< SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP > const &summary_full_sorted_list,
		core::Real const &output_fraction,
		bool const select_highest
	) const;

	/// @brief Write all the collected results from below to disk.
	/// @details Assumes silent output.
	void emperor_write_to_disk( std::string const &output ) const;

	/// ------------- Intermediate Master Methods --------

	/// @brief Is this an intermediate master node?
	/// @details The masters are responsible for distributing jobs to other masters lower in the hierarchy, and/or to slaves, and for
	/// collecting results from masters/slaves lower in the hierarchy and sending them up to the emperor.
	bool i_am_intermediate_master() const;

	/// @brief The jobs done by the intermediate masters during parallel execution.
	/// @details The masters are responsible for distributing jobs to other masters lower in the hierarchy, and/or to slaves, and for
	/// collecting results from masters/slaves lower in the hierarchy and sending them up to the emperor.
	void run_intermediate_master() const;

	/// @brief Relay the jobs received from below, held as a concatenated string, up the hierarchy.
	/// @details Transmission to be received with receive_pose_batch_as_string().
	void intermediate_master_send_poses_as_string_upward( std::string const &results, int const target_node ) const;

	/// ------------- Slave Methods ----------------------

	/// @brief Is this a slave (or worker) node?
	/// @details The slaves receive jobs from higher in the hierarchy, do them, and send results back up the hierarchy.
	bool i_am_slave() const;

	/// @brief The jobs done by the slaves during parallel execution.
	/// @details The slaves receive jobs from higher in the hierarchy, do them, and send results back up the hierarchy.
	/// @note In multi-threaded mode, if threads_per_slave_process_ is greater than 1, then the slaves launch threads
	/// to do the work.  Only the master thread does MPI calls.
	void run_slave() const;

	/// @brief Actually carry out the jobs.  This is where SimpleCycpepPredictApplication is created and invoked.
	///
	void slave_carry_out_njobs( core::Size &njobs_from_above, utility::vector1 < SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP > &jobsummaries, utility::vector1 < core::io::silent::SilentStructOP > &all_output ) const;

#ifdef MULTI_THREADED
	/// @brief Decrement the job counter.  Return true if the job counter was greater than zero.
	/// @details Does this with proper locking to prevent threads from stepping on one another.
	bool slave_decrement_jobcount_multithreaded( core::Size * available_job_count, core::Size &already_completed_job_count, core::Size const jobs_in_this_batch, core::Size const thread_index ) const;

	/// @brief Actually carry out the jobs.  This is where SimpleCycpepPredictApplication is created and invoked.
	/// @details This is the multi-threaded version, which locks the job count to decrement it, and locks the jobsummaries and all_output vectors to add to them.
	void slave_carry_out_njobs_in_thread(
		core::Size const thread_index,
		core::Size * njobs_from_above,
		core::Size const jobs_in_this_batch,
		utility::vector1 < SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP > * jobsummaries,
		utility::vector1 < core::io::silent::SilentStructOP > * all_output,
		core::scoring::ScoreFunctionOP sfxn,
		core::pose::PoseOP native,
		std::string const sequence, /*deliberately passed by copy*/
		core::Size const batch_index //Number of batches that have been sent out on this proc.
	) const;
#endif

	/// @brief Given a list of jobs that have been requested from above, send the corresponding poses up the hierarchy.
	/// @details Throws an error if any jbo was completed on a different node than this slave.
	void slave_send_poses_upward( utility::vector1< SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP > const &requested_jobs, utility::vector1 < core::io::silent::SilentStructOP > const &all_output ) const;

private:
	/// ------------- Data -------------------------------
	/// -------- When you add new data to this class, ----
	/// -------- you must update the copy constructor ----

	/// @brief The rank of THIS process.
	///
	int MPI_rank_;

	/// @brief The total number of processes.
	///
	int MPI_n_procs_;

	/// @brief The number of jobs that have been assigned to this process and completed, if it is a slave process.
	/// @details Must be mutable since it's incremented as jobs are assigned.
	mutable core::Size slave_job_count_;

	/// @brief The default scorefunction to use.
	/// @details The high-hbond version is constructed from this one. 
	/// If necessary, the aa_composition score term will be turned on
	/// in that one; it needn't be turned on in this one.
	core::scoring::ScoreFunctionOP scorefxn_;

	/// @brief The level in the communications hierarchy of this process.
	/// @details One-based: the emperor is level 1, and the slaves are level total_hierarchy_levels_.
	core::Size hierarchy_level_;

	/// @brief The total number of levels in the communication.
	/// @details Note that the slaves are total_hierarchy_levels_, and the emperor is level 1.
	core::Size total_hierarchy_levels_;

	/// @brief The number of processes at each level of the communications hierarchy.
	///
	utility::vector1 < core::Size > procs_per_hierarchy_level_;

	/// @brief The number of jobs sent per batch from this node to child nodes.
	///
	core::Size batchsize_;

	/// @brief The process indices of the children of the current process.
	/// @details Will be empty for slave processes.	
	utility::vector1 < int > my_children_;

	/// @brief The process index of the parent of the current process.
	/// @details Will be 0 for emperor.
	core::Size my_parent_;

	/// @brief The amino acid sequence.
	/// @details Read by emperor and transmitted to all other processes.
	std::string sequence_;

	/// @brief The native pose.
	/// @details Will be null if one is not provided.  Read by emperor and transmitted to all
	/// other processes.
	core::pose::PoseCOP native_;

	/// @brief What criterion should be used to sort solutions?
	/// @details Defaults to energy.  The SIMPLE_CYCPEP_PREDICT_MPI_SORT_TYPE enum is defined in protocols/cyclic_peptide_predict/util.hh.
	SIMPLE_CYCPEP_PREDICT_MPI_SORT_TYPE sort_type_;

	/// @brief When selecting the top N% to output, should we select the highest-scoring (true) or lowest-scoring (false) of the metric
	/// used to sort?
	/// @details Defalut false (lowest-scoring).
	bool select_highest_;

	/// @brief What fraction of total output structures should be written to disk?
	/// @details If less than 1.0, then the top N% of structures are written out based on the sort_type_ sort criterion.
	core::Real output_fraction_;

	/// @brief Filename for silent output.
	///
	std::string output_filename_;

	/// @brief Allowed canonical residues at each position.
	/// @details Map key 0 stores default settings applied to all positions not specified.	
	std::map< core::Size, utility::vector1 < std::string > > allowed_canonicals_;

	/// @brief Allowed noncanonical residues at each position.
	/// @details Map key 0 stores default settings applied to all positions not specified.	
	std::map< core::Size, utility::vector1 < std::string > > allowed_noncanonicals_;

	/// @brief Has an aa_composition setup file ben provided for residues in the L-alpha helix region of
	/// Ramachadran space?
	bool L_alpha_comp_file_exists_;

	/// @brief Has an aa_composition setup file ben provided for residues in the D-alpha helix region of
	/// Ramachadran space?
	bool D_alpha_comp_file_exists_;

	/// @brief Has an aa_composition setup file ben provided for residues in the L-beta strand region of
	/// Ramachadran space?
	bool L_beta_comp_file_exists_;

	/// @brief Has an aa_composition setup file ben provided for residues in the D-beta strand region of
	/// Ramachadran space?
	bool D_beta_comp_file_exists_;

	/// @brief Storage for the composition constraint setup for the L-alpha helix region of Ramachandran space.
	/// @details Cached to prevent repeated read from disk.
	std::string comp_file_contents_L_alpha_;

	/// @brief Storage for the composition constraint setup for the D-alpha helix region of Ramachandran space.
	/// @details Cached to prevent repeated read from disk.
	std::string comp_file_contents_D_alpha_;

	/// @brief Storage for the composition constraint setup for the L-beta strand region of Ramachandran space.
	/// @details Cached to prevent repeated read from disk.
	std::string comp_file_contents_L_beta_;

	/// @brief Storage for the composition constraint setup for the D-beta strand region of Ramachandran space.
	/// @details Cached to prevent repeated read from disk.
	std::string comp_file_contents_D_beta_;

	/// @brief Storage for a bin definition file.
	/// @details Cached to prevent repeated read from disk.
	std::string abba_bins_;

	/// @brief Lambda, the bredth of the Gaussian used for calculating funnel quality (PNear).
	/// @details Read from options system; default 0.5.
	core::Real lambda_;

	/// @brief The Boltzmann temperature, kB*T, used for calculating funnel quality (PNear).
	/// @details Read from options system; default 1.0.
	core::Real kbt_;

#ifdef MULTI_THREADED
// Private member variables only used in multi-threaded build.
private:

	/// @brief The number of threads per slave process.  Setting this to 1 (the default)
	/// reproduces the same behaviour as the non-threaded build.
	core::Size threads_per_slave_proc_;

	/// @brief Lock the results list for read or for write.
	mutable std::mutex results_mutex_;

	/// @brief Lock the available jobs list and next job counter for read or for write.
	mutable std::mutex joblist_mutex_;

	/// @brief Are we using a constant random seed?
	/// @details If false, time is used as seed.
	bool use_const_random_seed_;

	/// @brief The random seed offset to use.
	/// @details Different threads will be offset by 1, and different processes, by the
	/// number of threads.  This value will be added as well.
	int random_seed_offset_;

	/// @brief The random seed to use.
	int random_seed_;

	/// @brief What type of random generator to use?
	std::string rgtype_;

#endif //ifdef MULTI_THREADED

};

} //cyclic_peptide
} //protocols

#endif //INCLUDED_protocols_cyclic_peptide_predict_SimpleCycpepPredictApplication_MPI_hh

#endif //USEMPI
