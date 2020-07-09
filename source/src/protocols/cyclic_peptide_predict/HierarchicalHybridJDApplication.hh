// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide_predict/HierarchicalHybridJDApplication.hh
/// @brief Application-level code for the simple_cycpep_predict app, MPI version.
/// @details  This application predicts structures of simple backbone-cyclized peptides made of alpha-, beta-, or gamma-amino acids (of any chirality)
/// using generalized kinematic closure (GenKIC) for cyclization, and enforcing user-defined requiresments for numbers of mainchain hydrogen bonds.  This
/// version uses MPI for cross-communication with parallel processes.
/// On 4 Aug 2017, work commenced to make this appliction compatible with a job distribution scheme in which the last level of the hierarchy splits
/// jobs over many threads (hierarchical process-level parallelism plus thread-level parallelism).
/// On 29 Oct 2018, this code was moved to the HierarchicalHybridJDApplication base class, from which both the SimpleCycpepPredictApplication_MPI and
/// HelicalBundlePredictApplication_MPI classes derive.
/// On 16 June 2020, this code was updated to replace "emperor"/"master"/"slave", which some people found objectionable.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifdef USEMPI

#ifndef INCLUDED_protocols_cyclic_peptide_predict_HierarchicalHybridJDApplication_hh
#define INCLUDED_protocols_cyclic_peptide_predict_HierarchicalHybridJDApplication_hh

// Unit Headers
#include <protocols/cyclic_peptide_predict/HierarchicalHybridJDApplication.fwd.hh>
#include <protocols/cyclic_peptide_predict/HierarchicalHybridJD_JobResultsSummary.fwd.hh>
#include <protocols/cyclic_peptide_predict/HierarchicalHybridJD_RMSDToBestSummary.fwd.hh>
#include <protocols/cyclic_peptide_predict/HierarchicalHybridJD_SASASummary.fwd.hh>
#include <protocols/cyclic_peptide_predict/HierarchicalHybridJD_PNearToArbitraryStateSummary.fwd.hh>
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
#include <utility/VirtualBase.hh>
#include <utility/vector1.hh>

// Basic Headers
#include <basic/Tracer.hh>

// C++ headers
#include <stdio.h>
#include <time.h>

#ifdef MULTI_THREADED
#include <mutex>
#endif

namespace protocols {
namespace cyclic_peptide_predict {

enum HIERARCHICAL_MPI_COMMUNICATION_TYPE {
	NULL_MESSAGE = 1,
	REQUEST_NEW_JOBS_BATCH_UPWARD,
	REQUEST_TOP_POSE_BCAST, // The director is asking that the top pose be shared with all nodes.
	OFFER_NEW_JOBRESULTSSUMMARY_BATCH_UPWARD,
	OFFER_NEW_RMSD_TO_BEST_SUMMARY_BATCH_UPWARD,
	OFFER_NEW_SASA_SUMMARY_BATCH_UPWARD,
	OFFER_NEW_JOBS_ATTEMPTED_COUNT_UPWARD,
	OFFER_NEW_POSE_BATCH_UPWARD,
	SILENT_STRUCT_TRANSMISSION,
	GIVE_COMPLETION_SIGNAL_UPWARD,
	OFFER_NEW_JOBS_BATCH_DOWNWARD,
	REQUEST_NEW_POSE_BATCH_DOWNWARD,
	REQUEST_SASA_SUMMARIES_DOWNWARD, // The director is asking that all workers transmit SASA summaries upward.

	BEGIN_PNEAR_TO_LOWEST_FRACT_DOWNWARD, // The director is indicating that we're going to compute the PNear values to the lowest fraction of states found.
	SKIP_PNEAR_TO_LOWEST_FRACT_DOWNWARD, // The director is indicating that we're NOT going to compute the PNear values to the lowest fraction of states found (because no states were sampled).
	REQUEST_PNEAR_TO_PARTICULAR_SAMPLE_DOWNWARD, //The director is requesting data for a PNear computation to a particular sample, and is about to transmit the information for the sample.
	END_PNEAR_TO_LOWEST_FRACT_DOWNWARD, //The director is indicating that we're finished computing the PNear values to the lowest fraction of states found.

	HALT_SIGNAL
};

enum HIERARCHICAL_HYBRID_JD_MPI_TAG_TYPE {
	GENERAL_REQUEST=1,
	NEW_JOBS_DOWNWARD,
	JOB_NODE_AND_INDEX_DOWNWARD,
	RESULTS_SUMMARY_UPWARD,
	JOBS_ATTEMPTED_COUNT_UPWARD
};

/// @brief Application-level code for the simple_cycpep_predict app, MPI version.
/// @details  This application predicts structures of simple backbone-cyclized peptides made of alpha-, beta-, or gamma-amino acids (of any chirality)
/// using generalized kinematic closure (GenKIC) for cyclization, and enforcing user-defined requiresments for numbers of mainchain hydrogen bonds.  This
/// version uses MPI for cross-communication with parallel processes.
/// On 4 Aug 2017, work commenced to make this appliction compatible with a job distribution scheme in which the last level of the hierarchy splits
/// jobs over many threads (hierarchical process-level parallelism plus thread-level parallelism).
/// On 29 Oct 2018, this code was moved to the HierarchicalHybridJDApplication base class, from which both the SimpleCycpepPredictApplication_MPI and
/// HelicalBundlePredictApplication_MPI classes derive.
/// On 16 June 2020, this code was updated to replace "emperor"/"master"/"slave", which some people found objectionable.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
class HierarchicalHybridJDApplication : public utility::VirtualBase
{
public:
	/// @brief Constructor
	/// @details Default constructor is explicitly deleted.
	HierarchicalHybridJDApplication() = delete;

	/// @brief Tracer initialization constructor
	/// @details Base class must be initialized to point to derived class tracers.
	HierarchicalHybridJDApplication( basic::Tracer & tracer, basic::Tracer & summary_tracer );

	/// @brief Constructor with options
	///
	HierarchicalHybridJDApplication(
		basic::Tracer & tracer,
		basic::Tracer & summary_tracer,
		int const MPI_rank,
		int const MPI_n_procs,
		core::scoring::ScoreFunctionCOP sfxn_in,
		core::Size const total_hierarchy_levels,
		utility::vector1 < core::Size > const & procs_per_hierarchy_level,
		utility::vector1 < core::Size > const &batchsize_per_level,
		std::string const &sort_type,
		bool const select_highest,
		core::Real const output_fraction,
		std::string const &output_filename,
		core::Real const lambda,
		core::Real const kbt,
		bool const compute_rmsd_to_lowest,
		core::Real const compute_pnear_to_this_fract,
		bool const compute_sasa_metrics,
		core::Size const threads_per_worker_proc //Only used in multi-threaded build.
	);

	/// @brief Explicit virtual destructor.
	///
	virtual ~HierarchicalHybridJDApplication();

	/// @brief Explicit copy constructor.
	///
	HierarchicalHybridJDApplication( HierarchicalHybridJDApplication const &src );

	/// @brief Clone function: make a copy of this object and return an owning pointer to the copy.
	virtual HierarchicalHybridJDApplicationOP clone() const = 0;

#ifdef MULTI_THREADED
	/// @brief Initialize private member variables related to thread random seeds from the options
	/// system.  Does nothing if this isn't a multi-threaded compilation.
	void init_random_info_from_options();
#endif

	/// @brief Actually run the application.
	/// @details On worker nodes, this creates an instance of the relevant application and runs that.  Nodes higher
	/// in the communications hierarchy are just involved in sending out jobs and pulling in results.
	/// @note The run() function is nonconst due to some setup steps that it performs.  It then calls const run functions
	/// for director, manager, and worker nodes.
	void run();

public:
	/// ------------- Public Methods ---------------------

	/// @brief Set the sort type, by string.
	///	
	void set_sort_type( std::string const &sort_type );

	/// @brief Set the sort type, by enum.
	///	
	void set_sort_type( HIERARCHICAL_HYBRID_JD_MPI_SORT_TYPE const sort_type );

	/// @brief Set the ouput fraction.
	/// @details Checks that this is between 0 and 1.
	void set_output_fraction( core::Real const &val );

protected:

	/// @brief Get the protocol-specific settings.
	/// @details The director reads these from disk and broadcasts them to all other nodes.  This function should be called from all nodes;
	/// it figures out which behaviour it should be performing.
	/// @note Pure virtual.  Must be implemented by derived classes.
	virtual void get_protocol_specific_settings() = 0;

	/// @brief Create an instance of the appropriate app class, and carry out N jobs on a single process.
	/// @details This code is called in a single thread in multi-threaded mode, and is used in the single-threaded version too.
	/// @note Pure virutal function must be implemented by derived classes.
	virtual
	void derived_worker_carry_out_n_jobs(
		core::Size const njobs_from_above,
		utility::vector1 < HierarchicalHybridJD_JobResultsSummaryOP > &jobsummaries,
		utility::vector1 < core::io::silent::SilentStructOP > &all_output,
		core::scoring::ScoreFunctionOP sfxn,
		core::pose::PoseCOP native,
		std::string const &sequence	
	) const = 0;

	/// @brief Compute the RMSD between a pose and a reference pose.
	/// @details Must be implemented by derived classes, since this might be done differently for
	/// different classes of molecule.
	virtual core::Real
	derived_worker_compute_rmsd(
		core::pose::Pose const & pose,
		core::pose::Pose const & reference_pose,
		std::string const &sequence
	) const = 0;

	/// @brief Compute SASA metrics for this pose and bundle them in a SASASummary.
	/// @details Although derived classes MAY override this, they needn't.  The default behaviour
	/// is just to use SASA metrics, which should be pretty general.
	virtual HierarchicalHybridJD_SASASummaryOP
	generate_sasa_summary(
		core::pose::Pose const &pose,
		HierarchicalHybridJD_JobResultsSummaryCOP jobsummary
	) const;

	/// @brief Get a const owning pointer to the native pose.
	/// @details Returns nullptr if native_ == nullptr.
	inline core::pose::PoseCOP native() const { return native_; }

	/// @brief Get an owning pointer to the scorefunction.
	/// @details Will be nullptr if the scorefunction isn't defined.
	inline core::scoring::ScoreFunctionOP scorefxn() const { return scorefxn_; }

	/// @brief Get the sequence.
	inline std::string const &sequence() const { return sequence_; }

	/// @brief Get the MPI rank.
	inline int MPI_rank() const { return MPI_rank_; }

	/// @brief Get the already-completed job count.
	inline core::Size worker_job_count() const { return worker_job_count_; }

private:
	/// ------------- Methods ----------------------------

	/// @brief Reads a FASTA file and returns a string of space-separated full names for the sequence.
	/// @details TRIGGERS A READ FROM DISK.
	std::string sequence_from_fasta( std::string const &fasta_file ) const;

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
	/// @details The director reads this from disk and broadcasts it to all other nodes.  This function should be called from all nodes;
	/// it figures out which behaviour it should be performing.
	void get_sequence();

	/// @brief Get the native structure of the peptide we're going to predict (if the user has specified one with the -in:file:native flag).
	/// @details The director reads this from disk and broadcasts it to all other nodes.  This function should be called from all nodes;
	/// it figures out which behaviour it should be performing.
	void get_native();

protected:

	/// @brief Given a map of indicies to lists of residue names, broadcast it to all MPI ranks.
	void broadcast_res_list( std::map< core::Size, utility::vector1 < std::string > > &res_list) const;

private:

	/// @brief The director sends a message to everyone letting them know it's time to start.
	/// @details Following the go signal, workers send requests for jobs upward.
	void go_signal() const;

	/// @brief The director sends a message to everyone letting them know it's time to stop.
	/// @details Following the stop signal, HierarchicalHybridJDApplication::run() terminates.
	void stop_signal() const;

	/// @brief Any node can wait for a node, above or below in the hierarchy, to send it some sort of request.
	/// @details Only messags with tag GENERAL_REQUEST.
	/// @param[out] requesting_node The node from which the request came.
	/// @param[out] message The type of request received.
	void wait_for_request( int &requesting_node, HIERARCHICAL_MPI_COMMUNICATION_TYPE &message ) const;

	/// @brief Send a signal to stop job distribution to a set of intermediate managers.
	///
	void send_halt_signal( utility::vector1 < int > const &ranks_to_target ) const;

	/// @brief Any node can send some sort of request to a specific node in the hierarchy (parent or child).
	/// @details Sends messags with tag GENERAL_REQUEST.
	/// @param[in] target_node The node to which we're sending the request.
	/// @param[in] message The type of request to send.
	void send_request( int const target_node, HIERARCHICAL_MPI_COMMUNICATION_TYPE const message ) const;

	/// @brief Send an integer to a child indicating that the child is now holding this number of jobs.
	/// @details Sends messags with tag NEW_JOBS_DOWNWARD.  Sends zero if no jobs remain to send.
	/// @param[in] target_node The rank of the node to which we're sending the request.
	/// @param[in,out] structures_remaining_to_send The number of structures remaining in the send buffer.  Decremented by this operation to a minimum of zero.
	/// @param[in] number_to_send How many jobs should I send down?
	void send_new_jobs_downward( int const target_node, core::Size &structures_remaining_to_send, core::Size const number_to_send ) const;

	/// @brief Receive an integer from my parent indicating that I should add N jobs to my set to carry out or to pass to my children.
	/// @details Recieves message with tag NEW_JOBS_DOWNARD.  Receives zero if no jobs remain in parent to send.
	/// @param[in,out] njobs The number of jobs held on this node that are to be done.  Incremented by this function with however many are received from above.
	void receive_njobs_from_above( core::Size &njobs ) const;

	/// @brief Non-originating nodes must call this when the originating node calls broadcast_silent_struct_from_this_node.
	/// @details This will build a pose and return an owning pointer to it.
	/// @note If skip_pose_build is true, we don't bother to convert the structure into a pose; we just participate in the
	/// broadcast.  In that case, this function returns nullptr.
	core::pose::PoseCOP receive_broadcast_silent_struct_and_build_pose( int node_to_receive_from = 0, bool const skip_pose_build = false ) const;

	/// @brief Convert a vector of silent structs into a character string and send it to a node.
	/// @details Intended to be used with receive_pose_batch_as_string() to allow another node to receive the transmission.
	/// Message is tagged with SILENT_STRUCT_TRANSMISSION.
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
	void receive_and_sort_job_summaries( utility::vector1< HierarchicalHybridJD_JobResultsSummaryOP > &original_summary_list, int const originating_node, bool const append_to_handler_list ) const;

	/// @brief Send a list of job summaries.
	/// @details To be used in conjunction with receive_and_sort_job_summaries().  Sending and receiving procs must send messages to synchronize, first.
	void send_job_summaries( utility::vector1< HierarchicalHybridJD_JobResultsSummaryOP > const &summary_list, int const target_node ) const;

	/// @brief Given a short list of job summaries, split the list by the index of the node that I'd have to send the request to, and send requests for full poses down the hierarchy.
	/// @details Throws an error if any of the jobs in the list cannot be reached by sending a request by way of my_children_.  To be used with receive_pose_requests_from_above().
	/// @param[in] summary_shortlist The list of job summaries to split up and send downard.
	/// @param[out] children_receiving_nonzero_requests The number of children receiving a request for at least one pose.
	void request_poses_from_below( utility::vector1< HierarchicalHybridJD_JobResultsSummaryOP > const &summary_shortlist, core::Size &children_receiving_nonzero_requests ) const;

	/// @brief Recieve a short list of job summaries from above.
	/// @param[out] summary_shortlist The job summaries list, cleared and populated by this function.
	void receive_pose_requests_from_above( utility::vector1< HierarchicalHybridJD_JobResultsSummaryOP > &summary_shortlist ) const;

	/// @brief The director is sending out a request that the worker that produced the top pose broadcast
	/// that pose to all other nodes.  All nodes participate in the broadcast.  At the end of this operation,
	/// all nodes know:
	/// - The node index that produced the top structure.
	/// - The index of the top structure on that node index.
	/// - The top structure (the pose), as a binary silent structure.
	/// They can then compare the top structure to all of their structures and compute an RMSD for each.
	/// @note If broadcast_no_poses_found is true, the director announces to all processes that there is no best
	/// pose with which to compare.
	/// @returns TRUE for FAILURE (Director reports no poses found), FALSE for SUCCESS.
	bool
	top_pose_bcast(
		HierarchicalHybridJD_JobResultsSummaryOP top_summary=nullptr,
		core::io::silent::SilentStructOP top_pose_silentstruct=nullptr,
		bool const broadcast_no_poses_found=false
	) const;

protected:

	/// @brief Given a string on the director node, send it to all nodes.
	/// @details The "mystring" string is the input on the director node, and the output on all other nodes.
	/// @note Protected, not private.
	void broadcast_string_from_director( std::string &mystring ) const;

	/// @brief Given a string on a given node, send it to all nodes.
	/// @details The "mystring" string is the input on the originating node, and the output on all other nodes.
	/// @note Protected, not private.
	void broadcast_string_from_node( std::string &mystring, int const originating_node_index ) const;

	/// @brief Given a vector of RMSDs to best summaries, send them up the hierarchy.
	/// @details Made for use with receive_and_sort_all_rmsd_to_best_summaries(), which calls receive_rmsd_to_best_summaries_from_below().
	void
	send_rmsds_to_best_summaries_upward(
		utility::vector1< HierarchicalHybridJD_RMSDToBestSummaryOP > const &rmsds_to_best_summaries,
		int const receiving_node
	) const;

	/// @brief From one of the nodes below me in the next level down, receive one vector of RMSD-to-best summaries.
	/// @details Made for use with send_rmsds_to_best_summaries_upward().  THe new_rmsd_to_best_summaries vector is cleared
	/// and populated by this operation.
	void
	receive_rmsd_to_best_summaries_from_below(
		utility::vector1< HierarchicalHybridJD_RMSDToBestSummaryOP > &new_rmsd_to_best_summaries,
		int const requesting_node,
		bool const append_to_handler_list
	) const;

	/// @brief From all of the nodes below me in the next level down, receive RMSD-to-best summaries.  Then sort them
	/// based on the already-received job summaries.
	/// @details Made for use with send_rmsds_to_best_summaries_upward().  THe rmsds_to_best_summaries vector is cleared
	/// and populated by this operation.
	void
	receive_and_sort_all_rmsd_to_best_summaries (
		utility::vector1< HierarchicalHybridJD_RMSDToBestSummaryOP > &rmsds_to_best_summaries,
		utility::vector1< HierarchicalHybridJD_JobResultsSummaryOP > const &job_summary_list,
		bool const append_to_handler_list
	) const;

	/// @brief Given a vector of SASA summaries, send them to a target node.
	/// @details This function complements receive_sasa_summaries_from_below().
	void
	send_sasa_summaries_upward(
		utility::vector1< HierarchicalHybridJD_SASASummaryOP > const & sasa_summaries,
		int const receiving_node
	) const;

	/// @brief A child node has sent a batch of SASA summaries.  Receive them.
	/// @details This function complements send_sasa_summaries_upward().
	void
	receive_sasa_summaries_from_below(
		utility::vector1< HierarchicalHybridJD_SASASummaryOP > & sasa_summaries,
		int const requesting_node,
		bool const append_to_handler_list
	) const;

	/// @brief Recieve all SASA summaries from all children, and sort them into the order that jobsummaries is in.
	/// @details This is intended for use with send_sasa_summaries_upward().
	void
	receive_and_sort_all_sasa_summaries(
		utility::vector1< HierarchicalHybridJD_SASASummaryOP > & sasa_summaries,
		utility::vector1< HierarchicalHybridJD_JobResultsSummaryOP > const & job_summary_list,
		bool const append_to_handler_list
	) const;

	/// ------------- Director Methods --------------------

	/// @brief Is this an director (root) node?
	/// @details The director is responsible for sending out and retrieving all jobs, and for all file output.
	/// @note Protected, not private.
	bool i_am_director() const;

private:

	/// @brief The jobs done by the director during parallel execution.
	/// @details The director is responsible for sending out and retrieving all jobs, and for all file output.
	void run_director() const;

	/// @brief Convert a silent struct into a character string and broadcast it to all nodes.
	/// @details Intended to be used with receive_broadcast_silent_struct_and_build_pose() to allow all other nodes to receive the broadcast.
	void broadcast_silent_struct_from_this_node( core::io::silent::SilentStructOP ss ) const;

	/// @brief Write out a summary of the jobs completed (node, job index on node, total energy, rmsd, handler path) to the summary tracer.
	/// @details The RMSD to best pose vector will only be populated if the -compute_rmsd_to_lowest option is used.
	void
	director_write_summaries_to_tracer(
		utility::vector1< HierarchicalHybridJD_JobResultsSummaryOP > const &summary_list,
		utility::vector1< HierarchicalHybridJD_RMSDToBestSummaryOP > const &rmsds_to_best_pose,
		utility::vector1< HierarchicalHybridJD_SASASummaryOP > const &sasa_summaries,
		utility::vector1< HierarchicalHybridJD_PNearToArbitraryStateSummaryCOP > const & pnear_to_lowest_fract_summaries
	) const;

	/// @brief Based on the sorted list of summaries, populate a short list of jobs, the results of which will be collected from below for output to disk.
	/// @param[out] summary_shortlist The short list of job summaries populated by this function.
	/// @param[in] summary_full_sorted_list The full list of job summaries collected from below, sorted.
	/// @param[in] output_fraction The fraction of total jobs to collect.
	/// @param[in] select_highest Should we select from the top of the summary list (lowest values) or from the bottom (highest)?
	void director_select_best_summaries(
		utility::vector1< HierarchicalHybridJD_JobResultsSummaryOP > &summary_shortlist,
		utility::vector1< HierarchicalHybridJD_JobResultsSummaryOP > const &summary_full_sorted_list,
		core::Real const &output_fraction,
		bool const select_highest
	) const;

	/// @brief Write all the collected results from below to disk.
	/// @details Assumes silent output.
	void director_write_to_disk( std::string const &output ) const;

	/// @brief The director is asking for SASA metrics to be sent up the hierarchy.
	void director_send_request_for_sasa_summaries_downward() const;

	/// @brief The director is asking for data with which to compute PNear values from below.
	/// @details The steps are:
	///     - Sends a request for each state in turn from the director node to all worker nodes.
    ///     - Each worker node broadcasts that state to all other worker nodes.
    ///     - All worker nodes compute RMSD to that state.
    ///     - Director collects RMSDs up the hierarchy and carries out PNear calculation.
    ///     - Repeat for each relevant state.
	/// @note The pnears_to_lowest_fract vector is cleared and populated by this operation.
 	void
	director_compute_pnear_to_lowest_fract(
		core::Real const & fraction,
		utility::vector1< HierarchicalHybridJD_JobResultsSummaryOP > const & results_summary,
		utility::vector1< HierarchicalHybridJD_PNearToArbitraryStateSummaryCOP > & pnears_to_lowest_fract
	) const;

	/// ------------- Intermediate Manager Methods --------

protected:

	/// @brief Is this an intermediate manager node?
	/// @details The managers are responsible for distributing jobs to other managers lower in the hierarchy, and/or to workers, and for
	/// collecting results from managers/workers lower in the hierarchy and sending them up to the director.
	bool i_am_intermediate_manager() const;

private:

	/// @brief The jobs done by the intermediate managers during parallel execution.
	/// @details The managers are responsible for distributing jobs to other managers lower in the hierarchy, and/or to workers, and for
	/// collecting results from managers/workers lower in the hierarchy and sending them up to the director.
	void run_intermediate_manager() const;

	/// @brief Relay the jobs received from below, held as a concatenated string, up the hierarchy.
	/// @details Transmission to be received with receive_pose_batch_as_string().
	void intermediate_manager_send_poses_as_string_upward( std::string const &results, int const target_node ) const;

	/// @brief Receive a request for SASA summaries from above, and send it to all children.
	/// @details This function expects that the only possible message that can be received
	/// at this point is the request for SASA summaries!
	void intermediate_manager_relay_request_for_sasa_summaries_downward() const;

	/// @brief Send requests for PNear data for the lowest-energy N% of structures down
	/// the hierarchy, and facilitate results going up the hierarchy.
	void
	intermediate_manager_compute_pnear_to_lowest_fract(
		utility::vector1 < HierarchicalHybridJD_JobResultsSummaryOP > const & jobsummaries
	) const;

	/// ------------- Worker Methods ----------------------

protected:

	/// @brief Is this a worker (or worker) node?
	/// @details The workers receive jobs from higher in the hierarchy, do them, and send results back up the hierarchy.
	bool i_am_worker() const;

private:

	/// @brief The jobs done by the workers during parallel execution.
	/// @details The workers receive jobs from higher in the hierarchy, do them, and send results back up the hierarchy.
	/// @note In multi-threaded mode, if threads_per_worker_process_ is greater than 1, then the workers launch threads
	/// to do the work.  Only the manager thread does MPI calls.
	void run_worker() const;

	/// @brief Actually carry out the jobs.  This is where SimpleCycpepPredictApplication is created and invoked.
	///
	void worker_carry_out_njobs( core::Size &njobs_from_above, utility::vector1 < HierarchicalHybridJD_JobResultsSummaryOP > &jobsummaries, utility::vector1 < core::io::silent::SilentStructOP > &all_output ) const;

#ifdef MULTI_THREADED
	/// @brief Decrement the job counter.  Return true if the job counter was greater than zero.
	/// @details Does this with proper locking to prevent threads from stepping on one another.
	bool worker_decrement_jobcount_multithreaded( core::Size * available_job_count, core::Size &already_completed_job_count, core::Size const jobs_in_this_batch, core::Size const thread_index ) const;

	/// @brief Actually carry out the jobs.  This is where SimpleCycpepPredictApplication is created and invoked.
	/// @details This is the multi-threaded version, which locks the job count to decrement it, and locks the jobsummaries and all_output vectors to add to them.
	void worker_carry_out_njobs_in_thread(
		core::Size const thread_index,
		core::Size * njobs_from_above,
		core::Size const jobs_in_this_batch,
		utility::vector1 < HierarchicalHybridJD_JobResultsSummaryOP > * jobsummaries,
		utility::vector1 < core::io::silent::SilentStructOP > * all_output,
		core::scoring::ScoreFunctionOP sfxn,
		core::pose::PoseOP native,
		std::string const sequence, /*deliberately passed by copy*/
		core::Size const batch_index //Number of batches that have been sent out on this proc.
	) const;
#endif

	/// @brief If we're computing the RMSDs to the very best pose, do so.
	/// @details This function clears and populates the rmsds_to_best_summaries vector, ensuring that RMSDs to the
	/// pose represented by top_pose_silentstruct are computed in the order that matches jobsummaries.
	void
	worker_compute_sorted_rmsds_to_best(
		core::io::silent::SilentStruct const &top_pose_silentstruct,
		utility::vector1< HierarchicalHybridJD_RMSDToBestSummaryOP > & rmsds_to_best_summaries,
		utility::vector1< HierarchicalHybridJD_JobResultsSummaryOP > const &jobsummaries,
		utility::vector1< core::io::silent::SilentStructOP > const &poses_from_this_worker
	) const;

	/// @brief Wait for requests from the director for RMSDs to a given structure, participate in a broadcast of
	/// that structure, and compute RMSDs for all output to that structure to send up the hierarchy.
	void
	worker_compute_pnear_to_lowest_fract(
		utility::vector1< HierarchicalHybridJD_JobResultsSummaryCOP > const & jobsummaries,
		utility::vector1 < core::io::silent::SilentStructCOP > const & all_output 
	) const;

	/// @brief Given a list of jobs that have been requested from above, send the corresponding poses up the hierarchy.
	/// @details Throws an error if any jbo was completed on a different node than this worker.
	void worker_send_poses_upward( utility::vector1< HierarchicalHybridJD_JobResultsSummaryOP > const &requested_jobs, utility::vector1 < core::io::silent::SilentStructOP > const &all_output ) const;


	/// @brief Receive a request for SASA metrics from above.
	/// @details This must be the ONLY type of request that this worker can receive at this time!
	void worker_receive_request_for_sasa_summaries() const;

	/// @brief Generate SASA metrics, and sort these for transmission up the hierarchy.
	/// @details Sort order matches the order of jobsummaries.
	void
	worker_generate_and_sort_sasa_summaries(
		utility::vector1< HierarchicalHybridJD_SASASummaryOP > & sasa_summaries,
		utility::vector1< HierarchicalHybridJD_JobResultsSummaryOP > const &jobsummaries,
		utility::vector1< core::io::silent::SilentStructOP > const &poses_from_this_worker
	) const;

private:
	/// ------------- Data -------------------------------
	/// -------- When you add new data to this class, ----
	/// -------- you must update the copy constructor ----

	/// @brief The tracer from the derived class.
	basic::Tracer & derivedTR_;

	/// @brief The summary tracer from the derived class.
	basic::Tracer & derivedTR_summary_;

	/// @brief The rank of THIS process.
	///
	int MPI_rank_;

	/// @brief The total number of processes.
	///
	int MPI_n_procs_;

	/// @brief The number of jobs that have been assigned to this process and completed, if it is a worker process.
	/// @details Must be mutable since it's incremented as jobs are assigned.
	mutable core::Size worker_job_count_;

	/// @brief The default scorefunction to use.
	/// @details The high-hbond version is constructed from this one. 
	/// If necessary, the aa_composition score term will be turned on
	/// in that one; it needn't be turned on in this one.
	core::scoring::ScoreFunctionOP scorefxn_;

	/// @brief The level in the communications hierarchy of this process.
	/// @details One-based: the director is level 1, and the workers are level total_hierarchy_levels_.
	core::Size hierarchy_level_;

	/// @brief The total number of levels in the communication.
	/// @details Note that the workers are total_hierarchy_levels_, and the director is level 1.
	core::Size total_hierarchy_levels_;

	/// @brief The number of processes at each level of the communications hierarchy.
	///
	utility::vector1 < core::Size > procs_per_hierarchy_level_;

	/// @brief The number of jobs sent per batch from this node to child nodes.
	///
	core::Size batchsize_;

	/// @brief The process indices of the children of the current process.
	/// @details Will be empty for worker processes.	
	utility::vector1 < int > my_children_;

	/// @brief The process index of the parent of the current process.
	/// @details Will be 0 for director.
	core::Size my_parent_;

	/// @brief The amino acid sequence.
	/// @details Read by director and transmitted to all other processes.
	std::string sequence_;

	/// @brief The native pose.
	/// @details Will be null if one is not provided.  Read by director and transmitted to all
	/// other processes.
	core::pose::PoseCOP native_;

	/// @brief What criterion should be used to sort solutions?
	/// @details Defaults to energy.  The HIERARCHICAL_HYBRID_JD_MPI_SORT_TYPE enum is defined in protocols/cyclic_peptide_predict/util.hh.
	HIERARCHICAL_HYBRID_JD_MPI_SORT_TYPE sort_type_;

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

	/// @brief Lambda, the bredth of the Gaussian used for calculating funnel quality (PNear).
	/// @details Read from options system; default 0.5.
	core::Real lambda_;

	/// @brief The Boltzmann temperature, kB*T, used for calculating funnel quality (PNear).
	/// @details Read from options system; default 1.0.
	core::Real kbt_;

	/// @brief If true, the RMSD to the lowest-energy state found is computed.  False by default.
	bool compute_rmsd_to_lowest_ = false;

	/// @brief If nonzero, the PNear to the lowest energy N% of states is computed.  If set
	/// to zero, this doesn't happen.
	core::Real compute_pnear_to_this_fract_ = 0.0;

	/// @brief If true, sasa, polar sasa, and hydrophobic sasa are computed for each structure and for
	/// the ensemble.  False by default.
	bool compute_sasa_metrics_ = false;

#ifdef MULTI_THREADED
// Private member variables only used in multi-threaded build.
private:

	/// @brief The number of threads per worker process.  Setting this to 1 (the default)
	/// reproduces the same behaviour as the non-threaded build.
	core::Size threads_per_worker_proc_;

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

#endif //INCLUDED_protocols_cyclic_peptide_predict_HierarchicalHybridJDApplication_hh

#endif //USEMPI
