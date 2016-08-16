// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   devel/mmt_msd/MMTDriver.hh
/// @brief  declaration for class MMTDriver
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_devel_mmt_msd_MMTDriver_HH
#define INCLUDED_devel_mmt_msd_MMTDriver_HH

// Unit headers
#include <devel/mmt_msd/MMTDriver.fwd.hh>

// Package headers
#include <devel/mmt_msd/MMTReceiver.fwd.hh>

// Core headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/ResfileReader.fwd.hh>

// Protocols headers
#include <protocols/genetic_algorithm/GeneticAlgorithm.fwd.hh>
#include <protocols/genetic_algorithm/Entity.fwd.hh>
#include <protocols/genetic_algorithm/EntityRandomizer.fwd.hh>
#include <protocols/multistate_design/MultiStatePacker.fwd.hh>

#include <protocols/pack_daemon/DynamicAggregateFunction.fwd.hh>
#include <protocols/pack_daemon/MultistateFitnessFunction.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/io/FileContentsMap.fwd.hh>
#include <utility/vector1.hh>

// C++ headers
#include <algorithm>
#include <list>
#include <map>
#include <set>
#include <string>

namespace devel {
namespace mmt_msd {

class PackingJobRecord {
public:
	typedef std::pair< core::Size, core::Real > npd_property_and_value;
	typedef std::list< npd_property_and_value > npd_properties;

public:
	PackingJobRecord();
	void state_index( core::Size setting );
	void sequence_string( std::string setting );
	void node_run_on( core::Size setting );
	void energy( core::Real setting );
	void running_time( core::Real setting );
	void npd_props( npd_properties const & setting );

	core::Size state_index() const;
	std::string const & sequence_string() const;
	core::Size node_run_on() const;
	core::Real energy() const;
	core::Real running_time() const;
	npd_properties const & npd_props() const;

	//bool operator < ( PackingJobRecord const & ) const;

private:

	core::Size  state_index_;
	std::string sequence_string_;
	core::Size  node_run_on_;
	core::Real  energy_;
	core::Real  running_time_;
	npd_properties npd_props_;
};

struct WorkerNodeStats {
public:
	core::Size mpi_id_;
	core::Size max_n_threads_;
	core::Size curr_n_running_threads_;
};

class OneGenerationJobInfo : public utility::pointer::ReferenceCount
{
public:
	typedef std::pair< core::Size, core::Size > SeqAndStateInds;
	typedef SeqAndStateInds JobID;

public:
	OneGenerationJobInfo(
		protocols::genetic_algorithm::GeneticAlgorithmBase const & ga,
		protocols::pack_daemon::DynamicAggregateFunction const & daf
	);

	static JobID job_id( core::Size sequence_index, core::Size state_index );

	bool unassigned_jobs_remain() const;
	bool unfinished_jobs_outstanding() const;

	core::Size n_new_sequences() const { return n_new_sequences_; }

	utility::vector1< core::Size >  const & new_seq_inds() const { return new_seq_inds_; }
	utility::vector1< std::string > const & sequences() const { return sequences_; }
	utility::vector1< core::Size >  const & full_seqinds_2_newseqinds() const { return full_seqinds_2_newseqinds_; }
	std::map< std::string, core::Size > const & sequence_to_index_map() const { return sequence_to_index_map_; }
	utility::vector1< utility::vector1< core::Real > > job_completed() const { return job_completed_; }
	std::list< JobID > const & job_order() const { return job_order_; }
	std::set< JobID > const & work_outstanding() const { return work_outstanding_; }

	JobID pop_job();
	void push_outstanding( JobID const & job );
	void remove_outstanding( JobID const & job );

private:
	core::Size n_new_sequences_;

	utility::vector1< core::Size >  new_seq_inds_;
	utility::vector1< std::string > sequences_;
	utility::vector1< core::Size >  full_seqinds_2_newseqinds_;
	std::map< std::string, core::Size > sequence_to_index_map_;
	utility::vector1< utility::vector1< core::Real > > job_completed_;
	std::list< JobID > job_order_;
	std::set< JobID > work_outstanding_;
};

class JobsForSequence : public utility::pointer::ReferenceCount
{
public:
	JobsForSequence( core::Size n_states, core::Size n_npd_properties );
	~JobsForSequence();

	void entity( protocols::genetic_algorithm::EntityOP setting );
	protocols::genetic_algorithm::EntityOP entity();

	PackingJobRecord & job_record( core::Size state_index );
	PackingJobRecord const & job_record ( core::Size state_index ) const;

	void finalize_state_energies_and_npd_properties();

	utility::vector1< core::Real > const & state_energies() const;
	utility::vector1< core::Real > const & npd_properties() const;

private:
	utility::vector1< PackingJobRecord > jobs_;
	protocols::genetic_algorithm::EntityOP entity_;
	utility::vector1< core::Real > state_energies_;
	utility::vector1< core::Real > npd_properties_;
};

class MMTDriver : public utility::pointer::ReferenceCount
{
public:
	typedef utility::vector1< JobsForSequenceOP > JobsForGeneration;
	typedef utility::vector1< WorkerNodeStats > WorkerNodes;
	typedef std::map< std::string, JobsForSequenceOP > SavedJobsForSequence;
public:
	MMTDriver();
	~MMTDriver();

	void set_nworkers( core::Size setting );
	void set_ngenerations( core::Size setting );
	void set_pop_size( core::Size setting );
	void set_frac_by_recomb( core::Real setting );
	void set_randomizer( protocols::genetic_algorithm::PositionSpecificRandomizerOP setting );
	void set_sfxn( core::scoring::ScoreFunction const & sfxn );
	void set_file_contents( utility::io::FileContentsMapOP file_contents );
	void set_n_results_to_output( core::Size setting );
	void set_daf_fname( std::string const & setting );
	void set_entity_resfile_fname( std::string const & setting );

	void setup();
	void run();

	std::map< std::string, std::string > retrieve_optimal_solutions();
	void write_optimal_solutions_to_disk();

private:

	void initial_handshake();
	void main_optimization_loop();
	bool optimize_generation();
	void send_new_job_to_node( core::Size communicating_node );

	void send_sequence_to_node( core::Size communicating_node, std::string const & seq ) const;

	void send_state_info_to_node(
		core::Size node_index,
		core::Size state_index
	);

	void
	receive_completed_job(
		core::Size communicating_node
	);

	void evaluate_entity_fitnesses();

	void
	instruct_receivers_to_keep_job_data_for_entity(
		core::Size seq_index
	);

	void
	instruct_receivers_to_drop_old_job_data_for_entity(
		protocols::genetic_algorithm::EntityOP entity
	);

	std::string
	retrieve_pdb_from_node(
		core::Size node,
		core::Size state_index,
		std::string const & sequence
	);

private:
	protocols::genetic_algorithm::GeneticAlgorithmBaseOP ga_;
	protocols::genetic_algorithm::PositionSpecificRandomizerOP ga_randomizer_;

	core::pack::task::ResfileContentsOP entity_resfile_contents_;
	core::pack::task::PackerTaskOP entity_task_;
	core::Size num_entities_;

	protocols::pack_daemon::DynamicAggregateFunctionOP daf_;

	protocols::pack_daemon::TopEntitySetOP top_entities_;

	core::Size n_workers_;
	WorkerNodes node_stats_;

	utility::io::FileContentsMapOP file_contents_;

	core::scoring::ScoreFunctionOP sfxn_;

	OneGenerationJobInfoOP this_gen_work_;
	JobsForGeneration this_gen_results_;
	SavedJobsForSequence top_jobs_archive_;

	std::string entity_resfile_fname_;
	std::string daf_fname_;

	core::Size ga_pop_size_;
	core::Size ga_num_to_propagate_;
	core::Size ga_ngenerations_;
	core::Real ga_frac_by_recomb_;
	core::Size n_top_results_to_output_;


	// I imagine I might add quite a lot more data here to
	// enhance the decision making process on what jobs to
	// assign when.

};

}
}

#endif
