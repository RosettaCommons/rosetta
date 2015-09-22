// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/mmt_msd/MMTDriver.cc
/// @brief  Implementation for class MMTDriver
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <devel/mmt_msd/MMTDriver.hh>

// Package headers
#include <devel/mmt_msd/MMTReceiver.hh>

// Core headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/import_pose/import_pose.hh>
#include <core/import_pose/import_pose_options.hh>

// Protocols headers
#include <protocols/genetic_algorithm/GeneticAlgorithm.hh>
#include <protocols/genetic_algorithm/Entity.hh>
#include <protocols/genetic_algorithm/Mutate1Randomizer.hh>
#include <protocols/multistate_design/MultiStatePacker.hh>

#include <protocols/pack_daemon/DynamicAggregateFunction.hh>
#include <protocols/pack_daemon/MultistateFitnessFunction.hh>
#include <protocols/pack_daemon/util.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/mpi_util.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/FileContentsMap.hh>

#ifdef MULTI_THREADED
#ifdef CXX11
// C++11 headers
#include <chrono>
#endif
#endif

namespace devel {
namespace mmt_msd {

static THREAD_LOCAL basic::Tracer TR( "devel.mmt_msd.MMTDriver" );

std::string
sequence_from_entity(
	protocols::genetic_algorithm::Entity const & entity
)
{
	using protocols::multistate_design::PosType;
	using protocols::multistate_design::PosTypeCOP;
	core::Size length = entity.traits().size();
	std::string sequence;
	sequence.reserve( length );
	for ( core::Size ii = 1; ii <= length; ++ii ) {
		protocols::genetic_algorithm::EntityElementCOP iiee = entity.traits()[ ii ];
		PosTypeCOP ii_pos_type = utility::pointer::dynamic_pointer_cast< protocols::multistate_design::PosType const > ( iiee );
		if ( ! ii_pos_type ) {
			throw utility::excn::EXCN_Msg_Exception( "Position " + utility::to_string( ii ) + " of entity " + entity.to_string() + " could not be dynamically casted to PosType" );
		}
		sequence.push_back( core::chemical::oneletter_code_from_aa( ii_pos_type->type() ));
	}
	return sequence;
}

PackingJobRecord::PackingJobRecord()
{}

void PackingJobRecord::state_index( core::Size setting )
{
	state_index_ = setting;
}

void PackingJobRecord::sequence_string( std::string setting )
{
	sequence_string_ = setting;
}

void PackingJobRecord::node_run_on( core::Size setting )
{
	node_run_on_ = setting;
}

void PackingJobRecord::energy( core::Real setting )
{
	energy_ = setting;
}

void PackingJobRecord::running_time( core::Real setting )
{
	running_time_ = setting;
}

void PackingJobRecord::npd_props( npd_properties const & setting )
{
	npd_props_ = setting;
}

core::Size PackingJobRecord::state_index() const { return state_index_; }

std::string const & PackingJobRecord::sequence_string() const { return sequence_string_; }

core::Size PackingJobRecord::node_run_on() const { return node_run_on_; }

core::Real PackingJobRecord::energy() const { return energy_; }

core::Real PackingJobRecord::running_time() const { return running_time_; }

PackingJobRecord::npd_properties const &
PackingJobRecord::npd_props() const { return npd_props_; }

//bool operator < ( PackingJobRecord const & ) const{}


OneGenerationJobInfo::OneGenerationJobInfo(
	protocols::genetic_algorithm::GeneticAlgorithmBase const & ga,
	protocols::pack_daemon::DynamicAggregateFunction const & daf
)
{
	using protocols::genetic_algorithm::GeneticAlgorithmBase;

	full_seqinds_2_newseqinds_.resize( ga.max_population_size() );

	core::Size count( 0 ), count_new( 0 );
	for ( GeneticAlgorithmBase::pop_const_iter
			iter = ga.current_generation_begin(),
			iter_end = ga.current_generation_end();
			iter != iter_end; ++iter ) {
		++count;
		if ( ! (*iter)->fitness_valid() ) {
			++count_new;
			std::string seqstring = sequence_from_entity( **iter );
			new_seq_inds_.push_back( count );
			sequences_.push_back( seqstring );
			sequence_to_index_map_[ seqstring ] = count;
			full_seqinds_2_newseqinds_[ count ] = count_new;
		}
	}

	n_new_sequences_ = count_new;
	job_completed_.resize( count_new );

	for ( core::Size ii = 1; ii <= n_new_sequences_; ++ii ) {
		job_completed_[ ii ].resize( daf.num_states(), false );
		for ( core::Size jj = 1; jj <= daf.num_states(); ++jj ) {
			job_order_.push_back( job_id( ii, jj ) );
		}
	}

}

OneGenerationJobInfo::JobID
OneGenerationJobInfo::job_id( core::Size sequence_index, core::Size state_index )
{
	return std::make_pair( sequence_index, state_index );
}


bool OneGenerationJobInfo::unassigned_jobs_remain() const {
	return job_order_.begin() != job_order_.end();
}

bool OneGenerationJobInfo::unfinished_jobs_outstanding() const {
	return work_outstanding_.begin() != work_outstanding_.end();
}

OneGenerationJobInfo::JobID
OneGenerationJobInfo::pop_job()
{
	JobID job = job_order_.front();
	job_order_.pop_front();
	return job;
}

void OneGenerationJobInfo::push_outstanding( JobID const & job )
{
	assert( work_outstanding_.find( job ) == work_outstanding_.end() );
	work_outstanding_.insert( job );
}

void OneGenerationJobInfo::remove_outstanding( JobID const & job )
{
	std::set< JobID >::iterator iter = work_outstanding_.find( job );
	assert( iter != work_outstanding_.end() );
	work_outstanding_.erase( iter );
	job_completed_[ job.first ][ job.second ] = true;
}

JobsForSequence::JobsForSequence( core::Size n_states, core::Size n_npd_properties )
{
	jobs_.resize( n_states );
	state_energies_.resize( n_states );
	npd_properties_.resize( n_npd_properties );
}

JobsForSequence::~JobsForSequence() {}

void JobsForSequence::entity( protocols::genetic_algorithm::EntityOP setting )
{
	entity_ = setting;
}

protocols::genetic_algorithm::EntityOP JobsForSequence::entity()
{
	return entity_;
}

PackingJobRecord & JobsForSequence::job_record( core::Size state_index )
{
	return jobs_[ state_index ];
}

PackingJobRecord const & JobsForSequence::job_record ( core::Size state_index ) const
{
	return jobs_[ state_index ];
}

void JobsForSequence::finalize_state_energies_and_npd_properties()
{
	for ( core::Size ii = 1; ii <= jobs_.size(); ++ii ) {
		state_energies_[ ii ] = jobs_[ ii ].energy();
		for ( PackingJobRecord::npd_properties::const_iterator
				iter = jobs_[ ii ].npd_props().begin(),
				iter_end = jobs_[ ii ].npd_props().end();
				iter != iter_end; ++iter ) {
			npd_properties_[ iter->first ] = npd_properties_[ iter->second ];
		}
	}
}

utility::vector1< core::Real > const &
JobsForSequence::state_energies() const
{
	return state_energies_;
}

utility::vector1< core::Real > const &
JobsForSequence::npd_properties() const
{
	return npd_properties_;
}


MMTDriver::MMTDriver() :
	n_workers_( 0 ),
	ga_pop_size_( 0 ),
	ga_num_to_propagate_( 0 ),
	ga_ngenerations_( 0 ),
	ga_frac_by_recomb_( 0.0 ),
	n_top_results_to_output_( 1 )
{}

MMTDriver::~MMTDriver()
{}

void MMTDriver::set_nworkers( core::Size setting ) { n_workers_ = setting; node_stats_.resize( setting ); }
void MMTDriver::set_ngenerations( core::Size setting ) { ga_ngenerations_ = setting; }
void MMTDriver::set_pop_size( core::Size setting ) { ga_pop_size_ = setting; ga_num_to_propagate_ = setting / 2;}
void MMTDriver::set_frac_by_recomb( core::Real setting ) { ga_frac_by_recomb_ = setting; }
void MMTDriver::set_randomizer( protocols::genetic_algorithm::PositionSpecificRandomizerOP setting ) { ga_randomizer_ = setting; }
void MMTDriver::set_sfxn( core::scoring::ScoreFunction const & sfxn ) { sfxn_ = sfxn.clone(); }
void MMTDriver::set_file_contents( utility::io::FileContentsMapOP file_contents ) { file_contents_ = file_contents; }
void MMTDriver::set_n_results_to_output( core::Size setting )  { n_top_results_to_output_ = setting; }
void MMTDriver::set_daf_fname( std::string const & setting )  { daf_fname_ = setting; }
void MMTDriver::set_entity_resfile_fname( std::string const & setting )  { entity_resfile_fname_ = setting; }

/// @throws Throws a utility::excn::EXCN_Msg_Exception if the input entity resfile or the
/// input dynamic aggregate function file cannot be properly initialized.
void
MMTDriver::setup()
{
	if ( ! sfxn_ ) sfxn_ = core::scoring::get_score_function();
	if ( ! file_contents_ ) file_contents_ = utility::io::FileContentsMapOP( new utility::io::FileContentsMap );
	try {
		std::istringstream entity_resfile_stream( file_contents_->get_file_contents_ref( entity_resfile_fname_ ) );
		protocols::pack_daemon::create_entity_resfile_contents(
			entity_resfile_stream, entity_resfile_fname_,
			entity_resfile_contents_, entity_task_, num_entities_ );
	} catch ( utility::excn::EXCN_Msg_Exception & e ) {
		throw utility::excn::EXCN_Msg_Exception( "MMTDriver failed to initialize the entity resfile\n" + e.msg() );
	}

	daf_ = protocols::pack_daemon::DynamicAggregateFunctionOP( new protocols::pack_daemon::DynamicAggregateFunction );
	daf_->set_num_entity_elements( num_entities_ );
	daf_->set_score_function( *sfxn_ );
	try {
		std::istringstream daf_stream( file_contents_->get_file_contents_ref( daf_fname_ ) );
		daf_->read_all_variables_from_input_file( daf_stream );
	} catch ( utility::excn::EXCN_Msg_Exception & e ) {
		throw utility::excn::EXCN_Msg_Exception( "MMTDriver::setup failed to initialize DyanmicAggregateFunction\n" + e.msg() );
	}

	if ( ! ga_randomizer_ ) {
		ga_randomizer_ = protocols::genetic_algorithm::PositionSpecificRandomizerOP( new protocols::genetic_algorithm::Mutate1Randomizer );
		// reset the default value of 1.0, but the mutation-rate variable is not used by the Mutate1Randomizer!
		ga_randomizer_->set_mutation_rate( 0.0 /*option[ ms::mutate_rate ]()*/ );
		protocols::pack_daemon::initialize_ga_randomizer_from_entity_task( ga_randomizer_, entity_task_ );
	}

	ga_ = protocols::genetic_algorithm::GeneticAlgorithmBaseOP( new protocols::genetic_algorithm::GeneticAlgorithmBase );
	ga_->set_max_generations( ga_ngenerations_ );
	ga_->set_max_pop_size( ga_pop_size_ );
	ga_->set_num_to_propagate( ga_num_to_propagate_ );
	ga_->set_frac_by_recomb( ga_frac_by_recomb_ );
	ga_->set_rand( ga_randomizer_ );
	ga_->fill_with_random_entities();

	top_entities_ = protocols::pack_daemon::TopEntitySetOP( new protocols::pack_daemon::TopEntitySet );
	top_entities_->desired_entity_history_size( n_top_results_to_output_ );
}

void
MMTDriver::run()
{
	initial_handshake();
	main_optimization_loop();
}


std::map< std::string, std::string >
MMTDriver::retrieve_optimal_solutions()
{
	using protocols::pack_daemon::TopEntitySet;

	std::map< std::string, std::string > opt_sols;

	// there might be ties, so take the number of structures to output
	// from the top_entities_ data structure.
	core::Size n_entities_to_output = top_entities_->size();
	for ( core::Size ii = n_entities_to_output; ii > 0; --ii ) {
		TopEntitySet::EntityAndScore ent_and_score = top_entities_->pop();
		std::string ii_sequence = sequence_from_entity( *ent_and_score.first );
		JobsForSequenceOP ii_jobs = top_jobs_archive_[ ii_sequence ];

		protocols::pack_daemon::MultistateAggregateFunction::StateIndices relevant_states =
			daf_->select_relevant_states( ii_jobs->state_energies(), ii_jobs->npd_properties(), *ii_jobs->entity() );
		// now retrieve poses (pdbs) from the nodes on which these states were originally run

		for ( core::Size jj = 1; jj <= relevant_states.size(); ++jj ) {
			core::Size jj_state_index = relevant_states[ jj ];
			core::Size jj_node_run_on = ii_jobs->job_record( jj_state_index ).node_run_on();
			std::string pdb_string = retrieve_pdb_from_node( jj_node_run_on, jj_state_index, ii_sequence );

			core::pose::Pose jj_pose;

			std::string jj_out_name = "msd_output_" + utility::to_string( ii ) + "_" + daf_->state_name( jj_state_index ) + ".pdb";

			core::import_pose::ImportPoseOptions import_opts;
			core::import_pose::pose_from_pdbstring( jj_pose, pdb_string, import_opts, jj_out_name );

			TR << "Writing structure " << jj_out_name << " with score: " << (*sfxn_)( jj_pose ) << std::endl;

			std::ostringstream ostr;
			core::io::pdb::dump_pdb( jj_pose, ostr );
			core::io::pdb::extract_scores( jj_pose, jj_out_name, ostr );

			opt_sols[ jj_out_name ] = ostr.str();
		}

	}

	for ( core::Size ii = 1; ii <= n_workers_; ++ii ) {
		utility::send_integer_to_node( ii, spin_down );
	}

	return opt_sols;
}

void
MMTDriver::write_optimal_solutions_to_disk()
{
	std::map< std::string, std::string > opt_solutions = retrieve_optimal_solutions();
	for ( std::map< std::string, std::string >::const_iterator
			iter = opt_solutions.begin(),
			iter_end = opt_solutions.end();
			iter != iter_end; ++iter ) {
		utility::io::ozstream outfile( iter->first );
		outfile << iter->second;
	}

}
void
MMTDriver::initial_handshake() {
	TR << "Beginning handshake" << std::endl;
	for ( core::Size ii = 1; ii <= n_workers_; ++ii ) {
		TR << "Beginning handshake with node " << ii << std::endl;
		utility::send_integer_to_node( ii, handshake_begin );
		mmt_message message = mmt_message( utility::receive_integer_from_node( ii ) );
		if ( message == handshake_acknowledged ) {
			core::Size ii_capacity = utility::receive_integer_from_node( ii );
			node_stats_[ ii ].mpi_id_ = ii;
			node_stats_[ ii ].max_n_threads_ = ii_capacity;
			node_stats_[ ii ].curr_n_running_threads_ = 0;

			utility::send_string_to_node( ii, entity_resfile_fname_ );
			utility::send_string_to_node( ii, file_contents_->get_file_contents_ref( entity_resfile_fname_ ) );
		} else {
			// I don't know?! This is kind of unreachable.
		}
	}
	TR << "Handshake successfull" << std::endl;
}

void
MMTDriver::main_optimization_loop()
{
	while ( ! ga_->complete() ) {
		TR << "Beginning generation " << ga_->current_generation() << std::endl;
		if ( ga_->current_generation_complete() ) {
			// only call this after the first generation
			ga_->evolve_next_generation();
		}

		optimize_generation();
		evaluate_entity_fitnesses();


	}
}


bool
MMTDriver::optimize_generation()
{
	using namespace protocols::genetic_algorithm;

#ifdef MULTI_THREADED
#ifdef CXX11
	auto starttime = std::chrono::system_clock::now();
#endif
#else
	clock_t starttime = clock();
#endif


	this_gen_work_ = OneGenerationJobInfoOP( new OneGenerationJobInfo( *ga_, *daf_ ) );
	this_gen_results_.clear();
	this_gen_results_.resize( this_gen_work_->n_new_sequences() );
	for ( core::Size ii = 1; ii <= this_gen_work_->n_new_sequences(); ++ii ) {
		this_gen_results_[ ii ] = JobsForSequenceOP( new JobsForSequence( daf_->num_states(), daf_->num_npd_properties() ) );
	}

	// broadcast that we're starting a new generation
	for ( core::Size ii = 1; ii <= n_workers_; ++ii ) {
		utility::send_integer_to_node( ii, new_generation );
	}

	// main optimization loop
	// 1st approach: wait for each node to request work, and dole out the work as requested

	std::set< core::Size > nodes_still_working;
	for ( core::Size ii = 1; ii <= n_workers_; ++ii ) nodes_still_working.insert( ii );

	while ( ! nodes_still_working.empty() ||
			this_gen_work_->unfinished_jobs_outstanding() ||
			this_gen_work_->unassigned_jobs_remain() ) {

		// wait for a message
		core::Size communicating_node = utility::receive_integer_from_anyone();
		mmt_message mess = mmt_message( utility::receive_integer_from_node( communicating_node ) );

		if ( mess == requesting_new_job ) {
			if ( ! this_gen_work_->unassigned_jobs_remain() ) {
				assert( nodes_still_working.find( communicating_node ) != nodes_still_working.end() );
				utility::send_integer_to_node( communicating_node, generation_complete );
				nodes_still_working.erase( communicating_node );
			} else {
				send_new_job_to_node( communicating_node );
			}
		} else if ( mess == job_complete ) {
			receive_completed_job( communicating_node );
		}

	}

	assert( ! this_gen_work_->unassigned_jobs_remain() && ! this_gen_work_->unfinished_jobs_outstanding() );

	core::Real running_time;
#ifdef MULTI_THREADED
#ifdef CXX11
	auto stoptime = std::chrono::system_clock::now();
	running_time = std::chrono::duration_cast< std::chrono::seconds >( stoptime - starttime ).count();
#endif
#else
	clock_t stoptime = clock();
	running_time = ((double) stoptime - starttime ) / CLOCKS_PER_SEC;
#endif
	TR  << "Generation " << ga_->current_generation() << " took " << running_time << " seconds" << std::endl;

	return true;
}

void
MMTDriver::send_new_job_to_node(
	core::Size communicating_node
)
{
	utility::send_integer_to_node( communicating_node, new_job_ready );

	OneGenerationJobInfo::JobID job = this_gen_work_->pop_job();
	std::string const & job_seq = this_gen_work_->sequences()[ job.first ];

	TR << "Sending job " << job.first << " " << job.second << " to node " << communicating_node << std::endl;
	send_state_info_to_node( communicating_node, job.second );
	send_sequence_to_node( communicating_node, job_seq );

	PackingJobRecord & job_record = this_gen_results_[ job.first ]->job_record( job.second );
	job_record.state_index( job.second );
	job_record.sequence_string( job_seq );
	job_record.node_run_on( communicating_node );

	this_gen_work_->push_outstanding( job );
}

void MMTDriver::send_sequence_to_node(
	core::Size communicating_node,
	std::string const & seq
) const
{
	utility::send_string_to_node( communicating_node, seq );
}

void
MMTDriver::send_state_info_to_node(
	core::Size node_index,
	core::Size state_index
)
{
	// This function pairs with MMTReceiver::receive_state_input_data
	utility::send_integer_to_node( node_index, state_index );

	protocols::pack_daemon::StructureFileNames const & sfn = daf_->file_inputs_for_job( state_index );

	utility::send_string_to_node( node_index, sfn.pdb_name_ );
	utility::send_string_to_node( node_index, file_contents_->get_file_contents_ref( sfn.pdb_name_ ) );

	utility::send_string_to_node( node_index, sfn.correspondence_file_name_ );
	utility::send_string_to_node( node_index, file_contents_->get_file_contents_ref( sfn.correspondence_file_name_ ) );

	utility::send_string_to_node( node_index, sfn.resfile_name_ );
	utility::send_string_to_node( node_index, file_contents_->get_file_contents_ref( sfn.resfile_name_ ) );

	// now send the list of non-pairwise decomposable properties that need to be computed for the
	// final structure for this state.
	utility::send_integer_to_node( node_index, daf_->num_npd_properties_for_state( state_index ) );
	for ( std::list< std::pair< core::Size, std::string > >::const_iterator
			npditer     = daf_->npd_variable_indices_for_state_begin( state_index ),
			npditer_end = daf_->npd_variable_indices_for_state_end( state_index );
			npditer != npditer_end; ++npditer ) {
		utility::send_integer_to_node( node_index, npditer->first );
		utility::send_string_to_node( node_index, npditer->second );
	}


}

void
MMTDriver::receive_completed_job(
	core::Size communicating_node
)
{
	TR << "Job completed message from " << communicating_node << std::endl;
	core::Size  state_index  = utility::receive_integer_from_node( communicating_node );
	std::string sequence     = utility::receive_string_from_node( communicating_node );
	core::Real  final_energy = utility::receive_double_from_node( communicating_node );
	core::Real  running_time = utility::receive_double_from_node( communicating_node );

	core::Size n_npd_properties = utility::receive_integer_from_node( communicating_node );
	PackingJobRecord::npd_properties npds;
	for ( core::Size ii = 1; ii <= n_npd_properties; ++ii ) {
		core::Size prop_ind = utility::receive_integer_from_node( communicating_node );
		core::Real prop_val = utility::receive_double_from_node( communicating_node );
		npds.push_back( std::make_pair( prop_ind, prop_val ) );
	}

	core::Size seq_index = this_gen_work_->full_seqinds_2_newseqinds()[ this_gen_work_->sequence_to_index_map().find( sequence )->second ];

	TR << "Job completed message from " << communicating_node << "; job " << state_index << " " << seq_index << " took: " << running_time << " and produced an energy of " << final_energy << std::endl;

	PackingJobRecord & job_record = this_gen_results_[ seq_index ]->job_record( state_index );
	assert( job_record.node_run_on() == communicating_node );

	job_record.energy( final_energy );
	job_record.running_time( running_time );
	job_record.npd_props( npds );

	this_gen_work_->remove_outstanding( OneGenerationJobInfo::job_id( seq_index, state_index ) );
}

void
MMTDriver::evaluate_entity_fitnesses()
{
	utility::vector1< core::Real > state_energies( daf_->num_states() );
	utility::vector1< core::Real > npd_properties( daf_->num_npd_properties() );
	for ( core::Size ii = 1; ii <= this_gen_work_->n_new_sequences(); ++ii ) {
		core::Size ii_entity_index = this_gen_work_->new_seq_inds()[ ii ];
		protocols::genetic_algorithm::EntityOP ii_entity = ga_->curr_gen_entity( ii_entity_index );
		JobsForSequenceOP iijobs = this_gen_results_[ii];
		iijobs->finalize_state_energies_and_npd_properties();
		iijobs->entity( ii_entity );

		core::Real ii_fitness = daf_->evaluate( iijobs->state_energies(), iijobs->npd_properties(), *ii_entity );
		ii_entity->set_fitness( ii_fitness );

		bool keep( false );
		protocols::pack_daemon::TopEntitySet::StateEnergiesAndNPDs seanpds =
			std::make_pair( iijobs->state_energies(), iijobs->npd_properties() );
		std::list< protocols::genetic_algorithm::EntityOP > dropped = top_entities_->update_entity_history( *ii_entity, seanpds, keep );
		if ( keep ) {
			instruct_receivers_to_keep_job_data_for_entity( ii );
		}

		for ( std::list< protocols::genetic_algorithm::EntityOP >::iterator iter = dropped.begin(); iter != dropped.end(); ++iter ) {
			instruct_receivers_to_drop_old_job_data_for_entity( *iter );
		}

	}
}

void
MMTDriver::instruct_receivers_to_keep_job_data_for_entity(
	core::Size seq_index
)
{
	std::string const & sequence   = this_gen_work_->sequences()[ seq_index ];
	JobsForSequenceOP jobs_for_seq = this_gen_results_[ seq_index ];
	for ( core::Size ii = 1; ii <= daf_->num_states(); ++ii ) {
		core::Size ii_node = jobs_for_seq->job_record( ii ).node_run_on();
		utility::send_integer_to_node( ii_node, result_from_last_generation_needs_saving );
		utility::send_integer_to_node( ii_node, ii );
		utility::send_string_to_node(  ii_node, sequence );
	}
	top_jobs_archive_[ sequence ] = jobs_for_seq;
}

void
MMTDriver::instruct_receivers_to_drop_old_job_data_for_entity(
	protocols::genetic_algorithm::EntityOP entity
)
{
	std::string sequence = sequence_from_entity( *entity );

	SavedJobsForSequence::iterator iter = top_jobs_archive_.find(sequence );
	assert( iter != top_jobs_archive_.end() );
	JobsForSequence const & jobs_for_seq =  *iter->second;
	for ( core::Size ii = 1; ii <= daf_->num_states(); ++ii ) {
		core::Size ii_node = jobs_for_seq.job_record( ii ).node_run_on();
		utility::send_integer_to_node( ii_node, old_result_can_be_discarded );
		utility::send_integer_to_node( ii_node, ii );
		utility::send_string_to_node(  ii_node, sequence );
	}

	top_jobs_archive_.erase( iter );

}

std::string
MMTDriver::retrieve_pdb_from_node(
	core::Size node,
	core::Size state_index,
	std::string const & sequence
)
{
	utility::send_integer_to_node( node, recover_pose_for_result );
	send_state_info_to_node( node, state_index );
	send_sequence_to_node( node, sequence );

	mmt_message retrieval_status = mmt_message( utility::receive_integer_from_node( node ) );
	if ( retrieval_status == error ) {
		std::string emsg = utility::receive_string_from_node( node );
		// oh no! -- spin_down_everything_but_node( node );
		throw utility::excn::EXCN_Msg_Exception( emsg );
	}

	return utility::receive_string_from_node( node );
}


}
}
