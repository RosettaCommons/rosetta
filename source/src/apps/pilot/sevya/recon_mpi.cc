// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   apps/pilot/sevya/recon_mpi.cc
/// @brief  Multistate design through RECON run through MPI
/// @author Alex Sevy (alex.sevy@gmail.com)

#ifdef USEMPI
#include <mpi.h>
#endif

#include <devel/init.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/jd2.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/parser.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

#include <protocols/jd2/JobInputter.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobsContainer.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/VectorPoseMover.hh>
#include <protocols/filters/VectorPoseFilter.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/chemical/AA.hh>
#include <core/scoring/constraints/ResidueTypeConstraint.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/conformation/Residue.hh>

#include <protocols/rosetta_scripts/RosettaScriptsParser.hh>
#include <protocols/rosetta_scripts/ParsedProtocol.hh>

#include <utility/mpi_util.hh>
#include <utility/string_util.hh>

#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/MutateResidue.hh>

static basic::Tracer TR ("pilot_apps.sevya.recon_mpi");

namespace protocols {
namespace jd2 {

class RECONMPIJobDistributor : public protocols::jd2::JobDistributor {

public:
	// TODO add something to these pure virtual methods
	core::Size get_new_job_id() {
		return current_nstruct_ * utility::mpi_rank();
	}

	protocols::jd2::JobOP current_job() const { return this_nodes_job_; }
	void mark_current_job_id_for_repetition() {};
	void job_failed( core::pose::Pose & ,
		bool  ) {};
	void handle_interrupt() {};
	// TODO add something to these pure virtual methods


	utility::vector1< core::scoring::constraints::ConstraintCOP >
	apply_linking_constraints( core::pose::PoseOP & pose, utility::vector1< std::string > other_pose_sequences, core::Real constraint_weight ) {
		using namespace core::scoring::constraints;
		utility::vector1< ConstraintCOP > constraints;
		for ( core::Size ii = 1; ii <= other_pose_sequences.size(); ++ii ) {
			std::string other_pose_sequence = other_pose_sequences[ ii ];

			for ( core::Size jj = 1; jj <= my_designable_residues_.size(); ++jj ) {
				core::Size current_position = my_designable_residues_[ jj ];
				char constrained_AA = other_pose_sequence[ jj-1 ]; // string is zero indexed but designable_residues is one indexed
				std::string constrained_AA_3letter = core::chemical::name_from_aa( core::chemical::aa_from_oneletter_code( constrained_AA ) );
				TR << "adding linking constraint to position " <<
					current_position <<
					" for residue " <<
					constrained_AA_3letter <<
					" for weight " <<
					constraint_weight <<
					std::endl;
				ResidueTypeConstraintCOP temp_cst ( new ResidueTypeConstraint(
					*pose,
					current_position, //seqpos
					constrained_AA_3letter, //AAname
					constraint_weight //constraint weight
					) );
				constraints.push_back( pose->add_constraint( temp_cst ) );

			}
		}
		return constraints;
	}


	bool design_this_nodes_pose( core::pose::PoseOP & pose, std::string & this_nodes_resfile, core::Real constraint_weight ) {
		using namespace core::pack::task;

		/// Create design mover and add resfile behavior if not already initialized
		if ( packer_ ) TR << "packer exists on node " << utility::mpi_rank() << std::endl;
		else TR << "packer does not exist on node " << utility::mpi_rank() << std::endl;

		if ( !packer_ ) {
			packer_ = protocols::simple_moves::PackRotamersMoverOP( new protocols::simple_moves::PackRotamersMover );
			TaskFactoryOP factory ( new TaskFactory );
			operation::ReadResfileOP rrf ( new operation::ReadResfile( this_nodes_resfile ) );
			factory->push_back( rrf );

			operation::InitializeFromCommandlineOP ifcl ( new operation::InitializeFromCommandline );
			factory->push_back( ifcl );

			packer_->task_factory( factory );
		}
		/// END Create design mover

		/// Create score function properly weighted for residue type constraints
		if ( !sfxn_ ) {
			sfxn_ = core::scoring::get_score_function();
			sfxn_->set_weight( core::scoring::res_type_constraint, 1.0 );
			packer_->score_function( sfxn_ );
		}
		/// END Create score function

		/// Find out my designable residues from my pose and my resfile if they're not already defined
		if ( my_designable_residues_.size() == 0 ) {
			PackerTaskOP packer_task = TaskFactory::create_packer_task( *pose );
			parse_resfile( *pose, *packer_task, this_nodes_resfile );
			utility::vector1<bool> designing = packer_task->designing_residues();
			for ( core::Size i = 1; i <= designing.size(); ++i ) {
				if ( designing[ i ] ) my_designable_residues_.push_back( i );
			}

			/// Let rank 0 gather the number of designable residues from each rank and compare them
			if ( utility::mpi_rank() == 0 ) {
				for ( core::Size rank_to_receive = 1; rank_to_receive < number_jobs_; ++rank_to_receive ) {
					core::Size number_designable_positions = utility::receive_integer_from_node( rank_to_receive );
					if ( number_designable_positions != my_designable_residues_.size() ) {
						TR << "Rank " << rank_to_receive << " has " << number_designable_positions <<
							" and rank 0 has " << my_designable_residues_.size() << std::endl;
						utility_exit_with_message( "Error: all states must have the same number of designable residues" );
					}
				}
			} else {
				core::Size number_designable_positions = my_designable_residues_.size();
				utility::send_integer_to_node( 0, number_designable_positions  );
			}
			/// END Let rank 0 gather
		}
		/// END Find out my designable residues


		/// Make a string out of the AAs at my designable positions in the current state
		std::string my_sequence = "";
		for ( core::Size resno = 1; resno <= my_designable_residues_.size(); ++resno ) {
			my_sequence += pose->residue( my_designable_residues_[ resno ] ).name1();
		}
		/// END Make a string

		/// Get the AAs at designable positions of the other states I need to cooperate with
		utility::vector1<std::string> other_pose_sequences;
		for ( core::Size jobno = 1; jobno <= number_jobs_; ++jobno ) {
			core::Size rank_to_receive = jobno - 1;

			if ( (int)rank_to_receive == utility::mpi_rank() ) {
				for ( core::Size rank_to_send = 0; rank_to_send < number_jobs_; ++rank_to_send ) {
					if ( (int)rank_to_send != utility::mpi_rank() ) utility::send_string_to_node( rank_to_send, my_sequence );
				}
			} else {
				// Get all residues that are on the pose from a different processor
				std::string received_sequence = utility::receive_string_from_node( rank_to_receive );
				other_pose_sequences.push_back( received_sequence );
			}
		}
		/// END Get the AAs at designable positions

		/// Run my multistate design
		utility::vector1< core::scoring::constraints::ConstraintCOP > constraints = apply_linking_constraints( pose, other_pose_sequences, constraint_weight );
		packer_->apply( *pose );
		pose->remove_constraints( constraints, true );
		/// END Run my multistate design

		return 1;
	}

	void pick_consensus_AA( core::pose::PoseOP & this_nodes_pose, std::string my_sequence, utility::vector1<std::string> other_pose_sequences, core::Size zero_index ) {
		utility::vector1<char> candidate_AAs;
		candidate_AAs.push_back( my_sequence[ zero_index ] );
		for ( core::Size ii = 1; ii <= other_pose_sequences.size(); ++ii ) {
			char current_state_AA = other_pose_sequences[ ii ][ zero_index ];
			if ( std::find(candidate_AAs.begin(), candidate_AAs.end(), current_state_AA)
					== candidate_AAs.end() ) candidate_AAs.push_back( current_state_AA );
		}

		utility::vector1<core::Real> candidate_AA_scores;
		protocols::simple_moves::MutateResidueOP mutation_mover;
		/// Create packer to be used post-threading
		using namespace core::pack::task;
		protocols::simple_moves::PackRotamersMoverOP post_packer = protocols::simple_moves::PackRotamersMoverOP( new protocols::simple_moves::PackRotamersMover );
		core::scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function();
		sfxn->set_weight( core::scoring::res_type_constraint, 1.0 );
		post_packer->score_function( sfxn );

		operation::RestrictToRepackingOP rtr ( new operation::RestrictToRepacking );
		TaskFactoryOP factory = TaskFactoryOP( new TaskFactory );
		factory->push_back( rtr );
		post_packer->task_factory( factory );
		/// END Create packer

		for ( core::Size ii = 1; ii <= candidate_AAs.size(); ++ii ) {
			std::string trial_AA = core::chemical::name_from_aa( core::chemical::aa_from_oneletter_code( candidate_AAs[ ii ] ) );
			mutation_mover = protocols::simple_moves::MutateResidueOP ( new protocols::simple_moves::MutateResidue (
				my_designable_residues_[ zero_index+1 ], //position
				trial_AA //residue
				) );
			core::pose::PoseOP tmp_pose = this_nodes_pose->clone();
			mutation_mover->apply( *tmp_pose );
			post_packer->apply( *tmp_pose );
			candidate_AA_scores.push_back( sfxn->score( *tmp_pose ) );
		}

		/// From all my candidate AA scores compile them to make a fitness
		utility::vector1< core::Real > candidate_AA_fitnesses (candidate_AA_scores.size());
		for ( core::Size ii = 1; ii <= candidate_AA_scores.size(); ++ii ) {
			// Let the master do all the work
			if ( utility::mpi_rank() == 0 ) {
				candidate_AA_fitnesses[ii] = candidate_AA_scores[ ii ];
				for ( core::Size jj = 1; jj < number_jobs_; ++jj ) {
					core::Real other_AA_score = utility::receive_double_from_node( jj );
					candidate_AA_fitnesses[ ii ] += other_AA_score;
				}

			} else {
				utility::send_double_to_node( 0, candidate_AA_scores[ ii ] );
			}
		}
		/// END compile candidate AA scores

		/// Tell everyone what the best AA is
		utility::vector1<core::Real>::iterator best_fitness_iter = std::min_element(std::begin(candidate_AA_fitnesses), std::end(candidate_AA_fitnesses));
		core::Real best_fitness = *best_fitness_iter;
		core::Size best_candidate_index = candidate_AA_fitnesses.index( best_fitness );
		std::string best_AA;
		if ( utility::mpi_rank() == 0 ) {
			best_AA = core::chemical::name_from_aa( core::chemical::aa_from_oneletter_code( candidate_AAs[ best_candidate_index ] ) );
			for ( core::Size jj = 1; jj < number_jobs_; ++jj ) {
				utility::send_string_to_node( jj, best_AA );
			}
		} else {
			best_AA = utility::receive_string_from_node( 0 );
		}
		// END Tell everyone

		/// Thread best AA over my state
		mutation_mover = protocols::simple_moves::MutateResidueOP ( new protocols::simple_moves::MutateResidue (
			my_designable_residues_[ zero_index+1 ], //position
			best_AA //residue
			) );
		mutation_mover->apply( *this_nodes_pose );
		post_packer->apply( *this_nodes_pose );
		/// END thread best AA over state
	}

	bool find_consensus_sequence( core::pose::PoseOP & this_nodes_pose ) {

		/// Make a string out of the AAs at my designable positions in the current state
		std::string my_sequence = "";
		for ( core::Size resno = 1; resno <= my_designable_residues_.size(); ++resno ) {
			my_sequence += this_nodes_pose->residue( my_designable_residues_[ resno ] ).name1();
		}
		/// END Make a string

		/// Get the AAs at designable positions of the other states I need to cooperate with
		utility::vector1<std::string> other_pose_sequences;
		for ( core::Size jobno = 1; jobno <= number_jobs_; ++jobno ) {
			core::Size rank_to_receive = jobno - 1;

			if ( (int)rank_to_receive == utility::mpi_rank() ) {
				for ( core::Size rank_to_send = 0; rank_to_send < number_jobs_; ++rank_to_send ) {
					if ( (int)rank_to_send != utility::mpi_rank() ) utility::send_string_to_node( rank_to_send, my_sequence );
				}
			} else {
				// Get all residues that are on the pose from a different processor
				std::string received_sequence = utility::receive_string_from_node( rank_to_receive );
				other_pose_sequences.push_back( received_sequence );
			}
		}
		/// END Get the AAs at designable positions

		protocols::simple_moves::MutateResidueOP mutation_mover;
		for ( core::Size i = 0; i < my_sequence.size(); ++i ) {
			bool diff = false;
			for ( core::Size j = 1; j <= other_pose_sequences.size(); ++j ) {
				std::string other_pose_sequence = other_pose_sequences[ j ];

				if ( my_sequence[ i ] != other_pose_sequence[ i ] ) {
					diff = true;
					break;
				}
			}

			// if there is some positive state that is not uniform at this spot
			// move on with the protocol
			if ( diff ) pick_consensus_AA( this_nodes_pose, my_sequence, other_pose_sequences, i );
		}
		return 1;
	}

	void go_hardcoded( protocols::moves::MoverOP ) {

		time_t const allstarttime = time(NULL);
		JobInputterOP job_inputter = this->job_inputter();
		JobOutputterOP job_outputter = this->job_outputter();

		JobsContainer jobs = get_jobs();

		core::Size rank = utility::mpi_rank();
		core::Size n_procs = utility::mpi_nprocs();

		// Since each job is not its own independent task, this groups jobs by nstruct to create one super-job, so to speak
		for ( core::Size i = 1; i <= jobs.size(); ++i ) {
			job_map_[ jobs[i]->nstruct_index() ].push_back(jobs[ i ]);
		}

		current_nstruct_ = 1;

		// Find number of jobs -> if they are less than number of processors then give an error message
		number_jobs_ = job_map_.find( current_nstruct_ )->second.size();

		if ( number_jobs_ != n_procs ) {
			TR << "you passed " << number_jobs_ << " jobs and " << n_procs << " processors" << std::endl;
			utility_exit_with_message( "Error: to use RECON with MPI please use as many processors as you have states" );
		}

		if ( number_jobs_ < 2 ) {
			utility_exit_with_message( "At least two input structures expected for multistate design" );
		}

		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		std::string this_nodes_resfile;
		/// Get resfile list from command line
		if ( rank == 0 ) { // if I'm the master read them in
			utility::vector1< std::string> resfiles = option[ packing::resfile ]();

			if ( resfiles.size() != number_jobs_ ) {
				utility_exit_with_message( "Error - number of states must be the same as number of resfiles" );
			}

			for ( core::Size i = 1; i <= resfiles.size(); ++i ) {
				TR << "resfile " << i << " name: " << resfiles[ i ] << std::endl;
				if ( i == 1 ) {
					this_nodes_resfile = resfiles[ i ];
				} else {
					utility::send_string_to_node( i-1, resfiles[i] );
				}
			}

		} else { // if I'm the slave just get passed a resfile
			this_nodes_resfile = utility::receive_string_from_node( 0 );
		}

		while ( job_map_.find( current_nstruct_ ) != job_map_.end() ) {

			time_t job_time = time(NULL);

			Jobs current_jobs = job_map_.find( current_nstruct_ )->second;

			this_nodes_job_ = current_jobs[ rank+1 ];
			core::pose::PoseOP this_nodes_pose ( new core::pose::Pose );
			job_inputter->pose_from_job( *this_nodes_pose, this_nodes_job_ );
			core::pose::setPoseExtraScore( *this_nodes_pose, "msd_job_dist_index", rank+1 );


			using namespace protocols::rosetta_scripts;

			//    RosettaScriptsParser parser;
			//
			//    parser.generate_mover_from_pose( this_nodes_job, *this_nodes_pose, mover, 1,
			//     option[ parser::protocol ]() );
			//    ParsedProtocolOP this_nodes_protocol =
			//     utility::pointer::dynamic_pointer_cast< ParsedProtocol > ( mover );

			/// As of now, don't support XML input - just run a typical RECON run
			utility::vector1< core::Real > constraint_weights;
			constraint_weights.push_back( 0.5 );
			constraint_weights.push_back( 1.0 );
			constraint_weights.push_back( 1.5 );
			constraint_weights.push_back( 2.0 );

			bool passed;
			for ( core::Size i = 1; i <= constraint_weights.size(); ++i ) {
				TR << "Running MSD with constraint weight " << constraint_weights[ i ] << std::endl;
				passed = design_this_nodes_pose( this_nodes_pose, this_nodes_resfile, constraint_weights[ i ] );
				if ( !passed ) {
					utility_exit_with_message( "something went wrong..." );
				}
			}

			TR << "Running find consensus sequence " << std::endl;
			passed = find_consensus_sequence( this_nodes_pose );

			if ( passed ) {
				for ( core::Size i = 1; i <= current_jobs.size(); ++i ) {
					//      protocols[ i ]->final_score( *this_nodes_pose );
					job_outputter->final_pose( this_nodes_job_, *this_nodes_pose, "");
					job_outputter->flush();
				}
				TR << "Current job completed in " << (time(NULL) - job_time) << " seconds " << std::endl;
				break;
			}

			current_nstruct_++;
		}
		TR << "All jobs completed in " << (time(NULL) - allstarttime ) << " seconds " << std::endl;
	}

	bool apply_parsed_protocol( core::pose::PoseOP & pose,
		protocols::rosetta_scripts::ParsedProtocolOP & protocol ) {

		for ( core::Size mover_it = 1; mover_it <= protocol->size(); ++mover_it ) {

			protocols::moves::MoverOP current_mover = protocol->
				get_mover_filter_pair( mover_it ).first.first;
			protocols::filters::FilterOP current_filter = protocol->
				get_mover_filter_pair( mover_it ).second;

			TR << "=================running mover " << current_mover->get_name() << " - "
				<< protocol->get_mover_filter_pair( mover_it ).first.second
				<< "======================" << std::endl;

			// i refers to mover type -> j refers to each pose

			// Check if current mover is a VectorPoseMover
			// if so hand it the working poses
			moves::VectorPoseMoverOP msd_mover = utility::pointer::dynamic_pointer_cast< moves::VectorPoseMover >( current_mover );
			core::pose::PoseOP temp_pose = pose->clone();

			if ( msd_mover ) {
				msd_mover->apply_mpi( *temp_pose );
			} else {
				current_mover->apply( *temp_pose );
			}


			TR << "=================end mover " << "======================" << std::endl;
			TR << "=================running filter " << current_filter->get_type() << "======================" << std::endl;

			// Check if current filter is a VectorPoseFilter
			// if so hand it the working poses
			filters::VectorPoseFilterOP msd_filter = utility::pointer::dynamic_pointer_cast< filters::VectorPoseFilter >( current_filter );

			bool filter_passed;

			if ( msd_filter ) {
				filter_passed = msd_filter->apply_mpi( *temp_pose );
			} else {
				filter_passed = current_filter->apply( *temp_pose );
			}

			TR << "=================end filter " << current_filter->get_type() <<": " << (filter_passed ? " passed! " : " failed ;-(" ) << "======================" << std::endl;
			if ( filter_passed ) {
				pose = temp_pose;

			} else {
				return false;
			}
		}

		return true;
	}

	void go( protocols::moves::MoverOP mover ) {

		core::Size ntrials = 1;
		if ( basic::options::option[ basic::options::OptionKeys::jd2::ntrials ].user() ) {
			ntrials = basic::options::option[ basic::options::OptionKeys::jd2::ntrials ];
		}

		time_t const allstarttime = time(NULL);
		JobInputterOP job_inputter = this->job_inputter();
		JobOutputterOP job_outputter = this->job_outputter();

		JobsContainer jobs = get_jobs();

		core::Size rank = utility::mpi_rank();
		core::Size n_procs = utility::mpi_nprocs();

		// Since each job is not its own independent task, this groups jobs by nstruct to create one super-job, so to speak
		for ( core::Size i = 1; i <= jobs.size(); ++i ) {
			job_map_[ jobs[i]->nstruct_index() ].push_back(jobs[ i ]);
		}

		current_nstruct_ = 1;

		// Find number of jobs -> if they are less than number of processors then give an error message
		number_jobs_ = job_map_.find( current_nstruct_ )->second.size();

		if ( number_jobs_ != n_procs ) {
			TR << "you passed " << number_jobs_ << " jobs and " << n_procs << " processors" << std::endl;
			utility_exit_with_message( "Error: to use RECON with MPI please use as many processors as you have states" );
		}

		if ( number_jobs_ < 2 ) {
			utility_exit_with_message( "At least two input structures expected for multistate design" );
		}

		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		while ( job_map_.find( current_nstruct_ ) != job_map_.end() ) {

			time_t job_time = time(NULL);
			Jobs current_jobs = job_map_.find( current_nstruct_ )->second;
			this_nodes_job_ = current_jobs[ rank+1 ];
			core::pose::PoseOP this_nodes_pose ( new core::pose::Pose );
			job_inputter->pose_from_job( *this_nodes_pose, this_nodes_job_ );
			core::pose::setPoseExtraScore( *this_nodes_pose, "msd_job_dist_index", rank+1 );

			core::Size retry_counter = 0;

			while ( retry_counter < ntrials ) {

				using namespace protocols::rosetta_scripts;

				RosettaScriptsParser parser;

				parser.generate_mover_from_pose( this_nodes_job_, *this_nodes_pose, mover, 1,
					option[ parser::protocol ]() );
				ParsedProtocolOP this_nodes_protocol =
					utility::pointer::dynamic_pointer_cast< ParsedProtocol > ( mover );

				bool passed = apply_parsed_protocol( this_nodes_pose, this_nodes_protocol );

				if ( passed ) {
					for ( core::Size i = 1; i <= current_jobs.size(); ++i ) {
						this_nodes_protocol->final_score( *this_nodes_pose );
						job_outputter->final_pose( this_nodes_job_, *this_nodes_pose, "");
						job_outputter->flush();
					}
					TR << "Current job completed in " << (time(NULL) - job_time) << " seconds " << std::endl;
					break;
				}

				++retry_counter;
			}

			if ( retry_counter >= ntrials ) {
				TR.Error << "Error: failures exceeded limit of " << ntrials << std::endl;
			}

			++current_nstruct_;

		}
		TR << "All jobs completed in " << (time(NULL) - allstarttime ) << " seconds " << std::endl;
	}

private:
	protocols::jd2::JobOP this_nodes_job_;
	core::Size number_jobs_;
	std::map<core::Size, protocols::jd2::Jobs> job_map_;
	core::Size current_nstruct_;
	utility::vector1< core::Size > my_designable_residues_;
	protocols::simple_moves::PackRotamersMoverOP packer_;
	core::scoring::ScoreFunctionOP sfxn_;
};
} // jd2
} // protocols

int
main( int argc, char * argv [] )
{
	try {
		// setup random numbers and options
		devel::init(argc, argv);
		protocols::moves::MoverOP mover;
#ifdef USEMPI
  TR << "I'm using MPI ;-)" << std::endl;
#else
		TR << "I'm not using MPI ;-(" << std::endl;
#endif
		protocols::jd2::RECONMPIJobDistributor jobdist;
		core::Size rank = utility::mpi_rank();
		core::Size n_procs = utility::mpi_nprocs();
		jobdist.go( mover );
	} catch( utility::excn::EXCN_Base& excn ) {
		excn.display();
		std::exit( 1 );
	}
	return 0;
}

//bool apply_parsed_protocol( utility::vector1< core::pose::PoseOP > & working_poses,
// utility::vector1<protocols::rosetta_scripts::ParsedProtocolOP> & protocols,
// utility::vector1< core::Size > & pose_order ) {
//
// for ( core::Size mover_it = 1; mover_it <= protocols[ 1 ]->size(); ++mover_it ) {
//
//  // If flag is passed to randomize input pose order, randomize the list
//  if ( randomize_input_ ) {
//   numeric::random::random_permutation( pose_order, numeric::random::rg() );
//  }
//
//
//  for ( core::Size i = 1; i <= pose_order.size(); ++i ) {
//
//   current_pose_ = pose_order[ i ];
//
//   protocols::moves::MoverOP current_mover = protocols[ current_pose_ ]->
//    get_mover_filter_pair( mover_it ).first.first;
//   protocols::filters::FilterOP current_filter = protocols[ current_pose_ ]->
//    get_mover_filter_pair( mover_it ).second;
//   TR << "=================running mover " << current_mover->get_name() << " - "
//    << protocols[ current_pose_ ]->get_mover_filter_pair( mover_it ).first.second
//    << "======================" << std::endl;
//
//   // i refers to mover type -> j refers to each pose
//
//   // Check if current mover is a VectorPoseMover
//   // if so hand it the working poses
//   moves::VectorPoseMoverOP msd_mover = utility::pointer::dynamic_pointer_cast< moves::VectorPoseMover >( current_mover );
//   core::pose::PoseOP temp_pose = working_poses[ current_pose_ ]->clone();
//
//   if ( msd_mover ) {
//    msd_mover->set_poses( working_poses );
//    msd_mover->apply( *temp_pose );
//   } else {
//    current_mover->apply( *temp_pose );
//   }
//
//
//   TR << "=================end mover " << "======================" << std::endl;
//   TR << "=================running filter " << current_filter->get_type() << "======================" << std::endl;
//
//   // Check if current filter is a VectorPoseFilter
//   // if so hand it the working poses
//   filters::VectorPoseFilterOP msd_filter = utility::pointer::dynamic_pointer_cast< filters::VectorPoseFilter >( current_filter );
//
//   bool filter_passed;
//
//   if ( msd_filter ) {
//    msd_filter->set_poses( working_poses );
//    filter_passed = msd_filter->apply( *temp_pose );
//   } else {
//    filter_passed = current_filter->apply( *temp_pose );
//   }
//
//   TR << "=================end filter " << current_filter->get_type() <<": " << (filter_passed ? " passed! " : " failed ;-(" ) << "======================" << std::endl;
//   if ( filter_passed ) {
//    working_poses[ current_pose_ ] = temp_pose;
//
//   } else {
//    return false;
//   }
//
//  }
// }
//
// return true;
//}
