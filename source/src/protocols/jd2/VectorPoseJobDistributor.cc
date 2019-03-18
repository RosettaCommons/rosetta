// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/VectorPoseJobDistributor.cc
/// @brief  Job distributor subclass for running RECON multistate design
/// @brief  Takes in all input poses from the command line and passes them to any mover or filter that
/// @brief  derives from VectorPoseMover or VectorPoseFilter, meaning that it is able to receive and
/// @brief  operate on multiple poses simultaneously. Only accessible through recon application.
/// @author Alex Sevy (alex.sevy@gmail.com)

#ifdef USEMPI
#include <mpi.h>
#endif

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/jd2.OptionKeys.gen.hh>
#include <basic/options/keys/recon.OptionKeys.gen.hh>
#include <basic/options/keys/parser.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <basic/Tracer.hh>
#include <protocols/jd2/VectorPoseJobDistributor.hh>
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
#include <core/pose/extra_pose_info_util.hh>
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

// #include <protocols/minimization_packing/PackRotamersMover.hh>
// #include <protocols/simple_moves/MutateResidue.hh>

#include <numeric/random/random_permutation.hh>


namespace protocols {
namespace jd2 {

static basic::Tracer TR ("protocols.jd2.VectorPoseJobDistributor");

VectorPoseJobDistributor::VectorPoseJobDistributor() :
	JobDistributor() {}

VectorPoseJobDistributor::~VectorPoseJobDistributor() {}

/// @brief Get job ID of next job. Robust to both MPI and non-MPI input
core::Size VectorPoseJobDistributor::get_new_job_id() {
#ifdef USEMPI
	return current_nstruct_ * utility::mpi_rank();
#else
	return current_nstruct_ * current_pose_;
#endif
}

/// @brief Get current job. Robust to both MPI and non-MPI input
JobOP VectorPoseJobDistributor::current_job() const {
#ifdef USEMPI
	return this_nodes_job_;
#else
	return job_map_.at( current_nstruct_ )[ current_pose_ ];
#endif

}

void VectorPoseJobDistributor::mark_current_job_id_for_repetition() {}

void VectorPoseJobDistributor::job_failed( core::pose::Pose &, bool ) {}

void VectorPoseJobDistributor::handle_interrupt() {}


void VectorPoseJobDistributor::go( protocols::moves::MoverOP mover ) {
#ifdef USEMPI
	go_mpi( mover );
#else
	go_serial( mover );
#endif
}

void VectorPoseJobDistributor::go_mpi( protocols::moves::MoverOP mover ) {

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
	bool master = (rank==0);

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

			parser.generate_mover_from_pose(
				*this_nodes_pose,
				mover,
				1,
				option[ parser::protocol ]()
			);

			// parser.generate_mover_from_pose( this_nodes_job_, *this_nodes_pose, mover, 1,
			//  option[ parser::protocol ]() );
			ParsedProtocolOP this_nodes_protocol =
				utility::pointer::dynamic_pointer_cast< ParsedProtocol > ( mover );

			bool passed = apply_parsed_protocol_mpi( this_nodes_pose, this_nodes_protocol );

			if ( passed ) {
				this_nodes_protocol->final_score( *this_nodes_pose );
				// Synchronize to make sure that only one node is writing at a time
				if ( master ) {
					job_outputter->final_pose( this_nodes_job_, *this_nodes_pose, "");
					job_outputter->flush();

					// Go through other nodes and give them permission to go
					for ( core::Size ii = 1; ii < n_procs; ++ii ) {
						utility::send_integer_to_node( ii, 0 );
						utility::receive_integer_from_node( ii );
					}
				} else {
					utility::receive_integer_from_node( 0 );
					job_outputter->final_pose( this_nodes_job_, *this_nodes_pose, "");
					job_outputter->flush();
					utility::send_integer_to_node( 0, 0 );
				}
				//    job_outputter->final_pose( this_nodes_job_, *this_nodes_pose, "");
				//    job_outputter->flush();
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

	/// Synchronize mover status before moving on
#ifdef USEMPI
	MPI_Barrier( MPI_COMM_WORLD );

	MPI_Finalize();
#endif
}

/// Go method used for non-MPI jobs. First we gather all of the input poses, then we pass
/// them to the movers specified in the ParsedProtocol either one after another (if it's a
/// regular mover/filter) or all together (if it's a VectorPoseMover or VectorPoseFilter)
void
VectorPoseJobDistributor::go_serial( protocols::moves::MoverOP mover ) {

	time_t const allstarttime = time(NULL);

	// mover is not defined yet - can't instantiate until it's defined as a parsed protocol below
	current_nstruct_ = 1;
	job_inputter_ = job_inputter();
	job_outputter_ = job_outputter();
	JobsContainer jobs = get_jobs();

	randomize_input_ = basic::options::option[ basic::options::OptionKeys::recon::randomize ].user();
	core::Size ntrials = 1;
	if ( basic::options::option[ basic::options::OptionKeys::jd2::ntrials ].user() ) {
		ntrials = basic::options::option[ basic::options::OptionKeys::jd2::ntrials ];
	}

	// Since each job is not its own independent task, this groups jobs by nstruct to create one super-job, so to speak
	for ( core::Size i = 1; i <= jobs.size(); ++i ) {
		job_map_[ jobs[i]->nstruct_index() ].push_back(jobs[ i ]);
	}
	while ( job_map_.find( current_nstruct_ ) != job_map_.end() ) {
		Jobs current_jobs = job_map_.find( current_nstruct_ )->second;
		utility::vector1< core::pose::PoseOP > working_poses;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		core::Size retry_counter = 0;

		while ( retry_counter < ntrials ) {

			time_t job_time = time(NULL);

			// Keep a record of the order that poses should be handled for randomization of inputs
			utility::vector1< core::Size > pose_order;

			for ( core::Size i = 1; i <= current_jobs.size(); ++i ) {
				core::pose::PoseOP pose ( new core::pose::Pose );
				job_inputter_->pose_from_job( *pose, current_jobs[i] );
				core::pose::setPoseExtraScore( *pose, "msd_job_dist_index", i );
				working_poses.push_back( pose );
				pose_order.push_back( i );
			}

			if ( working_poses.size() < 2 ) {
				utility_exit_with_message( "At least two input structures expected for multi state design" );
			}
			protocols::rosetta_scripts::RosettaScriptsParser parser;
			utility::vector1<protocols::rosetta_scripts::ParsedProtocolOP> protocols;


			for ( core::Size ii = 1; ii <= working_poses.size(); ++ii ) {
				current_pose_ = ii;

				parser.generate_mover_from_pose(
					*working_poses[ ii ],
					mover,
					1,
					option[ parser::protocol ]()
				);

				protocols::rosetta_scripts::ParsedProtocolOP parsed_protocol =
					utility::pointer::dynamic_pointer_cast< protocols::rosetta_scripts::ParsedProtocol > ( mover );
				protocols.push_back( parsed_protocol );
			}

			bool passed = apply_parsed_protocol_serial( working_poses, protocols, pose_order );

			if ( passed ) {
				for ( core::Size i = 1; i <= current_jobs.size(); ++i ) {
					protocols[ i ]->final_score( *working_poses[ i ] );
					job_outputter_->final_pose(current_jobs[ i ], *working_poses[ i ], "");
					job_outputter_->flush();
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

bool VectorPoseJobDistributor::apply_parsed_protocol_mpi( core::pose::PoseOP & pose,
	protocols::rosetta_scripts::ParsedProtocolOP & protocol ) {

	for ( core::Size mover_it = 1; mover_it <= protocol->size(); ++mover_it ) {

		protocols::moves::MoverOP current_mover = protocol->
			get_mover_filter_pair( mover_it ).first.first;
		protocols::filters::FilterOP current_filter = protocol->
			get_mover_filter_pair( mover_it ).second;

		TR << "=================running mover " << current_mover->get_name() << " - "
			<< protocol->get_mover_filter_pair( mover_it ).first.second
			<< "======================" << std::endl;
		time_t mover_time = time(NULL);
		// i refers to mover type -> j refers to each pose

		// Check if current mover is a VectorPoseMover
		// if so hand it the working poses
		moves::VectorPoseMoverOP vpm = utility::pointer::dynamic_pointer_cast< moves::VectorPoseMover >( current_mover );
		core::pose::PoseOP temp_pose = pose->clone();

		if ( vpm ) {
			vpm->apply_mpi( *temp_pose );

			/// Synchronize mover status before moving on
#ifdef USEMPI
				MPI_Barrier( MPI_COMM_WORLD );
#endif

		} else {
			current_mover->apply( *temp_pose );
		}

		TR << "mover " << protocol->get_mover_filter_pair( mover_it ).first.second << " finished in " << time(NULL) - mover_time << " seconds" << std::endl;
		TR << "=================end mover " << "======================" << std::endl;
		TR << "=================running filter " << current_filter->get_type() << "======================" << std::endl;

		// Check if current filter is a VectorPoseFilter
		// if so hand it the working poses
		filters::VectorPoseFilterOP vpf = utility::pointer::dynamic_pointer_cast< filters::VectorPoseFilter >( current_filter );

		bool filter_passed;

		if ( vpf ) {
			filter_passed = vpf->apply_mpi( *temp_pose );

			/// Synchronize mover status before moving on
#ifdef USEMPI
				MPI_Barrier( MPI_COMM_WORLD );
#endif
			/// END Synchronize

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

bool VectorPoseJobDistributor::apply_parsed_protocol_serial( utility::vector1< core::pose::PoseOP > & working_poses,
	utility::vector1<protocols::rosetta_scripts::ParsedProtocolOP> & protocols,
	utility::vector1< core::Size > & pose_order ) {

	for ( core::Size mover_it = 1; mover_it <= protocols[ 1 ]->size(); ++mover_it ) {

		// If flag is passed to randomize input pose order, randomize the list
		if ( randomize_input_ ) {
			numeric::random::random_permutation( pose_order, numeric::random::rg() );
		}


		for ( core::Size i = 1; i <= pose_order.size(); ++i ) {

			current_pose_ = pose_order[ i ];

			protocols::moves::MoverOP current_mover = protocols[ current_pose_ ]->
				get_mover_filter_pair( mover_it ).first.first;
			protocols::filters::FilterOP current_filter = protocols[ current_pose_ ]->
				get_mover_filter_pair( mover_it ).second;
			TR << "=================running mover " << current_mover->get_name() << " - "
				<< protocols[ current_pose_ ]->get_mover_filter_pair( mover_it ).first.second
				<< "======================" << std::endl;

			// i refers to mover type -> j refers to each pose

			// Check if current mover is a VectorPoseMover
			// if so hand it the working poses
			moves::VectorPoseMoverOP vpm = utility::pointer::dynamic_pointer_cast< moves::VectorPoseMover >( current_mover );
			core::pose::PoseOP temp_pose = working_poses[ current_pose_ ]->clone();

			if ( vpm ) {
				vpm->set_poses( working_poses );
				vpm->apply( *temp_pose );
			} else {
				current_mover->apply( *temp_pose );
			}


			TR << "=================end mover " << "======================" << std::endl;
			TR << "=================running filter " << current_filter->get_type() << "======================" << std::endl;

			// Check if current filter is a VectorPoseFilter
			// if so hand it the working poses
			filters::VectorPoseFilterOP vpf = utility::pointer::dynamic_pointer_cast< filters::VectorPoseFilter >( current_filter );

			bool filter_passed;

			if ( vpf ) {
				vpf->set_poses( working_poses );
				filter_passed = vpf->apply( *temp_pose );
			} else {
				filter_passed = current_filter->apply( *temp_pose );
			}

			TR << "=================end filter " << current_filter->get_type() <<": " << (filter_passed ? " passed! " : " failed ;-(" ) << "======================" << std::endl;
			if ( filter_passed ) {
				working_poses[ current_pose_ ] = temp_pose;

			} else {
				return false;
			}

		}
	}

	return true;
}

} // jd2
} // protocols

