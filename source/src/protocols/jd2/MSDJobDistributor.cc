// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/MSDJobDistributor.cc
/// @brief  Job distributor subclass for running restrained multistate design
/// @brief  Takes in all input poses from the command line and passes them to any mover that
/// @brief  derives from VectorPoseMover, meaning that it is able to receive  and operate on multiple poses simultaneously
/// @author Alex Sevy (alex.sevy@gmail.com)

#include <basic/Tracer.hh>
#include <core/types.hh>
#include <utility/pointer/owning_ptr.hh>
#include <core/pose/Pose.hh>
#include <protocols/jd2/MSDJobDistributor.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobDistributorFactory.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/VectorPoseMover.hh>
//#include <protocols/simple_moves/MSDMover.hh>
#include <protocols/filters/Filter.hh>
#include <utility/excn/EXCN_Base.hh>
//option key includes
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/parser.OptionKeys.gen.hh>
#include <basic/options/keys/jd2.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/option.hh>
// protocol specific includes
//#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/rosetta_scripts/RosettaScriptsParser.hh>
#include <protocols/rosetta_scripts/ParsedProtocol.hh>
#include <core/pose/util.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <numeric/random/random_permutation.hh>

namespace protocols {
namespace jd2 {

static THREAD_LOCAL basic::Tracer TR ("protocols.jd2.MSDJobDistributor");

MSDJobDistributor::~MSDJobDistributor() = default;

core::Size
MSDJobDistributor::get_new_job_id() {
	return current_nstruct_ * current_pose_;
}

JobOP
MSDJobDistributor::current_job() const {
	return job_map_.at( current_nstruct_ )[ current_pose_ ];
}

void
MSDJobDistributor::mark_current_job_id_for_repetition() {
	// do nothing - increment nstruct and continue
	// current_nstruct_++;
}

void
MSDJobDistributor::job_failed( core::pose::Pose & /*pose*/,
	bool /*will_retry*/ ) {
	//  current_nstruct_++;
}

void
MSDJobDistributor::handle_interrupt() {
	// does nothing - no temporary files are being created
}

void
MSDJobDistributor::go( protocols::moves::MoverOP mover ) {

	time_t const allstarttime = time(nullptr);

	// mover is not defined yet - can't instantiate until it's defined as a parsed protocol below
	current_nstruct_ = 1;
	job_inputter_ = job_inputter();
	job_outputter_ = job_outputter();
	JobsContainer jobs = get_jobs();

	randomize_input_ = basic::options::option[ basic::options::OptionKeys::run::msd_randomize ].user();
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

			time_t job_time = time(nullptr);

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
				parser.generate_mover_from_pose( *working_poses[ ii ], mover, 1,
					option[ parser::protocol ]() );
				protocols::rosetta_scripts::ParsedProtocolOP parsed_protocol =
					utility::pointer::dynamic_pointer_cast< protocols::rosetta_scripts::ParsedProtocol > ( mover );
				protocols.push_back( parsed_protocol );
			}

			bool passed = apply_parsed_protocol( working_poses, protocols, pose_order );

			if ( passed ) {
				for ( core::Size i = 1; i <= current_jobs.size(); ++i ) {
					protocols[ i ]->final_score( *working_poses[ i ] );
					job_outputter_->final_pose(current_jobs[ i ], *working_poses[ i ], "");
					job_outputter_->flush();
				}
				TR << "Current job completed in " << (time(nullptr) - job_time) << " seconds " << std::endl;
				break;
			}

			++retry_counter;
		}

		if ( retry_counter >= ntrials ) {
			TR.Error << "Failures exceeded limit of " << ntrials << std::endl;
		}

		++current_nstruct_;

	}
	TR << "All jobs completed in " << (time(nullptr) - allstarttime ) << " seconds " << std::endl;

}

bool MSDJobDistributor::apply_parsed_protocol( utility::vector1< core::pose::PoseOP > & working_poses,
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

			moves::VectorPoseMoverOP msd_mover = utility::pointer::dynamic_pointer_cast< moves::VectorPoseMover >( current_mover );
			core::pose::PoseOP temp_pose = working_poses[ current_pose_ ]->clone();

			if ( msd_mover ) {
				msd_mover->set_poses( working_poses );
				msd_mover->apply( *temp_pose );
			} else {
				current_mover->apply( *temp_pose );
			}


			TR << "=================end mover " << "======================" << std::endl;
			TR << "=================running filter " << current_filter->get_type() << "======================" << std::endl;

			bool filter_passed = current_filter->apply( *temp_pose );
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


} //jd2
} //protocols
