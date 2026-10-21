// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/recon_design/FindConsensusSequence.cc
/// @brief Takes in multiple poses from the VectorPoseJobDistributor and finds the consensus
/// sequence that optimizes energy of all input poses. Used in conjuction with MSDMover
/// at the end of a protocol to make sure that you end up with one multistate solution.
/// Only accessible through recon application.
/// @author Alex Sevy (alex.sevy@gmail.com)

#include <protocols/recon_design/FindConsensusSequence.hh>
#include <protocols/recon_design/FindConsensusSequenceCreator.hh>
#include <protocols/moves/VectorPoseMover.hh>

// type headers
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/tag/Tag.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/scoring/xml_util.hh>
#include <core/pack/task/xml_util.hh>

#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/simple_moves/MutateResidue.hh>

#include <utility/mpi_util.hh>
#include <utility/string_util.hh>

#include <protocols/recon_design/recon_util.hh>

#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

#include <core/pack/task/TaskFactory.hh> // AUTO IWYU For TaskFactory

namespace protocols {
namespace recon_design {

static basic::Tracer TR("protocols.recon_design.FindConsensusSequence");

std::string
FindConsensusSequenceCreator::keyname() const {
	return FindConsensusSequenceCreator::mover_name();
}

moves::MoverOP
FindConsensusSequenceCreator::create_mover() const {
	return FindConsensusSequenceOP( new FindConsensusSequence );
}

std::string
FindConsensusSequenceCreator::mover_name() {
	return "FindConsensusSequence";
}

FindConsensusSequence::FindConsensusSequence() :
	moves::VectorPoseMover( FindConsensusSequenceCreator::mover_name() ),
	sfxn_( core::scoring::get_score_function() ),
	task_factory_( core::pack::task::TaskFactoryOP( new core::pack::task::TaskFactory )),
	debug_( false )
{

}

FindConsensusSequence::~FindConsensusSequence() {}

moves::MoverOP FindConsensusSequence::clone() const {
	return FindConsensusSequenceOP( new FindConsensusSequence( *this ) );
}

moves::MoverOP FindConsensusSequence::fresh_instance() const {
	return FindConsensusSequenceOP ( new FindConsensusSequence );
}

void FindConsensusSequence::parse_my_tag(
	utility::tag::TagCOP const tag,
	basic::datacache::DataMap & datamap
)  {

	sfxn_ = core::scoring::parse_score_function( tag, "scorefxn", datamap );

	task_factory_ = core::pack::task::parse_task_operations( tag, datamap );

	if ( tag->hasOption( "resfiles" ) ) {
		resfiles_ = utility::string_split( tag->getOption<std::string>( "resfiles" ), ',' );
	} else {
		utility_exit_with_message("Error: must input a resfile for FindConsensusSequence");
	}

	debug_ = tag->getOption<bool>( "debug", false );

#ifndef NDEBUG
	debug_ = true;
#endif
}


void FindConsensusSequence::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	core::scoring::attributes_for_parse_score_function( attlist );
	core::pack::task::attributes_for_parse_task_operations( attlist );

	attlist
		+ XMLSchemaAttribute( "resfiles", xs_string, "A list of resfiles to define designable and repackable residues for all states in multistate design. Multiple resfiles can be used for multiple states - in this case the first resfile in the tag will be applied to the first structure, etc. One single resfile used for all states is also supported.")
		+ XMLSchemaAttribute::attribute_w_default( "debug", xsct_rosetta_bool, "if true, outputs extra messages during protocol", "false");

	protocols::moves::xsd_type_definition_w_attributes( xsd, FindConsensusSequenceCreator::mover_name(), "Takes in multiple poses from the VectorPoseJobDistributor and finds the consensus sequence that optimizes energy of all input poses", attlist );
}

void FindConsensusSequenceCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	FindConsensusSequence::provide_xml_schema( xsd );
}

/// Populates designable_residues_ with a list of designable residues corresponding
/// to the poses in poses_, corresponding element-wise (designables_residues[0] matches
/// to poses_[0], etc)
void
FindConsensusSequence::parse_resfiles () {

	// This vector should be empty, but clear it just in case
	designable_residues_.clear();

	for ( core::Size ii = 1; ii <= poses_.size(); ++ii ) {
		designable_residues_.push_back( get_designable_residues( *poses_[ ii ], resfile_at( ii ) ));

		if ( designable_residues_[ ii ].size() != designable_residues_[ 1 ].size() ) {
			utility_exit_with_message( "All resfiles must have the same number of designable residues");
		}
	}
}

/// Get the resfile corresponding to the pose at index. If
/// only one resfile is present then it will be returned regardless
/// of the value of index
std::string
FindConsensusSequence::resfile_at ( core::Size index ) {
	if ( resfiles_.size() == 1 ) {
		return resfiles_[ 1 ];
	} else if ( resfiles_.size() == 0 ) {
		utility_exit_with_message("No resfiles initialized");
	} else {
		if ( index > resfiles_.size() ) {
			utility_exit_with_message("Number of resfiles does not match number of structures");
		} else {
			return resfiles_[ index ];
		}
	}
}

/// Set up the packer to be used when substituting different
/// candidate AAs and repacking before evaluating energy
void FindConsensusSequence::initialize_packer() {
	using namespace core::pack::task;

	// If it's the first time I run this, set up all the boilerplate for the packer
	if ( !packer_ ) {
		packer_ = minimization_packing::PackRotamersMoverOP( new minimization_packing::PackRotamersMover );
		packer_->score_function( sfxn_ );
		operation::RestrictToRepackingOP rtr ( new operation::RestrictToRepacking );
		if ( !task_factory_ ) {
			task_factory_ = TaskFactoryOP( new TaskFactory );
		}
		task_factory_->push_back( rtr );
		packer_->task_factory( task_factory_ );
	}
}

void FindConsensusSequence::apply( core::pose::Pose & pose ) {
	if ( poses_.size() == 0 ) {
		utility_exit_with_message("Error: no poses initialized. "
			"Did you try to run this from Rosetta Scripts? "
			"You need to pass the -run:msd_job_dist flag");
	}

	// Find index of your apply pose in the pose vector
	core::Size current_pose = find_pose_in_vector( pose, poses_ );

	// Parse resfiles to get a list of designable residues for each pose
	parse_resfiles();

	// Get the designable sequences of all poses
	utility::vector1< utility::vector1< std::string > > all_pose_sequences;
	for ( core::Size ii = 1; ii <= designable_residues_.size(); ++ii ) {
		all_pose_sequences.push_back(
			get_designable_sequence ( *poses_[ii], designable_residues_[ii] )
		);
	}

	// See which positions differ between my poses to force them to converge
	for ( core::Size str_position = 1; str_position <= all_pose_sequences[ 1 ].size(); ++str_position ) {
		utility::vector1< std::string > candidate_AAs = get_candidate_AAs( all_pose_sequences, str_position );
		if ( candidate_AAs.size() > 1 ) {
			pick_consensus_AA( str_position, candidate_AAs );
		}
	}

	// Magic that prevents mover from running multiple times
	pose = *poses_[ current_pose ];
}

void FindConsensusSequence::apply_mpi( core::pose::Pose & pose ) {

	using namespace core::pack::task;

	/// Get my node's information
	rank_ = utility::mpi_rank()+1; // I want to make the rank 1-indexed
	n_procs_ = utility::mpi_nprocs();
	master_ = (rank_==1);

	/// Get my resfile
	std::string this_nodes_resfile = resfile_at( rank_ );

	utility::vector1< core::Size > my_designable_residues = get_designable_residues( pose, this_nodes_resfile );

	/// Make a string out of the AAs at my designable positions in the current state
	utility::vector1< std::string > my_sequence = get_designable_sequence ( pose, my_designable_residues );

	std::string pass_seq;
	for ( std::string const & resi_base_name: my_sequence ) {
		pass_seq += resi_base_name + " ";
	}

	utility::vector1< utility::vector1<std::string> > all_pose_sequences (n_procs_);
	for ( core::Size ii = 1; ii <= n_procs_; ++ii ) {
		if ( rank_ == ii ) {
			for ( core::Size jj = 1; jj <= n_procs_; ++jj ) {
				if ( rank_!=jj ) utility::send_string_to_node( jj-1, pass_seq ); // node ranks are 0-indexed
				else all_pose_sequences[jj] = my_sequence;
			}
		} else {
			//all_pose_sequences[ii] = utility::receive_string_from_node( ii-1 ); // node ranks are 0-indexed
			std::string passed_seq;
			passed_seq = utility::receive_string_from_node( ii-1 ); // node ranks are 0-indexed
			//Need to split passed_seq by spaces
			std::istringstream iss(passed_seq);
			std::string resi_name;
			while ( iss >> resi_name ) {
				all_pose_sequences[ii].push_back(resi_name);
			}
		}
	}

	if ( debug_ ) {
		TR << "printing all AAs at designable positions to check indices: "<< std::endl;
		for ( core::Size i = 1; i <= all_pose_sequences.size(); ++i ) {
			TR << i << ":";
			for ( const std::string& resi : all_pose_sequences[i] ) {
				TR << resi << "-";
			}
			TR << ", ";
		}
		TR << std::endl;
	}

	utility::vector1< std::string > candidate_AAs;
	for ( core::Size jj = 1; jj <= all_pose_sequences[ 1 ].size(); ++jj ) {
		candidate_AAs.clear();
		candidate_AAs = get_candidate_AAs( all_pose_sequences, jj );

		// If there are more than one candidate amino acids then you have to choose a consensus
		if ( candidate_AAs.size() > 1 ) {
			pick_consensus_AA_mpi( pose, candidate_AAs, my_designable_residues[jj] );
		}
	}
}


/// Based on all the input poses, find the optimal AA at position
/// res_link_index. candidate_AAs specifies all of the AAs present in any
/// of poses_ at position res_link_index. Fitness is defined implicitly
/// as the sum of energy of all poses.
void FindConsensusSequence::pick_consensus_AA_mpi(
	core::pose::Pose & pose,
	utility::vector1<std::string> candidate_AAs,
	core::Size pose_position) {

	if ( debug_ ) {
		TR << "testing amino acids at position " << pose_position << std::endl;
	}

	utility::vector1<core::Real> candidate_AA_scores;
	simple_moves::MutateResidueOP mutation_mover;

	initialize_packer();

	for ( std::string const & trial_AA: candidate_AAs ) {
		candidate_AA_scores.push_back( test_AA( pose.clone(), trial_AA, pose_position) );
	}

	/// From all my candidate AA scores compile them to make a fitness
	utility::vector1< core::Real > candidate_AA_fitnesses (candidate_AA_scores.size());
	for ( core::Size ii = 1; ii <= candidate_AA_scores.size(); ++ii ) {
		// Let the master do all the work
		if ( master_ ) {
			candidate_AA_fitnesses[ ii ] = candidate_AA_scores[ ii ];
			for ( core::Size jj = 1; jj < n_procs_; ++jj ) {
				core::Real other_AA_score = utility::receive_double_from_node( jj );
				candidate_AA_fitnesses[ ii ] += other_AA_score;
			}

		} else {
			utility::send_double_to_node( 0, candidate_AA_scores[ ii ] );
		}
	}

	/// Tell everyone what the best AA is
	utility::vector1<core::Real>::iterator best_fitness_iter = std::min_element(candidate_AA_fitnesses.begin(), candidate_AA_fitnesses.end() );
	core::Real best_fitness = *best_fitness_iter;
	core::Size best_candidate_index = candidate_AA_fitnesses.index( best_fitness );
	std::string best_AA = candidate_AAs[ best_candidate_index ];

	if ( master_ ) {
		for ( core::Size ii = 2; ii <= n_procs_; ++ii ) {
			utility::send_string_to_node( ii-1, best_AA ); // nodes are 0-indexed
		}
	} else {
		best_AA = utility::receive_string_from_node( 0 ); // get the best AA from master
	}
	/// Thread best AA over my state
	mutation_mover = simple_moves::MutateResidueOP ( new simple_moves::MutateResidue (
		pose_position , //position
		best_AA //residue
		) );
	mutation_mover->apply( pose );
	packer_->apply( pose );
}

/// Place a candidate AA (trial_AA) onto a pose (pose_copy)
/// at position pose_position and repack and measure the energy.
/// Function returns the energy of that pose with the trial_AA
core::Real
FindConsensusSequence::test_AA( core::pose::PoseOP pose_copy,
	std::string trial_AA,
	core::Size pose_position ) {

	if ( debug_ ) {
		TR << "trying AA " << trial_AA << " at position "<< pose_position << std::endl;
	}
	simple_moves::MutateResidueOP mutation_mover = simple_moves::MutateResidueOP ( new simple_moves::MutateResidue (
		pose_position, //position
		trial_AA //residue
		) );

	mutation_mover->apply( *pose_copy );
	packer_->apply( *pose_copy );

	return sfxn_->score( *pose_copy );
}

/// Based on all the input poses, find the optimal AA at position
/// res_link_index. candidate_AAs specifies all of the AAs present in any
/// of poses_ at position res_link_index. Fitness is defined implicitly
/// as the sum of energy of all poses.
void
FindConsensusSequence::pick_consensus_AA( core::Size str_position,
	utility::vector1< std::string > candidate_AAs ) {

	if ( debug_ ) {
		TR << "picking the best residue at position ";
		for ( core::Size i = 1; i <= designable_residues_.size(); ++i ) {
			TR << designable_residues_[ i ][ str_position ];
			if ( i != designable_residues_.size() ) TR << " or ";
		}
		TR << std::endl;
	}

	utility::vector1< core::Real > candidate_AA_fitnesses (candidate_AAs.size());
	for ( core::Size ii = 1; ii <= candidate_AAs.size(); ++ii ) {
		std::string trial_AA = candidate_AAs[ ii ];
		core::Real trial_AA_fitness = 0;

		for ( core::Size pose_ind = 1; pose_ind <= poses_.size(); ++pose_ind ) {
			core::pose::PoseOP tmp_pose = poses_[ pose_ind ]->clone();
			core::Size pose_position = designable_residues_[ pose_ind ][ str_position ];

			initialize_packer();
			trial_AA_fitness += test_AA( tmp_pose, trial_AA, pose_position );
		}

		candidate_AA_fitnesses[ ii ] = trial_AA_fitness;
	}

	/// Find what the best AA is
	utility::vector1<core::Real>::iterator best_fitness_iter =
		std::min_element(candidate_AA_fitnesses.begin(), candidate_AA_fitnesses.end() );
	core::Real best_fitness = *best_fitness_iter;
	core::Size best_candidate_index = candidate_AA_fitnesses.index( best_fitness );
	std::string best_AA = candidate_AAs[ best_candidate_index ];

	// Now I've checked all my candidate AAs
	// thread the best one over my final poses
	if ( debug_ ) {
		TR << "best residue is " << best_AA << std::endl;
	}

	simple_moves::MutateResidueOP mutation_mover;
	for ( core::Size pose_ind = 1; pose_ind <= poses_.size(); ++pose_ind ) {

		// if it doesn't already have this residue thread it over=
		std::string current_AA = poses_[ pose_ind ]->residue( designable_residues_[ pose_ind ][ str_position ] ).name3();
		if (  current_AA != best_AA ) {
			mutation_mover = simple_moves::MutateResidueOP ( new simple_moves::MutateResidue (
				designable_residues_[ pose_ind ][ str_position ], //position
				best_AA //residue
				) );
			mutation_mover->apply( *poses_[ pose_ind ] );
			packer_->apply( *poses_[ pose_ind ] );
		}
	}
}

std::string FindConsensusSequence::get_name() const {
	return "FindConsensusSequence";
}

core::scoring::ScoreFunctionOP FindConsensusSequence::score_function() const {
	return sfxn_;
}

void FindConsensusSequence::score_function( core::scoring::ScoreFunctionOP sfxn) {
	sfxn_ = sfxn;
}

core::pack::task::TaskFactoryOP FindConsensusSequence::task_factory() const {
	return task_factory_;
}

void FindConsensusSequence::task_factory( core::pack::task::TaskFactoryOP task_factory ) {
	task_factory_ = task_factory;
}

void FindConsensusSequence::resfiles( utility::vector1< std::string > & resfiles ) {
	resfiles_ = resfiles;
}

utility::vector1< utility::vector1< core::Size > > FindConsensusSequence::designable_residues () { return designable_residues_; }

} //recon_design
} //protocols
