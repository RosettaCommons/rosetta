// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/MSDMover.cc
/// @brief Multistate design mover used for restrained multistate design
/// @brief Takes in multiple poses from MSDJobDistributor, applies residue linking constraints based on
/// @brief sequence of all input poses and runs a design submover that has been specified in the tag
/// @author Alex Sevy (alex.sevy@gmail.com)

// unit headers
#include <protocols/simple_moves/MSDMover.hh>
#include <protocols/simple_moves/MSDMoverCreator.hh>

// type headers
#include <core/types.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/VectorPoseMover.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/tag/Tag.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/ResfileReader.hh>
//#include <core/scoring/constraints/ConstraintSet.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <core/scoring/constraints/ResidueTypeConstraint.hh>
#include <core/chemical/ResidueType.hh>
#include <utility/string_util.hh>
#include <core/pose/util.hh>
//#include <protocols/simple_moves/MutateResidue.hh>
//#include <core/chemical/AA.hh>


namespace protocols {
namespace simple_moves {

static basic::Tracer TR("protocols.simple_moves.MSDMover");

std::string
MSDMoverCreator::keyname() const
{
	return MSDMoverCreator::mover_name();
}

moves::MoverOP
MSDMoverCreator::create_mover() const {
	return MSDMoverOP( new MSDMover );
}

std::string
MSDMoverCreator::mover_name()
{
	return "MSDMover";
}

MSDMover::MSDMover() :
	moves::VectorPoseMover( MSDMoverCreator::mover_name() ),
	weight_(0.5),
	current_pose_( 1 ),
	debug_( false )
{}

MSDMover::MSDMover( protocols::moves::MoverOP mover, utility::vector1< std::string > resfiles ) :
	design_mover_( mover ),
	resfiles_( resfiles ),
	weight_(0.5),
	current_pose_( 1 ),
	debug_( false )
{}

MSDMover::~MSDMover() {}

moves::MoverOP MSDMover::clone() const {
	return MSDMoverOP ( new MSDMover( *this ) );
}
moves::MoverOP MSDMover::fresh_instance() const {
	return MSDMoverOP ( new MSDMover );
}

std::string
MSDMover::get_name() const {
	return MSDMoverCreator::mover_name();
}

void
MSDMover::apply( core::pose::Pose & pose ) {

	// Set up class variables
	setup_mover( pose );

	// Add linking constraints to pose
	utility::vector1< core::scoring::constraints::ConstraintCOP > constraints = apply_linked_constraints( pose );

	// Apply design mover
	design_mover_->apply( pose );

	if ( pose.energies().weights()[core::scoring::res_type_constraint] == 0 ) {
		TR.Warning << "Warning: res type constraint is set at zero. You need to reweight it for multistate design to work properly" << std::endl;
	}
	// Remove constraints from pose
	pose.remove_constraints( constraints, true );

	if ( debug_ ) {
		TR << "fitness of pose " << current_pose_ << ": " << pose.energies().total_energy() << std::endl;
	}

}

void
MSDMover::setup_mover ( core::pose::Pose & pose ) {
	// Set up movers, resfiles, packer tasks, throw exceptions, etc.
	if ( poses_.size() == 0 ) {
		throw utility::excn::EXCN_RosettaScriptsOption("Error: no poses initialized. If you run MSDMover from RosettaScripts you need to pass the -run:msd_job_dist flag");
	}

	if ( !design_mover_ ) {
		throw utility::excn::EXCN_RosettaScriptsOption("Error: MSD Mover must have design mover");
	}

	// Find index of your apply pose in the pose vector
	current_pose_ = 0;
	core::Real ref_pose_index;
	bool return_value = core::pose::getPoseExtraScore(  pose, "msd_job_dist_index", ref_pose_index );
	if ( !return_value ) {
		throw utility::excn::EXCN_RosettaScriptsOption("Error: poses are not indexed correctly. If you run MSDMover from RosettaScripts you need to pass the -run:msd_job_dist flag");
	}
	for ( core::Size i = 1; i <= poses_.size(); ++i ) {
		core::Real current_pose_index;
		core::pose::getPoseExtraScore( *poses_[ i ], "msd_job_dist_index", current_pose_index );

		if ( ref_pose_index == current_pose_index ) {
			current_pose_ = i;
			break;
		}
	}

	if ( current_pose_ == 0 ) {
		throw utility::excn::EXCN_RosettaScriptsOption("Error: pose not found in vector. Strange...");
	}
	update_packer_task();
	parse_resfiles();

}


void MSDMover::parse_my_tag(
	utility::tag::TagCOP const tag,
	basic::datacache::DataMap & ,//datamap
	filters::Filters_map const &,
	moves::Movers_map const & movers,
	core::pose::Pose const &
) {

	if ( tag->hasOption("design_mover") ) {
		std::string const design_mover_key = tag->getOption<std::string>("design_mover");
		moves::Movers_map::const_iterator  find_mover = movers.find( design_mover_key );

		if ( find_mover == movers.end() && design_mover_key != "" ) {
			throw utility::excn::EXCN_RosettaScriptsOption("Mover " + design_mover_key + " not found in data map");
		}

		design_mover_ =find_mover->second;

	} else { // throw exception if there's no design mover specified
		throw utility::excn::EXCN_RosettaScriptsOption("Error: MSD Mover must have design mover");
	}

	if ( tag->hasOption( "resfiles" ) ) {
		resfiles_ = utility::string_split( tag->getOption<std::string>( "resfiles" ), ',' );
	} else {
		throw utility::excn::EXCN_RosettaScriptsOption( "Error: must provide resfile in MSD Mover tag" );
	}

	if ( resfiles_.size() > 1 && design_mover_->get_name() != "PackRotamersMover" ) {
		throw utility::excn::EXCN_RosettaScriptsOption( "Error: multiple resfiles only supported for PackRotamersMover" );
	}


	weight_ = tag->getOption<core::Real>( "constraint_weight", 1.0 );

	debug_ = tag->getOption<bool>( "debug", false );

#ifndef NDEBUG
	debug_ = true;
#endif

}

utility::vector1< core::Size >
MSDMover::parse_resfile ( core::pack::task::PackerTaskCOP design_task )
{
	utility::vector1< core::Size > vector;
	utility::vector1<bool> designing = design_task->designing_residues();
	for ( core::Size i = 1; i <= designing.size(); ++i ) {
		if ( designing[ i ] ) {
			vector.push_back( i );
		}
	}
	return vector;
}

utility::vector1< core::scoring::constraints::ConstraintCOP >
MSDMover::apply_linked_constraints( core::pose::Pose & pose )
{
	using namespace core::scoring::constraints;
	utility::vector1< ConstraintCOP > constraints;
	if ( weight_ != 0 ) { // If we are applying a constraint weight
		/* Iterate through poses in the outer loop */
		for ( core::Size i = 1; i <= poses_.size(); ++i ) {
			if ( i == current_pose_ ) continue;

			/* Iterate through residues to be linked */
			for ( core::Size k = 1; k <= res_links_[ current_pose_ ].size(); ++k ) {

				if ( debug_ ) {
					TR << "adding linking constraint to position " <<
						res_links_[ current_pose_ ][ k ] <<
						" for residue " <<
						poses_[ i ]->residue_type( res_links_[ i ][ k ] ).name3() <<
						" of magnitude " << weight_ << std::endl;
				}
				ResidueTypeConstraintCOP temp_cst ( new ResidueTypeConstraint(
					pose,
					res_links_[ current_pose_ ][ k ], //seqpos
					poses_[ i ]->residue_type( res_links_[ i ][ k ] ).name3(), //AAname
					weight_ //favor native bonus
					) );
				// When I add a constraint to the pose it returns a const copy - I return this vector so that later
				// I can remove these specific constraints
				constraints.push_back( pose.add_constraint( temp_cst ) );

			}
		}
	}
	return constraints;
}

void
MSDMover::parse_resfiles ()
{
	res_links_.clear();
	utility::vector1< core::pack::task::PackerTaskOP > design_tasks ( poses_.size() );
	for ( core::Size ii = 1; ii <= poses_.size(); ++ii ) {
		design_tasks[ ii ] = core::pack::task::TaskFactory::create_packer_task( *poses_[ ii ] );
		core::pack::task::parse_resfile( *poses_[ ii ], *design_tasks[ ii ], resfile_at( ii ) );

	}
	core::Size designable_residues = -1;
	for ( core::Size i = 1; i <= design_tasks.size(); ++i ) {
		res_links_.push_back( parse_resfile( design_tasks[ i ]) );
		if ( i == 1 )  designable_residues = res_links_[ 1 ].size();
		else {
			if ( res_links_[ i ].size() != designable_residues ) {
				utility_exit_with_message( "All resfiles must have the same number of designable residues");
			}
		}
	}
}

void
MSDMover::update_packer_task () {
	// if my submover is pack rotamers, allow it to take in multiple resfiles from the tag
	if ( design_mover_->get_name() == "PackRotamersMover" ) {
		simple_moves::PackRotamersMoverOP design_packer = utility::pointer::dynamic_pointer_cast< simple_moves::PackRotamersMover >( design_mover_ );

		// If mover already has task factory modify it -> otherwise create a new one
		core::pack::task::TaskFactoryOP factory;
		if ( design_packer->task_factory() == 0 ) {
			factory = core::pack::task::TaskFactoryOP ( new core::pack::task::TaskFactory );
		} else {
			factory =  design_packer->task_factory()->clone();
		}
		core::pack::task::operation::ReadResfileOP rrf ( new core::pack::task::operation::ReadResfile( resfile_at ( current_pose_ ) ) );
		factory->push_back( rrf );
		design_packer->task_factory( factory );
		//design_mover_ = design_packer;
	}
}

std::string
MSDMover::resfile_at ( core::Size index ) {
	if ( resfiles_.size() == 1 ) {
		return resfiles_[ 1 ];
	} else if ( resfiles_.size() == 0 ) {
		utility_exit_with_message("No resfiles set - make sure to set resfiles before updating packer task");
	} else {
		if ( index > resfiles_.size() ) {
			utility_exit_with_message("Number of resfiles does not match number of structures");
		} else {
			return resfiles_[ index ];
		}
	}
}

moves::MoverOP MSDMover::design_mover() { return design_mover_; }

void MSDMover::design_mover( moves::MoverOP design_mover ) { design_mover_ = design_mover; }

utility::vector1< std::string > MSDMover::resfiles() { return resfiles_; }

void MSDMover::resfiles ( utility::vector1< std::string > resfiles ) { resfiles_ = resfiles; }

core::Real MSDMover::weight () { return weight_; }

void MSDMover::weight( core::Real weight ) { weight_ = weight; }

void MSDMover::set_current_pose( core::Size current_pose ) { current_pose_ = current_pose; }

utility::vector1< utility::vector1< core::Size > > MSDMover::res_links () { return res_links_; }


} // simple_moves
} // protocols
