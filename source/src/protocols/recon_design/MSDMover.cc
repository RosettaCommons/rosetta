// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/recon_design/MSDMover.cc
/// @brief Multistate design mover used for RECON multistate design
/// Takes in multiple poses, applies residue linking constraints based on
/// sequence of all input poses and runs a design submover that has been specified in the tag
/// Only accessible through recon application.
/// @author Alex Sevy (alex.sevy@gmail.com)

// unit headers
#include <protocols/recon_design/MSDMover.hh>
#include <protocols/recon_design/MSDMoverCreator.hh>

// type headers
//#include <core/types.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/VectorPoseMover.hh>
#include <protocols/rosetta_scripts/util.hh>

#include <core/pose/Pose.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/constraints/ResidueTypeConstraint.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>


#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/tag/Tag.hh>

#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>

#include <protocols/minimization_packing/PackRotamersMover.hh>

#include <sstream>
#include <utility/string_util.hh>
#include <utility/mpi_util.hh>

#include <protocols/recon_design/recon_util.hh>

#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

namespace protocols {
namespace recon_design {

static basic::Tracer TR("protocols.recon_design.MSDMover");

MSDMover::MSDMover() :
	moves::VectorPoseMover( MSDMoverCreator::mover_name() ),
	weight_(0.5),
	current_pose_( 1 ),
	debug_( false )
{}

MSDMover::MSDMover( protocols::moves::MoverOP mover, utility::vector1< std::string > const & resfiles ) :
	design_mover_( mover ),
	resfiles_( resfiles ),
	weight_(0.5),
	current_pose_( 1 ),
	debug_( false )
{}

MSDMover::~MSDMover() {}

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
MSDMoverCreator::mover_name() {
	return "MSDMover";
}

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

void MSDMover::parse_my_tag(
	utility::tag::TagCOP const tag,
	basic::datacache::DataMap & data
) {

	if ( tag->hasOption("design_mover") ) {
		std::string const design_mover_key = tag->getOption<std::string>("design_mover");
		design_mover_ = protocols::rosetta_scripts::parse_mover( design_mover_key, data );
	} else { // throw exception if there's no design mover specified
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,"Error: MSD Mover mustohave design mover");
	}

	if ( tag->hasOption("post_mover") ) {
		std::string const post_mover_key = tag->getOption<std::string>("post_mover");
		post_mover_ = protocols::rosetta_scripts::parse_mover( post_mover_key, data );
	}

	if ( tag->hasOption( "resfiles" ) ) {
		resfiles_ = utility::string_split( tag->getOption<std::string>( "resfiles" ), ',' );
	} else {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,"Error: must provide resfile in MSD Mover tag" );
	}

	if ( resfiles_.size() > 1 && design_mover_->get_name() != "PackRotamersMover" ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "Error: multiple resfiles only supported for PackRotamersMover" );
	}


	weight_ = tag->getOption<core::Real>( "constraint_weight", 1.0 );

	debug_ = tag->getOption<bool>( "debug", false );

#ifndef NDEBUG
	debug_ = true;
#endif

}

void MSDMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	attlist
		+ XMLSchemaAttribute( "design_mover", xs_string, "A previously defined mover that will design the input states. Note that since resfiles are applied as a TaskOperation within the MSDMover, there's no need to specify ReadResfile behavior for the design_mover - however all other desired TaskOperations (InitializeFromCommandLine, etc) should be specified in the design_mover tag")
		+ XMLSchemaAttribute( "post_mover", xs_string, "A previously defined mover that will act on the input states after the design step.")
		+ XMLSchemaAttribute( "resfiles", xs_string, "A comma-separated list of resfiles to define designable and repackable residues for all states in multistate design. Multiple resfiles can be used for multiple states - in this case the first resfile in the tag will be applied to the first structure, etc. One single resfile used for all states is also supported.")
		+ XMLSchemaAttribute::attribute_w_default( "constraint_weight", xsct_real, "The weight of amino acid linking constraints during the RECON protocol. Generally weights will be ramped from 0.5 to 2.0, to allow searching of more sequence space in early rounds.", "1.0")
		+ XMLSchemaAttribute::attribute_w_default( "debug", xsct_rosetta_bool, "if true, outputs extra messages during protocol", "false");

	protocols::moves::xsd_type_definition_w_attributes( xsd, MSDMoverCreator::mover_name(), "Multistate design mover used for RECON protocol", attlist );
}

void MSDMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	MSDMover::provide_xml_schema( xsd );
}

void
MSDMover::apply( core::pose::Pose & pose ) {

	using namespace core::scoring::constraints;

	// Set up class variables
	setup_mover( pose );

	// Find the sequence of the other poses
	utility::vector1< utility::vector1< std::string > > other_pose_sequences;
	for  ( core::Size ii = 1; ii <= poses_.size(); ++ii ) {
		if ( ii != current_pose_ ) {
			other_pose_sequences.push_back(
				get_designable_sequence ( *poses_[ii], designable_residues_[ ii ] )
			);
		}
	}

	// Add linking constraints to pose
	// utility::vector1< ConstraintCOP > constraints = apply_linked_constraints( pose );
	utility::vector1< ConstraintCOP > constraints = apply_linked_constraints(
		pose, other_pose_sequences, designable_residues_[ current_pose_ ]
	);



	if ( debug_ ) {
		utility::vector1< ConstraintCOP > constraint_set = pose.constraint_set()->get_all_constraints();
		TR << "this pose has " << constraint_set.size() << " constraints." << std::endl;
		for ( ConstraintCOP constraint: constraint_set ) {
			TR << constraint->to_string() << std::endl;
		}
	}

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

	// If there's a post-relaxation mover, apply it
	run_post_design_mover( pose );

}

/// Based on the sequence of the other poses, apply a residue type constraint to
/// encourage poses to adopt same sequence. Returns a COP of the constraints that were
/// applied so they can be removed later
utility::vector1< core::scoring::constraints::ConstraintCOP >
MSDMover::apply_linked_constraints( core::pose::Pose & pose,
	utility::vector1< utility::vector1< std::string > > other_pose_sequences,
	utility::vector1< core::Size > my_designable_residues ) {

	using namespace core::scoring::constraints;
	utility::vector1< ConstraintCOP > constraints;

	// If there is no weight then there's no point in adding constraints
	if ( weight_ == 0 ) return constraints;

	// Build up a list of candidate AAs at each position to be constrained
	utility::vector1< std::string > candidate_AAs;

	// Iterate through each position that's being designed
	for ( core::Size jj = 1; jj <= other_pose_sequences[ 1 ].size(); ++jj ) {
		core::Size pose_position = my_designable_residues[ jj ];
		candidate_AAs.clear(); // TODO check if this is unnecessary
		candidate_AAs = get_candidate_AAs( other_pose_sequences, jj );

		// Now I have all my candidate AAs
		// constrain my pose to these AAs
		for ( core::Size ii = 1; ii <= candidate_AAs.size(); ++ii ) {
			std::string constrained_AA = candidate_AAs[ ii ];
			if ( debug_ ) {
				TR << "adding linking constraint to position " <<
					pose_position <<
					" for residue " <<
					constrained_AA <<
					" for weight " <<
					weight_ <<
					std::endl;
			}
			ResidueTypeConstraintCOP temp_cst ( new ResidueTypeConstraint(
				pose,
				pose_position, //seqpos
				constrained_AA, //AAname
				weight_, //constraint weight
				true // using base_name
				) );
			constraints.push_back( pose.add_constraint( temp_cst ) );
		}
	}
	return constraints;
}

void
MSDMover::apply_mpi( core::pose::Pose & pose ) {
	using namespace core::pack::task;
	using namespace core::scoring::constraints;

	/// Get my node status first
	core::Size rank = utility::mpi_rank()+1; // I want to make the rank 1-indexed
	core::Size n_procs = utility::mpi_nprocs();
	bool master = (rank==1);

	current_pose_ = rank;

	/// Get my resfile
	std::string this_nodes_resfile = resfile_at( rank );

	update_packer_task();

	/// Find out my designable residues from my pose and my resfile
	utility::vector1< core::Size > my_designable_residues = get_designable_residues( pose, this_nodes_resfile );

	/// Make a string vector out of the AAs at my designable positions in the current state
	utility::vector1< std::string > my_sequence = get_designable_sequence ( pose, my_designable_residues );

	std::string pass_seq;
	for ( const std::string& resi_base_name: my_sequence ) {
		pass_seq += resi_base_name + " ";
	}

	/// Get the AAs at designable positions of the other states I need to cooperate with
	utility::vector1< utility::vector1<std::string> > other_pose_sequences( n_procs );
	for ( core::Size ii = 1; ii <= n_procs; ++ii ) {
		if ( rank == ii ) {
			for ( core::Size jj = 1; jj <= n_procs; ++jj ) {
				if ( rank!=jj ) utility::send_string_to_node( jj-1, pass_seq ); // node ranks are 0-indexed
				else other_pose_sequences[jj] = my_sequence;
			}
		} else {
			std::string passed_seq = utility::receive_string_from_node( ii-1 ); // node ranks are 0-indexed
			//Need to split passed_seq by spaces
			utility::vector1<std::string> split_seq = utility::split_whitespace(passed_seq);
			other_pose_sequences[ii].append(split_seq);
		}
	}

	/// Let the master make sure all the sequences are the same length,
	/// i.e. all the states have same number of designable residues
	if ( master ) {
		for ( utility::vector1< std::string > const & sequence: other_pose_sequences ) {
			if ( sequence.size() != my_sequence.size() ) {
				utility_exit_with_message( "Error: all states must have the same number of designable residues" );
			}
		}
	}

	if ( debug_ ) {
		TR << "my sequence is ";
		for ( core::Size i = 1; i <= my_sequence.size(); ++i ) {
			TR << my_sequence[i] << " ";
		}
		TR << std::endl;
	}

	/// Run my multistate design
	utility::vector1< ConstraintCOP > constraints =
		apply_linked_constraints( pose, other_pose_sequences, my_designable_residues );
	design_mover_->apply( pose );

	if ( pose.energies().weights()[core::scoring::res_type_constraint] == 0 ) {
		TR.Warning << "Warning: res type constraint is set at zero. You need to reweight it for multistate design to work properly" << std::endl;
	}

	pose.remove_constraints( constraints, true );

	if ( debug_ ) {
		TR << "fitness of pose " << pose << ": " << pose.energies().total_energy() << std::endl;
	}


	run_post_design_mover( pose );

}

/// Helper functions

/// Initialize mover by checking that input poses were passed correctly,
/// a design mover was specified, and finding the pose given to apply() in the
/// poses_ vector
void MSDMover::setup_mover ( core::pose::Pose & pose ) {
	// Set up movers, resfiles, packer tasks, throw exceptions, etc.
	if ( poses_.size() == 0 ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,"Error: no poses initialized. If you run MSDMover from RosettaScripts you need to pass the -run:msd_job_dist flag");
	}

	if ( !design_mover_ ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,"Error: MSD Mover must have design mover");
	}

	// Find index of your apply pose in the pose vector
	current_pose_ = find_pose_in_vector( pose, poses_ );
	update_packer_task();
	parse_resfiles();
}

/// Populates designable_residues_ with a list of designable residues corresponding
/// to the poses in poses_, corresponding element-wise (designables_residues[0] matches
/// to poses_[0], etc)
void MSDMover::parse_resfiles () {
	designable_residues_.clear();

	for ( core::Size ii = 1; ii <= poses_.size(); ++ii ) {
		designable_residues_.push_back( get_designable_residues( *poses_[ii], resfile_at(ii) ) );
		if ( designable_residues_[ ii ].size() != designable_residues_[1].size() ) {
			utility_exit_with_message( "All resfiles must have the same number of designable residues");
		}
	}
}

/// The design mover can't be given its task operations on a global
/// level in the RosettaScript, bc each of the input poses can have
/// different designable residues. This function assigns the design mover
/// its tasks based on the designable residues of the current pose
void MSDMover::update_packer_task () {
	// if my submover is pack rotamers, allow it to take in multiple resfiles from the tag
	if ( design_mover_->get_name() == "PackRotamersMover" ) {
		minimization_packing::PackRotamersMoverOP design_packer = utility::pointer::dynamic_pointer_cast< minimization_packing::PackRotamersMover >( design_mover_ );

		// If mover already has task factory modify it -> otherwise create a new one
		core::pack::task::TaskFactoryOP factory;
		if ( design_packer->task_factory() == nullptr ) {
			factory = core::pack::task::TaskFactoryOP ( new core::pack::task::TaskFactory );
		} else {
			factory =  design_packer->task_factory()->clone();
		}

		using namespace core::pack::task::operation;
		ReadResfileOP rrf ( new ReadResfile( resfile_at ( current_pose_ ) ) );
		factory->push_back( rrf );
		design_packer->task_factory( factory );
		//design_mover_ = design_packer;
	}
}

/// Get the resfile corresponding to the pose at index. If
/// only one resfile is present then it will be returned regardless
/// of the value of index.
std::string MSDMover::resfile_at ( core::Size index ) {
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

void MSDMover::run_post_design_mover ( core::pose::Pose & pose ) {
	if ( post_mover_ ) {
		TR << "Design finished - applying post design mover " << post_mover_->get_name() << std::endl;
		post_mover_->apply( pose );
	}
}

// Getters and setters
moves::MoverOP MSDMover::design_mover() { return design_mover_; }

void MSDMover::design_mover( moves::MoverOP design_mover ) { design_mover_ = design_mover; }

utility::vector1< std::string > MSDMover::resfiles() { return resfiles_; }

void MSDMover::resfiles ( utility::vector1< std::string > resfiles ) { resfiles_ = resfiles; }

core::Real MSDMover::weight () { return weight_; }

void MSDMover::weight( core::Real weight ) { weight_ = weight; }

void MSDMover::set_current_pose( core::Size current_pose ) { current_pose_ = current_pose; }

utility::vector1< utility::vector1< core::Size > > MSDMover::designable_residues () { return designable_residues_; }


} // recon_design
} // protocols
