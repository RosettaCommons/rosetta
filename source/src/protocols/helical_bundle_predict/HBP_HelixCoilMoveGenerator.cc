// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/helical_bundle_predict/HBP_HelixCoilMoveGenerator.cc
/// @brief A class for a module to generate ParsedProtocols for the next move in a Monte Carlo trajectory, based on the
/// current state of the pose.  This version uses helix-coil transition theory to nucleate and extend helices.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

// Project headers:
#include <protocols/helical_bundle_predict/HBP_HelixCoilMoveGenerator.hh>

// Core headers:
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/conformation/Residue.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>

// Protocols headers:
#include <protocols/rosetta_scripts/ParsedProtocol.hh>
#include <protocols/simple_moves/bb_sampler/SmallBBSampler.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/TorsionSetMover.hh>
#include <protocols/backbone_moves/RandomizeBBByRamaPrePro.hh>

// Basic headers:
#include <basic/Tracer.hh>

// Utility headers:
#include <utility/pointer/memory.hh>

static basic::Tracer TR( "protocols.helical_bundle_predict.HBP_HelixCoilMoveGenerator" );


namespace protocols {
namespace helical_bundle_predict {

/// @brief Constructor.
HBP_HelixCoilMoveGenerator::HBP_HelixCoilMoveGenerator():
	HBP_MoveGenerator(),
	user_helix_assignments_(),
	current_helix_assignments_(),
	candidate_helix_assignments_(),
	ramaprepro_sfxn_( utility::pointer::make_shared< core::scoring::ScoreFunction >() )
{
	ramaprepro_sfxn_->set_weight( core::scoring::rama_prepro, 1.0 );
}

/// @brief Copy constructor.
HBP_HelixCoilMoveGenerator::HBP_HelixCoilMoveGenerator( HBP_HelixCoilMoveGenerator const &src ) :
	HBP_MoveGenerator(src),
	user_helix_assignments_(src.user_helix_assignments_),
	current_helix_assignments_(src.current_helix_assignments_),
	candidate_helix_assignments_(src.candidate_helix_assignments_),
	ramaprepro_sfxn_( src.ramaprepro_sfxn_ == nullptr ? nullptr : src.ramaprepro_sfxn_->clone() )
{}

/// @brief Destructor.
HBP_HelixCoilMoveGenerator::~HBP_HelixCoilMoveGenerator() {}

/// @brief Clone this object: that is, make a copy and return an owning pointer to the copy.
HBP_MoveGeneratorOP
HBP_HelixCoilMoveGenerator::clone() const {
	return HBP_MoveGeneratorOP( utility::pointer::make_shared< HBP_HelixCoilMoveGenerator >( *this ) );
}

//////////////////////////////////////////// PUBLIC METHODS ////////////////////////////////////////////

/// @brief Given the current step index, the total number of steps in the trajectory, and the pose for analysis,
/// construct a ParsedProtocol of things to do to this pose for this move in the Monte Carlo trajectory.
/// @details This override uses helix-coil transition theory to nucleate, extend, or remove helix segments.  Non-helix
/// segments are simply subjected to small moves.
/// @returns A ParsedProtocolOP with the new protocol for success, or nullptr for failure.
protocols::rosetta_scripts::ParsedProtocolOP
HBP_HelixCoilMoveGenerator::generate_monte_carlo_move(
	core::Size const ,//current_step,
	core::Size const ,//num_steps,
	core::pose::Pose const &pose
) const {
	candidate_helix_assignments_ = current_helix_assignments_;

	protocols::rosetta_scripts::ParsedProtocolOP protocol( utility::pointer::make_shared< protocols::rosetta_scripts::ParsedProtocol >() );

	//Add helix-shrinking changes to the helix assignments:
	add_helix_retraction_moves( pose );

	//Add helix elongation changes to the helix assignments:
	add_helix_elongation_moves( pose );

	//Add helix nucleation changes to the helix assignments:
	add_helix_nucleation_moves( pose );

	//Add global parameter perturbations:
	add_global_parameter_perturbations( pose );

	//Add local parameter perturbations:
	add_local_parameter_perturbations( pose );

	//Add moves to the protocol to update helices based on helix assignments:
	if ( add_helix_update_moves( pose, *protocol ) ) return nullptr; //Return nullptr if we couldn't generate a helix due to bad parameters.

	//Small moves on non-helix stretches.  (TODO: ensure that this only affects non-helix stretches.)
	add_nonhelix_small_moves( pose, *protocol );

	return protocol;
}

/// @brief When a move is accepted, set the current helix assignments to the candidate helix assignments.
void
HBP_HelixCoilMoveGenerator::mark_move_accepted() const {
	current_helix_assignments_ = candidate_helix_assignments_;
}


/// @brief Given a helix assignment file's contents, set up the helix assignments.
/// @details Involves no read from disk.
void
HBP_HelixCoilMoveGenerator::set_up_user_helix_assignments(
	std::string const &file_contents
) {
	user_helix_assignments_.clear();
	user_helix_assignments_.initialize_from_file_contents( file_contents );
}

//////////////////////////////////////////// PRIVATE METHODS ////////////////////////////////////////////


/// @brief Update the helix assignments to nucleate helices.
void
HBP_HelixCoilMoveGenerator::add_helix_nucleation_moves(
	core::pose::Pose const &pose
) const {
	for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) {
		candidate_helix_assignments_.try_nucleating_helix(i, pose, user_helix_assignments_);
	}
}

/// @brief Update the helix assignments to elongate helices.
void
HBP_HelixCoilMoveGenerator::add_helix_elongation_moves(
	core::pose::Pose const &pose
) const {
	for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) {
		candidate_helix_assignments_.try_elongating_helix(i, pose, user_helix_assignments_);
	}
}

/// @brief Update the helix assignments to shrink helices.
void
HBP_HelixCoilMoveGenerator::add_helix_retraction_moves(
	core::pose::Pose const &pose
) const {
	for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) {
		candidate_helix_assignments_.try_retracting_helix(i, pose, user_helix_assignments_);
	}
}

/// @brief Add global parameter perturbations.
void
HBP_HelixCoilMoveGenerator::add_global_parameter_perturbations(
	core::pose::Pose const &//pose
) const {
	candidate_helix_assignments_.try_perturbing_global_helix_parameters(); //Updates locals from globals.
}

/// @brief Add local parameter perturbations.
void
HBP_HelixCoilMoveGenerator::add_local_parameter_perturbations(
	core::pose::Pose const &//pose
) const {
	candidate_helix_assignments_.try_perturbing_local_helix_parameters();
}

/// @brief Add movers to the move list to update the helices based on the helix assignments.
/// @details Returns true for FAILURE (due to inability to generate helix from parameters), false for success.
bool
HBP_HelixCoilMoveGenerator::add_helix_update_moves(
	core::pose::Pose const & pose,
	protocols::rosetta_scripts::ParsedProtocol &protocol
) const {
	std::stringstream helix_residues;
	core::select::residue_selector::ResidueIndexSelectorOP selector( utility::pointer::make_shared< core::select::residue_selector::ResidueIndexSelector >() );

	using namespace protocols::simple_moves;
	using namespace protocols::backbone_moves;

	// Indices in the following vector correspond to helix ids:
	utility::vector1< std::map< core::id::TorsionID, core::Real > > torsion_values_for_helices;
	if ( candidate_helix_assignments_.generate_torsion_values_for_helices( torsion_values_for_helices, pose ) ) return true; //Abort with failure signal if we can't generate a helix.

	utility::vector1< core::id::TorsionID > ids;
	utility::vector1< core::Real > values;
	bool at_least_one_frayed_position(false);
	for ( core::Size ir(1), irmax(pose.total_residue()); ir<=irmax; ++ir ) {
		if ( candidate_helix_assignments_.is_in_helix(ir) ) {
			core::Size const curhelix( candidate_helix_assignments_.get_containing_helix_index(ir) );
			core::conformation::Residue const & curres( pose.residue(ir) );
			std::map< core::id::TorsionID, core::Real > const & curmap( torsion_values_for_helices[curhelix] );

			for ( core::Size itors(1), itorsmax( curres.mainchain_torsions().size() ); itors<=itorsmax; ++itors ) {
				core::id::TorsionID const curtors( ir, core::id::BB, itors );
				ids.push_back( curtors );
				runtime_assert_string_msg( curmap.count(curtors) != 0, "Residue torsion mismatch in helix " + std::to_string( curhelix ) + " (at residue " + std::to_string(ir) + ")." );
				values.push_back( curmap.at(curtors) );
			}

			if ( TR.visible() ) helix_residues << "H";
		} else {
			if ( TR.visible() ) helix_residues << "-";
			if ( current_helix_assignments_.is_in_helix(ir) ) {
				//Randomize positions that were helix positions, which have frayed.
				selector->append_index(ir);
				at_least_one_frayed_position = true;
			}
		}
	}
	debug_assert( ids.size() == values.size() );

	// Set torsions for helix positions
	if ( ids.size() != 0 ) {
		TorsionSetMoverOP torsset( utility::pointer::make_shared< TorsionSetMover >(ids, values) );
		protocol.add_mover_filter_pair( torsset, "Update_Helices", nullptr, false);
	}

	// Randomize positions that used to be helical, but have frayed:
	if ( at_least_one_frayed_position ) {
		RandomizeBBByRamaPreProOP randomize( utility::pointer::make_shared< RandomizeBBByRamaPrePro >() );
		randomize->set_residue_selector(selector);
		protocol.add_mover_filter_pair( randomize, "Randomize_Frayed_Positions", nullptr, false );
	}

	if ( TR.visible() ) { TR << "Helix pattern: " << helix_residues.str() << std::endl; }

	return false; //Success.
}

/// @brief Add the SmallMover moves to the protocol, at non-helix positions.
void
HBP_HelixCoilMoveGenerator::add_nonhelix_small_moves(
	core::pose::Pose const &pose,
	protocols::rosetta_scripts::ParsedProtocol & protocol
) const {
	core::kinematics::MoveMapOP movemap( utility::pointer::make_shared< core::kinematics::MoveMap >() );
	movemap->set_bb(true); //Update this.
	for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) {
		if ( candidate_helix_assignments_.is_in_helix(i) ) movemap->set_bb(i, false); //Disable at helix positions.
	}
	protocols::simple_moves::SmallMoverOP smallmover( utility::pointer::make_shared< protocols::simple_moves::SmallMover >( movemap, 5.0, 500 ) );
	smallmover->scorefxn( ramaprepro_sfxn_ );
	smallmover->angle_max( 15 );
	protocol.add_mover_filter_pair( smallmover, "Small", nullptr, false );
}


} //protocols
} //helical_bundle_predict






