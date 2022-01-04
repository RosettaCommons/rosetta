// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/CycpepRigidBodyPermutationMover.cc
/// @brief Implementations of functions for a mover that takes a cyclic peptide and alters its position (rigid-body transform), superimposing
/// it on a cyclic permutation or reverse cyclic permutation of itself.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

// Unit headers
#include <protocols/cyclic_peptide/CycpepRigidBodyPermutationMover.hh>
#include <protocols/cyclic_peptide/CycpepRigidBodyPermutationMoverCreator.hh>

// Core headers
#include <core/select/residue_selector/util.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/pose/Pose.hh>
#include <core/id/AtomID.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/scoring/rms_util.hh>

// Protocols headers
#include <protocols/moves/mover_schemas.hh>
#include <protocols/rosetta_scripts/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/pointer/memory.hh>
#include <utility/vector1.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

// Numeric headers
#include <numeric/random/random.hh>
#include <numeric/random/random_xyz.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/xyz.functions.hh>

// Citation Manager
#include <basic/citation_manager/UnpublishedModuleInfo.hh>

// STL headers
#include <map>

static basic::Tracer TR( "protocols.cyclic_peptide.CycpepRigidBodyPermutationMover" );

namespace protocols {
namespace cyclic_peptide {

////////////////////////////////////////////////////////////////////////////////
// ENUM AND UTILITY FUNCTIONS
////////////////////////////////////////////////////////////////////////////////

/// @brief Given a mover mode string, get the mode enum.
/// @details Returns CycpepRigidBodyPermutationMoverMode::UNKNOWN_MODE if string not recognized.
CycpepRigidBodyPermutationMoverMode
cycpep_rigid_body_permutation_mover_mode_from_string(
	std::string const & modestring
) {
	for ( core::Size i(1); i <= static_cast< core::Size >(CycpepRigidBodyPermutationMoverMode::N_MODES); ++i ) {
		if ( string_from_cycpep_rigid_body_permutation_mover_mode( static_cast< CycpepRigidBodyPermutationMoverMode >(i) ) == modestring ) {
			return static_cast< CycpepRigidBodyPermutationMoverMode >(i);
		}
	}
	return CycpepRigidBodyPermutationMoverMode::UNKNOWN_MODE;
}

/// @brief Given a mover mode enum, get the mode string.
std::string
string_from_cycpep_rigid_body_permutation_mover_mode(
	CycpepRigidBodyPermutationMoverMode const mode
) {
	runtime_assert( static_cast< core::Size >(mode) > 0 && mode <= CycpepRigidBodyPermutationMoverMode::N_MODES );
	switch(mode) {
	case CycpepRigidBodyPermutationMoverMode::SET_PERMUTATION :
		return "set_permutation";
	case CycpepRigidBodyPermutationMoverMode::RANDOM_PERMUTATION :
		return "randomized_permutation";
	default :
		utility_exit_with_message( "Error in protocols::cyclic_peptide::string_from_cycpep_rigid_body_permutation_mover_mode(): Unknown mode!" );
	}
	return "ERROR"; //Should never reach here, but keeps older compilers happy.
}

/// @brief Get a comma-separated list of valid mover modes.
std::string
list_cycpep_rigid_body_permutation_mover_modes() {
	std::stringstream ss;
	for ( core::Size i(1); i <= static_cast< core::Size >(CycpepRigidBodyPermutationMoverMode::N_MODES); ++i ) {
		if ( static_cast< core::Size >(CycpepRigidBodyPermutationMoverMode::N_MODES) == 2 && i > 1 ) {
			ss << " and ";
		} else if ( i > 1 ) {
			ss << ", ";
			if ( i == static_cast< core::Size >(CycpepRigidBodyPermutationMoverMode::N_MODES) ) {
				ss << " and ";
			}
		}
		ss << string_from_cycpep_rigid_body_permutation_mover_mode( static_cast< CycpepRigidBodyPermutationMoverMode >(i) );
	}
	return ss.str();
}

////////////////////////////////////////////////////////////////////////////////
// CONSTRUCTORS AND DESTRUCTORS
////////////////////////////////////////////////////////////////////////////////

/// @brief Default constructor
CycpepRigidBodyPermutationMover::CycpepRigidBodyPermutationMover():
	protocols::moves::Mover( CycpepRigidBodyPermutationMover::mover_name() )
{}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
CycpepRigidBodyPermutationMover::~CycpepRigidBodyPermutationMover() = default;

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Apply the mover
void
CycpepRigidBodyPermutationMover::apply(
	core::pose::Pose & pose
) {
	// Generate vector of whether residues are selected.
	utility::vector1< bool > selected_residues_boolvec(
		residue_selector_ == nullptr ?
		utility::vector1< bool >( pose.total_residue(), true ) :
		residue_selector_->apply(pose)
	);

	// Generate vector of selected residue indices.
	utility::vector1< core::Size > const selected_residues(
		core::select::residue_selector::selection_positions( selected_residues_boolvec )
	);

	// Check that the pose or selection is a cyclic peptide.
	confirm_is_cyclic_peptide( pose, selected_residues );

	// Copy cyclic peptide part of pose.
	core::pose::PoseOP cycpep_pose( generate_cycpep_pose_copy( pose, selected_residues ) );

	// Figure out the mapping of cyclic peptide pose residues to original pose residues, based on mode and either
	// the set offset and inversion or a random offset and inversion.
	utility::vector1< core::Size > const residue_index_map( generate_residue_index_map( selected_residues ) ); //Contains indices of *original pose*, where key is index in *cycpep_pose*.

	// Align the cyclic peptide pose to the original.
	align_cycpep_pose_to_original_pose( *cycpep_pose, pose, residue_index_map );

	// Perturb the cyclic peptide pose position.
	perturb_cycpep_pose_position( *cycpep_pose );

	// Perturb the cyclic peptide pose orientation.
	perturb_cycpep_pose_orientation( *cycpep_pose );

	// Copy the cyclic peptide pose coordinates to the original.
	copy_cycpep_pose_coordinates_to_original_pose( *cycpep_pose, pose, selected_residues );

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
CycpepRigidBodyPermutationMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief Parse XML tag (to use this Mover in RosettaScripts).
void
CycpepRigidBodyPermutationMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap
) {
	std::string const errmsg( "Error in CycpepRigidBodyPermutationMover::parse_my_tag(): " );

	if ( tag->hasOption("mode") ) {
		set_mover_mode( tag->getOption<std::string>("mode") );
	}

	if ( mode_ == CycpepRigidBodyPermutationMoverMode::SET_PERMUTATION ) {
		runtime_assert_string_msg( !tag->hasOption("allow_random_inversion"),
			errmsg + "The \"allow_random_inversion\" option is not allowed in \"set_permutation\" mode."
		);
	} else if ( mode_ == CycpepRigidBodyPermutationMoverMode::RANDOM_PERMUTATION ) {
		runtime_assert_string_msg( !tag->hasOption("set_inverse_alignment"),
			errmsg + "The \"set_inverse_alignment\" option is not allowed in \"randomized_permutation\" mode."
		);
		runtime_assert_string_msg( !tag->hasOption("set_permutation_offset"),
			errmsg + "The \"set_permutation_offset\" option is not allowed in \"randomized_permutation\" mode."
		);
	}

	if ( tag->hasOption("set_permutation_offset") ) {
		set_permutation_setting( tag->getOption<signed long int>( "set_permutation_offset" ) );
	}
	if ( tag->hasOption("set_inverse_alignment") ) {
		set_inversion_setting( tag->getOption<bool>( "set_inverse_alignment" ) );
	}
	if ( tag->hasOption("allow_random_inversion") ) {
		set_inversion_random( tag->getOption<bool>( "allow_random_inversion" ) );
	}

	if ( tag->hasOption("random_position_offset") ) {
		set_random_position_offset( tag->getOption<core::Real>( "random_position_offset" ) );
	}
	if ( tag->hasOption("random_orientation_perturbation") ) {
		set_random_orientation_perturbation( tag->getOption<core::Real>( "random_orientation_perturbation" ) );
	}

	if ( tag->hasOption("residue_selector") ) {
		set_residue_selector( protocols::rosetta_scripts::parse_residue_selector( tag, datamap, "residue_selector" ) );
	}
}

/// @brief Describe the XML interface for this mover, in a machine-readable way that permits
/// automatic checking of user inputs.
void
CycpepRigidBodyPermutationMover::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) {

	using namespace utility::tag;
	AttributeList attlist;

	attlist
		+ XMLSchemaAttribute( "mode", xs_string,
		"The mode for this mover, used to set whether we're setting a particular "
		"permutation or whether we're drawing randomly from possible permutations.  "
		"Allowed settings are: " + list_cycpep_rigid_body_permutation_mover_modes() + "."
		)
		+ XMLSchemaAttribute( "allow_random_inversion", xsct_rosetta_bool,
		"If true (the default), then in randomized_permutation mode, allowed permutations "
		"include those that align the peptide to a reversed permutation of itself.  If false, "
		"then only forward permutations are allowed.  Only used in randomized_permutation mode.  "
		"See set_inverse_alignment for the related setting for set_permutation mode."
		)
		+ XMLSchemaAttribute( "set_permutation_offset", xs_integer,
		"The cyclic permutation offset.  If aligning to the forward sequence, this is the offset "
		"in the forward direction along the sequence.  If aligning to the reversed sequence, this is "
		"the offset in the backward direction.  Only used in set_permutation mode."
		)
		+ XMLSchemaAttribute( "set_inverse_alignment", xsct_rosetta_bool,
		"Sets whether we align the peptide to its forward sequence or to its reversed sequence.  Only "
		"used in set_permutation mode.  See allow_random_inversion for the related setting for "
		"randomized_permutation mode."
		)
		+ XMLSchemaAttribute( "random_position_offset", xsct_real,
		"The magnitude of a random offset, in Angstroms, that will be applied to the peptide position after every cyclic permutation. "
		"The direction is wholly random.  Set this to 0.0 to disable.  Used in all modes."
		)
		+ XMLSchemaAttribute( "random_orientation_perturbation", xsct_real,
		"The magnitude of a random perturbation, in degrees, that will be applied to the peptide orientation after every cyclic permutation. "
		"The rotation axis is wholly random, while the center of rotation is the peptide centroid.  Set this to 0.0 to disable.  Used in all modes."
	)
		;

	protocols::rosetta_scripts::attributes_for_parse_residue_selector(
		attlist,
		"An optional residue selector, used to select the cyclic peptide part of a pose."
	);

	protocols::moves::xsd_type_definition_w_attributes(
		xsd, mover_name(),
		"A mover that takes a cyclic peptide and alters its position (rigid-body transform), "
		"superimposing it on a cyclic permutation or reverse cyclic permutation of itself.  Created "
		"22 November 2021 by Vikram K. Mulligan and Stephan Kudlacek.",
		attlist
	);
}


////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
CycpepRigidBodyPermutationMover::fresh_instance() const
{
	return utility::pointer::make_shared< CycpepRigidBodyPermutationMover >();
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
CycpepRigidBodyPermutationMover::clone() const
{
	return utility::pointer::make_shared< CycpepRigidBodyPermutationMover >( *this );
}

std::string CycpepRigidBodyPermutationMover::get_name() const {
	return mover_name();
}

std::string CycpepRigidBodyPermutationMover::mover_name() {
	return "CycpepRigidBodyPermutationMover";
}



/////////////// Creator ///////////////

protocols::moves::MoverOP
CycpepRigidBodyPermutationMoverCreator::create_mover() const
{
	return utility::pointer::make_shared< CycpepRigidBodyPermutationMover >();
}

std::string
CycpepRigidBodyPermutationMoverCreator::keyname() const
{
	return CycpepRigidBodyPermutationMover::mover_name();
}

void CycpepRigidBodyPermutationMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	CycpepRigidBodyPermutationMover::provide_xml_schema( xsd );
}

/// @brief This mover is unpublished.  It returns Vikram K. Mulligan as its author.
void
CycpepRigidBodyPermutationMover::provide_citation_info(
	basic::citation_manager::CitationCollectionList & citations
) const {
	basic::citation_manager::UnpublishedModuleInfoOP unpub_module_info(
		utility::pointer::make_shared< basic::citation_manager::UnpublishedModuleInfo >(
		"CycpepRigidBodyPermutationMover", basic::citation_manager::CitedModuleType::Mover,
		"Vikram K. Mulligan",
		"Systems Biology Group, Center for Computational Biology, Flatiron Institute",
		"vmulligan@flatironinstitute.org",
		"Wrote the CycpepRigidBodyPermutationMover."
		)
	);
	unpub_module_info->add_author(
		"Stephan Kudlacek",
		"Menten AI",
		"stephan.kudlacek@menten.ai",
		"Initial testing of mover."
	);
	citations.add(
		unpub_module_info
	);
}

////////////////////////////////////////////////////////////////////////////////
// PUBLIC SETTERS
////////////////////////////////////////////////////////////////////////////////

/// @brief Set how we will permute the peptide, by string: setting its permutation or randomizing it.
/// @details Throws if string not recognized.
void
CycpepRigidBodyPermutationMover::set_mover_mode(
	std::string const & setting
) {
	CycpepRigidBodyPermutationMoverMode const setting_enum( cycpep_rigid_body_permutation_mover_mode_from_string( setting ) );
	runtime_assert_string_msg(
		setting_enum != CycpepRigidBodyPermutationMoverMode::UNKNOWN_MODE,
		"Error in CycpepRigidBodyPermutationMover::set_mover_mode(): Could not parse \"" + setting + "\".  Valid "
		"modes are: " + list_cycpep_rigid_body_permutation_mover_modes() + "."
	);
	set_mover_mode(setting_enum);
}

/// @brief Set how we will permute the peptide: setting its permutation or randomizing it.
void
CycpepRigidBodyPermutationMover::set_mover_mode(
	CycpepRigidBodyPermutationMoverMode const setting
) {
	runtime_assert( static_cast< core::Size >(setting) > 0 && setting <= CycpepRigidBodyPermutationMoverMode::N_MODES );
	mode_ = setting;
	TR << "Set CycpepRigidBodyPermutationMover mode to " << string_from_cycpep_rigid_body_permutation_mover_mode(mode_) << "." << std::endl;
}

/// @brief Set how much we are permuting the peptide in setting mode.  If inversion_setting is
/// false, this is residues in the forward direction; if it is true, this is residues in the reverse
/// direction.
void
CycpepRigidBodyPermutationMover::set_permutation_setting(
	signed long int const offset
) {
	permutation_setting_ = offset;
	TR << "Set CycpepRigidBodyPermutationMover offset setting to " << permutation_setting_ << ".  (Note that this assumes set_permutation mode.)" << std::endl;
}

/// @brief Set whether we are overlaying the peptide on a permuted version of itself (false) or on
/// a permuted version of the reverse sequence (true).  Only used in setting mode.
void
CycpepRigidBodyPermutationMover::set_inversion_setting(
	bool const inversion_setting
) {
	inversion_setting_ = inversion_setting;
	TR << "Set CycpepRigidBodyPermutationMover inversion setting to " << (inversion_setting_ ? "TRUE" : "FALSE")
		<< ".  (Note that this assumes set_permutation mode.)" << std::endl;
}

/// @brief Set whether we are only overlaying the peptide on a permuted version of itself (false) or
/// whether we are allowed to consider overlays on a permuted version of the reverse sequence (true).
/// Only used in randomization mode.
void
CycpepRigidBodyPermutationMover::set_inversion_random(
	bool const inversion_random
) {
	inversion_random_ = inversion_random;
	TR << "Set CycpepRigidBodyPermutationMover inversion setting to " << (inversion_random_ ? "TRUE" : "FALSE")
		<< ".  (Note that this assumes randomize_permutation mode.)" << std::endl;
}

/// @brief Set whether we should add a small random offset to the peptide's position.
/// @details The default, 0.0, means no offset.  Larger values set the magnitude of the
/// random displacement, in Angstroms.
void
CycpepRigidBodyPermutationMover::set_random_position_offset(
	core::Real const offset
) {
	runtime_assert_string_msg(
		offset >= 0.0,
		"Error in CycpepRigidBodyPermutationMover::set_random_position_offset(): "
		"The random position offset magnitude must be greater than or equal to zero Angstroms."
	);
	random_position_offset_ = offset;
}

/// @brief Set whether we should add a small random perturbation to the peptide's orientation.
/// @details The default, 0.0, means no perturbation.  Larger values set the magnitude of the
/// random roll, in degrees.  The roll is about the peptide centroid.
void
CycpepRigidBodyPermutationMover::set_random_orientation_perturbation(
	core::Real const perturbation
) {
	runtime_assert_string_msg(
		perturbation >= 0.0,
		"Error in CycpepRigidBodyPermutationMover::set_random_orientation_perturbation(): "
		"The random orientation pertubation magnitude must be greater than or equal to zero degrees."
	);
	random_orientation_perturbation_ = perturbation;
}

/// @brief Set a residue selector to select a cyclic peptide that is a part of a pose.
/// @details Input selector is used directly, not cloned.  Pass nullptr to disable.
void
CycpepRigidBodyPermutationMover::set_residue_selector(
	core::select::residue_selector::ResidueSelectorCOP selector_in
) {
	residue_selector_ = selector_in;
	if ( residue_selector_ != nullptr ) {
		TR << "Set a " << residue_selector_->get_name() << " residue selector to select the cyclic peptide." << std::endl;
	} else {
		TR << "Clearing residue selector.  The whole pose will be assumed to be a cyclic peptide." << std::endl;
	}
}

////////////////////////////////////////////////////////////////////////////////
// PUBLIC GETTERS
////////////////////////////////////////////////////////////////////////////////

/// @brief Get how we will permute the peptide, by enum: setting its permutation or randomizing it.
CycpepRigidBodyPermutationMoverMode CycpepRigidBodyPermutationMover::mover_mode() const { return mode_; }

/// @brief Get how much we are permuting the peptide in setting mode.  If inversion_setting is
/// false, this is residues in the forward direction; if it is true, this is residues in the reverse
/// direction.
/// @note Deliberately signed.
signed long int CycpepRigidBodyPermutationMover::permutation_setting() const { return permutation_setting_; }

/// @brief Set whether we are overlaying the peptide on a permuted version of itself (false) or on
/// a permuted version of the reverse sequence (true).  Only used in setting mode.
bool CycpepRigidBodyPermutationMover::inversion_setting() const { return inversion_setting_; }

/// @brief Get whether we are only overlaying the peptide on a permuted version of itself (false) or
/// whether we are allowed to consider overlays on a permuted version of the reverse sequence (true).
/// Only used in randomization mode.
bool CycpepRigidBodyPermutationMover::inversion_random() const { return inversion_random_; }

/// @brief Get whether we should add a small random offset to the peptide's position.
/// @details The default, 0.0, means no offset.  Larger values set the magnitude of the
/// random displacement, in Angstroms.
core::Real CycpepRigidBodyPermutationMover::random_position_offset() const { return random_position_offset_; }

/// @brief Get whether we should add a small random perturbation to the peptide's orientation.
/// @details The default, 0.0, means no perturbation.  Larger values set the magnitude of the
/// random roll, in degrees.  The roll is about the peptide centroid.
core::Real CycpepRigidBodyPermutationMover::random_orientation_perturbation() const { return random_orientation_perturbation_; }

/// @brief Get the residue selector to select a cyclic peptide that is a part of a pose.
/// @details Could be nullptr if no selector is set.
core::select::residue_selector::ResidueSelectorCOP CycpepRigidBodyPermutationMover::residue_selector() const { return residue_selector_; }

////////////////////////////////////////////////////////////////////////////////
// PRIVATE METHODS
////////////////////////////////////////////////////////////////////////////////

/// @brief Check that the pose or selection is a cyclic peptide, and throw with an
/// informative error message if it is not.
/// @details The checks are:
/// - Does each residue have a lower connect and an upper connect?
/// - Is the upper of each connected to the lower of the next?
/// - Is there a connection from the last to the first?
/// @note This is general enough for peptides with other cyclization chemistries,
/// but NOT general enough for lariats or peptides with tails.  (Though this function
/// and generate_residue_index_map could be generalized to allow that in the future.)
void
CycpepRigidBodyPermutationMover::confirm_is_cyclic_peptide(
	core::pose::Pose const & pose,
	utility::vector1< core::Size > const & selected_residues
) const {
	std::string const errmsg( "Error in CycpepRigidBodyPermutationMover::confirm_is_cyclic_peptide(): ");
	runtime_assert_string_msg( selected_residues.size() >= 3, errmsg + "The pose or selected resiudes are " + std::to_string(selected_residues.size()) + ", but at least 3 residues are needed for a macrocycle." );

	for ( core::Size i(1), imax(selected_residues.size()); i<=imax; ++i ) {
		core::Size const curres( selected_residues[i] );
		core::chemical::ResidueType const & restype( pose.residue_type(curres) );
		core::conformation::Residue const & rsd( pose.residue(curres) );

		if ( i < imax ) {
			runtime_assert_string_msg( restype.upper_connect_id() != 0, errmsg + "Residue " + restype.base_name() + std::to_string(curres) + " has no upper connection." );
		}
		if ( i > 1 ) {
			runtime_assert_string_msg( restype.lower_connect_id() != 0, errmsg + "Residue " + restype.base_name() + std::to_string(curres) + " has no lower connection." );
			core::Size const prevres( selected_residues[i-1] );
			runtime_assert_string_msg( rsd.connected_residue_at_lower() == prevres, errmsg + "Expected residue " + restype.base_name() + std::to_string( curres ) + " to be connected to residue " + std::to_string( prevres ) + " at its lower connection, but instead it is connected to " + std::to_string(rsd.connected_residue_at_lower()) + "!" );
			runtime_assert_string_msg( pose.residue(prevres).connected_residue_at_upper() == curres, errmsg + "Expected residue " + pose.residue_type(prevres).base_name() + std::to_string( prevres ) + " to be connected to residue " + std::to_string( curres ) + " at its upper connection, but instead it is connected to " + std::to_string(pose.residue(prevres).connected_residue_at_upper()) + "!" );
		}
	}
	core::Size const firstres( selected_residues[1] ), lastres( selected_residues[selected_residues.size()] );
	runtime_assert_string_msg( pose.residue(lastres).is_bonded( firstres ), errmsg + "Expected residue " + pose.residue_type(lastres).base_name() + std::to_string(lastres) + " to be bonded to " + pose.residue_type(firstres).base_name() + std::to_string(firstres) + ", but they share no bond!" );
}

/// @brief Given a selection of residues, plus the current configuration of this mover,
/// generate a vector of residue indices that has been suitably circularly permuted or
/// inverse circular permuted.
utility::vector1< core::Size >
CycpepRigidBodyPermutationMover::generate_residue_index_map(
	utility::vector1< core::Size > const & selected_residues
) const {

	TR << "Generating residue index map for " << string_from_cycpep_rigid_body_permutation_mover_mode( mode_ ) << " mode." << std::endl;

	core::Size const nresidue( selected_residues.size() );
	runtime_assert( nresidue > 0 );

	// Determine the cyclic offset:
	signed long int const cyclic_offset(
		mode_ == CycpepRigidBodyPermutationMoverMode::SET_PERMUTATION ?
		permutation_setting_ :
		numeric::random::random_range( 0, nresidue )
	);
	TR << "\tCyclic offset = " << cyclic_offset << std::endl;

	bool const inverted(
		mode_ == CycpepRigidBodyPermutationMoverMode::SET_PERMUTATION ?
		inversion_setting_ :
		(numeric::random::uniform() > 0.5)
	);
	TR << "\tInversion = " << (inverted ? "TRUE" : "FALSE") << std::endl;

	utility::vector1< core::Size > returnvec( nresidue, 0 );

	core::Size startval( static_cast< core::Size >( numeric::modulo( cyclic_offset, static_cast<signed long int>(nresidue) ) ) + 1 );
	core::Size curval( startval ), counter(1);
	do {
		debug_assert( counter <= nresidue && counter > 0 );
		debug_assert( curval <= nresidue && curval > 0 );
		returnvec[counter] = selected_residues[curval];
		++counter;
		if ( inverted ) {
			--curval;
			if ( curval == 0 ) {
				curval = nresidue;
			}
		} else {
			++curval;
			if ( curval > nresidue ) {
				curval = 1;
			}
		}
	} while( curval != startval );

	TR << "\tPermuted indices: ";
	for ( core::Size i(1); i<=nresidue; ++i ) {
		TR << returnvec[i];
		if ( i < nresidue ) {
			TR << ", ";
		}
	}
	TR << std::endl;

	return returnvec;
}

/// @brief Given a pose, cut out selected residues.
/// @details Assumes that we've got a cyclic peptide.
core::pose::PoseOP
CycpepRigidBodyPermutationMover::generate_cycpep_pose_copy(
	core::pose::Pose const & pose,
	utility::vector1< core::Size > const & selected_residues
) const {
	runtime_assert( !selected_residues.empty() );
	runtime_assert( !pose.empty() );
	TR << "Building peptide pose." << std::endl;

	core::pose::PoseOP outpose( utility::pointer::make_shared< core::pose::Pose >() );
	for ( core::Size i(1), imax(selected_residues.size()); i<=imax; ++i ) {
		core::Size const curres( selected_residues[i] );
		runtime_assert( curres > 0 && curres <= pose.total_residue() );

		TR << "\tAdding " << pose.residue_type(curres).base_name() << curres << " as residue " << i << "." << std::endl;

		if ( i == 1 ) {
			outpose->append_residue_by_jump( pose.residue(curres), 0 );
		} else {
			core::Size const prevres( selected_residues[i - 1] );
			runtime_assert( pose.residue_type(curres).lower_connect_id() != 0 );
			runtime_assert( outpose->residue_type(i-1).upper_connect_id() != 0 );
			runtime_assert( pose.residue(curres).connected_residue_at_lower() == prevres );
			runtime_assert( pose.residue(prevres).connected_residue_at_upper() == curres );
			outpose->append_residue_by_bond( pose.residue(curres) );

			//Check for additional chemical bonds and add these:
			if ( i > 2 ) {
				for ( core::Size j(1), jmax(i-2); j<=jmax; ++j ) {
					core::Size const otherres( selected_residues[j] );
					for ( core::Size connid(1), connidmax( pose.residue_type(curres).n_possible_residue_connections() ); connid <= connidmax; ++connid ) {
						if ( pose.residue(curres).connected_residue_at_resconn(connid) == otherres ) {
							//Found a connection.  Add it.
							std::string const curatm( pose.residue_type(curres).atom_name(
								pose.residue_type(curres).residue_connection(connid).atomno() )
							);

							core::Size const other_connid( pose.residue(curres).residue_connection_conn_id( connid ) );

							std::string const otheratm( pose.residue_type(otherres).atom_name(
								pose.residue_type(otherres).residue_connection(other_connid).atomno() )
							);

							TR << "\t\tAdding bond between " << pose.residue_type(curres).base_name() << curres
								<< " (residue " << i << "), atom " << curatm << " and " << pose.residue_type(otherres).base_name()
								<< otherres << " (residue " << j << ") atom " << otheratm << "." << std::endl;
							outpose->conformation().declare_chemical_bond( i, curatm, j, otheratm );
						}
					}
				}
			}
		}
	}

	outpose->update_residue_neighbors();
	return outpose;
}

/// @brief Align the cyclic peptide pose to the original.
/// @details The offset in the residue index map is used.  N, CA, C, O, CB, CM (beta-aas), and CA1 (peptoids) atoms
/// are used, provided that pairs of residues have both.
/// @param[inout] cycpep_pose The pose whose coordinates are updated by this operation.
/// @param[in] pose The pose that we're aligning TO.  This is unaltered by this operation.
/// @param[in] residue_index_map A vector of residue indices in pose, where the index in
/// the vector is the residue index in cycpep_pose.
void
CycpepRigidBodyPermutationMover::align_cycpep_pose_to_original_pose(
	core::pose::Pose & cycpep_pose,
	core::pose::Pose const & pose,
	utility::vector1 < core::Size > const & residue_index_map
) const {
	std::map< core::id::AtomID, core::id::AtomID > atom_map;
	runtime_assert_string_msg( residue_index_map.size() >= 3, "Error in CycpepRigidBodyPermutationMover::align_cycpep_pose_to_original_pose(): At least 3 residues must be selected that define the cyclic peptide." );

	utility::fixedsizearray1< std::string, 7 > const atnames{ "N", "CA", "C", "O", "CB", "CM", "CA1" };

	for ( core::Size cycpep_res(1), cycpep_res_max(residue_index_map.size()); cycpep_res<=cycpep_res_max; ++cycpep_res ) {
		core::Size const orig_res( residue_index_map[cycpep_res] );
		for ( std::string const & atname : atnames ) {
			if ( atname == "CM" ) {
				// Only align CM is this is a beta-amino acid.
				if ( (!pose.residue_type(orig_res).is_beta_aa()) || (!cycpep_pose.residue_type(cycpep_res).is_beta_aa()) ) {
					continue;
				}
			} else if ( atname == "CA1" ) {
				// Only align CA1 if this is a peptoid.
				if ( (!pose.residue_type(orig_res).is_peptoid()) || (!cycpep_pose.residue_type(cycpep_res).is_peptoid()) ) {
					continue;
				}
			}
			if ( cycpep_pose.residue_type(cycpep_res).has(atname) && pose.residue_type(orig_res).has(atname) ) {
				atom_map[ core::id::AtomID( cycpep_pose.residue_type(cycpep_res).atom_index(atname), cycpep_res ) ] =
					core::id::AtomID( pose.residue_type(orig_res).atom_index(atname), orig_res );
			}
		}
	}

	if ( TR.visible() ) {
		TR << "Aligning the following:" << std::endl;
		TR << "ORIGINAL_RES\tORIGINAL_ATOM\tCYCPEP_RES\tCYCPEP_ATOM" << std::endl;
		TR << "------------\t-------------\t----------\t-----------" << std::endl;
		for ( std::map< core::id::AtomID, core::id::AtomID >::const_iterator it( atom_map.begin() ); it != atom_map.end(); ++it ) {
			TR << pose.residue_type(it->second.rsd()).base_name() << it->second.rsd() << "\t"
				<< pose.residue_type(it->second.rsd()).atom_name( it->second.atomno() ) << "(" << it->second.atomno() << ")\t"
				<< cycpep_pose.residue_type(it->first.rsd()).base_name() << it->first.rsd() << "\t"
				<< cycpep_pose.residue_type(it->first.rsd()).atom_name( it->first.atomno() ) << "(" << it->first.atomno() << ")" << std::endl;
		}
	}

	core::Real const rmsd( core::scoring::superimpose_pose( cycpep_pose, pose, atom_map ) );
	TR << "Superimposed cyclic peptide with offset.  RMSD = " << rmsd << " A." << std::endl;
}

/// @brief Perturb the cyclic peptide pose position.
/// @details Does nothing if random_position_offset_ == 0.0.  Otherwise,
/// generates a random vector from a Gaussian distribution with mean length
/// random_position_offset_ and adds this to the coordinates of the pose.
/// @param[inout] pose The pose whose position we're offsetting.  Altered by
/// this operation.
void
CycpepRigidBodyPermutationMover::perturb_cycpep_pose_position(
	core::pose::Pose & pose
) const {
	if ( random_position_offset_ == 0.0 ) return;

	utility::vector1< core::id::AtomID > all_ids;
	all_ids.reserve( pose.total_atoms() );
	for ( core::Size ir(1), irmax(pose.total_residue()); ir<=irmax; ++ir ) {
		for ( core::Size ia(1), iamax(pose.residue_type(ir).natoms()); ia<=iamax; ++ia ) {
			all_ids.emplace_back( core::id::AtomID( ia, ir ) );
		}
	}
	utility::vector1< numeric::xyzVector< core::Real > > xyzvals;
	pose.batch_get_xyz( all_ids, xyzvals );

	numeric::xyzVector< core::Real > const offset_vec( numeric::random::random_vector_spherical() * random_position_offset_ );
	for ( numeric::xyzVector< core::Real > & xyzval : xyzvals ) {
		xyzval += offset_vec;
	}

	pose.batch_set_xyz( all_ids, xyzvals );
}

/// @brief Perturb the cyclic peptide pose orientation.
/// @details Does nothing if random_orientation_perturbation_ == 0.0.  Otherwise,
/// generates a random unit vector from a spherically-symmetric distribution, and rotates
/// by X degrees about this vector, where X is drawn from a Gaussian distribution of
/// standard deviation random_orientation_perturbation_ degrees.  The centre of rotation
/// is the centroid of the N, C, CA, O, CB (if present), CM (if beta-aa), and CA1 (if peptoid)
/// atoms.
/// @param[inout] pose The pose whose orientation we're perturbing.  Altered by
/// this operation.
void
CycpepRigidBodyPermutationMover::perturb_cycpep_pose_orientation(
	core::pose::Pose & pose
) const {

	if ( random_orientation_perturbation_ == 0 ) return; //Do nothing in this case.

	// Find centroid:
	numeric::xyzVector< core::Real > centroid(0.0, 0.0, 0.0);
	core::Size atom_counter(0);
	utility::fixedsizearray1< std::string, 7 > const atnames{ "N", "CA", "C", "O", "CB", "CM", "CA1" };

	for ( core::Size ir(1), irmax(pose.total_residue()); ir<=irmax; ++ir ) {
		for ( std::string const & atname : atnames ) {
			if ( atname == "CM" ) {
				// Only use CM is this is a beta-amino acid.
				if ( !pose.residue_type(ir).is_beta_aa() ) {
					continue;
				}
			} else if ( atname == "CA1" ) {
				// Only use CA1 if this is a peptoid.
				if ( !pose.residue_type(ir).is_peptoid() ) {
					continue;
				}
			}
			if ( pose.residue_type(ir).has( atname ) ) {
				centroid += pose.xyz( core::id::AtomID( pose.residue_type(ir).atom_index(atname), ir ) );
				++atom_counter;
			}
		}
	}
	runtime_assert( atom_counter > 0 );
	centroid /= static_cast< core::Real >( atom_counter );

	//Generate a random rotation axis, drawn from a spherically-symmetric distribution.
	numeric::xyzVector< core::Real > const axis( numeric::random::uniform_vector_sphere() );

	//Generate the random rotation, and the rotation matrix:
	core::Real const random_rotation_degrees( numeric::random::gaussian() * random_orientation_perturbation_ );
	numeric::xyzMatrix< core::Real > rotation_matrix( numeric::rotation_matrix_degrees( axis, random_rotation_degrees ) );

	//Get the XYZ coordinates from the pose.
	utility::vector1< numeric::xyzVector< core::Real > > coords;
	utility::vector1< core::id::AtomID > allids;
	allids.reserve( pose.total_atoms() );
	for ( core::Size ir(1), irmax(pose.total_residue()); ir<=irmax; ++ir ) {
		for ( core::Size ia(1), iamax(pose.residue_type(ir).natoms()); ia<=iamax; ++ia ) {
			allids.emplace_back( core::id::AtomID( ia, ir ) );
		}
	}
	pose.batch_get_xyz( allids, coords );

	//Transform the coordinates:
	for ( numeric::xyzVector< core::Real > & coord : coords ) {
		coord = (rotation_matrix * (coord - centroid)) + centroid;
	}

	//Set XYZ coordinates for the pose:
	pose.batch_set_xyz( allids, coords );
}

/// @brief Copy the cyclic peptide pose coordinates to the original.
/// @param[in] cycpep_pose The cyclic peptide pose, the pose that we are copying FROM.
/// @param[inout] pose The original pose, the pose whose coordinates we are UPDATING.
/// @param[in] selected_residues A vector of indices of the residues of the cyclic peptide.
void
CycpepRigidBodyPermutationMover::copy_cycpep_pose_coordinates_to_original_pose(
	core::pose::Pose const & cycpep_pose,
	core::pose::Pose & pose,
	utility::vector1< core::Size > const & selected_residues
) const {
	std::string const errmsg( "Error in protocols::cyclic_peptide::CycpepRigidBodyPermutationMover::copy_cycpep_pose_coordinates_to_original_pose():  " );
	core::Size const nresidue( selected_residues.size() );

	runtime_assert( cycpep_pose.total_residue() == nresidue );
	runtime_assert( pose.total_residue() >= nresidue );

	utility::vector1< core::id::AtomID > original_pose_atomids, cycpep_pose_atomids;
	for ( core::Size ir(1); ir<=nresidue; ++ir ) {
		core::Size const original_ir( selected_residues[ir] );
		core::Size const natom( cycpep_pose.residue_type(ir).natoms() );

		runtime_assert_string_msg(
			pose.residue_type( original_ir ).natoms() == natom,
			errmsg + "Original pose residue " + pose.residue_type( original_ir ).base_name() + std::to_string(original_ir) + " has "
			+ std::to_string( pose.residue_type( original_ir ).natoms() ) + " atoms, but " + std::to_string( natom )
			+ " were expected for cyclic peptide pose residue " + cycpep_pose.residue_type(ir).base_name()
			+ std::to_string(ir) + " !"
		);

		for ( core::Size ia(1); ia<=natom; ++ia ) {
			original_pose_atomids.emplace_back( core::id::AtomID( ia, original_ir ) );
			cycpep_pose_atomids.emplace_back( core::id::AtomID( ia, ir ) );
		}
	}

	utility::vector1< core::PointPosition > xyz_coords;
	cycpep_pose.batch_get_xyz( cycpep_pose_atomids, xyz_coords );
	pose.batch_set_xyz( original_pose_atomids, xyz_coords );
	pose.update_residue_neighbors();

	TR << "Copied coordinates of perturbed peptide to original." << std::endl;

}

std::ostream &
operator<<( std::ostream & os, CycpepRigidBodyPermutationMover const & mover )
{
	mover.show(os);
	return os;
}


} //cyclic_peptide
} //protocols
