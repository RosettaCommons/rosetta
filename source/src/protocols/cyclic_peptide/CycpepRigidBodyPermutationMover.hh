// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/CycpepRigidBodyPermutationMover.hh
/// @brief Headers for a mover that takes a cyclic peptide and alters its position (rigid-body transform), superimposing
/// it on a cyclic permutation or reverse cyclic permutation of itself.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_protocols_cyclic_peptide_CycpepRigidBodyPermutationMover_HH
#define INCLUDED_protocols_cyclic_peptide_CycpepRigidBodyPermutationMover_HH

// Unit headers
#include <protocols/cyclic_peptide/CycpepRigidBodyPermutationMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.fwd.hh>

// STL headers:
#include <string>

namespace protocols {
namespace cyclic_peptide {

////////////////////////////////////////////////////////////////////////////////
// ENUM AND UTILITY FUNCTIONS
////////////////////////////////////////////////////////////////////////////////

/// @brief An enum for the modes for this mover.
/// @details If you add to this list, update the
/// string_from_cycpep_rigid_body_permutation_mover_mode()
/// function.
enum class CycpepRigidBodyPermutationMoverMode {
	SET_PERMUTATION = 1,
	RANDOM_PERMUTATION, //Keep third-to-last.
	N_MODES=RANDOM_PERMUTATION, //Keep second-to-last.
	UNKNOWN_MODE //Keep last.
};


/// @brief Given a mover mode string, get the mode enum.
/// @details Returns CycpepRigidBodyPermutationMoverMode::UNKNOWN_MODE if string not recognized.
CycpepRigidBodyPermutationMoverMode
cycpep_rigid_body_permutation_mover_mode_from_string(
	std::string const & modestring
);

/// @brief Given a mover mode enum, get the mode string.
std::string
string_from_cycpep_rigid_body_permutation_mover_mode(
	CycpepRigidBodyPermutationMoverMode const mode
);

/// @brief Get a comma-separated list of valid mover modes.
std::string
list_cycpep_rigid_body_permutation_mover_modes();


/// @brief A mover that takes a cyclic peptide and alters its position (rigid-body transform), superimposing it on
/// a cyclic permutation or reverse cyclic permutation of itself.
class CycpepRigidBodyPermutationMover : public protocols::moves::Mover {

public:

	////////////////////////////////////////////////////////////////////////////////
	// CONSTRUCTORS AND DESTRUCTORS
	////////////////////////////////////////////////////////////////////////////////

	/// @brief Default constructor
	CycpepRigidBodyPermutationMover();

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~CycpepRigidBodyPermutationMover() override;


public:

	////////////////////////////////////////////////////////////////////////////////
	// MOVER METHODS
	////////////////////////////////////////////////////////////////////////////////


	/// @brief Apply the mover
	void
	apply( core::pose::Pose & pose ) override;

	void
	show( std::ostream & output = std::cout ) const override;


public:

	////////////////////////////////////////////////////////////////////////////////
	// ROSETTASCRIPTS SUPPORT
	////////////////////////////////////////////////////////////////////////////////

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data ) override;

	//CycpepRigidBodyPermutationMover & operator=( CycpepRigidBodyPermutationMover const & src );

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	clone() const override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	/// @brief Describe the XML interface for this mover, in a machine-readable way that permits
	/// automatic checking of user inputs.
	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

public: //Function overrides needed for the citation manager:

	/// @brief This mover is unpublished.  It returns Vikram K. Mulligan as its author.
	void provide_citation_info(basic::citation_manager::CitationCollectionList & citations) const override;

public:

	////////////////////////////////////////////////////////////////////////////////
	// PUBLIC SETTERS
	////////////////////////////////////////////////////////////////////////////////

	/// @brief Set how we will permute the peptide, by string: setting its permutation or randomizing it.
	/// @details Throws if string not recognized.
	void set_mover_mode( std::string const & setting );

	/// @brief Set how we will permute the peptide, by enum: setting its permutation or randomizing it.
	void set_mover_mode( CycpepRigidBodyPermutationMoverMode const setting );

	/// @brief Set how much we are permuting the peptide in setting mode.  If inversion_setting is
	/// false, this is residues in the forward direction; if it is true, this is residues in the reverse
	/// direction.
	/// @note Deliberately signed.
	void set_permutation_setting( signed long int const offset );

	/// @brief Set whether we are overlaying the peptide on a permuted version of itself (false) or on
	/// a permuted version of the reverse sequence (true).  Only used in setting mode.
	void set_inversion_setting( bool const inversion_setting );

	/// @brief Set whether we are only overlaying the peptide on a permuted version of itself (false) or
	/// whether we are allowed to consider overlays on a permuted version of the reverse sequence (true).
	/// Only used in randomization mode.
	void set_inversion_random( bool const inversion_random );

	/// @brief Set whether we should add a small random offset to the peptide's position.
	/// @details The default, 0.0, means no offset.  Larger values set the magnitude of the
	/// random displacement, in Angstroms.
	void set_random_position_offset( core::Real const offset );

	/// @brief Set whether we should add a small random perturbation to the peptide's orientation.
	/// @details The default, 0.0, means no perturbation.  Larger values set the magnitude of the
	/// random roll, in degrees.  The roll is about the peptide centroid.
	void set_random_orientation_perturbation( core::Real const perturbation );

	/// @brief Set a residue selector to select a cyclic peptide that is a part of a pose.
	/// @details Input selector is used directly, not cloned.  Pass nullptr to disable.
	void set_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector_in );

public:

	////////////////////////////////////////////////////////////////////////////////
	// PUBLIC GETTERS
	////////////////////////////////////////////////////////////////////////////////

	/// @brief Get how we will permute the peptide, by enum: setting its permutation or randomizing it.
	CycpepRigidBodyPermutationMoverMode mover_mode() const;

	/// @brief Get how much we are permuting the peptide in setting mode.  If inversion_setting is
	/// false, this is residues in the forward direction; if it is true, this is residues in the reverse
	/// direction.
	/// @note Deliberately signed.
	signed long int permutation_setting() const;

	/// @brief Set whether we are overlaying the peptide on a permuted version of itself (false) or on
	/// a permuted version of the reverse sequence (true).  Only used in setting mode.
	bool inversion_setting() const;

	/// @brief Get whether we are only overlaying the peptide on a permuted version of itself (false) or
	/// whether we are allowed to consider overlays on a permuted version of the reverse sequence (true).
	/// Only used in randomization mode.
	bool inversion_random() const;

	/// @brief Get whether we should add a small random offset to the peptide's position.
	/// @details The default, 0.0, means no offset.  Larger values set the magnitude of the
	/// random displacement, in Angstroms.
	core::Real random_position_offset() const;

	/// @brief Get whether we should add a small random perturbation to the peptide's orientation.
	/// @details The default, 0.0, means no perturbation.  Larger values set the magnitude of the
	/// random roll, in degrees.  The roll is about the peptide centroid.
	core::Real random_orientation_perturbation() const;

	/// @brief Get the residue selector to select a cyclic peptide that is a part of a pose.
	/// @details Could be nullptr if no selector is set.
	core::select::residue_selector::ResidueSelectorCOP residue_selector() const;

private:

	////////////////////////////////////////////////////////////////////////////////
	// PRIVATE FUNCTIONS
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
	confirm_is_cyclic_peptide(
		core::pose::Pose const & pose,
		utility::vector1< core::Size > const & selected_residues
	) const;


	/// @brief Given a selection of residues, plus the current configuration of this mover,
	/// generate a vector of residue indices that has been suitably circularly permuted or
	/// inverse circular permuted.
	utility::vector1< core::Size >
	generate_residue_index_map(
		utility::vector1< core::Size > const & selected_residues
	) const;

	/// @brief Given a pose, cut out selected residues.
	/// @details Assumes that we've got a cyclic peptide.
	core::pose::PoseOP
	generate_cycpep_pose_copy(
		core::pose::Pose const & pose,
		utility::vector1< core::Size > const & selected_residues
	) const;

	/// @brief Align the cyclic peptide pose to the original.
	/// @details The offset in the residue index map is used.  N, CA, C, O, CB, CM (beta-aas), and CA1 (peptoids) atoms
	/// are used, provided that pairs of residues have both.
	/// @param[inout] cycpep_pose The pose whose coordinates are updated by this operation.
	/// @param[in] pose The pose that we're aligning TO.  This is unaltered by this operation.
	/// @param[in] residue_index_map A vector of residue indices in pose, where the index in
	/// the vector is the residue index in cycpep_pose.
	void
	align_cycpep_pose_to_original_pose(
		core::pose::Pose & cycpep_pose,
		core::pose::Pose const & pose,
		utility::vector1 < core::Size > const & residue_index_map
	) const;

	/// @brief Perturb the cyclic peptide pose position.
	/// @details Does nothing if random_position_offset_ == 0.0.  Otherwise,
	/// generates a random vector from a Gaussian distribution with mean length
	/// random_position_offset_ and adds this to the coordinates of the pose.
	/// @param[inout] pose The pose whose position we're offsetting.  Altered by
	/// this operation.
	void
	perturb_cycpep_pose_position(
		core::pose::Pose & pose
	) const;

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
	perturb_cycpep_pose_orientation(
		core::pose::Pose & pose
	) const;

	/// @brief Copy the cyclic peptide pose coordinates to the original.
	/// @param[in] cycpep_pose The cyclic peptide pose, the pose that we are copying FROM.
	/// @param[inout] pose The original pose, the pose whose coordinates we are UPDATING.
	/// @param[in] selected_residues A vector of indices of the residues of the cyclic peptide.
	void
	copy_cycpep_pose_coordinates_to_original_pose(
		core::pose::Pose const & cycpep_pose,
		core::pose::Pose & pose,
		utility::vector1< core::Size > const & selected_residues
	) const;

private:

	////////////////////////////////////////////////////////////////////////////////
	// PRIVATE MEMBER DATA
	////////////////////////////////////////////////////////////////////////////////

	/// @brief How we will permute the peptide: setting its permutation or randomizing it.
	CycpepRigidBodyPermutationMoverMode mode_ = CycpepRigidBodyPermutationMoverMode::RANDOM_PERMUTATION;

	/// @brief If we are setting the permutation, what is the setting?  Zero by default (no permutation).
	/// @note Deliberately signed.
	signed long int permutation_setting_ = 0;

	/// @brief If we are setting the permutation, are we setting this to overlay on the forward or backwards
	/// sequence?  False by default (overlaying on forward sequence).
	bool inversion_setting_ = false;

	/// @brief If we are randomizing the permutation, do we also allow overlay on the backwards sequence?
	/// True by default.
	bool inversion_random_ = true;

	/// @brief Should we add a small random offset to the peptide's position?
	/// @details The default, 0.0, means no offset.  Larger values set the magnitude of the
	/// random displacement, in Angstroms.
	core::Real random_position_offset_ = 0.0;

	/// @brief Should we add a small random perturbation to the peptide's orientation?
	/// @details The default, 0.0, means no perturbation.  Larger values set the magnitude of the
	/// random roll, in degrees.  The roll is about the peptide centroid.
	core::Real random_orientation_perturbation_ = 0.0;

	/// @brief An optional residue selector, used to select a cyclic peptide in a pose.
	/// @details Not used if nullptr.
	core::select::residue_selector::ResidueSelectorCOP residue_selector_;

};

std::ostream &
operator<<( std::ostream & os, CycpepRigidBodyPermutationMover const & mover );

} //cyclic_peptide
} //protocols

#endif //protocols_cyclic_peptide_CycpepRigidBodyPermutationMover_HH
