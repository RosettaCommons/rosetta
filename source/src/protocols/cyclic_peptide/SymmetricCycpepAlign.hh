// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/SymmetricCycpepAlign.hh
/// @brief Given a quasi-symmetric cyclic peptide, this mover aligns the peptide so that the cyclic symmetry axis lies along the Z-axis and the centre of mass is at the origin.
/// It then optionally removes all but one symmetry repeat, so that true symmetry may be set up with the SetupForSymmetry mover.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

#ifndef INCLUDED_protocols_cyclic_peptide_SymmetricCycpepAlign_HH
#define INCLUDED_protocols_cyclic_peptide_SymmetricCycpepAlign_HH

// Unit headers
#include <protocols/cyclic_peptide/SymmetricCycpepAlign.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
//#include <utility/tag/XMLSchemaGeneration.fwd.hh> //transcluded from Mover

namespace protocols {
namespace cyclic_peptide {

///@brief Given a quasi-symmetric cyclic peptide, this mover aligns the peptide so that the cyclic symmetry axis lies along the Z-axis and the centre of mass is at the origin.  It then optionally removes all but one symmetry repeat, so that true symmetry may be set up with the SetupForSymmetry mover.
class SymmetricCycpepAlign : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	SymmetricCycpepAlign();

	/// @brief Copy constructor (not needed unless you need deep copies)
	SymmetricCycpepAlign( SymmetricCycpepAlign const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~SymmetricCycpepAlign() override;

	/////////////////////
	/// Mover Methods ///
	/////////////////////

public:
	/// @brief Apply the mover
	void
	apply( core::pose::Pose & pose ) override;

	void
	show( std::ostream & output = std::cout ) const override;

	///////////////////////////////
	/// Rosetta Scripts Support ///
	///////////////////////////////

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;

	//SymmetricCycpepAlign & operator=( SymmetricCycpepAlign const & src );

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

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

public: // setters and getters

	/// @brief Get the number of symmetry repeats.  For example, for c4 or c4/m symmetry,
	/// this would be "4".  Defaults to "2".
	inline core::Size symmetry_repeats() const { return symmetry_repeats_; }

	/// @brief Get whether we're using mirror symmetry (e.g. c2/m symmetry vs. c2 symmetry).
	///
	inline bool mirror_symmetry() const { return mirror_symmetry_; }

	/// @brief Get whether we're auto-detecting symmetry.
	///
	inline bool auto_detect_symmetry() const { return auto_detect_symmetry_; }

	/// @brief Get the angle threshold for counting poses as symmetric.
	/// @details Two dihedral values from different residues must fall within this cutoff, in degrees, for them to
	/// be considered the "same".  Default 10 degrees.
	inline core::Real const & angle_threshold( ) const { return angle_threshold_; }

	/// @brief Get whether all geometry that is not the protein backbone of a single symmetry repeat will be
	/// deleted.  (This includes any crosslinkers.)
	inline bool trim_to_single_repeat() const { return trim_to_single_repeat_; }

	/// @brief If trim_to_single_repeat_ is true, this is the symmetry repeat to preserve.
	///
	inline core::Size repeat_to_preserve() const { return repeat_to_preserve_; }

	/// @brief Report the symmetry of the last auto-detected peptide.
	/// @details Returns 0 if auto-detection is not enabled.
	inline core::Size last_symmetry_repeats() const { return auto_detect_symmetry_ ? last_symmetry_repeats_ : 0; }

	/// @brief Report whether the symmetry of the last auto-detected peptide was mirror symmetry.
	/// @details Returns false always if auto-detection is not enabled.
	inline bool last_symmetry_mirror() const { return auto_detect_symmetry_ ? last_symmetry_mirror_ : false; }

	/// @brief Are we aligning the peptide normal with the Z-axis (false) or the inverse Z-axis (true)?
	///
	inline bool invert() const {return invert_;}

	/// @brief Set the number of symmetry repeats and whether we're using mirror symmetry.  For example, for c4
	/// symmetry, inputs are "4", "false".  For c4/m, they'd be "4", "true".
	void set_symmetry( core::Size const repeats_in, bool const mirror_in );

	/// @brief Set whether we're auto-detecting symmetry.
	///
	void set_auto_detect_symmetry( bool const setting );

	/// @brief Set the angle threshold for counting poses as symmetric.
	/// @details Two dihedral values from different residues must fall within this cutoff, in degrees, for them to
	/// be considered the "same".  Default 10 degrees.
	void set_angle_threshold( core::Real const &setting );

	/// @brief Set whether we're going to delete all geometry that isn't the protein backbone of a single
	/// symmetry repeat, and which repeat to preserve.
	/// @param[in] do_trim If set to true, the peptide will be trimmed down to a single repeat.  False by default (no trimming).
	/// @param[in] repeat_to_preserve If do_trim is set to true, this is the repeat that should be preserved.  1 by default.
	void set_trim_info( bool const do_trim, core::Size const repeat_to_preserve=1 );

	/// @brief Set whether we're aligning the peptide normal with the Z-axis (false) or the inverse Z-axis (true).
	///
	inline void set_invert( bool const setting ) { invert_ = setting; }

private: // methods

	/// @brief Set the symmetry of the last auto-detected peptide symmetry.
	///
	void set_last_symmetry( core::Size const last_symm_in, bool const last_symm_mirror_in );

	/// @brief Given a quasi-symmetric pose, figure out its symmetry.
	/// @details Starts with the maximum possible, and tries every possible symmetry type, favouring mirror symmetry over non-mirror symmetry.
	/// @param[in] pose The quasi-symmetric pose.
	/// @param[out] symmrepeats The number of symmetry repeats.
	/// @param[out] mirrorsymm Is this a mirror symmetry type (cN/m) or not (cN)?
	bool do_auto_detection_of_symmetry( core::pose::Pose const &pose, core::Size &symmrepeats, bool &mirrorsymm ) const;

	/// @brief Given a quasi-symmetric pose, confirm that it has the specified symmetry.
	///
	bool do_symmetry_checks( core::pose::Pose const &pose, core::Size const symmrepeats, bool const mirrorsymm ) const;

	/// @brief Given a pose, count the protein residues in it.
	///
	core::Size count_protein_residues( core::pose::Pose const &pose ) const;

	/// @brief Given a pose, return a ResidueSelector that selects all the protein residues in the pose.
	///
	core::select::residue_selector::ResidueSelectorCOP select_protein_residues( core::pose::Pose const &pose ) const;

	/// @brief Given a pose, center it on the origin.
	///
	void align_to_origin( core::pose::Pose &pose ) const;

	/// @brief Given a pose centered on the origin, align it to the z-axis based on its symmetry.
	///
	void align_to_zaxis( core::pose::Pose &pose, core::Size const symmrepeats, bool const mirrorsymm ) const;

	/// @brief Given a pose that has been properly centered on the origin and aligned to the z-axis, delete all but
	/// a single symmetry repeat.
	void do_trim_to_single_repeat( core::pose::Pose &pose, core::Size const symmrepeats, core::Size const repeat_to_preserve) const;

private: // data

	/// @brief The number of symmetry repeats.  For example, for c4 or c4/m symmetry,
	/// this would be "4".  Defaults to "2".
	core::Size symmetry_repeats_;

	/// @brief Does this peptide have mirror symmetry?  Default "false".
	///
	bool mirror_symmetry_;

	/// @brief Should the mover auto-detect symmetry?  Default "false".
	/// @details Starts with the highets-fold allowed and works down, testing mirror symmetry before non-mirror
	/// symmetry.  Only based on backbone.
	bool auto_detect_symmetry_;

	/// @brief The angle threshold for detecting symmetry.
	/// @details Two dihedral values from different residues must fall within this cutoff, in degrees, for them to
	/// be considered the "same".  Default 10 degrees.
	core::Real angle_threshold_;

	/// @brief If true, all geometry that is not the protein backbone of a single symmetry repeat will be
	/// deleted.  (This includes any crosslinkers.)  False by default.
	bool trim_to_single_repeat_;

	/// @brief If trim_to_single_repeat_ is true, this is the symmetry repeat to preserve.  Default 1.
	/// @details Does nothing if trim_to_single_repeat_ is false.
	core::Size repeat_to_preserve_;

	/// @brief The symmetry of the last auto-detected peptide.
	/// @details Initialized to 0; only set if auto-detection is used.
	core::Size last_symmetry_repeats_;

	/// @brief The mirror status of the symmetry of the last auto-detected peptide.
	/// @details Initialized to false; only set if auto-detection is used.
	bool last_symmetry_mirror_;

	/// @brief Aligns the peptide normal with the inverse Z-axis instead of with the Z-axis.
	///
	bool invert_;

};

std::ostream &
operator<<( std::ostream & os, SymmetricCycpepAlign const & mover );

} //protocols
} //cyclic_peptide

#endif //protocols_cyclic_peptide_SymmetricCycpepAlign_HH
