// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/select/residue_selector/GlycanPositionSelector.hh
/// @brief  A Residue Selector for selecting specific parts of arbitrary glycans.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_select_residue_selector_GlycanPositionSelector_HH
#define INCLUDED_core_select_residue_selector_GlycanPositionSelector_HH

// Unit headers
#include <core/select/residue_selector/GlycanPositionSelector.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/ResidueRanges.hh>

#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.hh>

// C++ headers
#include <set>

namespace core {
namespace select {
namespace residue_selector {

/// @brief A Residue Selector for selecting specific parts of arbitrary glycans by 'position'.
///
/// @details.
///
///  Background
///
///    Lets assume that the ASN that a particular N-linked Glycan is attached to, starts from Residue 0.
///    The residues off this continue with 1 being the next residue, 2, and so on.  Each branch corresponds to a number.
///
///    This allows you to choose parts of the glycan, without knowing the actual glycan residue numbers.  For example, maybe you want to
///    select the outer part of all glycans or between specific positions.
///
///  Tips For use
///
///		This Selector works on all glycans of the pose at once.
///
///     Settings are:
///
///       range (start, end)
///       positions (specific resnums, such as the carbohydrate at position 4)

///       from_residue (all glycan foliage from this and including this residue.)
///       to_residue (all glycan foliage up to and including this residue.)
///
///     Use the 'glycan_info' application to determine the glycan position numbers from a pose.
///
///     Combine with the GlycanResidueSelector to get unions of specific glycans
///      (such as the leaf of all Man5 residues or the stem of the glycan that starts at ASN85.)
///

///
class GlycanPositionSelector : public core::select::residue_selector::ResidueSelector {
public:

	/// @brief Constructor.
	GlycanPositionSelector();

    /// @brief Copy Constructor.  Use if you have non-basic private variables (classes, OPs, etc.)
	GlycanPositionSelector(GlycanPositionSelector const & src);
	
public:
	
	/// @brief Add a range to the list.
	void
	add_range(ResidueRange const & range);
	
	/// @brief Set the range of glycan positions to select from.
	void
	set_range(utility::vector1<ResidueRange> const & ranges);
	
	/// @brief Set a specific set of positions to select on.
	void
	set_positions(utility::vector1< Size > const & positions);
	
	/// @brief Set the position from which to select all outer foliage from, (and includeing) this position.
	void
	set_select_from_residue_position( Size const select_from_residue_position);
	
	/// @brief Set the position where we will select all glycan residues up to this specific glycan position.
	void
	set_select_to_residue_position( Size const select_to_residue_position);
	
	
public:


	/// @brief Destructor.
	virtual
	~GlycanPositionSelector();

	/// @brief Clone operator.
	/// @details Copy the current object (creating the copy on the heap) and return an owning pointer
	/// to the copy.  All ResidueSelectors must implement this.
	virtual
	core::select::residue_selector::ResidueSelectorOP clone() const;

	/// @brief "Apply" function.
	/// @details Given the pose, generate a vector of bools with entries for every residue in the pose
	/// indicating whether each residue is selected ("true") or not ("false").
	virtual
	core::select::residue_selector::ResidueSubset apply( core::pose::Pose const & pose ) const;

	/// @brief XML parse.
	/// @details Parse RosettaScripts tags and set up this mover.
	virtual void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap
	);

	/// @brief Get the mover class name.
	virtual
	std::string
	get_name() const;

	/// @brief Get the mover class name.
	static std::string
	class_name();

	/// @brief Provide XSD information, enabling mechanical validation of input XML.
	static void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	
	// All Mutually exclusive options to choose parts of the glycan.
	utility::vector1< ResidueRange > ranges_;
	utility::vector1< Size > positions_;
	Size from_residue_;
	Size to_residue_;
	
};


} //core
} //select
} //residue_selector


#endif //INCLUDEDcore/select/residue_selector_GlycanPositionSelector_hh
