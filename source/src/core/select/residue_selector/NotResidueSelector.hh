// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/select/residue_selector/NotResidueSelector.hh
/// @brief  The NotResidueSelector negates the logic of its loaded ResidueSelector
/// @author Robert Lindner (rlindner@mpimf-heidelberg.mpg.de)

#ifndef INCLUDED_core_select_residue_selector_NotResidueSelector_HH
#define INCLUDED_core_select_residue_selector_NotResidueSelector_HH

// Unit headers
#include <core/select/residue_selector/NotResidueSelector.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/pose/Pose.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>

namespace core {
namespace select {
namespace residue_selector {

/// @brief The NotResidueSelector negates the input of one loaded ResidueSelector, i.e., it is a logical NOT -
/// it selects all unselected residues and deselects the selected ones.  True becomes false, false becomes true.
/// The ResidueSelector to be negated can be pulled in through RosettaScipt using the selector option, subtags for
/// ResidueSelectors known to the ResidueSelectorFactory or programmatically using set_residue_selector.
/// Note that since most ResidueSelectors clear the input ResidueSubset, NOT can be thought of as simply selecting
/// the opposite of the passed in selector.
class NotResidueSelector : public ResidueSelector {
public:
	// derived from base class
	NotResidueSelector();

	/// @brief Copy constructor
	///
	NotResidueSelector( NotResidueSelector const &src);

	/// @brief Clone operator.
	/// @details Copy this object and return an owning pointer to the new object.
	virtual ResidueSelectorOP clone() const;

	NotResidueSelector( ResidueSelectorCOP selector );
	virtual ~NotResidueSelector();

	virtual ResidueSubset apply( core::pose::Pose const & pose ) const;
	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap
	);

	virtual
	std::string
	get_name() const;

	static std::string class_name();

	//unit-specific
	/**
	* @brief sets a ResidueSelector
	*/
	void set_residue_selector(ResidueSelectorCOP selector);

private: // data members
	ResidueSelectorCOP selector_;

};


} //namespace residue_selector
} //namespace select
} //namespace core


#endif
