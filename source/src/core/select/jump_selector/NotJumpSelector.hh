// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/jump_selector/NotJumpSelector.hh
/// @brief  The NotJumpSelector negates the logic of its loaded JumpSelector
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_select_jump_selector_NotJumpSelector_HH
#define INCLUDED_core_select_jump_selector_NotJumpSelector_HH

// Unit headers
#include <core/select/jump_selector/NotJumpSelector.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/select/jump_selector/JumpSelector.hh>
#include <core/pose/Pose.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace select {
namespace jump_selector {

/// @brief The NotJumpSelector negates the input of one loaded JumpSelector, i.e., it is a logical NOT -
/// it selects all unselected jumps and deselects the selected ones.  True becomes false, false becomes true.
/// The JumpSelector to be negated can be pulled in through RosettaScipt using the selector option, subtags for
/// JumpSelectors known to the JumpSelectorFactory or programmatically using set_jump_selector.
/// Note that since most JumpSelectors clear the input JumpSubset, NOT can be thought of as simply selecting
/// the opposite of the passed in selector.
class NotJumpSelector : public JumpSelector {
public:
	// derived from base class
	NotJumpSelector();

	/// @brief Copy constructor
	///
	NotJumpSelector( NotJumpSelector const &src);

	/// @brief Clone operator.
	/// @details Copy this object and return an owning pointer to the new object.
	JumpSelectorOP clone() const override;

	NotJumpSelector( JumpSelectorCOP selector );
	virtual ~NotJumpSelector();

	JumpSubset apply( core::pose::Pose const & pose ) const override;
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap
	) override;

	std::string
	get_name() const override;

	static std::string class_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	//unit-specific
	/**
	* @brief sets a JumpSelector
	*/
	void set_jump_selector(JumpSelectorCOP selector);

private: // data members
	JumpSelectorCOP selector_;

};


} //namespace jump_selector
} //namespace select
} //namespace core


#endif
