// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/jump_selector/OrJumpSelector.hh
/// @brief  The OrJumpSelector combines logic from multiple JumpSelectors
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_select_jump_selector_OrJumpSelector_HH
#define INCLUDED_core_select_jump_selector_OrJumpSelector_HH

// Unit headers
#include <core/select/jump_selector/OrJumpSelector.fwd.hh>
#include <core/select/jump_selector/JumpSelector.hh>

// Package headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/select/jump_selector/JumpSelectorCreator.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.hh>

// C++ headers
#include <list>

namespace core {
namespace select {
namespace jump_selector {

/// @brief The OrJumpSelector combines the output of multiple JumpSelectors using OR
/// logic, i.e., jumps selected by ANY of the contained JumpSelectors will be selected.
/// JumpSelectors can be pulled in from a DataMap, from subtags (for JumpSelectors
/// known to the JumpSelectorFactory) or programmatically through %add_jump_selector.
class OrJumpSelector : public JumpSelector {
public:
	// derived from base class
	OrJumpSelector();

	/// @brief Copy constructor
	///
	OrJumpSelector( OrJumpSelector const &src);

	OrJumpSelector( JumpSelectorCOP selector1, JumpSelectorCOP selector2 );
	virtual ~OrJumpSelector();

	/// @brief Clone operator.
	/// @details Copy this object and return an owning pointer to the new object.
	JumpSelectorOP clone() const override;

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
	* @brief adds a JumpSelector
	*/
	void add_jump_selector(JumpSelectorCOP selector);

	Size num_selectors() const;

	/**
	* @brief applies newSubset to existingSubset and thereby modifies the latter
	*/
	void apply_or_to_subset(JumpSubset const & newSubset, JumpSubset & existingSubset) const;

private: // data members

	std::list< JumpSelectorCOP > selectors_;

};


} //namespace jump_selector
} //namespace select
} //namespace core

#endif
