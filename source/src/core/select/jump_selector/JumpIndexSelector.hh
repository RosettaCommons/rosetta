// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/jump_selector/JumpUpstreamSelector.hh
/// @brief  The JumpIndexSelector selects a specific jump
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_select_jump_selector_JumpIndexSelector_HH
#define INCLUDED_core_select_jump_selector_JumpIndexSelector_HH

// Unit headers
#include <core/select/jump_selector/JumpIndexSelector.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/select/jump_selector/JumpSelector.hh>
#include <core/pose/Pose.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.hh>

// C++ headers
#include <set>

namespace core {
namespace select {
namespace jump_selector {

/// @brief The JumpIndexSelector returns a JumpSubset, i.e. a utility::vector1< bool > containing
/// 'true' for jump positions which lie upstream of a given jump in the FoldTree. The jump is
/// specified by its integer index.
class JumpIndexSelector : public JumpSelector {
public:
	// derived from base class
	JumpIndexSelector();

	/// @brief Copy constructor
	///
	JumpIndexSelector( JumpIndexSelector const & src );

	/// @brief Clone operator.
	/// @details Copy this object and return an owning pointer to the new object.
	JumpSelectorOP clone() const override;

	JumpIndexSelector( int jump );
	virtual ~JumpIndexSelector();

	JumpSubset apply( core::pose::Pose const & pose ) const override;
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &
	) override;

	std::string
	get_name() const override;

	static std::string class_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	int  jump() const;
	void jump( int jump );

private: // data members
	int jump_;

};

} //namespace jump_selector
} //namespace select
} //namespace core

#endif
