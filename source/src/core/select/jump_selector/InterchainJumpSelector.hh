// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/jump_selector/InterchainJumpSelector.hh
/// @brief  The InterchainJumpSelector selects all jumps that land on a different chain from the take off
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_select_jump_selector_InterchainJumpSelector_HH
#define INCLUDED_core_select_jump_selector_InterchainJumpSelector_HH

// Unit headers
#include <core/select/jump_selector/InterchainJumpSelector.fwd.hh>

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

/// @brief The InterchainJumpSelector returns a JumpSubset, i.e. a utility::vector1< bool > containing
/// 'true' for jump where the take-off residue is on a different chain from the landing residue.
class InterchainJumpSelector : public JumpSelector {
public:
	// derived from base class
	InterchainJumpSelector();

	virtual ~InterchainJumpSelector();

	/// @brief Copy constructor
	///
	InterchainJumpSelector( InterchainJumpSelector const &src);

	/// @brief Clone operator.
	/// @details Copy this object and return an owning pointer to the new object.
	JumpSelectorOP clone() const override;
	JumpSubset apply( core::pose::Pose const & pose ) const override;
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &
	) override;

	std::string
	get_name() const override;

	static std::string class_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

};

} //namespace jump_selector
} //namespace select
} //namespace core


#endif
