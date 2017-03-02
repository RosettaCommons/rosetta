// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/jump_selector/AndJumpSelectorCreator.hh
/// @brief  Creator for Boolean "AND" jump selector
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_select_jump_selector_AndJumpSelectorCreator_HH
#define INCLUDED_core_select_jump_selector_AndJumpSelectorCreator_HH

// Package headers
#include <core/select/jump_selector/JumpSelector.fwd.hh>
#include <core/select/jump_selector/JumpSelectorCreator.hh>

namespace core {
namespace select {
namespace jump_selector {

class AndJumpSelectorCreator : public core::select::jump_selector::JumpSelectorCreator {
public:
	core::select::jump_selector::JumpSelectorOP create_jump_selector() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const override;
};

} //namespace jump_selectors
} //namespace select
} //namespace core


#endif
