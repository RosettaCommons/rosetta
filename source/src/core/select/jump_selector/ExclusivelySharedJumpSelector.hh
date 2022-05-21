// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/jump_selector/ExclusivelySharedJumpSelector.hh
/// @brief  ExclusivelySharedJumpSelector selects the jump that builds ALL and ONLY the residues passed in
/// @details unit tested in EnsureExclusivelySharedJumpMover.cxxtest.hh
/// @author Jack Maguire

#ifndef INCLUDED_core_select_jump_selector_ExclusivelySharedJumpSelector_HH
#define INCLUDED_core_select_jump_selector_ExclusivelySharedJumpSelector_HH

// Unit headers
#include <core/select/jump_selector/ExclusivelySharedJumpSelector.fwd.hh>

// Package headers
#include <core/select/jump_selector/JumpSelector.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers

namespace core {
namespace select {
namespace jump_selector {

class ExclusivelySharedJumpSelector : public JumpSelector {
public:
	// derived from base class
	ExclusivelySharedJumpSelector() = default;
	~ExclusivelySharedJumpSelector() override = default;

	ExclusivelySharedJumpSelector( ExclusivelySharedJumpSelector const & ) = default;
	ExclusivelySharedJumpSelector & operator=( ExclusivelySharedJumpSelector const & ) = default;

	ExclusivelySharedJumpSelector( ExclusivelySharedJumpSelector && ) = default;
	ExclusivelySharedJumpSelector & operator=( ExclusivelySharedJumpSelector && ) = default;

	///@brief helpful ctor
	ExclusivelySharedJumpSelector( residue_selector::ResidueSelectorCOP selector );

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

	///@brief We will select the jumps that build the residues selected by this residue selector
	void
	set_residue_selector( residue_selector::ResidueSelectorCOP selector ){
		selector_ = selector;
	}

private:
	residue_selector::ResidueSelectorCOP selector_ = nullptr;
};

} //namespace jump_selector
} //namespace select
} //namespace core


#endif
