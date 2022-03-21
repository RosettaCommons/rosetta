// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/jump_selector/JumpForResidue.hh
/// @brief  JumpForResidue selects the jump that builds the residues passed in
/// @author Jack Maguire, jackmaguire1444@gmail.com

#ifndef INCLUDED_core_select_jump_selector_JumpForResidue_HH
#define INCLUDED_core_select_jump_selector_JumpForResidue_HH

// Unit headers
#include <core/select/jump_selector/JumpForResidue.fwd.hh>

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

class JumpForResidue : public JumpSelector {
public:
	// derived from base class
	JumpForResidue() = default;
	~JumpForResidue() override = default;

	JumpForResidue( JumpForResidue const & ) = default;
	JumpForResidue & operator=( JumpForResidue const & ) = default;

	JumpForResidue( JumpForResidue && ) = default;
	JumpForResidue & operator=( JumpForResidue && ) = default;

	// helpful ctor
	JumpForResidue( core::Size const resid );

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
	set_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector ){
		selector_ = selector;
	}

	///@brief If false, we will assert that every residue passed in is built by the same jump
	void
	set_allow_multiple_results( bool const setting ){
		allow_multiple_results_ = setting;
	}

private:
	core::select::residue_selector::ResidueSelectorCOP selector_ = nullptr;
	bool allow_multiple_results_ = false;
};

} //namespace jump_selector
} //namespace select
} //namespace core


#endif
