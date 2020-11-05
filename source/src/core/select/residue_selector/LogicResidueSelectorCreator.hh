// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/select/residue_selector/LogicResidueSelectorCreator.hh
/// @brief A residue selector that takes arbitrarily many residue selectors and performs and/or/not boolean logic within the "residue_selectors" option.
/// @author francesc45 (francesc345@gmail.com)

#ifndef INCLUDED_core_select_residue_selector_LogicResidueSelectorCreator_HH
#define INCLUDED_core_select_residue_selector_LogicResidueSelectorCreator_HH

// Package headers
#include <core/select/residue_selector/ResidueSelectorCreator.hh>

namespace core {
namespace select {
namespace residue_selector {

class LogicResidueSelectorCreator : public  core::select::residue_selector::ResidueSelectorCreator {
public:
	core::select::residue_selector::ResidueSelectorOP create_residue_selector() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const override;
};

} //residue_selector
} //select
} //core

#endif //INCLUDED_core_select_residue_selector_LogicResidueSelectorCreator_HH
