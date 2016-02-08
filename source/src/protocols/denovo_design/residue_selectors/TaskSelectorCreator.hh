// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/denovo_design/residue_selectors/TaskCreators.hh
/// @brief  Creator for denovo design residue selectors
/// @author Tom Linsky (tlinsky@gmail.com)

#ifndef INCLUDED_protocols_denovo_design_residue_selectors_TaskSelectorCreator_HH
#define INCLUDED_protocols_denovo_design_residue_selectors_TaskSelectorCreator_HH

// Package headers
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/select/residue_selector/ResidueSelectorCreator.hh>

namespace protocols {
namespace denovo_design {
namespace residue_selectors {

class TaskSelectorCreator : public core::select::residue_selector::ResidueSelectorCreator {
public:
	virtual core::select::residue_selector::ResidueSelectorOP create_residue_selector() const;
	virtual std::string keyname() const;
	virtual void provide_selector_xsd( utility::tag::XMLSchemaDefinition & ) const;
};

} //namespace residue_selectors
} //namespace denovo_design
} //namespace protocols


#endif
