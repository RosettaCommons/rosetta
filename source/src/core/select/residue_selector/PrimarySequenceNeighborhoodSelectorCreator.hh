// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/PrimarySequenceNeighborhoodCreators.hh
/// @brief  Creator for denovo design residue selectors
/// @author Tom Linsky (tlinsky@gmail.com)

#ifndef INCLUDED_core_select_residue_selector_PrimarySequenceNeighborhoodSelectorCreator_HH
#define INCLUDED_core_select_residue_selector_PrimarySequenceNeighborhoodSelectorCreator_HH

// Package headers
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/select/residue_selector/ResidueSelectorCreator.hh>

namespace core {
namespace select {
namespace residue_selector {

class PrimarySequenceNeighborhoodSelectorCreator : public ResidueSelectorCreator {
public:
	virtual ResidueSelectorOP create_residue_selector() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const;
};

} //namespace residue_selectors
} //namespace select
} //namespace core


#endif
