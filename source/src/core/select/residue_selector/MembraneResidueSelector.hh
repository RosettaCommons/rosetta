// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/select/residue_selector/MembraneResidueSelector.hh
/// @brief  Select the memrbane residue in a membrane pose
/// @author Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_select_residue_selector_MembraneResidueSelector_HH
#define INCLUDED_core_select_residue_selector_MembraneResidueSelector_HH

// Unit headers
#include <core/select/residue_selector/MembraneResidueSelector.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pose/Pose.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>

// C++ headers
#include <set>

namespace core {
namespace select {
namespace residue_selector {

/// @brief Select the memrbane residue in a membrane pose
class MembraneResidueSelector : public core::select::residue_selector::ResidueSelector {
public:

	MembraneResidueSelector();
	virtual ~MembraneResidueSelector();

	virtual core::select::residue_selector::ResidueSubset apply( core::pose::Pose const & pose ) const;
	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap
	);

	virtual std::string get_name() const;
	virtual ResidueSelectorOP clone() const;

	static std::string class_name();
	static void provide_selector_xsd( utility::tag::XMLSchemaDefinition & xsd ); 

};

} // residue_selector
} // select
} // core

#endif // INCLUDED_core_select_residue_selector_MembraneResidueSelector_hh
