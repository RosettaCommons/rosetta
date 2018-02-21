// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/residue_selectors/StoredResidueSubsetSelector.hh
/// @brief  The StoredResidueSubsetSelector selects residues using a previously stored residue subset
/// @author Tom Linsky (tlinsky@uw.edu))

#ifndef INCLUDED_protocols_residue_selectors_StoredResidueSubsetSelector_HH
#define INCLUDED_protocols_residue_selectors_StoredResidueSubsetSelector_HH

// Unit headers
#include <protocols/residue_selectors/StoredResidueSubsetSelector.fwd.hh>

// Package headers
#include <core/pose/Pose.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.hh>

// Utility Headers
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers

namespace protocols {
namespace residue_selectors {

/// @brief The StoredResidueSubsetSelector returns a ResidueSubset, i.e. a utility::vector1< bool > containing
/// 'true' for residue positions which are located near the given selected residues in primary sequence space
class StoredResidueSubsetSelector : public core::select::residue_selector::ResidueSelector {
public:
	StoredResidueSubsetSelector();
	StoredResidueSubsetSelector( std::string const & subset_name ); // move-constructed

	~StoredResidueSubsetSelector() override;

	/// @brief Clone operator.
	/// @details Copy this object and return an owning pointer to the new object.
	core::select::residue_selector::ResidueSelectorOP
	clone() const override;

	core::select::residue_selector::ResidueSubset
	apply( core::pose::Pose const & pose ) const override;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & ) override;


	std::string
	get_name() const override;

	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & );

	static std::string class_name();

	std::string const &
	get_subset_name() const {
		return subset_name_;
	}

private: // data members
	std::string subset_name_;
};

} //namespace residue_selectors
} //namespace protocols


#endif
