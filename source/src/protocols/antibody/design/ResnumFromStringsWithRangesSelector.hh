// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/ResnumFromStringsWithRangesSelector.hh
/// @brief Select residues with the protocols::antibody::design::get_resnums_from_strings_with_ranges() approach.
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_core_select_residue_selector_ResnumFromStringsWithRangesSelector_HH
#define INCLUDED_core_select_residue_selector_ResnumFromStringsWithRangesSelector_HH

// Unit headers

// Package headers
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.hh>

// C++ headers

namespace protocols {
namespace antibody {
namespace design {

/// @brief A selector which uses the protocols::antibody::design::get_resnums_from_strings_with_ranges() approach.
///
class ResnumFromStringsWithRangesSelector : public core::select::residue_selector::ResidueSelector {
public:
	typedef core::select::residue_selector::ResidueSelectorOP ResidueSelectorOP;
	typedef core::select::residue_selector::ResidueSubset ResidueSubset;

public:

	/// @brief Constructor.
	ResnumFromStringsWithRangesSelector();

	ResnumFromStringsWithRangesSelector( utility::vector1<std::string> const & pdb_residues );

	/// @brief Clone operator.
	/// @details Copy the current object (creating the copy on the heap) and return an owning pointer
	/// to the copy.  All ResidueSelectors must implement this.
	ResidueSelectorOP clone() const override;

public:

	/// @brief "Apply" function.
	/// @details Return the set subset.
	ResidueSubset apply( core::pose::Pose const & pose ) const override;

	/// @brief XML parse.
	/// @details Parse RosettaScripts tags and set up this mover.
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap
	) override;

	/// @brief Get the mover class name.
	std::string
	get_name() const override;

	/// @brief Get the mover class name.
	static std::string
	class_name();

	/// @brief Provide XSD information, enabling mechanical validation of input XML.
	static void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:

	utility::vector1<std::string> pdb_residues_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION


};


} //design
} //antibody
} //protocols


#endif //INCLUDEDcore/select/residue_selector_ResnumFromStringsWithRangesSelector_hh
