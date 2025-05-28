// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/drug_design/ReactionGrow.hh
/// @brief apply RDKit's reaction mechanism to add a fragment to a ResidueType
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_protocols_drug_design_ReactionGrow_hh
#define INCLUDED_protocols_drug_design_ReactionGrow_hh

#include <protocols/drug_design/ReactionGrow.fwd.hh>
#include <protocols/drug_design/ReactionChemistry.hh>

#include <core/chemical/rdkit/RDKit.fwd.hh>

#include <core/chemical/AtomRefMapping.hh>

#include <core/chemical/MutableResidueType.fwd.hh>

#include <utility/tag/XMLSchemaGeneration.fwd.hh>

#include <string>

namespace protocols {
namespace drug_design {

class ReactionGrow  : public ReactionChemistry {
public:
	ReactionGrow();

	void apply( core::chemical::MutableResidueType & ) override;

	/// @brief The file which contains the fragments to add to input residue type.
	void fragment_database( std::string filename, bool append=false );

	/// @brief Reduce the fragment set to those which are compatible with the reactions.
	/// @details Discarded fragments are discarded permanently.
	/// Only call after all reactions are finalized for this Chemistry
	void prefilter_fragments();

	/// @brief Filter reaction list for those compatible with the given fragments
	/// @details Discarded reactions will be discarded permanently.
	/// Only call with the finalized fragment list.
	void
	prefilter_reactions() { ReactionChemistry::prefilter_reactions( fragments_, true ); }

	/// @brief If not empty, use property weighting based on the given property.
	void weight_by_property( std::string const & setting ) { property_name_ = setting; }
	std::string const & weight_by_property() const { return property_name_; }

	/// @brief Initialize any data members of this instance from an input tag
	/// and a DataMap object
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datacache
	) override;

	core::chemical::VDVDMapping
	get_mapping() const override;

	/// @brief Add a reaction to the list of reactions to use.
	/// Reaction should be written in the synthetic direction with a single product
	void
	add_reaction( ::RDKit::ChemicalReactionOP rxn, core::Real weight ) override;

	static std::string class_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:

	/// @brief The fragments to apply
	utility::vector1< ::RDKit::ROMolOP > fragments_;

	/// @brief If not empty, pick fragments based on the weighting by the given property.
	std::string property_name_;

	core::chemical::VDVDMapping mapping_;

};

} // namespace drug_design
} // namespace protocols

#endif
