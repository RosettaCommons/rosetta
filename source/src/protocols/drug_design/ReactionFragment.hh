// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/drug_design/ReactionFragment.hh
/// @brief apply RDKit's reaction mechanism to split a ResidueType
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_protocols_drug_design_ReactionFragment_hh
#define INCLUDED_protocols_drug_design_ReactionFragment_hh

#include <protocols/drug_design/ReactionFragment.fwd.hh>
#include <protocols/drug_design/ReactionChemistry.hh>

#include <core/chemical/rdkit/RDKit.fwd.hh>

#include <core/chemical/AtomRefMapping.hh>

#include <core/chemical/MutableResidueType.fwd.hh>

#include <utility/tag/XMLSchemaGeneration.fwd.hh>

#include <string>

namespace protocols {
namespace drug_design {

class ReactionFragment  : public ReactionChemistry {
public:
	ReactionFragment();

	void keep_bigger( bool setting ) { keep_bigger_ = setting; }
	void keep_random( bool setting ) { keep_random_ = setting; }
	void keep_atom( std::string const & keep_atom );

	bool keep_bigger() const { return keep_bigger_; }
	bool keep_random() const { return keep_random_; }
	std::string keep_atom() const { return keep_atom_; }

	void apply( core::chemical::MutableResidueType & ) override;

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

	/// @brief Keep the fragment which is bigger.
	bool keep_bigger_;

	/// @brief Keep a random fragment.
	bool keep_random_;

	/// @brief Keep the fragment which contains the given atom.
	std::string keep_atom_;

	core::chemical::VDVDMapping mapping_;

};

} // namespace drug_design
} // namespace protocols

#endif
