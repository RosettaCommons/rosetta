// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/drug_design/ReactionChemistry.hh
/// @brief An abstract base class for Chemistries which use RDKit reactions
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_protocols_drug_design_ReactionChemistry_hh
#define INCLUDED_protocols_drug_design_ReactionChemistry_hh

#include <core/chemical/rdkit/RDKit.fwd.hh>

#include <protocols/chemistries/Chemistry.hh>
#include <core/chemical/AtomRefMapping.hh>

#include <core/chemical/MutableResidueType.fwd.hh>

#include <numeric/random/WeightedSampler.hh>

#include <string>

#include <rdkit/GraphMol/Substruct/SubstructMatch.h>

namespace protocols {
namespace drug_design {

class ReactionChemistry  : public protocols::chemistries::Chemistry {

private:
	ReactionChemistry(); // Have to set the type string

public:
	ReactionChemistry( std::string const & type );

	void apply( core::chemical::MutableResidueType & ) override = 0;

	/// @brief The file which contains the reactions which to use.
	virtual void
	reaction_file( std::string filename, bool append=false );

	/// @brief Add a reaction to the list of reactions known by this Chemistry.
	virtual void
	add_reaction( ::RDKit::ChemicalReactionOP rxn, core::Real weight);

	/// @brief Filter reaction list for those compatible with the given reactants
	/// @details Discarded reactions will be discarded permanently.
	/// Only call with the finalized fragment list.
	void
	prefilter_reactions( utility::vector1< ::RDKit::ROMolOP >  const & reactants, bool exclude_first = true);

	/// @brief Initialize any data members of this instance from an input tag
	/// and a DataMap object
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datacache
	) override = 0;

	core::chemical::VDVDMapping
	get_mapping() const override = 0;

protected: // methods

	/// @brief Filter the reactions to just the ones which can be applied to rdmol
	void
	filter_reactions(
		::RDKit::ROMol const & rdmol,
		utility::vector1< ::RDKit::ChemicalReactionOP > & rxns, // Return by reference.
		numeric::random::WeightedSampler & rxn_sampler ) const;

	/// @brief attempt to clean up the product of an RDKit reaction
	/// Modifies the RWMol in place - if modification unsuccessful, it will return true.
	bool
	cleanup_product( ::RDKit::RWMol & prod ) const;

	utility::vector1< std::pair< ::RDKit::ChemicalReactionOP, core::Real > > const &
	get_reactions() const;
	

private: // data

	/// @brief The chemical reactions to use
	utility::vector1< std::pair< ::RDKit::ChemicalReactionOP, core::Real > > reactions_;

};

} // namespace drug_design
} // namespace protocols

#endif
