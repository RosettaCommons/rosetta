// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/drug_design/ChemicalReaction.hh
/// @brief A reaction object that keeps track of its fragment sets.
/// @author Yidan Tang (yidan.tang@vanderbilt.edu)

#ifndef INCLUDED_protocols_drug_design_ChemicalReaction_hh
#define INCLUDED_protocols_drug_design_ChemicalReaction_hh

#include <protocols/drug_design/ChemicalReaction.fwd.hh>

#include <core/chemical/rdkit/RDKit.fwd.hh>
#include <core/types.hh>

#include <utility/VirtualBase.hh>
#include <utility/vector1.hh>
#include <utility/vector0.hh>

#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/ChemReactions/Reaction.h>
#include <rdkit/GraphMol/ChemReactions/ReactionParser.h>

namespace protocols {
namespace drug_design {

class ChemicalReaction : public utility::VirtualBase {
public:

	//Constructor
	ChemicalReaction(
		std::string const & reaction_dir,
		std::string const & reaction_name,
		std::string const & reaction_smirks
	);

	ChemicalReaction(
			std::string const & reaction_name,
			std::string const & reaction_smirks
	);

	~ChemicalReaction();

	/// @brief Is the reaction valid, or was there an error loading the reaction?
	/// (Doesn't look at reactants.)
	bool
	reaction_valid() const;

	/// @brief is the reaction valid, with usable reaction list (deprecated)
	bool
	is_reaction_usable();

	std::string const &
	reaction_name() const { return name_; }


	// The number of components for this reaction
	core::Size
	nreagents() const;

	// The number of loaded reagents for the given reagent number (deprecated)
	core::Size
	n_availible_reagents(core::Size reag_no);

	// Get a specific reagent (deprecated)
	::RDKit::ROMolOP
	reagent(core::Size reag_no, core::Size reag_index);

	/// Return a representative product structure
	::RDKit::ROMolOP
	representative_prod();

	// Lazy loading of reagents -- no-op if already loaded (deprecated)
	void load_reagents();

	/// Do the corresponding reaction, and return the result.
	/// The passed parameter is the selection for each reagent.
	/// Can return a nullptr if the reaction doesn't work.
	::RDKit::RWMolOP
	reaction_result( ::RDKit::MOL_SPTR_VECT const & reagents ) const;

	::RDKit::ChemicalReactionOP const &
	 get_reaction() const;

private:

	bool cleanup_product( ::RDKit::RWMol & prod ) const;

	// Returns the line-by-line contents of the relevant file. (deprecated)
	static
	utility::vector1< std::string >
	load_file( std::string const & reaction_dir, std::string const & filename );

private:

	std::string dir_;

	// The designation for this reaction
	std::string name_;

	// The reaction itself.
	::RDKit::ChemicalReactionOP rxn_;

	// The reagents.
	// The outer list is indexed reaction reactant number.
	// The inner list is the set of reagents
	// This is dynamically loaded, as needed.
	utility::vector0< utility::vector1< std::string > > reagent_smiles_;
	utility::vector0< utility::vector1< ::RDKit::RWMolOP > > reagents_;

	// A way to short-circuit re-trying the loading of reagents.
	bool reagents_bad_ = false;
};


}
}

#endif


