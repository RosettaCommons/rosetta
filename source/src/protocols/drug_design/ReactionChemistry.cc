// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/drug_design/ReactionChemistry.hh
/// @brief apply RDKit's reaction-based fragment addition to a ResidueType
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <protocols/drug_design/ReactionChemistry.hh>

#include <core/chemical/rdkit/RDMolToRestype.hh>
#include <core/chemical/rdkit/RestypeToRDMol.hh>
#include <core/chemical/rdkit/util.hh>
#include <core/chemical/MutableResidueType.hh>

#include <utility/string_util.hh>
#include <utility/io/izstream.hh>
#include <utility/tag/Tag.hh>
#include <basic/database/open.hh>
#include <basic/Tracer.hh>

#include <numeric/random/WeightedSampler.hh>
#include <numeric/random/random.hh>

#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/ChemReactions/Reaction.h>
#include <rdkit/GraphMol/ChemReactions/ReactionParser.h>
#include <rdkit/GraphMol/FileParsers/MolSupplier.h>
#include <rdkit/GraphMol/MolOps.h>
// For 3D conformation building
#include <rdkit/GraphMol/DistGeomHelpers/Embedder.h>
#include <rdkit/GraphMol/Conformer.h>
#include <rdkit/ForceField/ForceField.h>
// For substructure mapping
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/SmilesParse/SmilesWrite.h> // For MolToSmiles
#include <rdkit/GraphMol/FileParsers/FileParsers.h> // For MolToMolBlock


namespace protocols {
namespace drug_design {

static basic::Tracer TR("protocols.drug_design.ReactionChemistry");

ReactionChemistry::ReactionChemistry( std::string const & type ) :
	Chemistry(type)
{}

/// @brief The file which contains the reactions which to use.
void
ReactionChemistry::reaction_file( std::string filename, bool append /*=false*/ ) {

	if ( ! append ) {
		reactions_.clear();
	}

	if ( filename.size() == 0 ) {
		utility_exit_with_message("Cannot open empty reaction file in protocols::drug_design::ReactionChemistry.");
	}

	// Use a local reaction file if present, else use the one in the database
	std::string actual_filename( filename );
	utility::io::izstream data( actual_filename );
	if ( ! data.good() ) {
		actual_filename = basic::database::full_name("protocol_data/drug_design/"+filename);
		utility::io::izstream data( actual_filename );
		if ( ! data.good() ) {
			TR.Error << "Can't find reaction file, either (./)" << filename << " or " << actual_filename << std::endl;
			utility_exit_with_message("ERROR: cannot find reaction file.");
		}
	}

	if ( utility::endswith(actual_filename, ".sdf" ) ) {
		// This should be possible, but I don't know how yet.
		TR.Error << "Cannot currently read reactions in sdf format." << std::endl;
		utility_exit_with_message("Error trying to read an sdf-format reaction file.");
	}
	// SMARTS format file
	std::string line;
	while ( getline( data, line ) ) {
		utility::vector1< std::string > parts( utility::split_whitespace(line) );
		if ( parts.size() == 0 ) { continue; }
		if ( parts[1].size() < 1 || parts[1][0] == '#' ) { continue; }
		core::Real weight( 1.0 );
		if ( parts.size() > 1 && parts[2].size() >=1 && parts[2][0] != '#' ) {
			weight = utility::string2Size( parts[2] );
		}
		TR.Debug << "Parsing Reaction '" << parts[1] << "' " << std::endl;
		::RDKit::ChemicalReactionOP rxn( ::RDKit::RxnSmartsToChemicalReaction( parts[1] ) );
		if ( ! rxn ) {
			TR.Error << "Cannot parse as chemical reaction: " << parts[1] << std::endl;
			continue;
		}
		rxn->initReactantMatchers();
		TR.Debug << "Reaction as parsed: " << ::RDKit::ChemicalReactionToRxnSmarts( *rxn ) << std::endl;
		add_reaction( rxn, weight );
	}

	if ( reactions_.size() == 0 ) {
		TR.Warning << "No reactions found in reaction file " << actual_filename << std::endl;
	}
}

/// @brief Default reaction addition - just add as is;
void
ReactionChemistry::add_reaction( ::RDKit::ChemicalReactionOP rxn, core::Real weight) {
	if ( rxn->getNumProductTemplates() < 1 ) {
		TR.Error << "Reaction does not have listed products: " << ::RDKit::ChemicalReactionToRxnSmarts( *rxn ) << std::endl;
		return;
	}
	if ( rxn->getNumReactantTemplates() < 1 ) {
		TR.Error << "Reaction does not have listed reactants: " << ::RDKit::ChemicalReactionToRxnSmarts( *rxn ) << std::endl;
		return;
	}
	std::pair< ::RDKit::ChemicalReactionOP, core::Real > rxn_pair( rxn, weight );
	reactions_.push_back( rxn_pair );
}

void
ReactionChemistry::prefilter_reactions( utility::vector1< ::RDKit::ROMolOP > const & reactants, bool exclude_first ) {
	if ( reactions_.empty() ) {
		utility_exit_with_message("No reactions found when trying to pre-filter!");
	}
	utility::vector1< std::pair< ::RDKit::ChemicalReactionOP, core::Real > > filtered_reactions;

	for ( core::Size rr(1); rr <= reactions_.size(); ++rr ) {
		bool valid = true;
		::RDKit::ChemicalReactionOP rxn( reactions_[rr].first );
		for ( ::RDKit::MOL_SPTR_VECT::const_iterator itr( rxn->beginReactantTemplates() ),
				itr_end( rxn->endReactantTemplates() ); itr != itr_end; ++itr ) {
			if ( exclude_first && itr == rxn->beginReactantTemplates() ) { continue; }
			bool reactant_found = false;
			for ( core::Size ff(1); ff <= reactants.size(); ++ff ) {
				::RDKit::MatchVectType tvect;
				if ( SubstructMatch(*reactants[ff],**itr,tvect) ) {
					reactant_found = true;
					break;
				}
			}
			if ( !reactant_found ) {
				valid = false;
				break;
			}
		}
		if ( valid ) {
			filtered_reactions.push_back( reactions_[rr] );
		}
	}

	TR << "Prefiltered reactions from " << reactions_.size() << " to " << filtered_reactions.size() << " that are compatible with the reactants." << std::endl;
	swap( reactions_, filtered_reactions ); // Replace reactions_, and discard the old one.
	if ( reactions_.empty() ) {
		utility_exit_with_message("After filtering, no usable reactions remain!");
	}
}

/// @brief Filter the reactions to just the ones which can be applied to rdmol
void
ReactionChemistry::filter_reactions(
	::RDKit::ROMol const & rdmol,
	utility::vector1< ::RDKit::ChemicalReactionOP > & rxns,
	numeric::random::WeightedSampler & rxn_sampler ) const
{
	::RDKit::MatchVectType tvect;
	for ( core::Size rr(1); rr <= reactions_.size(); ++rr ) {
		::RDKit::ROMOL_SPTR first_reactant( *reactions_[rr].first->beginReactantTemplates() );
		if ( ::RDKit::SubstructMatch(rdmol,*first_reactant,tvect) ) {
			rxns.push_back( reactions_[rr].first );
			rxn_sampler.add_weight( reactions_[rr].second );
		}
	}
}

/// @brief attempt to clean up the product of an RDKit reaction
/// Modifies the RWMol in place - if modification unsuccessful, it will return true.
bool
ReactionChemistry::cleanup_product( ::RDKit::RWMol & prod ) const {
	// Clean up product molecule
	try {
		::RDKit::MolOps::sanitizeMol(prod);
	} catch (::RDKit::MolSanitizeException &se){
		TR.Error << "Cannot Sanitize product with RDKit: " << se.what() << std::endl;
		if ( TR.Debug.visible() ) {
			TR.Debug << "\n" << ::RDKit::MolToMolBlock( prod, /*includeStero*/ true, /*confId*/ -1, /*kekulize*/ false ) << std::endl;
		}
		return true;
	}

	// The geometry of the post-reaction compounds is typically bad (especially if we added fragments)
	// Re-embed the molecule to fix the conformation
	::RDKit::MolOps::addHs(prod,/*explicitOnly=*/false,/*addCoords=*/false); // Embedding (& min) works better if we have hydrogens.

	// This long call is needed such that we can put "useExpTorsionAnglePrefs=true" and "useBasicKnowledge=true" so embedding doesn't mess up rings, etc.
	// We also use a constant seed for embedding, to make the generated conformations consistent.
	int conf_num = ::RDKit::DGeomHelpers::EmbedMolecule(prod, /*maxIterations*/ 0, /*seed*/ 111111, /*clearConfs*/ true,
		/*useRandomCoords*/ false, /*boxSizeMult*/ 2.0, /*randNegEig*/ true, /*numZeroFail*/ 1, /*coordMap*/ nullptr, /*optimizerForceTol*/ 1e-3,
		/*ignoreSmoothingFailures*/ false, /*enforceChirality*/ true, /*useExpTorsionAnglePrefs*/ true, /*useBasicKnowledge*/ true);
	if ( conf_num == -1 ) {
		::RDKit::MolOps::removeHs(prod,false,true);
		TR.Warning << "Could not find 3D coordinates for reaction product: " << ::RDKit::MolToSmiles( prod ) << std::endl;
		if ( TR.Debug.visible() ) {
			TR.Debug << "\n" << ::RDKit::MolToMolBlock( prod, /*includeStero*/ true, /*confId*/ -1, /*kekulize*/ false ) << std::endl;
		}
		return true;
	}

	::RDKit::ForceFieldOP ff( core::chemical::rdkit::get_forcefield(prod, conf_num) );
	if ( ! ff ) {
		::RDKit::MolOps::removeHs(prod,false,true);
		return true; // More detailed error messages should be from the get_forcefield() function
	}

	TR.Debug << "Starting Energy: " << ff->calcEnergy() << std::endl;
	if ( ff->minimize(200) == 1 ) {
		TR.Debug << "Try with looser tolerances" << std::endl;
		TR.Debug << "Energy: " << ff->calcEnergy() << std::endl;
		// Try with looser tolerances
		if ( ff->minimize(200, 1e-2, 1e-4) == 1 ) {
			TR.Debug << "Okay, really loose tolerances" << std::endl;
			TR.Debug << "Energy: " << ff->calcEnergy() << std::endl;
			if ( ff->minimize(200, 1, 1) == 1 ) {
				TR.Warning << "In reaction cleanup, minimization did not converge." << std::endl;
				TR.Debug << "Energy: " << ff->calcEnergy() << std::endl;
				TR << "PostMin: \n" << ::RDKit::MolToMolBlock( prod ) << std::endl;
				::RDKit::MolOps::removeHs(prod);
				TR << "Product: " << ::RDKit::MolToSmiles( prod ) << std::endl;
				return true;
			}
		}
	}
	TR.Debug << "Final Energy: " << ff->calcEnergy() << std::endl;

	::RDKit::MolOps::removeHs(prod,false,true); // Renormalize the molecule into the "explicit" hydrogen form.
	return false;
}

utility::vector1< std::pair< ::RDKit::ChemicalReactionOP, core::Real > > const &
ReactionChemistry::get_reactions() const {
	return reactions_;
}

}
}
