// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/drug_design/ReactionFragment.hh
/// @brief apply RDKit's reaction-based fragment addition to a ResidueType
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <protocols/drug_design/ReactionFragment.hh>
#include <protocols/drug_design/ReactionFragmentCreator.hh>

#include <core/chemical/rdkit/RDMolToRestype.hh>
#include <core/chemical/rdkit/RestypeToRDMol.hh>
#include <core/chemical/rdkit/util.hh>
#include <core/chemical/MutableResidueType.hh>
#include <protocols/chemistries/util.hh>

#include <utility/string_util.hh>
#include <utility/io/izstream.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <basic/database/open.hh>
#include <basic/Tracer.hh>

#include <numeric/random/WeightedSampler.hh>
#include <numeric/random/random.hh>

#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/ChemReactions/Reaction.h>
#include <rdkit/GraphMol/ChemReactions/ReactionParser.h>
#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/SmilesParse/SmilesWrite.h> // For MolToSmiles
#include <rdkit/GraphMol/SmilesParse/SmartsWrite.h> // For MolToSmarts

namespace protocols {
namespace drug_design {

static basic::Tracer TR("protocols.drug_design.ReactionFragment");

//------------------------- Creator -----------------------------

protocols::chemistries::ChemistryOP
ReactionFragmentCreator::create_chemistry() const {
	return protocols::chemistries::ChemistryOP( new ReactionFragment );
}

std::string
ReactionFragmentCreator::keyname() const {
	return ReactionFragment::class_name();
}

void
ReactionFragmentCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	ReactionFragment::provide_xml_schema( xsd );
}


//------------------------- Chemistry -----------------------------

ReactionFragment::ReactionFragment():
	ReactionChemistry(class_name()),
	keep_bigger_( false ),
	keep_random_( false )
{}

void
ReactionFragment::keep_atom( std::string const & keep_atom ) {
	keep_atom_ = keep_atom;
	utility::strip_whitespace( keep_atom_ );
}

// Overloaded function which reverses the reaction
void
ReactionFragment::add_reaction( ::RDKit::ChemicalReactionOP rxn, core::Real weight ) {
	if ( rxn->getNumProductTemplates() != 1 ) {
		TR.Error << "Input reaction should have one and only one listed product: " << ::RDKit::ChemicalReactionToRxnSmarts( *rxn ) << std::endl;
		return;
	}
	if ( rxn->getNumReactantTemplates() < 1 ) {
		TR.Error << "Input reaction does not have listed reactants: " << ::RDKit::ChemicalReactionToRxnSmarts( *rxn ) << std::endl;
		return;
	}

	TR.Debug << "Reaction before reversing: " << ::RDKit::ChemicalReactionToRxnSmarts( *rxn ) << std::endl;
	::RDKit::ChemicalReactionOP rev_rxn( new ::RDKit::ChemicalReaction );
	rev_rxn->addReactantTemplate( *rxn->beginProductTemplates()  );
	for ( ::RDKit::MOL_SPTR_VECT::const_iterator itr( rxn->beginReactantTemplates()),
			itr_end(rxn->endReactantTemplates()); itr != itr_end; ++itr ) {
		rev_rxn->addProductTemplate( *itr );
	}
	rev_rxn->initReactantMatchers();

	TR.Debug << "Reaction after reversing: " << ::RDKit::ChemicalReactionToRxnSmarts( *rev_rxn ) << std::endl;
	ReactionChemistry::add_reaction( rev_rxn, weight );

}

/// @brief Check if the produced products are valid given the reaction.
/// (The ability to run the reaction comes from matching the *reacting* templates
/// the product templates may have additional restrictions which could invalidate a product.)
bool
invalid_products( ::RDKit::MOL_SPTR_VECT const & products, ::RDKit::ChemicalReactionCOP rxn ) {
	if ( products.size() == 0 ) {
		TR.Debug << "Skipping reaction because it doesn't produce any products! " << ::RDKit::ChemicalReactionToRxnSmarts( *rxn ) << std::endl;
		return true;
	}

	::RDKit::MOL_SPTR_VECT::const_iterator templ_itr( rxn->beginProductTemplates() ), templ_itr_end( rxn->endProductTemplates() );
	::RDKit::MOL_SPTR_VECT::const_iterator prod_itr( products.begin() ), prod_itr_end( products.end() );
	for ( /*none*/ ; templ_itr != templ_itr_end && prod_itr != prod_itr_end; ++templ_itr, ++prod_itr ) {
		// We need to tidy the products a bit before we can use them in a query.
		::RDKit::RWMOL_SPTR cleaned( new ::RDKit::RWMol(**prod_itr) );
		try {
			::RDKit::MolOps::sanitizeMol(*cleaned);
		} catch (::RDKit::MolSanitizeException &se){
			TR.Debug << "Cannot Sanitize fragment with RDKit: " << se.what() << std::endl;
			return true;
		}

		::RDKit::MatchVectType tvect;
		if ( ! SubstructMatch(*cleaned,**templ_itr,tvect) ) {
			// Not really an error - the reactions are generally going to be more lax in the reverse direction
			TR << "When fragmenting, fragment: " << ::RDKit::MolToSmiles( *cleaned ) << " doesn't match " << ::RDKit::MolToSmarts( **templ_itr ) << std::endl;
			TR << "   For reaction: " << ::RDKit::ChemicalReactionToRxnSmarts( *rxn ) << std::endl;
			TR << "   -- Ignoring potential fragmentation." << std::endl;
			return true;
		}
	}
	TR.Debug << "Products are valid!" << std::endl;
	return false;
}

void
ReactionFragment::apply( core::chemical::MutableResidueType & rsdtype )
{

	if ( keep_bigger_ && keep_random_ ) {
		utility_exit_with_message("Cannot set both keep_bigger and keep_random for ReactionFragment.");
	}
	if ( keep_bigger_ && keep_atom_ != "" ) {
		utility_exit_with_message("Cannot set both keep_bigger and keep_atom for ReactionFragment.");
	}
	if ( keep_random_ && keep_atom_ != "" ) {
		utility_exit_with_message("Cannot set both keep_random and keep_atom for ReactionFragment.");
	}

	core::chemical::rdkit::RestypeToRDMol to_converter(rsdtype);
	::RDKit::RWMOL_SPTR rdmol( to_converter.Mol() ); // Convert
	core::chemical::rdkit::label_with_index(*rdmol, "Original_Index"); // (Re)Label with "Original_Index" to use quick mapping behavior.
	// Properties carry through the reaction, except for a "bug" with the explicitly mentioned atoms.

	TR.Debug << "On molecule: " << ::RDKit::MolToSmiles( *rdmol ) << std::endl;

	// Select Reactions
	utility::vector1< ::RDKit::ChemicalReactionOP > rxns;
	numeric::random::WeightedSampler rxn_sampler;
	filter_reactions(*rdmol, rxns, rxn_sampler);

	if ( rxns.size() == 0 ) {
		TR.Warning << "No suitable reactions found. Doing nothing." << std::endl;
		mapping_.clear();
		mapping_.identity( true );
		set_last_status( core::chemical::modifications::FAIL_DO_NOT_RETRY );
		return;
	}

	// Keep trying until we get a reaction which works
	// but don't get into an infinite loop
	core::Size iterations_to_try( rxns.size()*10 );
	for ( core::Size cc(1); rxn_sampler.update_cumulative_distribution() && cc <= iterations_to_try; ++cc ) {

		// Pick a random reaction
		core::Size rxn_num( rxn_sampler.random_sample() );
		::RDKit::ChemicalReactionOP rxn( rxns[rxn_num] );

		TR.Debug << "   trying fragmenting with reaction: " << ::RDKit::ChemicalReactionToRxnSmarts( *rxn ) << std::endl;

		// Run the reaction
		::RDKit::MOL_SPTR_VECT reactants;
		reactants.push_back( rdmol );
		std::vector< ::RDKit::MOL_SPTR_VECT > all_products( rxn->runReactants( reactants ) );
		if ( all_products.empty() ) {
			// This shouldn't happen, because we should have pre-checked for it, but ...
			utility_exit_with_message("Error in pre-checking reaction: " + ::RDKit::ChemicalReactionToRxnSmarts( *rxn ) );
		}

		// Filter out the reactions which don't give valid products
		// (Since we reversed the reaction, they might be more specific than RDKit would generate.)
		// use the erase-remove idiom
		all_products.erase(
			std::remove_if( all_products.begin(), all_products.end(),
			std::bind( invalid_products, std::placeholders::_1, rxn ) ),
			all_products.end() );

		if ( all_products.empty() ) {
			// No proper products - remove this reaction from consideration and move on to the next.
			TR.Debug << "      => Reaction produced only products that didn't match the reactant conditions." << std::endl;
			rxn_sampler.set_weight(rxn_num, 0.0 );
			continue;
		}
		// Pick a random product set
		core::Size prod_num( numeric::random::random_range(0,all_products.size()-1) );
		::RDKit::MOL_SPTR_VECT const & products( all_products[ prod_num ] );

		if ( products.size() == 0 ) {
			TR.Warning << "Skipping reaction because it doesn't produce any products! " << ::RDKit::ChemicalReactionToRxnSmarts( *rxn ) << std::endl;
			// Ignore this reaction -- the arity of products should be the same for all variants
			rxn_sampler.set_weight(rxn_num, 0.0 );
			continue;
		}

		// Pick which product to use
		core::Size selected_product( -1 );
		if ( keep_atom_.size() ) {
			if ( ! rsdtype.has( keep_atom_ ) ) {
				utility_exit_with_message("Cannot find atom '" + keep_atom_ + "' in residue " + rsdtype.name() );
			}
			core::Size reactant_atom_index( to_converter.vd_to_index()[ rsdtype.atom_vertex( keep_atom_ ) ] );
			// Pick the fragment with a given atom
			for ( core::Size ii(0); ii < products.size(); ++ii ) {
				// ?? Do we need to tidy the product before finding a mapping?
				// Hopefully, the "Original_Index" property carries through, and can be used as a quick shortcut.
				core::chemical::IndexIndexMapping rxn_map( core::chemical::rdkit::find_mapping( rdmol, products[ii], "Original_Index") ); // -1 is invalid
				core::Size keep_index = rxn_map[ reactant_atom_index ];
				if ( keep_index != core::Size( -1 ) ) {
					selected_product = ii;
					break;
				}
			}
			TR.Debug << "Keeping product " << selected_product << " with atom " << keep_atom_ << std::endl;
		} else if ( keep_bigger_ ) {
			core::Size max_size(0);
			for ( core::Size ii(0); ii < products.size(); ++ii ) {
				if ( products[ii]->getNumHeavyAtoms() > max_size ) {
					max_size = products[ii]->getNumHeavyAtoms();
					selected_product = ii;
				}
			}
			TR.Debug << "Keeping biggest product " << selected_product << " of size " << max_size << std::endl;
		} else if ( keep_random_ ) {
			// Pick a random product
			selected_product = numeric::random::random_range(0,products.size());
			TR.Debug << "Keeping random product " << selected_product << std::endl;
		} else {
			// Pick the first product
			// This makes a nice symmetry with ReactionGrow, which goes the other direction and assumes you start with the first
			selected_product = 0;
		}

		if ( selected_product == core::Size(-1) ) {
			TR << "Unable to find the desired product! Using the first one instead." << std::endl;
			selected_product = 0;
		}

		::RDKit::RWMolOP prod( new ::RDKit::RWMol( *products[ selected_product ] ) );

		if ( cleanup_product( *prod ) ) {
			TR.Debug << "      => Reaction product failed in cleanup." << std::endl;
			continue; // Failed at cleanup - try the reaction again.
		}

		// We need to find a mapping from the pre-reaction molecule to the post-reaction molecule
		// Hopefully, the "Original_Index" property carries through, and can be used as a quick shortcut.
		core::chemical::IndexIndexMapping rxn_map( core::chemical::rdkit::find_mapping( rdmol, prod, "Original_Index") ); // -1 is invalid

		if ( rxn_map.empty() ) {
			TR.Warning << "Cannot find substructure mapping of fragment to parent." << std::endl;
			TR.Debug << "Fragment: " << ::RDKit::MolToSmiles( *prod ) << std::endl;
			TR.Debug << "Parent: " << ::RDKit::MolToSmiles( *rdmol ) << std::endl;
			TR.Debug << "Reaction: " << ::RDKit::ChemicalReactionToRxnSmarts( *rxn ) << std::endl;
			continue;
		}

		TR.Debug << "Selected fragmentation product: " << ::RDKit::MolToSmiles( *prod ) << std::endl;

		// Should be good. Now convert the residue into a Rosetta residue type.
		core::chemical::VDIndexMapping restype_prod_map( combine( to_converter.vd_to_index(), rxn_map ) );

		core::chemical::rdkit::RDMolToRestype from_converter(*prod);
		from_converter.set_nbr( restype_prod_map[rsdtype.nbr_vertex()] );

		core::chemical::MutableResidueTypeOP new_resop( from_converter.generate_restype(rsdtype,restype_prod_map) );

		TR << "Removed " << rsdtype.nheavyatoms() - new_resop->nheavyatoms() << " heavyatoms from " << rsdtype.name() << std::endl;

		mapping_ = combine( restype_prod_map, from_converter.index_to_vd() );
		rsdtype = *new_resop;
		mapping_ = combine( mapping_, combine( core::chemical::VDStringMapping(*new_resop), core::chemical::StringVDMapping(rsdtype)) );

		// That's it, we're good
		set_last_status( core::chemical::modifications::SUCCESS );
		return;
	}

	TR.Warning << "Tried to run the fragmenting reactions " << iterations_to_try << " times, but none of them worked. Giving up." << std::endl;
	set_last_status( core::chemical::modifications::FAIL_RETRY );

}

core::chemical::VDVDMapping
ReactionFragment::get_mapping() const {
	return mapping_;
}

std::string
ReactionFragment::class_name() {
	return "ReactionFragment";
}

void
ReactionFragment::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute("reactions", xs_string,
		"The name of the file containing the SMARTS-based reactions to use, written in the synthetic direction.")
		+ XMLSchemaAttribute::attribute_w_default("keep_random", xsct_rosetta_bool,
		"If true, keep a random reaction product (instead of the first).", "0")
		+ XMLSchemaAttribute::attribute_w_default("keep_bigger", xsct_rosetta_bool,
		"If true, keep the biggest reaction product", "0")
		+ XMLSchemaAttribute("keep_atom", xs_string,
		"If set, keep the reaction product with the given object");

	protocols::chemistries::xsd_type_definition_w_attributes(
		xsd, class_name(),
		"Split a molecule in two based on a list of reactions, keeping on one portion.",
		attlist );
}

/// @brief Initialize any data members of this instance from an input tag
/// and a DataMap object
void
ReactionFragment::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &) {

	reaction_file( tag->getOption<std::string>("reactions") );

	keep_bigger( tag->getOption<bool>("keep_bigger", false) );
	keep_random( tag->getOption<bool>("keep_random", false) );
	keep_atom( tag->getOption<std::string>("keep_atom", "") );

	if ( keep_bigger_ && keep_random_ ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "Cannot set both keep_bigger and keep_random for RandomFragmentLigand.");
	}
	if ( keep_bigger_ && keep_atom_ != "" ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "Cannnot set both keep_bigger and keep_atom for RandomFragmentLigand.");
	}
	if ( keep_random_ && keep_atom_ != "" ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "Cannot set both keep_random and keep_atom for RandomFragmentLigand.");
	}

}


}
}
