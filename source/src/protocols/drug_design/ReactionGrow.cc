// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/drug_design/ReactionGrow.hh
/// @brief apply RDKit's reaction-based fragment addition to a ResidueType
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <protocols/drug_design/ReactionGrow.hh>
#include <protocols/drug_design/ReactionGrowCreator.hh>
#include <protocols/drug_design/ReactionChemistry.hh>

#include <core/chemical/rdkit/RDMolToRestype.hh>
#include <core/chemical/rdkit/RestypeToRDMol.hh>
#include <core/chemical/rdkit/util.hh>
#include <core/chemical/MutableResidueType.hh>
#include <protocols/chemistries/util.hh>

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
#include <rdkit/GraphMol/FileParsers/MolSupplier.h>
#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/SmilesParse/SmilesWrite.h> // For MolToSmiles
#include <rdkit/GraphMol/SmilesParse/SmartsWrite.h> // For MolToSmarts

namespace protocols {
namespace drug_design {

static basic::Tracer TR("protocols.drug_design.ReactionGrow");

//------------------------- Creator -----------------------------

protocols::chemistries::ChemistryOP
ReactionGrowCreator::create_chemistry() const {
	return protocols::chemistries::ChemistryOP( new ReactionGrow );
}

std::string
ReactionGrowCreator::keyname() const {
	return ReactionGrow::class_name();
}

void
ReactionGrowCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	ReactionGrow::provide_xml_schema( xsd );
}

//------------------------- Chemistry -----------------------------

ReactionGrow::ReactionGrow():
	ReactionChemistry(class_name())
{}

void
ReactionGrow::add_reaction( ::RDKit::ChemicalReactionOP rxn, core::Real weight ) {
	if ( rxn->getNumProductTemplates() != 1 ) {
		TR.Error << "Input reaction should have one and only one listed product: " << ::RDKit::ChemicalReactionToRxnSmarts( *rxn ) << std::endl;
		return;
	}
	if ( rxn->getNumReactantTemplates() < 1 ) {
		TR.Error << "Input reaction does not have listed reactants: " << ::RDKit::ChemicalReactionToRxnSmarts( *rxn ) << std::endl;
		return;
	}

	ReactionChemistry::add_reaction( rxn, weight );
}

/// @brief The file which contains the fragments to add to input residue type.
void
ReactionGrow::fragment_database( std::string filename, bool append /*=false*/ ) {
	if ( ! append ) {
		fragments_.clear();
	}

	core::chemical::rdkit::load_sdf( filename, fragments_, /*removeHs=*/ true );

	if ( fragments_.size() == 0 ) {
		TR.Warning << "No molecule fragments found in file " << filename << std::endl;
	}
}

/// @brief Reduce the fragment set to those which are compatible with the reactions.
void
ReactionGrow::prefilter_fragments() {

	utility::vector1< ::RDKit::ROMolOP > filtered_fragments;

	utility::vector1< std::pair< ::RDKit::ChemicalReactionOP, core::Real > > const & reactions( get_reactions() );

	for ( core::Size ff(1); ff <= fragments_.size(); ++ff ) {
		bool found = false;
		for ( core::Size rr(1); rr <= reactions.size(); ++rr ) {
			::RDKit::ChemicalReactionOP rxn( reactions[rr].first );
			for ( ::RDKit::MOL_SPTR_VECT::const_iterator itr( ++(rxn->beginReactantTemplates()) ), // Skip the first reactant
					itr_end( rxn->endReactantTemplates() ); itr != itr_end; ++itr ) {
				::RDKit::MatchVectType tvect;
				if ( SubstructMatch(*fragments_[ff],**itr,tvect) ) {
					filtered_fragments.push_back( fragments_[ff] );
					found = true;
					break;
				}
			}
			if ( found ) { break; }
		}
	}

	TR << "Prefiltered fragments from " << fragments_.size() << " to " << filtered_fragments.size() << " that are compatible with the reactions." << std::endl;
	swap( fragments_, filtered_fragments ); // Replace fragments_, and discard the old one.
	if ( fragments_.empty() ) {
		// Fragment set may be empty, if someone is using a reaction set which is 1->1
		TR.Warning << "After filtering, no fragments remain!" << std::endl;
	}
}

void
ReactionGrow::apply( core::chemical::MutableResidueType & rsdtype )
{
	// Fragment set may be empty, if someone is using a reaction set which is 1->1

	core::chemical::rdkit::RestypeToRDMol to_converter(rsdtype); // Neutralize and remove hydrogens.
	::RDKit::RWMOL_SPTR rdmol( to_converter.Mol() ); // Convert
	core::chemical::rdkit::label_with_index(*rdmol, "Original_Index"); // (Re)Label with "Original_Index" to use quick mapping behavior.
	// Properties carry through the reaction, except for a "bug" with the explicitly mentioned atoms.

	// Select Reactions
	utility::vector1< ::RDKit::ChemicalReactionOP > rxns;
	numeric::random::WeightedSampler rxn_sampler;
	filter_reactions(*rdmol, rxns, rxn_sampler);

	TR.Debug << "Found " << rxns.size() << " reactions compatible with current molecule." << std::endl;

	if ( rxns.size() == 0 ) {
		TR.Warning << "No suitable reactions found. Doing nothing." << std::endl;
		mapping_.clear();
		mapping_.identity( true );
		set_last_status( core::chemical::modifications::FAIL_DO_NOT_RETRY );
		return;
	}

	// Keep trying until we get a reaction which works
	// but don't get into an infinite loop
	for ( core::Size cc(1); cc <= rxns.size()*10; ++cc ) {

		// Pick a random reaction
		core::Size rxn_num( rxn_sampler.random_sample() );
		::RDKit::ChemicalReactionOP rxn( rxns[rxn_num] );

		// Pick fragments which work as reactants
		::RDKit::MOL_SPTR_VECT reactants;
		for ( ::RDKit::MOL_SPTR_VECT::const_iterator itr( rxn->beginReactantTemplates() ),
				itr_end( rxn->endReactantTemplates() );
				itr != itr_end; ++itr ) {
			if ( reactants.size() == 0 ) {
				// First reactant is the molecule we're building.
				reactants.push_back( rdmol );
				continue;
			} else {
				// Fill out the rest of the reactions with fragments.
				utility::vector1< ::RDKit::ROMolOP > fragments;
				numeric::random::WeightedSampler fragment_sampler;
				for ( core::Size ff(1); ff <= fragments_.size(); ++ff ) {
					::RDKit::MatchVectType tvect;
					if ( SubstructMatch(*fragments_[ff],**itr,tvect) ) {
						core::Real weight( 1.0 );
						if ( property_name_.size() && fragments_[ff]->hasProp(property_name_) ) {
							weight = fragments_[ff]->getProp<core::Real>(property_name_);
						}
						fragments.push_back( fragments_[ff] );
						fragment_sampler.add_weight( weight );
					} else if ( TR.Trace.visible() ) {
						TR.Trace << "Fragment " << ::RDKit::MolToSmiles( *fragments_[ff] ) << " doesn't match " << ::RDKit::MolToSmarts( **itr ) << std::endl;
					}
				}
				TR.Debug << "From " << fragments_.size() << " total fragments, " << fragments.size() << " are compatible with reactant# " << reactants.size()+1 << std::endl;
				if ( fragments.size() == 0 ) {
					TR.Warning << "Reaction cannot be fufilled with current fragment set." << std::endl;
					TR.Warning << ::RDKit::ChemicalReactionToRxnSmarts( *rxn ) << " :: reactant #" << reactants.size()+1 << std::endl;
					break;
				}
				// Pick one of the applicable fragments to use.
				core::Size frag_num( fragment_sampler.random_sample() );
				reactants.push_back( fragments[ frag_num ] );
			}
		} // end pick reactants for loop

		if ( reactants.size() != rxn->getNumReactantTemplates() ) {
			// Couldn't find suitable products remove reaction from sampler.
			rxn_sampler.set_weight(rxn_num, 0.0 );
			if ( !rxn_sampler.update_cumulative_distribution() ) {
				// Removed all the reactions from sampler.
				TR.Warning << "No possible reactions with given fragment set. Doing nothing." << std::endl;
				mapping_.clear();
				mapping_.identity( true );
				set_last_status( core::chemical::modifications::FAIL_DO_NOT_RETRY );
				return;
			}
			TR.Debug << "Couldn't find fragments for selected reaction." << std::endl;
			continue;
		}

		if ( TR.Debug.visible() ) {
			TR.Debug << "Choosing the reaction " << ::RDKit::ChemicalReactionToRxnSmarts( *rxn ) << std::endl;
			for ( core::Size jj(0); jj < reactants.size(); ++jj ) {
				TR.Debug << "    Reactant " << jj+1 << ":  " << ::RDKit::MolToSmiles( *reactants[jj] ) << std::endl;
			}
		}

		// Actually run the reaction
		std::vector< ::RDKit::MOL_SPTR_VECT > products( rxn->runReactants( reactants ) );
		if ( products.size() == 0 ) {
			TR.Debug << "No products found for reaction = retrying." << std::endl;
			continue;
		}
		core::Size prod_num( numeric::random::random_range(0,products.size()-1) );
		if ( products[ prod_num].size() != 1 ) {
			TR.Error << "Bad reaction in ReactionGrow! Too many products (" << products[ prod_num].size() << "): "
				<< ::RDKit::ChemicalReactionToRxnSmarts( *rxn ) << std::endl;
			rxn_sampler.set_weight(rxn_num, 0.0 );
			if ( !rxn_sampler.update_cumulative_distribution () ) {
				// Removed all the reactions from sampler.
				TR.Warning << "No possible reactions with given fragment set. Doing nothing." << std::endl;
				mapping_.clear();
				mapping_.identity( true );
				set_last_status( core::chemical::modifications::FAIL_DO_NOT_RETRY );
				return;
			}
			continue;
		}

		::RDKit::RWMolOP prod( new ::RDKit::RWMol( *products[ prod_num ][0]) );

		if ( cleanup_product( *prod ) ) {
			TR.Debug << "Failed cleanup." << std::endl; // More extensive debugging information should be from the cleanup
			continue; // Failed at cleanup - try the reaction again.
		}

		TR.Debug << "    Product " << ::RDKit::MolToSmiles( *prod ) << std::endl;

		// We need to find a mapping from the pre-reaction molecule to the post-reaction molecule
		// Hopefully, the "Original_Index" property carries through, and can be used as a quick shortcut.
		core::chemical::IndexIndexMapping rxn_map( core::chemical::rdkit::find_mapping( rdmol, prod, "Original_Index") ); // -1 is invalid

		if ( rxn_map.empty() ) {
			TR.Warning << "Cannot find substructure mapping of reactant to product." << std::endl;
			continue;
		}

		// Should be good. Now convert the residue into a Rosetta residue type.

		core::chemical::VDIndexMapping restype_prod_map( combine( to_converter.vd_to_index(), rxn_map ) );

		core::chemical::rdkit::RDMolToRestype from_converter(*prod);
		from_converter.set_nbr( restype_prod_map[rsdtype.nbr_vertex()] );

		core::chemical::MutableResidueTypeOP new_resop( from_converter.generate_restype(rsdtype,restype_prod_map) );

		TR << "Added " << new_resop->nheavyatoms() - rsdtype.nheavyatoms() << " heavyatoms to " << rsdtype.name() << std::endl;

		mapping_ = combine( restype_prod_map, from_converter.index_to_vd() );
		rsdtype = *new_resop;
		mapping_ = combine( mapping_, combine( core::chemical::VDStringMapping(*new_resop), core::chemical::StringVDMapping(rsdtype)) );

		// That's it, we're successful
		set_last_status( core::chemical::modifications::SUCCESS );
		return;
	}

	TR.Warning << "Tried " << rxns.size()*10 << " reactions, but none of them worked. Giving up." << std::endl;
	set_last_status( core::chemical::modifications::FAIL_RETRY );

}

core::chemical::VDVDMapping
ReactionGrow::get_mapping() const {
	return mapping_;
}

std::string
ReactionGrow::class_name() {
	return "ReactionGrow";
}

void
ReactionGrow::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute("fragments", xs_string,
		"The name of the SDF file which contains the other reactants (fragments to add).")
		+ XMLSchemaAttribute::attribute_w_default("weight_by_property", xs_string,
		"When randomly picking the other reactants, weight by the given property from the SDF", "")
		+ XMLSchemaAttribute::required_attribute("reactions", xs_string,
		"The name of the file containing the SMARTS-based reactions to use.");

	protocols::chemistries::xsd_type_definition_w_attributes(
		xsd, class_name(),
		"Grow a molecule by adding another molecule to the structure based on a given list of reactions.",
		attlist );
}

/// @brief Initialize any data members of this instance from an input tag
/// and a DataMap object
void
ReactionGrow::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &) {

	fragment_database( tag->getOption<std::string>("fragments") );
	reaction_file( tag->getOption<std::string>("reactions") );
	weight_by_property( tag->getOption<std::string>("weight_by_property", "") );

	// Reduce size of the fragments and reactions, so we avoid picking any that are bad.
	prefilter_reactions();
	prefilter_fragments();
}


}
}
