// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/drug_design/ReactionMultiTransform.hh
/// @brief apply RDKit's reaction-based chemistry to transform a ResidueType
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <protocols/drug_design/ReactionMultiTransform.hh>
#include <protocols/drug_design/ReactionMultiTransformCreator.hh>
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
#include <numeric/random/random_permutation.hh>

#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/ChemReactions/Reaction.h>
#include <rdkit/GraphMol/ChemReactions/ReactionParser.h>
#include <rdkit/GraphMol/FileParsers/MolSupplier.h>
#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/SmilesParse/SmilesWrite.h> // For MolToSmiles

namespace protocols {
namespace drug_design {

static basic::Tracer TR("protocols.drug_design.ReactionMultiTransform");

//------------------------- Creator -----------------------------

protocols::chemistries::ChemistryOP
ReactionMultiTransformCreator::create_chemistry() const {
	return protocols::chemistries::ChemistryOP( new ReactionMultiTransform );
}

std::string
ReactionMultiTransformCreator::keyname() const {
	return ReactionMultiTransform::class_name();
}

void
ReactionMultiTransformCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	ReactionMultiTransform::provide_xml_schema( xsd );
}

//------------------------- Chemistry -----------------------------

ReactionMultiTransform::ReactionMultiTransform():
	ReactionChemistry(class_name()),
	last_product_(0)
{}

void
ReactionMultiTransform::add_reaction( ::RDKit::ChemicalReactionOP rxn, core::Real weight ) {
	if ( rxn->getNumProductTemplates() != 1 ) {
		TR.Error << "Input reaction should have one and only one listed product: " << ::RDKit::ChemicalReactionToRxnSmarts( *rxn ) << std::endl;
		return;
	}
	if ( rxn->getNumReactantTemplates() != 1 ) {
		TR.Error << "Input reaction should have one and only one listed reactant: " << ::RDKit::ChemicalReactionToRxnSmarts( *rxn ) << std::endl;
		return;
	}

	ReactionChemistry::add_reaction( rxn, weight );
}

void
ReactionMultiTransform::apply( core::chemical::MutableResidueType & rsdtype )
{

	core::chemical::rdkit::RestypeToRDMol to_converter(rsdtype);
	rdmol_ = to_converter.Mol();
	core::chemical::rdkit::label_with_index(*rdmol_, "Original_Index"); // (Re)Label with "Original_Index" to use quick mapping behavior.
	// Properties carry through the reaction, except for a "bug" with the explicitly mentioned atoms.
	input_map_ = to_converter.vd_to_index();

	::RDKit::RWMOL_SPTR rdmol( to_converter.Mol() ); // Convert
	::RDKit::MOL_SPTR_VECT reactants;
	reactants.push_back( rdmol );

	products_.clear();
	last_product_ = 0;

	// Select Reactions
	utility::vector1< ::RDKit::ChemicalReactionOP > rxns;
	numeric::random::WeightedSampler rxn_sampler;
	filter_reactions(*rdmol, rxns, rxn_sampler);

	// Try to apply each reaction, and assemble the different products
	for ( core::Size rxn_num(1); rxn_num <= rxns.size(); ++rxn_num ) {
		::RDKit::ChemicalReactionOP rxn( rxns[rxn_num] );
		std::vector< ::RDKit::MOL_SPTR_VECT > products( rxn->runReactants( reactants ) );
		if ( products.size() == 0 ) { continue; }
		for ( core::Size prod_num(0); prod_num < products.size(); ++ prod_num ) {
			if ( products[ prod_num].size() != 1 ) { // Shouldn't happen
				TR.Error << "Bad reaction in ReactionMultiTransform! Too many products (" << products[ prod_num].size() << "): "
					<< ::RDKit::ChemicalReactionToRxnSmarts( *rxn ) << std::endl;
				continue;
			}

			::RDKit::RWMolOP prod( new ::RDKit::RWMol( *products[ prod_num ][0]) );

			if ( cleanup_product( *prod ) ) {
				// Cleanup has more extensive error information
				TR.Debug << "Failed cleanup for product of " << ::RDKit::ChemicalReactionToRxnSmarts( *rxn ) << std::endl;
				continue;
			}

			products_.push_back( prod );
		}
	}

	if ( products_.size() == 0 ) {
		mapping_.clear();
		// Not having any reactions isn't actually an issue for this chemistry.
		TR << "[NOTICE]: No products found for reaction set. - Returning as no-op" << std::endl;
		set_last_status( core::chemical::modifications::FAIL_DO_NOT_RETRY );
		return;
	}

	// Randomize order in case we only want one reaction.
	numeric::random::random_permutation( products_ );
	last_product_ = 0;

	// For later conversion
	ref_restype_ = core::chemical::MutableResidueTypeOP( new core::chemical::MutableResidueType( rsdtype ) );

	// ( reuse the conversion code )
	core::chemical::MutableResidueTypeOP picked( get_additional_output() );

	rsdtype = *picked;

	//Need to massage mapping for reassignment
	core::chemical::VDStringMapping vd_name( *picked );
	core::chemical::StringVDMapping name_vd( rsdtype );

	mapping_ = combine( mapping_, combine( vd_name, name_vd ) );

	set_last_status( core::chemical::modifications::SUCCESS );

}

/// @brief Are there alternate ResidueTypes which are availible from the last time we called apply?
/// (That is, will get_addtional_output() return non-null?)
bool
ReactionMultiTransform::has_additional_output() const {
	if ( last_product_ >= products_.size() || products_.size() == 0 ) {
		return false;
	} else {
		return true;
	}
}

/// @brief Get additional generated ResidueTypes, if any.
core::chemical::MutableResidueTypeOP
ReactionMultiTransform::get_additional_output() {
	if ( ! has_additional_output() ) {
		if ( get_last_status() == core::chemical::modifications::SUCCESS ) {
			set_last_status( core::chemical::modifications::FAIL_RETRY );
		}
		return nullptr;
	}
	last_product_ += 1;
	debug_assert( last_product_ <= products_.size() );

	::RDKit::RWMolOP prod( products_[ last_product_ ] );

	// We need to find a mapping from the pre-reaction molecule to the post-reaction molecule
	// Hopefully, the "Original_Index" property carries through, and can be used as a quick shortcut.
	core::chemical::IndexIndexMapping rxn_map( core::chemical::rdkit::find_mapping( rdmol_, prod, "Original_Index") ); // -1 is invalid

	// Should be good. Now convert the residue into a Rosetta residue type.
	core::chemical::VDIndexMapping restype_prod_map( combine( input_map_, rxn_map ) );

	core::chemical::rdkit::RDMolToRestype from_converter(*prod);
	from_converter.set_nbr( restype_prod_map[ref_restype_->nbr_vertex()] );

	core::chemical::MutableResidueTypeOP new_resop( from_converter.generate_restype(*ref_restype_,restype_prod_map) );

	mapping_ = combine( restype_prod_map, from_converter.index_to_vd() );

	set_last_status( core::chemical::modifications::SUCCESS );
	return new_resop;
}


core::chemical::VDVDMapping
ReactionMultiTransform::get_mapping() const {
	return mapping_;
}

std::string
ReactionMultiTransform::class_name() {
	return "ReactionMultiTransform";
}

void
ReactionMultiTransform::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute("reactions", xs_string,
		"The name of the file containing the SMARTS-based reactions to use.");

	protocols::chemistries::xsd_type_definition_w_attributes(
		xsd, class_name(),
		"Apply all of the reactions listed in the file, and randomly pick one of the products.",
		attlist );
}

/// @brief Initialize any data members of this instance from an input tag
/// and a DataMap object
void
ReactionMultiTransform::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &)
{
	reaction_file( tag->getOption<std::string>("reactions") );
}


}
}
