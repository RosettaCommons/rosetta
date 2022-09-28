// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/drug_design/ApplyChemistryMover.cc
/// @brief Apply a given Chemistry modifier to the ResidueType at a given position, then replace the ResidueType at that position.
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Unit headers
#include <protocols/drug_design/ApplyChemistryMover.hh>
#include <protocols/drug_design/ApplyChemistryMoverCreator.hh>

// Core headers
#include <protocols/drug_design/util.hh>

#include <core/pose/Pose.hh>
#include <core/pose/selection.hh>
#include <core/conformation/Conformation.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <protocols/chemistries/Chemistry.hh>
#include <protocols/chemistries/ChemistryFactory.hh>
#include <core/chemical/AtomRefMapping.hh>
#include <core/chemical/CacheableResidueTypeSets.hh>
#include <core/chemical/PoseResidueTypeSet.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.drug_design.ApplyChemistryMover" );

namespace protocols {
namespace drug_design {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
ApplyChemistryMover::ApplyChemistryMover():
	protocols::moves::Mover( ApplyChemistryMover::mover_name() )
{}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
ApplyChemistryMover::ApplyChemistryMover( ApplyChemistryMover const & ) = default;

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
ApplyChemistryMover::~ApplyChemistryMover() = default;

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

void
ApplyChemistryMover::add_chemistry( protocols::chemistries::ChemistryOP setting ) {
	chemistries_.push_back( setting );
}

/// @brief Apply the mover
void
ApplyChemistryMover::apply( core::pose::Pose & pose ){
	using namespace core::chemical;

	if ( chemistries_.empty() ) {
		utility_exit_with_message( "ApplyChemistryMover called with no sub-chemistries defined!" );
		return;
	}

	utility::vector1<bool> residues;
	if ( residue_selector_ ) {
		residues = residue_selector_->apply(pose);
	} else if ( ! residues_.empty() ) {
		residues.resize( pose.size(), false );
		for ( core::Size ii: core::pose::get_resnum_list( residues_, pose ) ) {
			residues[ ii ] = true;
		}
	} else if ( residue_id_ >= 1 && residue_id_ <= pose.size() ) {
		residues.resize( pose.size(), false );
		residues[ residue_id_ ] = true;
	} else {
		utility_exit_with_message("In ApplyChemistryMover, no suitable residue designations found.");
	}

	// Calculate the new ResidueTypes
	std::map< core::chemical::ResidueTypeCOP, std::pair< core::chemical::MutableResidueTypeOP, core::chemical::IndexVDMapping > > type_mapping;
	for ( core::Size ii(1); ii <= pose.size(); ++ii ) {
		if ( ! residues[ ii ] ) { continue; }
		core::chemical::ResidueTypeCOP old_type( pose.residue_type_ptr(ii) );
		if ( type_mapping.count( old_type ) == 0 ) {
			// Only do the conversion once per residue type, to avoid wasted effort
			// (This mover isn't intended to generate multiple different types for repeated occurances.)

			core::chemical::MutableResidueTypeOP new_type( new core::chemical::MutableResidueType( *old_type ) );
			IndexVDMapping mapping( combine( IndexNameMapping(*old_type), NameVDMapping(*new_type) ) );
			for ( protocols::chemistries::ChemistryOP const & chemistry: chemistries_ ) {
				chemistry->apply( *new_type, pose ); // Allow context-sensitive chemistries
				mapping = combine( mapping, chemistry->get_mapping() );
			}
			type_mapping[ old_type ] = std::make_pair( new_type, mapping );
		}

	}

	// Stick the modified residues into the Pose's ResidueTypeSet.
	// We do this before substituting, so we can fiddle with names, if needed.
	core::chemical::CacheableResidueTypeSets restype_sets; // In case we're a mixed type case
	for ( auto const & outer_pair: type_mapping ) {
		core::chemical::MutableResidueTypeOP new_type( outer_pair.second.first );
		if ( ! restype_sets.has_res_type_set( new_type->mode() ) ) {
			restype_sets.set_res_type_set( pose.conformation().modifiable_residue_type_set_for_conf( new_type->mode() ) );
		}
		core::chemical::PoseResidueTypeSetOP type_set( restype_sets.get_res_type_set( new_type->mode() ) );
		runtime_assert( type_set != nullptr );
		std::string new_name( new_type->name() );
		if ( type_mapping.size() == 1 ) {
			new_name = new_name_;
		}
		while ( type_set->has_name( new_name ) ) {
			new_name += "_";
			new_name += tag_;
		}
		new_type->name( new_name );
		type_set->add_base_residue_type( new_type );
	}
	for ( int mode(1); mode <= core::chemical::TYPE_SET_MODES_LENGTH; ++mode ) {
		if ( restype_sets.has_res_type_set( core::chemical::TypeSetMode(mode) ) ) {
			pose.conformation().reset_residue_type_set_for_conf( restype_sets.get_res_type_set( core::chemical::TypeSetMode(mode) ) );
		}
	}

	// Do the substitution
	for ( core::Size ii(1); ii <= pose.size(); ++ii ) {
		if ( ! residues[ ii ] ) { continue; }
		// Do the substitution
		core::chemical::ResidueTypeCOP old_type( pose.residue_type_ptr(ii) );
		debug_assert( type_mapping.count( old_type ) != 0 );
		core::chemical::MutableResidueTypeCOP const & new_muttype( type_mapping[ old_type ].first );
		core::chemical::ResidueTypeCOP new_type( pose.residue_type_set_for_pose( new_muttype->mode() )->name_mapOP( new_muttype->name() ) );
		IndexIndexMapping const & mapping( type_mapping[ old_type ].second.downstream_combine( VDNameMapping(*new_muttype) ).downstream_combine( NameIndexMapping(*new_type) ) );

		place_new_restype(pose, ii, *new_type, mapping);
	}
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
ApplyChemistryMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
ApplyChemistryMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data
) {
	residue_selector_ = core::select::residue_selector::parse_residue_selector(tag, data);
	if ( residue_selector_ == nullptr ) {
		if ( !tag->hasOption("residue") ) {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "The ApplyChemistryMover must be given either a residue_selector or residue option." );
		}
		residues( tag->getOption<std::string>("residue") );
	}

	new_name( tag->getOption<std::string>("new_name", "") );
	this->tag( tag->getOption<std::string>("tag", this->tag()) );


	std::string const chemistry_name( tag->getOption< std::string >( "chemistry", "" ) );
	if ( !chemistry_name.empty() ) {
		if ( ! data.has( "chemistry", chemistry_name ) ) {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "Cannot find chemistry named " + chemistry_name + ". Was it defined in the <CHEMISTRY> tag?" );
		}
		add_chemistry( data.get_ptr< protocols::chemistries::Chemistry >( "chemistry", chemistry_name ) );
	}

	for ( auto const & subtag: tag->getTags() ) {
		protocols::chemistries::ChemistryOP chemistry = protocols::chemistries::ChemistryFactory::get_instance()->new_chemistry( subtag, data );
		add_chemistry( chemistry );
	}

}

////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
moves::MoverOP
ApplyChemistryMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new ApplyChemistryMover );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
ApplyChemistryMover::clone() const
{
	return protocols::moves::MoverOP( new ApplyChemistryMover( *this ) );
}

std::string ApplyChemistryMover::get_name() const {
	return mover_name();
}

std::string ApplyChemistryMover::mover_name() {
	return "ApplyChemistryMover";
}

void ApplyChemistryMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	attlist
		+ XMLSchemaAttribute( "chemistry", xs_string, "The name of the chemistry to use" )
		+ XMLSchemaAttribute( "residue", xs_string, "Which residues (as a comma separated list) to apply the chemistry to. Ignored if residue selector is set." )
		+ XMLSchemaAttribute( "new_name", xs_string, "What to call the new ResidueType. (Only if a single type, and the name isn't already chosen.)" )
		+ XMLSchemaAttribute::attribute_w_default( "tag", xs_string, "A suffix to use to distinguish residues which have the same name", "mod" );

	core::select::residue_selector::attributes_for_parse_residue_selector_default_option_name( attlist, "Which residues (as a residue selector) to apply the chemistry to.");

	XMLSchemaSimpleSubelementList ssl;
	ssl.add_group_subelement( & protocols::chemistries::ChemistryFactory::chemistry_xml_schema_group_name );

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements(
		xsd, mover_name(),
		"Apply the given chemistry to the given residues, replacing the residues in the pose with the new residue type.",
		attlist,
		ssl );
}

/////////////// Creator ///////////////

protocols::moves::MoverOP
ApplyChemistryMoverCreator::create_mover() const
{
	return protocols::moves::MoverOP( new ApplyChemistryMover );
}

std::string
ApplyChemistryMoverCreator::keyname() const
{
	return ApplyChemistryMover::mover_name();
}

void ApplyChemistryMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ApplyChemistryMover::provide_xml_schema( xsd );
}

////////////////////////////////////////////////////////////////////////////////
/// private methods ///
///////////////////////


std::ostream &
operator<<( std::ostream & os, ApplyChemistryMover const & mover )
{
	mover.show(os);
	return os;
}

} //protocols
} //drug_design
