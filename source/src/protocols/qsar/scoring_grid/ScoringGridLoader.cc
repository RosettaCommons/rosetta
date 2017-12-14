// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/qsar/scoring_grid/DataLoader.cc
/// @brief  Implementation of the XML parser's DataLoader base class (ctor & dstor)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <protocols/qsar/scoring_grid/ScoringGridLoader.hh>
#include <protocols/qsar/scoring_grid/ScoringGridLoaderCreator.hh>

// Project Headers
#include <protocols/qsar/scoring_grid/GridSet.hh>
#include <protocols/qsar/scoring_grid/GridFactory.hh>
//#include <protocols/qsar/qsarTypeManager.hh>

#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>

// Boost Headers

namespace protocols {
namespace qsar {
namespace scoring_grid {

static basic::Tracer TR( "protocols.qsar.scoring_grid.ScoringGridLoader" );

ScoringGridLoader::ScoringGridLoader() = default;
ScoringGridLoader::~ScoringGridLoader() = default;

void ScoringGridLoader::load_data(
	core::pose::Pose const &,
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data
) const
{
	using namespace utility::tag;
	using TagCOPs = utility::vector0<TagCOP>;

	TagCOPs const sub_tags( tag->getTags() );

	TagCOPs gridset_tags;
	TagCOPs default_gridset_tags;

	for ( TagCOP subtag : sub_tags ) {
		if ( subtag->getName() == "GridSet" ) {
			// Note that we'll never actually get here, as I haven't actually enabled this in the XSD
			// For now, just repeat the SCORINGGRIDS tags (with the `name` attribute) if you need multiple Scoring Grids.
			gridset_tags.push_back( subtag );
		} else {
			default_gridset_tags.push_back( subtag );
		}
	}

	if ( !gridset_tags.empty() && !default_gridset_tags.empty() ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "Cannot use 'GridSet' tags simultaneously with bare non-GridSet tags!");
	}

	if ( gridset_tags.empty() ) {
		parse_gridset_tag( tag, tag, data );
	} else {
		// Note that we'll never actually get here, as I haven't actually enabled this in the XSD
		// For now, just repeat the SCORINGGRIDS tags (with the `name` attribute) if you need multiple Scoring Grids.
		for ( TagCOP gridset_tag : gridset_tags ) {
			parse_gridset_tag( gridset_tag, tag, data );
		}
	}

}

void ScoringGridLoader::parse_gridset_tag( utility::tag::TagCOP tag, utility::tag::TagCOP parent, basic::datacache::DataMap & data ) const
{
	using namespace utility::tag;
	using TagCOPs = utility::vector0<TagCOP>;

	// Initialize the prototype grid set
	qsar::scoring_grid::GridSetOP grid_set( new qsar::scoring_grid::GridSet );

	std::string name = "default";
	if ( tag->hasOption("name") ) {
		name = tag->getOption<std::string>("name");
	}

	// We defer to parent, such that we can set things in the SCORINGGRIDS tag that applies to all sub tags
	if ( tag->hasOption("width") ) {
		grid_set->width(tag->getOption<core::Real>("width"));
	} else if ( parent->hasOption("width") ) {
		grid_set->width(parent->getOption<core::Real>("width"));
	}

	if ( tag->hasOption("resolution") ) {
		grid_set->resolution(tag->getOption<core::Real>("resolution"));
	} else if ( parent->hasOption("resolution") ) {
		grid_set->resolution(parent->getOption<core::Real>("resolution"));
	}

	if ( tag->hasOption("ligand_chain") ) {
		grid_set->chain(tag->getOption<char>("ligand_chain"));
	} else if ( parent->hasOption("ligand_chain") ) {
		grid_set->chain(parent->getOption<char>("ligand_chain"));
	}

	if ( tag->hasOption("normalize_mode") ) {
		std::string normalize_mode = tag->getOption<std::string>("normalize_mode");
		grid_set->set_normalization_function(normalize_mode);
	} else if ( parent->hasOption("normalize_mode") ) {
		std::string normalize_mode = parent->getOption<std::string>("normalize_mode");
		grid_set->set_normalization_function(normalize_mode);
	}

	// Add grids to the GridSet

	TagCOPs const grid_tags( tag->getTags() );
	if ( grid_tags.size()==0 ) {
		TR.Warning << "No Scoring grids specified! Scoring Grids for the " << name << " GridSet will be empty!" << std::endl;
	}

	for ( TagCOP tag : grid_tags ) {
		grid_set->make_new_grid(tag);
	}

	// Put the grid in the "default" slot for the DataMap.
	bool data_add_status = data.add( "scoring_grids", name, grid_set );

	if ( !data_add_status ) {
		utility_exit_with_message( "Scoring Grid " + name + " already exists in the DataMap. Please rename or combine into a single entry." );
	}
}

protocols::parser::DataLoaderOP
ScoringGridLoaderCreator::create_loader() const { return protocols::parser::DataLoaderOP( new ScoringGridLoader ); }

std::string
ScoringGridLoaderCreator::keyname() const { return ScoringGridLoader::loader_name(); }

ScoringGridLoaderCreator::DerivedNameFunction
ScoringGridLoaderCreator::schema_ct_naming_function() const
{
	return & ScoringGridLoader::scoring_grid_loader_ct_namer;
}

void ScoringGridLoaderCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ScoringGridLoader::provide_xml_schema( xsd );
}


std::string
ScoringGridLoader::loader_name()
{
	return "SCORINGGRIDS";
}

std::string
ScoringGridLoader::scoring_grid_loader_ct_namer( std::string const & element_name )
{
	return "scoring_grid_loader_" + element_name + "_type";
}

void ScoringGridLoader::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	GridFactory::get_instance()->provide_xml_schema( xsd );
	XMLSchemaSimpleSubelementList subelements;
	subelements.add_group_subelement( & GridFactory::scoring_grid_xml_schema_group_name );

	XMLSchemaRestriction grid_score_norm_enumeration;
	grid_score_norm_enumeration.name( "grid_score_normalization" );
	grid_score_norm_enumeration.base_type( xs_string );
	grid_score_norm_enumeration.add_restriction( xsr_enumeration, "HeavyAtomNormalization" );
	grid_score_norm_enumeration.add_restriction( xsr_enumeration, "AllAtomNormalization" );
	grid_score_norm_enumeration.add_restriction( xsr_enumeration, "ChiAngleNormalization" );
	grid_score_norm_enumeration.add_restriction( xsr_enumeration, "MolecularWeightNormalization" );
	xsd.add_top_level_element( grid_score_norm_enumeration );

	AttributeList attributes;
	attributes + XMLSchemaAttribute::attribute_w_default( "name", xs_string, "The name of the GridSet for the listed set of grids", "default" )
		+ XMLSchemaAttribute( "width", xsct_real, "The reach of the grid, in Angstroms" )
		+ XMLSchemaAttribute( "resolution", xsct_real, "The distance between grid points, in Angstroms" )
		+ XMLSchemaAttribute( "ligand_chain", xsct_char, "The chain in the input Pose that the ligand will be" )
		+ XMLSchemaAttribute( "normalize_mode", "grid_score_normalization", "The normalization function to use to"
		" aid the comparison of scores between different ligands." );

	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.element_name( loader_name() )
		.complex_type_naming_func( & ScoringGridLoader::scoring_grid_loader_ct_namer )
		.description( "The ScoringGridLoader will populate a GridSet with the set of"
		" ScoringGrids that are given as subelements of the " + loader_name() + " element."
		" By default, the grids with be loaded into the 'default' GridSet." )
		.set_subelements_repeatable( subelements )
		.add_attributes( attributes )
		.write_complex_type_to_schema( xsd );
}



} //namespace scoring_grid
} //namespace qsar
} //namespace protocols
