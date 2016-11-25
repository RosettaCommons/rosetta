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
#include <protocols/qsar/scoring_grid/GridManager.hh>
#include <protocols/qsar/scoring_grid/GridFactory.hh>
//#include <protocols/qsar/qsarTypeManager.hh>

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

static THREAD_LOCAL basic::Tracer TR( "protocols.qsar.scoring_grid.ScoringGridLoader" );

ScoringGridLoader::ScoringGridLoader() {}
ScoringGridLoader::~ScoringGridLoader() {}

void ScoringGridLoader::load_data(
	core::pose::Pose const &,
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &
) const
{
	using namespace utility::tag;
	typedef utility::vector0< TagCOP > TagCOPs;

	/// Setup the scoring grid_manager

	//core::Real width = 40.0;
	//core::Real resolution = 0.25;

	qsar::scoring_grid::GridManager* grid_manager(qsar::scoring_grid::GridManager::get_instance());

	if ( tag->hasOption("width") ) {
		grid_manager->set_width(tag->getOption<core::Real>("width"));

	}
	if ( tag->hasOption("resolution") ) {
		grid_manager->set_resolution(tag->getOption<core::Real>("resolution"));
	}

	if ( tag->hasOption("ligand_chain") ) {
		grid_manager->set_chain(tag->getOption<char>("ligand_chain"));
	}

	if ( tag->hasOption("normalize_mode") ) {
		std::string normalize_mode = tag->getOption<std::string>("normalize_mode");
		grid_manager->set_normalization_function(normalize_mode);
	}

	/// Add grids to the scoring grid manager

	TagCOPs const grid_tags( tag->getTags() );
	if ( grid_tags.size()==0 ) {
		TR <<"WARNING WARNING grid manager will be empty" <<std::endl;
	}

	for ( TagCOP tag : grid_tags ) {
		grid_manager->make_new_grid(tag);
	}

	TR.flush();
}

protocols::jd2::parser::DataLoaderOP
ScoringGridLoaderCreator::create_loader() const { return protocols::jd2::parser::DataLoaderOP( new ScoringGridLoader ); }

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
	attributes + XMLSchemaAttribute( "width", xsct_real, "The reach of the grid, in Angstroms" )
		+ XMLSchemaAttribute( "resolution", xsct_real, "The distance between grid points, in Angstroms" )
		+ XMLSchemaAttribute( "ligand_chain", xsct_char, "The chain in the input Pose that the ligand will be" )
		+ XMLSchemaAttribute( "normalize_mode", "grid_score_normalization", "The normalization function to use to"
		" aid the comparison of scores between different ligands." );

	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.element_name( loader_name() )
		.complex_type_naming_func( & ScoringGridLoader::scoring_grid_loader_ct_namer )
		.description( "The ScoringGridLoader will populate the singleton (?!!?!?!) GridManager with the set of"
		" ScoringGrids that are given as subelements of the " + loader_name() + " element. Instead, it ought to"
		" create an object as a collection of Grids and load that object into the datamap for other Movers and"
		" filters to retrieve. If someone starts using ScoringGrids again, please make that refactor! Until that"
		" happens, any protocol that relies on ScoringGrids cannot be run in a multi-threaded environment." )
		.set_subelements_repeatable( subelements )
		.add_attributes( attributes )
		.write_complex_type_to_schema( xsd );
}



} //namespace scoring_grid
} //namespace qsar
} //namespace protocols
