// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/parser/DataLoader.cc
/// @brief  Implementation of the XML parser's DataLoader base class (ctor & dstor)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <protocols/parser/MonteCarloLoader.hh>
#include <protocols/parser/StandardLoaderCreators.hh>

// Project Headers
#include <core/scoring/ScoreFunction.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

#include <basic/datacache/DataMap.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace parser {

static basic::Tracer TR( "protocols.jd2.parser.MonteCarloLoader" );

MonteCarloLoader::MonteCarloLoader() = default;
MonteCarloLoader::~MonteCarloLoader() = default;

void MonteCarloLoader::load_data(
	core::pose::Pose const &,
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data
) const
{
	using namespace utility::tag;
	using TagCOPs = utility::vector0<TagCOP>;

	TagCOPs const montecarlo_tags( tag->getTags() );

	for ( TagCOP montecarlo_tag : montecarlo_tags ) {
		runtime_assert( montecarlo_tag->getName() == "MonteCarlo" );
		//std::string const mc_name( montecarlo_tag->getName() );
		std::string const mc_name( montecarlo_tag->getOption< std::string >( "name" ));
		auto const mctemp( montecarlo_tag->getOption< core::Real >( "temperature", 2.0 ));
		core::scoring::ScoreFunctionOP scorefxn =
			rosetta_scripts::parse_score_function( montecarlo_tag, data )->clone();

		protocols::moves::MonteCarloOP mc( new protocols::moves::MonteCarlo( *scorefxn, mctemp ) );
		// add more options for the MonteCarlo object here, e.g.
		// 1. autotemp / quenchtemp
		// 2. heat after cycles
		data.add( "montecarlos" , mc_name, mc );
	}
	TR.flush();
}

std::string
MonteCarloLoader::loader_name()
{
	return "MONTECARLOS";
}


std::string
MonteCarloLoader::monte_carlo_loader_ct_namer( std::string const & element_name )
{
	return "monte_carlo_loader_" + element_name + "_type";
}

void
MonteCarloLoader::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	// attributs for each MonteCarlo tag; so far, there's only a temperature and a scorefxn attribue.
	AttributeList mc_attributes;
	mc_attributes + required_name_attribute( "The name for the MonteCarlo object that will be used for lookup in the DataMap" )
		+ XMLSchemaAttribute::attribute_w_default( "temperature", xsct_real, "The temperature for the MonteCarlo object to operate at, interpretted in units of kT (e.g. 0.6 = room temperature)", "2.0" );
	rosetta_scripts::attributes_for_parse_score_function( mc_attributes );

	XMLSchemaSimpleSubelementList monte_carlo_subelements;
	monte_carlo_subelements.add_simple_subelement( "MonteCarlo", mc_attributes,
		"Each MonteCarlo object will be loaded into the data map so that it can be retrieved by Movers and Filters."
		" Note that a MonteCarlo object creates a deep copy of the ScoreFunction it is given, so changes to the ScoreFunction "
		" object the MonteCarlo object was initialized with after its creation will not be noticed by it." );

	XMLSchemaComplexTypeGenerator mc_loader_type;
	mc_loader_type.element_name( loader_name() )
		.complex_type_naming_func( & monte_carlo_loader_ct_namer )
		.description( "Data loader for MonteCarlo objects that can be created and then added to the"
		" data map so that they could be shared between multiple movers, if desired; superceded by"
		" the GenericMonteCarlo mover, and I believe, currently unused by any mover" )
		.set_subelements_repeatable( monte_carlo_subelements )
		.write_complex_type_to_schema( xsd );
}


DataLoaderOP
MonteCarloLoaderCreator::create_loader() const { return DataLoaderOP( new MonteCarloLoader ); }

std::string
MonteCarloLoaderCreator::keyname() const { return MonteCarloLoader::loader_name(); }

MonteCarloLoaderCreator::DerivedNameFunction
MonteCarloLoaderCreator::schema_ct_naming_function() const
{
	return & MonteCarloLoader::monte_carlo_loader_ct_namer;
}

void MonteCarloLoaderCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	MonteCarloLoader::provide_xml_schema( xsd );
}


} //namespace parser
} //namespace protocols
