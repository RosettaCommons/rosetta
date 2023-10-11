// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/per_residue_metrics/LoadedProbabilitiesMetric.cc
/// @brief A class to load a probabilities weights file into a PerResidueProbabilitiesMetric
/// @author Moritz Ertelt (moritz.ertelt@gmail.com)

// Unit headers
#include <core/simple_metrics/per_residue_metrics/LoadedProbabilitiesMetric.hh>
#include <core/simple_metrics/per_residue_metrics/LoadedProbabilitiesMetricCreator.hh>

// Core headers
#include <core/simple_metrics/util.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <utility/pointer/memory.hh>
#include <utility/io/util.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/map.hpp>
#endif // SERIALIZATION

static basic::Tracer TR( "core.simple_metrics.per_residue_metrics.LoadedProbabilitiesMetric" );


namespace core {
namespace simple_metrics {
namespace per_residue_metrics {

using namespace core::select;
using namespace core::select::residue_selector;

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
LoadedProbabilitiesMetric::LoadedProbabilitiesMetric():
	core::simple_metrics::PerResidueProbabilitiesMetric()
{}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
LoadedProbabilitiesMetric::~LoadedProbabilitiesMetric(){}

core::simple_metrics::SimpleMetricOP
LoadedProbabilitiesMetric::clone() const {
	return utility::pointer::make_shared< LoadedProbabilitiesMetric >( *this );
}

std::string
LoadedProbabilitiesMetric::name() const {
	return name_static();
}

std::string
LoadedProbabilitiesMetric::name_static() {
	return "LoadedProbabilitiesMetric";

}
std::string
LoadedProbabilitiesMetric::metric() const {

	return "LoadedProbabilitiesMetric";
}

void
LoadedProbabilitiesMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data)
{

	SimpleMetric::parse_base_tag( tag );

	(void)data; // avoid un-used parameter error

	// get filename
	set_filename( tag->getOption< std::string >( "filename") );
	// load the probabilities from the file
	utility::vector1< std::string > lines = utility::io::get_lines_from_file_data(filename_ );
	set_probabilities_from_lines( lines );

}

void
LoadedProbabilitiesMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	AttributeList attlist;
	attlist + XMLSchemaAttribute::required_attribute( "filename", xs_string, "The name of the output weights file storing the probabilities." );

	core::simple_metrics::xsd_per_residue_probabilities_metric_type_definition_w_attributes(xsd, name_static(),
		"A class to load a probabilities weights file into a PerResidueProbabilitiesMetric. Does not involve the current pose for the calculate function, but stores the loaded values in the pose cache.", attlist);
}

std::map< core::Size, std::map< core::chemical::AA, core::Real >>
LoadedProbabilitiesMetric::calculate(const core::pose::Pose & pose) const {

	(void)pose; // keep compiler happy
	return probabilities_;
}

void
LoadedProbabilitiesMetric::set_filename( std::string const & filename ) {
	filename_ = filename;
}

void
LoadedProbabilitiesMetric::set_probabilities_from_lines( utility::vector1< std::string > const & lines ) {
	// Clear existing data
	probabilities_.clear();

	for ( std::string const& line : lines ) {
		// Skip comments and empty lines
		if ( line.empty() || line[0] == '#' ) continue;

		// Split the line into columns
		auto columns = utility::string_split(line);
		if ( columns.size() < 3 ) {
			utility_exit_with_message("Weights have to be specified in the following format: POSNUM RESIDUETYPE WEIGHT");
		}
		// Extract data from columns
		core::Size resi = utility::string2Size(columns[1]);
		if ( resi < 1 ) {
			utility_exit_with_message("Invalid value for POSNUM: " + columns[1] );
		}
		auto resn = columns[2];
		core::Real weight = utility::string2Real(columns[3]);

		// Check if weight is in a valid range
		if ( weight < 0.0 || weight > 1.0 ) {
			utility_exit_with_message("Probability value out of range [0,1]: " + std::to_string(weight));
		}
		// Convert the string residue name to AA enum
		core::chemical::AA aa = core::chemical::aa_from_one_or_three(resn);
		if ( aa == core::chemical::AA::aa_unk ) {
			utility_exit_with_message("Unknown amino acid type " + resn + " at position " + columns[1] + ". Did you misspell it?");
		}

		// Populate the map
		probabilities_[resi][aa] = weight;
	}
}



/// @brief This simple metric is unpublished.  It returns Moritz Ertelt as its author.
void
LoadedProbabilitiesMetric::provide_citation_info( basic::citation_manager::CitationCollectionList & citations ) const {
	citations.add(
		utility::pointer::make_shared< basic::citation_manager::UnpublishedModuleInfo >(
		"LoadedProbabilitiesMetric", basic::citation_manager::CitedModuleType::SimpleMetric,
		"Moritz Ertelt",
		"University of Leipzig",
		"moritz.ertelt@gmail.com",
		"Wrote the LoadedProbabilitiesMetric."
		)
	);
}

void
LoadedProbabilitiesMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	LoadedProbabilitiesMetric::provide_xml_schema( xsd );
}

std::string
LoadedProbabilitiesMetricCreator::keyname() const {
	return LoadedProbabilitiesMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
LoadedProbabilitiesMetricCreator::create_simple_metric() const {
	return utility::pointer::make_shared< LoadedProbabilitiesMetric >();
}

} //per_residue_metrics
} //simple_metrics
} //core


#ifdef    SERIALIZATION


template< class Archive >
void
core::simple_metrics::per_residue_metrics::LoadedProbabilitiesMetric::save( Archive & arc ) const {
	arc( cereal::base_class< core::simple_metrics::PerResidueProbabilitiesMetric>( this ) );
	arc( CEREAL_NVP( filename_ ) );
	arc( CEREAL_NVP( probabilities_ ) );

}

template< class Archive >
void
core::simple_metrics::per_residue_metrics::LoadedProbabilitiesMetric::load( Archive & arc ) {
	arc( cereal::base_class< core::simple_metrics::PerResidueProbabilitiesMetric >( this ) );
	arc( filename_ );
	arc( probabilities_ );
}

SAVE_AND_LOAD_SERIALIZABLE( core::simple_metrics::per_residue_metrics::LoadedProbabilitiesMetric );
CEREAL_REGISTER_TYPE( core::simple_metrics::per_residue_metrics::LoadedProbabilitiesMetric )

CEREAL_REGISTER_DYNAMIC_INIT( core_simple_metrics_per_residue_metrics_LoadedProbabilitiesMetric )
#endif // SERIALIZATION




