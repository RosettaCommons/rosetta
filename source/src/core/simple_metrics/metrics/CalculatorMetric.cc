// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/metrics/CalculatorMetric.cc
/// @brief A metric which can combine other metrics in a (semi) abitrary calculation
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Unit headers
#include <core/simple_metrics/metrics/CalculatorMetric.hh>
#include <core/simple_metrics/metrics/CalculatorMetricCreator.hh>

// Core headers
#include <core/simple_metrics/RealMetric.hh>
#include <core/simple_metrics/util.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/simple_metrics/util.hh>

#include <core/pose/extra_pose_info_util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <utility/pointer/memory.hh>
#include <numeric/Calculator.hh>
#include <numeric/random/random.hh>

// XSD Includes
#include <core/simple_metrics/SimpleMetricFactory.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.hh>
#include <basic/citation_manager/CitationCollection.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/map.hpp>
#endif // SERIALIZATION

static basic::Tracer TR( "core.simple_metrics.metrics.CalculatorMetric" );


namespace core {
namespace simple_metrics {
namespace metrics {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
CalculatorMetric::CalculatorMetric():
	core::simple_metrics::RealMetric()
{}

CalculatorMetric::CalculatorMetric(std::string const & equation):
	core::simple_metrics::RealMetric(),
	calc_( utility::pointer::make_shared< numeric::Calculator >(equation) )
{}

////////////////////////////////////////////////////////////////////////////////

core::simple_metrics::SimpleMetricOP
CalculatorMetric::clone() const {
	return utility::pointer::make_shared< CalculatorMetric >( *this );
}

std::string
CalculatorMetric::name() const {
	return name_static();
}

std::string
CalculatorMetric::name_static() {
	return "CalculatorMetric";

}

void
CalculatorMetric::add_simple_metric( std::string const & name, core::simple_metrics::SimpleMetricCOP metric ) {
	if ( ! metric ) {
		utility_exit_with_message("Calculator simple metric can't use non-existant (null pointer) metric with name "+ name);
	}
	core::simple_metrics::RealMetricCOP rmetric = utility::pointer::dynamic_pointer_cast< core::simple_metrics::RealMetric const >(metric);
	if ( ! rmetric ) {
		utility_exit_with_message("Metric " + metric->name() + " being added under name " + name + " is not a Real-valued metric.");
	}
	metrics_[name] = rmetric;
}

void
CalculatorMetric::add_reported_value( std::string const & name, std::string const & report_key ) {
	reported_values_[name] = report_key;
}

void
CalculatorMetric::add_constant( std::string const & name, core::Real value ) {
	values_[name] = value;
}

bool
CalculatorMetric::check_equation() {
	std::map< std::string, core::Real > vars(values_);
	for ( auto & metric : metrics_ ) {
		vars[ metric.first ] = 1.0 + 0.00001 * numeric::random::uniform(); // Additional random to avoid "1/(alpha - beta)" type situations.
	}
	for ( auto & report : reported_values_ ) {
		vars[ report.first ] = 1.0 + 0.00001 * numeric::random::uniform(); // Additional random to avoid "1/(alpha - beta)" type situations.
	}
	numeric::Real dummy;
	if ( calc_->compute(vars, dummy) ) {
		return false; // Issue with doing the calculator.
	}
	return true; // We were successfull
}

std::string
CalculatorMetric::metric() const {
	return "calculator";
}

void
CalculatorMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data )
{
	SimpleMetric::parse_base_tag( tag );

	std::string equation = tag->getOption< std::string >( "equation" );
	calc_ = utility::pointer::make_shared< numeric::Calculator >( equation );

	for ( utility::tag::TagCOP sub_tag : tag->getTags() ) {
		std::string const & varname( sub_tag->getOption<std::string>( "name" ) );

		if ( sub_tag->getName() == "VAR" ||  sub_tag->getName() == "Var" ||  sub_tag->getName() == "var" ) {
			core::Size count = 0;
			if ( sub_tag->hasOption("reported") ) {
				add_reported_value( varname, sub_tag->getOption<std::string>( "reported" ) );
				++count;
			}
			if ( sub_tag->hasOption("value") ) {
				add_constant( varname, sub_tag->getOption<core::Real>( "value" ) );
				++count;
			}
			if ( sub_tag->hasOption("metric") ) {
				add_simple_metric( varname, core::simple_metrics::get_metric_from_datamap_and_subtags( sub_tag, data ) );
				++count;
			}
			if ( count != 1 ) {
				throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "CalculatorMetric VAR subtag must have one and only one of 'metric', 'value', or 'reported'.");
			}
		} else { // Should be a simple metric tag.
			add_simple_metric( varname, SimpleMetricFactory::get_instance()->new_simple_metric( sub_tag->getName(), sub_tag, data ) );
		}
	}

	if ( ! check_equation() ) {
		utility_exit_with_message("Bad equation in CalculatorMetric: " + equation);
	}
}

void
CalculatorMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;

	AttributeList subelement_attlist;
	subelement_attlist
		+ XMLSchemaAttribute::required_attribute( "name", xs_string, "Unique name to identify value for use in equation." )
		+ XMLSchemaAttribute( "value", xsct_real, "Specifiy a constant value." )
		+ XMLSchemaAttribute( "metric", xs_string, "Evaluate given filter at calculator evaluation time." )
		+ XMLSchemaAttribute( "reported", xs_string, "Retrieve reported value. See 'report_at_end=false' documentation in ParsedProtocol.");

	XMLSchemaSimpleSubelementList subelements;
	subelements
		.add_simple_subelement( "VAR", subelement_attlist, "Specify values to be available in equation." )
		.add_simple_subelement( "Var", subelement_attlist, "Specify values to be available in equation." )
		.add_simple_subelement( "var", subelement_attlist, "Specify values to be available in equation." );

	subelements.add_group_subelement( & SimpleMetricFactory::get_instance()->simple_metric_xml_schema_group_name );

	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute( "equation", xs_string, "Equation to evaluate filter value." );

	core::simple_metrics::xsd_simple_metric_type_definition_w_attributes_and_repeatable_subelements(xsd, name_static(),
		"Calculate a value based on an equation and the results of other simple metrics.", attlist, subelements);
}

core::Real
CalculatorMetric::calculate(const core::pose::Pose & pose ) const {
	debug_assert(calc_);

	// Assemble the current set of values from the various sources.
	std::map< std::string, core::Real > vars(values_);
	for ( auto & metric : metrics_ ) {
		debug_assert(metric.second);
		vars[ metric.first ] = (metric.second)->calculate( pose );
	}
	for ( auto & reported : reported_values_ ) {
		if ( !getPoseExtraScore(pose, reported.second, vars[reported.first]) ) {
			utility_exit_with_message("CalculatorFilter required reported value not yet present in pose.");
		}
	}

	numeric::Real value(999999);
	if ( calc_->compute(vars, value) ) {
		TR.Error << "Problem calculating equation in CalculatorFilter - resultant value likely garbage." << std::endl;
	}
	return value;
}

bool
CalculatorMetric::simple_metric_is_unpublished() const {
	return true;
}

utility::vector1< basic::citation_manager::UnpublishedModuleInfoCOP >
CalculatorMetric::provide_authorship_info_for_unpublished() const {
	utility::vector1< basic::citation_manager::UnpublishedModuleInfoCOP > retlist;
	return retlist;
}


void
CalculatorMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	CalculatorMetric::provide_xml_schema( xsd );
}

std::string
CalculatorMetricCreator::keyname() const {
	return CalculatorMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
CalculatorMetricCreator::create_simple_metric() const {
	return utility::pointer::make_shared< CalculatorMetric >();
}

} //metrics
} //simple_metrics
} //core


#ifdef    SERIALIZATION

template< class Archive >
void
core::simple_metrics::metrics::CalculatorMetric::save( Archive & arc ) const {
	arc( cereal::base_class< core::simple_metrics::RealMetric>( this ) );
	std::string const & equation = calc_ ? calc_->equation() : "";
	arc( CEREAL_NVP( equation ) ); // EXEMPT calc_
	arc( CEREAL_NVP( values_ ) );
	arc( CEREAL_NVP( metrics_ ) );
	arc( CEREAL_NVP( reported_values_ ) );
}

template< class Archive >
void
core::simple_metrics::metrics::CalculatorMetric::load( Archive & arc ) {
	arc( cereal::base_class< core::simple_metrics::RealMetric >( this ) );
	std::string equation;
	arc( equation ); // EXEMPT calc_
	if ( equation.empty() ) {
		calc_ = utility::pointer::make_shared< numeric::Calculator >(equation);
	} else {
		calc_ = nullptr;
	}
	arc( values_ );
	arc( metrics_ );
	arc( reported_values_ );
}

SAVE_AND_LOAD_SERIALIZABLE( core::simple_metrics::metrics::CalculatorMetric );
CEREAL_REGISTER_TYPE( core::simple_metrics::metrics::CalculatorMetric )

CEREAL_REGISTER_DYNAMIC_INIT( core_simple_metrics_metrics_CalculatorMetric )
#endif // SERIALIZATION




