// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/composite_metrics/ProtocolSettingsMetric.cc
/// @brief This Metric reports options that have been set in the command line and splits script_vars.  Each option name is the type and the setting is the value in the map.  This is primarily aimed at benchmarking and record-keeping for large-scale rosetta runs or experiments.  It works with both the global and local OptionsCollection to enable its use in JD3.  It is analogous to the ProtocolFeatures reporter, with more options for xml-based variables.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Unit headers
#include <core/simple_metrics/composite_metrics/ProtocolSettingsMetric.hh>
#include <core/simple_metrics/simple_metric_creators.hh>

// Core headers
#include <core/simple_metrics/CompositeStringMetric.hh>
#include <core/simple_metrics/util.hh>
#include <core/simple_metrics/metrics/TimingProfileMetric.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>

// Basic/Utility headers
#include <basic/options/option.hh>
#include <utility/options/keys/OptionKeys.hh>
#include <utility/options/keys/OptionKey.hh>

#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <utility/options/OptionCollection.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>


static basic::Tracer TR( "core.simple_metrics.composite_metrics.ProtocolSettingsMetric" );


namespace core {
namespace simple_metrics {
namespace composite_metrics {

using namespace utility::options;

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
ProtocolSettingsMetric::ProtocolSettingsMetric(bool base_name_option_only, bool get_script_vars, bool get_user_options, bool skip_corrections):
	core::simple_metrics::CompositeStringMetric()
{
	parse_options( basic::options::option, base_name_option_only, get_script_vars, get_user_options, skip_corrections);
}

///@brief Parse the Local Options Collection
ProtocolSettingsMetric::ProtocolSettingsMetric(
	utility::options::OptionCollection const & options,
	bool base_name_option_only,
	bool get_script_vars ,
	bool get_user_options,
	bool skip_corrections  ):

	core::simple_metrics::CompositeStringMetric()
{
	parse_options( options, base_name_option_only, get_script_vars, get_user_options, skip_corrections);
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
ProtocolSettingsMetric::~ProtocolSettingsMetric(){}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
ProtocolSettingsMetric::ProtocolSettingsMetric( ProtocolSettingsMetric const & ) = default;

core::simple_metrics::SimpleMetricOP
ProtocolSettingsMetric::clone() const {
	return core::simple_metrics::SimpleMetricOP(new ProtocolSettingsMetric( *this ) );

}

std::string
ProtocolSettingsMetric::name() const {
	return name_static();
}

std::string
ProtocolSettingsMetric::name_static() {
	return "ProtocolSettingsMetric";

}
std::string
ProtocolSettingsMetric::metric() const {
	return "opt";
}

void
ProtocolSettingsMetric::set_only_report_these_options(const utility::vector1<std::string> &select_opts){
	limit_reporting_to_these_options_ = select_opts;
}

utility::vector1< std::string >
ProtocolSettingsMetric::get_metric_names() const {
	utility::vector1< std::string > opt_names;
	for ( auto & pair : options_values_ ) {
		if ( limit_reporting_to_these_options_.size() == 0 ) {
			opt_names.push_back(pair.first);
		} else if ( limit_reporting_to_these_options_.contains(pair.first) ) {
			opt_names.push_back(pair.first);
		} else {
			continue;
		}
	}
	return opt_names;
}

void
ProtocolSettingsMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap )
{
	SimpleMetric::parse_base_tag( tag );

	TR << tag->to_string() << std::endl;
	bool base_name_only = tag->getOption< bool >("base_name_only", true);
	bool get_user_options = tag->getOption< bool >("get_user_options", true);
	bool get_script_vars = tag->getOption< bool >("get_script_vars", true);
	bool skip_corrections = tag->getOption< bool >("skip_corrections", true);

	if ( tag->hasOption("limit_to_options") ) {
		std::string option_list_restriction = tag->getOption<std::string>("limit_to_options");
		utility::vector1< std::string > limit_to = utility::string_split(option_list_restriction, ',');
		set_only_report_these_options(limit_to);
	}
	if ( datamap.has_resource("options") ) {
		parse_options(*datamap.get_resource<OptionCollection const>("options"), base_name_only, get_user_options, get_script_vars, skip_corrections);
	} else {
		parse_options(basic::options::option, base_name_only, get_script_vars, get_user_options, skip_corrections);
	}


}


void
ProtocolSettingsMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default("base_name_only", xsct_rosetta_bool, "Use only the base option name instead of the whole path", "true");
	attlist + XMLSchemaAttribute::attribute_w_default("get_user_options", xsct_rosetta_bool, "Report all set cmd-line options", "true");
	attlist + XMLSchemaAttribute::attribute_w_default("get_script_vars", xsct_rosetta_bool, "Split script_vars and report", "true");
	attlist + XMLSchemaAttribute("limit_to_options", xsct_string_cslist, "Limit reporting to these options (comma-separated list).  Can be user-set or script_vars, but works with the get_user_options and get_script_vars options.");
	attlist + XMLSchemaAttribute::attribute_w_default("skip_corrections", xsct_rosetta_bool, "Skip ScoreFunction Corrections, which are set in-code at the beginning of a run. " ,"true");

	//attributes_for_parse_residue_selector( attlist, "residue_selector",
	// "Selector specifying residues." );

	std::string description = "This Metric reports options that have been set in the command line and splits script_vars."\
		"  Each option name is the type and the setting is the value in the map.  "\
		"This is primarily aimed at benchmarking and record-keeping for large-scale rosetta runs or experiments.\n"\
		"  It works with both the global and local OptionsCollection to enable its use in JD3.  \n" \
		"It is analogous to the ProtocolFeatures reporter, with more options for xml-based variables.";


	core::simple_metrics::xsd_simple_metric_type_definition_w_attributes(xsd, name_static(),
		description, attlist);
}

void
ProtocolSettingsMetric::split_script_vars(std::string const & script_vars_option_string, std::map< std::string, std::string> & options_values) const{
	using namespace utility;

	vector1< std::string > split1 = string_split(script_vars_option_string, ' ');

	for ( std::string const & s : split1 ) {
		vector1< std::string > opt_value = string_split( s,'=');
		if ( opt_value.size() != 2 ) {
			utility_exit_with_message("Script Var option is not long enough! "+s);
		}
		options_values[opt_value[1]] = opt_value[2];
	}
}


void
ProtocolSettingsMetric::parse_options(utility::options::OptionCollection const & options, bool base_name_option_only, bool get_script_vars, bool get_user_options, bool skip_corrections){

	//metrics::TimingProfileMetric timer = metrics::TimingProfileMetric();

	options_values_.clear();

	for ( auto i = OptionKeys::begin(), e = OptionKeys::end(); i != e; ++i ) {
		OptionKey const & key( *i );
		if ( options.has( key ) && options[ key ].user() ) { // Active option
			Option const & opt( options( key ) );
			std::string const opt_group( OptionCollection::prefix( opt.id() ) );
			std::string const option_name( OptionCollection::suffix( opt.id() ));

			std::string opt_string = utility::strip(opt.equals_string(), "=");

			if ( opt_string.empty() ) {
				opt_string = "true";
			}

			if ( skip_corrections && (utility::contains( opt_group, "corrections")) ) {
				TR << "Skipping "<< opt.id() << std::endl;
				continue;
			}
			if ( get_user_options ) {
				if ( base_name_option_only ) {
					options_values_[option_name] = opt_string;
				} else {
					options_values_[opt.id()] = opt_string;
				}
			}

			if ( get_script_vars ) {
				if ( option_name == "script_vars" ) {
					split_script_vars(opt_string, options_values_);
				}
			}
		}
	}
	//std::cout << "Timer: " << timer.calc_time() << std::endl;

}


std::map< std::string, std::string >
ProtocolSettingsMetric::calculate(const pose::Pose & ) const {

	std::map< std::string, std::string > final_options;
	if ( limit_reporting_to_these_options_.size() != 0 ) {
		for ( auto & pair: options_values_ ) {
			if ( limit_reporting_to_these_options_.contains(pair.first) ) {
				final_options[pair.first] = pair.second;
			}
		}
		return final_options;
	} else {
		return options_values_;
	}
}

void
ProtocolSettingsMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	ProtocolSettingsMetric::provide_xml_schema( xsd );
}

std::string
ProtocolSettingsMetricCreator::keyname() const {
	return ProtocolSettingsMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
ProtocolSettingsMetricCreator::create_simple_metric() const {
	return core::simple_metrics::SimpleMetricOP( new ProtocolSettingsMetric );

}

} //core
} //simple_metrics
} //composite_metrics






