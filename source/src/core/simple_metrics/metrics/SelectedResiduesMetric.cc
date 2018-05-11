// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/metrics/SelectedResiduesMetric.cc
/// @brief Output residue-selected residues to a score file as rosetta resnums or pdbnums.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Unit headers
#include <core/simple_metrics/metrics/SelectedResiduesMetric.hh>
#include <core/simple_metrics/metrics/SelectedResiduesMetricCreator.hh>

// Core headers
#include <core/simple_metrics/StringMetric.hh>
#include <core/simple_metrics/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>


static basic::Tracer TR( "core.simple_metrics.metrics.SelectedResiduesMetric" );


namespace core {
namespace simple_metrics {
namespace metrics {

using namespace core::select::residue_selector;

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
SelectedResiduesMetric::SelectedResiduesMetric():
	core::simple_metrics::StringMetric()
{}

SelectedResiduesMetric::SelectedResiduesMetric( select::residue_selector::ResidueSelectorCOP selector ):
	core::simple_metrics::StringMetric()
{
	set_residue_selector( selector );
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
SelectedResiduesMetric::~SelectedResiduesMetric(){}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
SelectedResiduesMetric::SelectedResiduesMetric( SelectedResiduesMetric const & src ):
	core::simple_metrics::StringMetric( src ),
	rosetta_nums_(src.rosetta_nums_)
{
	selector_ = src.selector_;
}

core::simple_metrics::SimpleMetricOP
SelectedResiduesMetric::clone() const {
	return core::simple_metrics::SimpleMetricOP(new SelectedResiduesMetric( *this ) );

}

std::string
SelectedResiduesMetric::name() const {
	return name_static();
}

std::string
SelectedResiduesMetric::name_static() {
	return "SelectedResiduesMetric";

}
std::string
SelectedResiduesMetric::metric() const {
	return "selection";
}

void
SelectedResiduesMetric::set_residue_selector( select::residue_selector::ResidueSelectorCOP selector ){
	selector_ = selector;
}

void
SelectedResiduesMetric::set_output_in_rosetta_num(bool output_in_rosetta_num ){
	rosetta_nums_ = output_in_rosetta_num;
}

void
SelectedResiduesMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &  datamap)
{

	if ( tag->hasOption("residue_selector") ) {
		set_residue_selector(select::residue_selector::parse_residue_selector( tag, datamap ));
	} else {
		utility_exit_with_message("SelectedResiduesMetric: Residue Selector is required.");
	}

	set_output_in_rosetta_num( tag->getOption< bool >("rosetta_numbering", false));
}

void
SelectedResiduesMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default("rosetta_numbering", xsct_rosetta_bool, "Set to output in Rosetta numbering instead of PDB numbering", "false");

	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "residue_selector",
		"Required.  Output those residues selected. " );

	core::simple_metrics::xsd_simple_metric_type_definition_w_attributes(xsd, name_static(),
		"A utility metric that outputs the residue selection in Pose or PDB numbering.  Comma-Separated.", attlist);
}

std::string
SelectedResiduesMetric::calculate(const pose::Pose & pose) const {

	if ( ! selector_ ) {
		utility_exit_with_message("SelectedResiduesMetric: Residue Selector required!");
	}

	utility::vector1< bool > subset = selector_->apply( pose );
	utility::vector1< core::Size > selected = core::select::get_residues_from_subset( subset );

	std::string output = "";

	//Explicit or Fall-back.
	if ( ( rosetta_nums_ ) || ( pose.pdb_info() == nullptr ) || ( pose.pdb_info()->obsolete() ) ) {

		for ( core::Size i = 1; i <= selected.size(); ++i ) {
			if ( i == selected.size() ) {
				output += utility::to_string(selected[i]);
			} else {
				output += utility::to_string(selected[i])+",";
			}

		}
	} else {
		for ( core::Size i = 1; i <= selected.size(); ++i ) {
			core::Size resnum = selected[i];
			std::string chain = utility::to_string( pose.pdb_info()->chain(resnum));
			std::string num = utility::to_string( pose.pdb_info()->number(resnum));
			std::string icode = utility::to_string( pose.pdb_info()->icode(resnum));

			if ( icode == utility::to_string(' ') ) {
				output += num+chain;
			} else {
				output += num+chain+":"+icode;
			}

			if ( i != selected.size() ) {
				output += ",";
			}

		}
	}
	return output;
}

void
SelectedResiduesMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	SelectedResiduesMetric::provide_xml_schema( xsd );
}

std::string
SelectedResiduesMetricCreator::keyname() const {
	return SelectedResiduesMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
SelectedResiduesMetricCreator::create_simple_metric() const {
	return core::simple_metrics::SimpleMetricOP( new SelectedResiduesMetric );

}

} //core
} //simple_metrics
} //metrics






