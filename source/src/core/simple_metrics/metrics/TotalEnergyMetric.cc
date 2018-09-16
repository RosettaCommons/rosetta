// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/metrics/TotalEnergyMetric.cc
/// @brief A metric to report the total energy of the system or the delta total energy between another input pose or the set native.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Unit headers
#include <core/simple_metrics/metrics/TotalEnergyMetric.hh>
#include <core/simple_metrics/simple_metric_creators.hh>
#include <core/simple_metrics/per_residue_metrics/PerResidueEnergyMetric.hh>

// Core headers
#include <core/simple_metrics/RealMetric.hh>
#include <core/simple_metrics/util.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/util.hh>
#include <core/scoring/ScoreFunction.hh>

#include <core/scoring/xml_util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/ref_pose.hh>


// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>

#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>


static basic::Tracer TR( "core.simple_metrics.metrics.TotalEnergyMetric" );


namespace core {
namespace simple_metrics {
namespace metrics {

using namespace core::scoring;
using namespace core::pose;
using namespace core::select::residue_selector;
using namespace core::simple_metrics;

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
TotalEnergyMetric::TotalEnergyMetric():
	core::simple_metrics::RealMetric()
{

}

TotalEnergyMetric::TotalEnergyMetric( ResidueSelectorCOP selector ):
	core::simple_metrics::RealMetric()
{

	set_residue_selector( selector );
}

TotalEnergyMetric::TotalEnergyMetric( ResidueSelectorCOP selector, ScoreFunctionCOP scorefxn ):
	core::simple_metrics::RealMetric()
{

	set_residue_selector( selector );
	set_scorefunction( scorefxn );
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
TotalEnergyMetric::~TotalEnergyMetric(){}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
TotalEnergyMetric::TotalEnergyMetric( TotalEnergyMetric const & ) = default;


core::simple_metrics::SimpleMetricOP
TotalEnergyMetric::clone() const {
	return SimpleMetricOP(new TotalEnergyMetric( *this ) );

}

void
TotalEnergyMetric::set_scorefunction(core::scoring::ScoreFunctionCOP scorefxn){
	scorefxn_ = scorefxn;
}

void
TotalEnergyMetric::set_residue_selector(core::select::residue_selector::ResidueSelectorCOP residue_selector){
	residue_selector_ = residue_selector;
}

void
TotalEnergyMetric::set_comparison_pose(PoseCOP pose){
	ref_pose_ = pose;
}


std::string
TotalEnergyMetric::name() const {
	return name_static();
}

std::string
TotalEnergyMetric::name_static() {
	return "TotalEnergyMetric";

}
std::string
TotalEnergyMetric::metric() const {
	return "total_energy";
}

void
TotalEnergyMetric::set_scoretype( scoring::ScoreType scoretype ){
	scoretype_ = scoretype;
}

/*
void
TotalEnergyMetric::load_native_pose_as_reference(){
TR << "Loading native as reference pose " << std::endl;
if ( option[ OptionKeys::in::file::native ].user() ) {
ref_pose_ = core::import_pose::pose_from_file( option[ OptionKeys::in::file::native ].value() , core::import_pose::PDB_file);
}
}
*/


void
TotalEnergyMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap  )
{
	SimpleMetric::parse_base_tag( tag );
	if ( tag->hasOption("residue_selector") ) {
		set_residue_selector(parse_residue_selector( tag, datamap ));
	}

	if ( tag->hasOption("scorefxn") ) {
		set_scorefunction(parse_score_function( tag, datamap ));
	}
	set_scoretype(core::scoring::score_type_from_name( tag->getOption<std::string>("scoretype", "total_score")));

	//Comparison pose.
	if ( tag->hasOption("reference_name") ) {
		ref_pose_ = saved_reference_pose(tag, datamap, "reference_name");
		TR<<"Loaded reference pose: "<<tag->getOption< std::string >( "reference_name" )<<" with "<<ref_pose_->size()<<" residues"<<std::endl;
	} else if ( tag->getOption< bool > ("use_native", false) && datamap.has_resource("native_pose") ) {
		ref_pose_ = saved_native_pose(datamap);
	}


}

void
TotalEnergyMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	AttributeList attlist;

	attributes_for_saved_reference_pose( attlist );

	attlist + XMLSchemaAttribute::attribute_w_default(
		"use_native", xsct_rosetta_bool, "Use the native as reference for DELTA score if present on the cmd-line.", "false"
	);
	attlist + XMLSchemaAttribute::attribute_w_default("scoretype", "scoretypes", "ScoreType to calculate.", "total_score");

	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "residue_selector",
		"If a residue selector is present, we calculate the total energy of only these residues." );

	attributes_for_get_score_function_name( attlist );

	std::string description = "A metric to report the total energy of the system or the delta total energy between another input pose or the cmd-line set native.  ";

	core::simple_metrics::xsd_simple_metric_type_definition_w_attributes(xsd, name_static(),
		description, attlist);
}

core::Real
TotalEnergyMetric::calculate(const core::pose::Pose & pose) const {
	using namespace core::simple_metrics::per_residue_metrics;



	PerResidueEnergyMetric e_metric = PerResidueEnergyMetric();

	e_metric.set_scoretype(scoretype_);

	if ( residue_selector_ ) {
		e_metric.set_residue_selector(residue_selector_);
	}

	e_metric.set_scorefunction(scorefxn_);

	if ( ref_pose_ ) {
		e_metric.set_comparison_pose(ref_pose_);
	}

	std::map< core::Size, core::Real > energies = e_metric.calculate( pose );
	core::Real total_energy = 0;
	for ( auto res_e_pair : energies ) {
		total_energy+= res_e_pair.second;
	}
	return total_energy;
}

void
TotalEnergyMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	TotalEnergyMetric::provide_xml_schema( xsd );
}

std::string
TotalEnergyMetricCreator::keyname() const {
	return TotalEnergyMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
TotalEnergyMetricCreator::create_simple_metric() const {
	return core::simple_metrics::SimpleMetricOP( new TotalEnergyMetric );

}

} //core
} //simple_metrics
} //metrics






