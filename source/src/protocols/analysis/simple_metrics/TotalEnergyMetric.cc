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
#include <protocols/analysis/simple_metrics/TotalEnergyMetric.hh>
#include <protocols/analysis/simple_metrics/TotalEnergyMetricCreator.hh>

// Core headers
#include <core/simple_metrics/RealMetric.hh>
#include <core/simple_metrics/util.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/Energies.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>

// Protocol headers
#include <protocols/rosetta_scripts/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <utility/options/OptionCollection.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>


static basic::Tracer TR( "core.simple_metrics.metrics.TotalEnergyMetric" );


namespace protocols {
namespace analysis {
namespace simple_metrics {

using namespace core::scoring;
using namespace core::pose;
using namespace core::select::residue_selector;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::simple_metrics;
using namespace protocols::rosetta_scripts;

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
TotalEnergyMetric::TotalEnergyMetric( TotalEnergyMetric const & src ):
	RealMetric( src )
{
	scorefxn_ = src.scorefxn_;
	residue_selector_ = src.residue_selector_;
	ref_pose_ = src.ref_pose_;
}

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
TotalEnergyMetric::set_comparison_pose(const core::pose::Pose &pose){
	ref_pose_ = PoseOP( new Pose( pose ));
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
TotalEnergyMetric::load_native_pose_as_reference(){
	TR << "Loading native as reference pose " << std::endl;
	if ( option[ OptionKeys::in::file::native ].user() ) {
		ref_pose_ = core::import_pose::pose_from_file( option[ OptionKeys::in::file::native ].value() , core::import_pose::PDB_file);
	}
}

void
TotalEnergyMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap  )
{
	SimpleMetric::parse_base_tag( tag );
	if ( tag->hasOption("residue_selector") ) {
		set_residue_selector(protocols::rosetta_scripts::parse_residue_selector( tag, datamap ));
	}

	if ( tag->hasOption("scorefxn") ) {
		set_scorefunction(protocols::rosetta_scripts::parse_score_function( tag, datamap ));
	}

	//Comparison pose.
	if ( tag->hasOption("reference_name") ) {
		ref_pose_ = protocols::rosetta_scripts::saved_reference_pose(tag, datamap, "reference_name");
		TR<<"Loaded reference pose: "<<tag->getOption< std::string >( "reference_name" )<<" with "<<ref_pose_->size()<<" residues"<<std::endl;
	} else if ( tag->getOption< bool >("use_native", false) ) {
		load_native_pose_as_reference();
	}


}

void
TotalEnergyMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;
	using namespace protocols::rosetta_scripts;

	AttributeList attlist;

	attributes_for_saved_reference_pose( attlist );

	attlist + XMLSchemaAttribute::attribute_w_default(
		"use_native", xsct_rosetta_bool, "Use the native as reference for DELTA score if present on the cmd-line.", "false"
	);

	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "residue_selector",
		"If a residue selector is present, we calculate the total energy of only these residues." );

	attributes_for_get_score_function_name( attlist );

	std::string description = "A metric to report the total energy of the system or the delta total energy between another input pose or the cmd-line set native.  ";

	core::simple_metrics::xsd_simple_metric_type_definition_w_attributes(xsd, name_static(),
		description, attlist);
}

core::Real
TotalEnergyMetric::calculate(const core::pose::Pose & pose) const {

	Pose local_pose = pose;
	ScoreFunctionOP local_scorefxn;
	if ( ! scorefxn_ ) {
		local_scorefxn = get_score_function();
	} else {
		local_scorefxn = scorefxn_->clone();
	}

	core::Real total = 0;
	core::Real ref_total = 0;

	if ( ref_pose_ ) {
		ref_total = local_scorefxn->score( *ref_pose_);
	}

	if ( residue_selector_ ) {

		//Sum residue energies, making sure that we have pair-energies to do the summation properly
		core::scoring::methods::EnergyMethodOptionsOP emopts( new core::scoring::methods::EnergyMethodOptions( local_scorefxn->energy_method_options() ) );
		emopts->hbond_options().decompose_bb_hb_into_pair_energies( true );
		local_scorefxn->set_energy_method_options( *emopts );

		utility::vector1< bool > subset = residue_selector_->apply( local_pose );

		local_scorefxn->score( local_pose );
		total = local_scorefxn->get_sub_score(local_pose, subset);

		if ( ref_pose_ ) {
			ref_total = local_scorefxn->get_sub_score(*ref_pose_, subset);
		}

	} else {
		total = local_scorefxn->score( local_pose );
	}


	if ( ref_pose_ ) {
		TR <<"Calculating energy DELTA between ref_pose" << std::endl;
		return total - ref_total;

	} else {
		return total;
	}
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






