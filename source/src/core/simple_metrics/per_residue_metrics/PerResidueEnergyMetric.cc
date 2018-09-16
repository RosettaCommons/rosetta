// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/per_residue_metrics/PerResidueEnergyMetric.cc
/// @brief A per-residue metric that will calculate/output per residue total energies or a specific score component.  Correctly decomposes energies to per-residue.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Unit headers
#include <core/simple_metrics/per_residue_metrics/PerResidueEnergyMetric.hh>
#include <core/simple_metrics/simple_metric_creators.hh>

// Core headers
#include <core/simple_metrics/PerResidueRealMetric.hh>
#include <core/simple_metrics/util.hh>
#include <core/pose/ref_pose.hh>
#include <core/pose/Pose.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/xml_util.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/Energies.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tag/util.hh>

static basic::Tracer TR( "core.simple_metrics.per_residue_metrics.PerResidueEnergyMetric" );


namespace core {
namespace simple_metrics {
namespace per_residue_metrics {

using namespace core::select;
using namespace core::select::residue_selector;
using namespace core::pose;
using namespace core::scoring;

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
PerResidueEnergyMetric::PerResidueEnergyMetric():
	core::simple_metrics::PerResidueRealMetric()
{}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
PerResidueEnergyMetric::~PerResidueEnergyMetric(){}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
PerResidueEnergyMetric::PerResidueEnergyMetric( PerResidueEnergyMetric const &  ) = default;

core::simple_metrics::SimpleMetricOP
PerResidueEnergyMetric::clone() const {
	return core::simple_metrics::SimpleMetricOP(new PerResidueEnergyMetric( *this ) );

}

std::string
PerResidueEnergyMetric::name() const {
	return name_static();
}

std::string
PerResidueEnergyMetric::name_static() {
	return "PerResidueEnergyMetric";

}
std::string
PerResidueEnergyMetric::metric() const {
	return "res_energy";
}

void
PerResidueEnergyMetric::set_scorefunction(core::scoring::ScoreFunctionCOP scorefxn){
	scorefxn_ = scorefxn;
}

void
PerResidueEnergyMetric::set_comparison_pose(core::pose::PoseCOP pose){
	ref_pose_ = pose;
}

void
PerResidueEnergyMetric::set_scoretype( scoring::ScoreType scoretype ){
	scoretype_ = scoretype;
}

void
PerResidueEnergyMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap)
{

	SimpleMetric::parse_base_tag( tag );
	PerResidueRealMetric::parse_per_residue_tag( tag, datamap );

	set_scoretype(core::scoring::score_type_from_name( tag->getOption<std::string>("scoretype", "total_score")));

	if ( tag->hasOption("scorefxn") ) {
		set_scorefunction(parse_score_function( tag, datamap ));
	}

	//Comparison pose.
	if ( tag->hasOption("reference_name") ) {
		ref_pose_ = saved_reference_pose(tag, datamap, "reference_name");
		TR<<"Loaded reference pose: "<<tag->getOption< std::string >( "reference_name" )<<" with "<<ref_pose_->size()<<" residues"<<std::endl;
	} else if ( tag->getOption< bool > ("use_native", false) && datamap.has_resource("native_pose") ) {
		ref_pose_ = saved_native_pose(datamap);
	}

	if ( (tag->getOption<std::string>("scoretype", "total_score") != "total_score") && get_custom_type() == "" ) {
		TR.Warning << "Specific scoretype set but custom_type not set." << std::endl;
	}
}

void
PerResidueEnergyMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::scoring;
	using namespace core::select::residue_selector;

	AttributeList attlist;
	attributes_for_saved_reference_pose( attlist );

	attlist + XMLSchemaAttribute::attribute_w_default(
		"use_native", xsct_rosetta_bool, "Use the native if present on the cmd-line.", "false"
	);

	utility::vector1< std::string > score_names;
	for ( core::Size ii = 1; ii <= n_score_types; ++ii ) {
		ScoreType score_type = static_cast<ScoreType>(ii);
		std::string score_name = name_from_score_type( score_type );
		score_names.push_back( score_name);
	}

	utility::tag::add_schema_restrictions_for_strings( xsd, "scoretypes", score_names);

	attlist + XMLSchemaAttribute::attribute_w_default("scoretype", "scoretypes", "ScoreType to calculate.", "total_score");

	attributes_for_get_score_function_name( attlist );

	std::string description = "A per-residue metric that will calculate/output per residue total energies or a specific score component. WEIGHTED."
		"Correctly decomposes energies to per-residue.";

	core::simple_metrics::xsd_per_residue_real_metric_type_definition_w_attributes(xsd, name_static(),
		description, attlist);
}

std::map< core::Size, core::Real >
PerResidueEnergyMetric::calculate(const pose::Pose & pose) const {
	core::pose::Pose local_pose = pose;
	core::pose::PoseOP local_ref_pose;
	if ( ref_pose_ ) {
		local_ref_pose = core::pose::PoseOP( new core::pose::Pose( * ref_pose_));
	}
	ScoreFunctionOP local_scorefxn;
	if ( scorefxn_ == nullptr ) {
		local_scorefxn = get_score_function();
	} else {
		local_scorefxn = scorefxn_->clone();
	}

	std::map< core::Size, core::Real > rsd_sub_energies;
	std::map< core::Size, core::Real > ref_sub_energies;

	EnergyMap weights( local_scorefxn->weights() );

	//Sum residue energies, making sure that we have pair-energies to do the summation properly
	core::scoring::methods::EnergyMethodOptionsOP emopts( new core::scoring::methods::EnergyMethodOptions( local_scorefxn->energy_method_options() ) );
	emopts->hbond_options().decompose_bb_hb_into_pair_energies( true );
	local_scorefxn->set_energy_method_options( *emopts );

	utility::vector1< bool > subset = get_selector()->apply( local_pose );
	local_scorefxn->score( local_pose );

	//Score refpose before summing
	if ( ref_pose_ ) {
		local_scorefxn->score(*local_ref_pose);
	}

	//Get either per-residue or total weights (sum over all residues or some quick way.).
	for ( core::Size res : get_residues_from_subset( subset ) ) {

		EnergyMap rsd_energies(weights * local_pose.energies().residue_total_energies(res));
		EnergyMap ref_energies;
		if ( ref_pose_ ) {
			ref_energies = weights * local_ref_pose->energies().residue_total_energies(res);
		}

		if ( scoretype_ == total_score ) {
			core::Real rsd_total_energy = 0;
			core::Real ref_total_energy = 0;
			for ( core::Size ii = 1; ii <= n_score_types; ++ii ) {
				ScoreType score_type = static_cast<ScoreType>(ii);

				if ( local_scorefxn->has_nonzero_weight( score_type ) ) {
					core::Real const value( rsd_energies[ score_type ] );
					rsd_total_energy+=value;

					if ( ref_pose_ ) {
						core::Real const ref_value( ref_energies[ score_type ] );
						ref_total_energy+=ref_value;
					}// if ref_pose_
				}
			}

			rsd_sub_energies[ res ] = rsd_total_energy;
			ref_sub_energies[ res ] = ref_total_energy;

		} else { //total_score
			rsd_sub_energies[ res ] = rsd_energies[ scoretype_ ];

			if ( ref_pose_ ) {
				ref_sub_energies[ res ] = ref_energies[ scoretype_ ];
			}
		}

	}


	if ( ref_pose_ != nullptr ) {
		std::map< core::Size, core::Real > result;
		for ( auto res_e_pair : rsd_sub_energies ) {
			result[ res_e_pair.first] = res_e_pair.second - ref_sub_energies[ res_e_pair.first];
		}
		return result;
	} else {
		return rsd_sub_energies;
	}
}

void
PerResidueEnergyMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	PerResidueEnergyMetric::provide_xml_schema( xsd );
}

std::string
PerResidueEnergyMetricCreator::keyname() const {
	return PerResidueEnergyMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
PerResidueEnergyMetricCreator::create_simple_metric() const {
	return core::simple_metrics::SimpleMetricOP( new PerResidueEnergyMetric );

}

} //core
} //simple_metrics
} //per_residue_metrics






