// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/metrics/CompositeEnergyMetric.cc
/// @brief A metric to report/calculate all of the energies of a scorefunction that are not 0 or the delta of each energy between another input pose or the set native.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Unit headers
#include <protocols/analysis/simple_metrics/CompositeEnergyMetric.hh>
#include <protocols/analysis/simple_metrics/CompositeEnergyMetricCreator.hh>

// Core headers
#include <core/simple_metrics/CompositeRealMetric.hh>
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
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>


// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>


static basic::Tracer TR( "core.simple_metrics.metrics.CompositeEnergyMetric" );


namespace protocols {
namespace analysis {
namespace simple_metrics {

using namespace core::scoring;
using namespace core::pose;
using namespace core::select;
using namespace core::select::residue_selector;
using namespace core::simple_metrics;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace protocols::rosetta_scripts;

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
CompositeEnergyMetric::CompositeEnergyMetric():
	core::simple_metrics::CompositeRealMetric()
{
}

CompositeEnergyMetric::CompositeEnergyMetric( ResidueSelectorCOP selector):
	core::simple_metrics::CompositeRealMetric()
{
	set_residue_selector( selector );
}

CompositeEnergyMetric::CompositeEnergyMetric( ResidueSelectorCOP selector, ScoreFunctionCOP scorefxn):
	core::simple_metrics::CompositeRealMetric()
{
	set_residue_selector( selector );
	set_scorefunction( scorefxn );
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
CompositeEnergyMetric::~CompositeEnergyMetric(){}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
CompositeEnergyMetric::CompositeEnergyMetric( CompositeEnergyMetric const & src ):
	core::simple_metrics::CompositeRealMetric( src )
{
	residue_selector_ = src.residue_selector_;
	scorefxn_ = src.scorefxn_;
	ref_pose_ = src.ref_pose_;
}

std::string
CompositeEnergyMetric::name() const {
	return name_static();
}

std::string
CompositeEnergyMetric::name_static() {
	return "CompositeEnergyMetric";

}
std::string
CompositeEnergyMetric::metric() const {
	return "composite_energy";
}

utility::vector1< std::string >
CompositeEnergyMetric::get_metric_names() const {
	ScoreFunctionOP local_scorefxn;
	if ( scorefxn_ == nullptr ) {
		local_scorefxn = get_score_function();
	} else {
		local_scorefxn = scorefxn_->clone();
	}
	utility::vector1< ScoreType > types = local_scorefxn->get_nonzero_weighted_scoretypes();

	utility::vector1< std::string > s_types;
	for ( auto score_type : types ) {
		s_types.push_back(name_from_score_type( score_type  ));
	}
	return s_types;
}

void
CompositeEnergyMetric::set_scorefunction(core::scoring::ScoreFunctionCOP scorefxn){
	scorefxn_ = scorefxn;
}

void
CompositeEnergyMetric::set_residue_selector(core::select::residue_selector::ResidueSelectorCOP residue_selector){
	residue_selector_ = residue_selector;
}

void
CompositeEnergyMetric::set_comparison_pose(const core::pose::Pose &pose){
	ref_pose_ = core::pose::PoseOP( new core::pose::Pose( pose ));
}

void
CompositeEnergyMetric::load_native_pose_as_reference(){
	TR << "Loading refpose as native " << std::endl;
	if ( option[ OptionKeys::in::file::native ].user() ) {
		ref_pose_ = core::import_pose::pose_from_file( option[ OptionKeys::in::file::native ].value() , core::import_pose::PDB_file);
	}
}

void
CompositeEnergyMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &  datamap)
{

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
CompositeEnergyMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;
	using namespace protocols::rosetta_scripts;

	AttributeList attlist;
	attributes_for_saved_reference_pose( attlist );

	attlist + XMLSchemaAttribute::attribute_w_default(
		"use_native", xsct_rosetta_bool, "Use the native if present on the cmd-line.", "false"
	);

	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "residue_selector",
		"If a residue selector is present, we calculate each energy component for only these residues." );

	attributes_for_get_score_function_name( attlist );

	std::string description = "A metric to report/calculate all of the energies of a scorefunction that are not 0 or the delta of each energy between another input pose or the set native.";

	core::simple_metrics::xsd_simple_metric_type_definition_w_attributes(xsd, name_static(),
		description, attlist);
}

std::map< std::string, core::Real >
CompositeEnergyMetric::calculate(const core::pose::Pose & ext_pose) const {

	core::pose::Pose local_pose = ext_pose;
	ScoreFunctionOP local_scorefxn;
	if ( scorefxn_ == nullptr ) {
		local_scorefxn = get_score_function();
	} else {
		local_scorefxn = scorefxn_->clone();
	}

	std::map< std::string, core::Real > energy_totals;
	std::map< std::string, core::Real > ref_totals;

	EnergyMap weights( local_scorefxn->weights() );

	if ( residue_selector_ ) {

		//Sum residue energies, making sure that we have pair-energies to do the summation properly
		core::scoring::methods::EnergyMethodOptionsOP emopts( new core::scoring::methods::EnergyMethodOptions( local_scorefxn->energy_method_options() ) );
		emopts->hbond_options().decompose_bb_hb_into_pair_energies( true );
		local_scorefxn->set_energy_method_options( *emopts );

		utility::vector1< bool > subset = residue_selector_->apply( local_pose );
		local_scorefxn->score( local_pose );

		//Score refpose before summing
		if ( ref_pose_ ) {
			local_scorefxn->score(*ref_pose_);
		}

		//Get either per-residue or total weights (sum over all residues or some quick way.).
		for ( core::Size res : get_residues_from_subset( subset ) ) {

			EnergyMap rsd_energies(weights * local_pose.energies().residue_total_energies(res));
			EnergyMap ref_energies;
			if ( ref_pose_ ) {
				ref_energies = weights * ref_pose_->energies().residue_total_energies(res);
			}

			for ( core::Size ii = 1; ii <= n_score_types; ++ii ) {
				ScoreType score_type = static_cast<ScoreType>(ii);

				if ( local_scorefxn->has_nonzero_weight( score_type ) ) {
					std::string keyname= name_from_score_type( score_type  );
					core::Real const value( rsd_energies[ score_type ] );
					if ( ! energy_totals.count( keyname ) ) {
						energy_totals[keyname] = value;
					} else {
						energy_totals[keyname] += value;
					}
					if ( ref_pose_ ) {
						core::Real const ref_value( ref_energies[ score_type ] );
						if ( ! ref_totals.count( keyname ) ) {
							ref_totals[keyname] = ref_value;
						} else {
							ref_totals[keyname] += ref_value;
						}
					}// if ref_pose_
				}
			}
		}
	} else { // if residue_selector_

		//No Residue Selector - no need to break down energies to component residues.
		local_scorefxn->score(local_pose);
		EnergyMap all_energies = ( weights * local_pose.energies().total_energies() );
		EnergyMap ref_energies;
		if ( ref_pose_ ) {
			ref_energies = ( weights * ref_pose_->energies().total_energies() );
		}
		for ( core::Size ii = 1; ii <= n_score_types; ++ii ) {
			ScoreType score_type = static_cast<ScoreType>(ii);

			if ( local_scorefxn->has_nonzero_weight( score_type ) ) {
				std::string keyname= name_from_score_type( score_type  );
				energy_totals[ keyname ] = all_energies[score_type ];
				if ( ref_pose_ ) {
					ref_totals[ keyname ] = ref_energies[ score_type ];
				}
			}
		}
	}

	if ( ref_pose_ != nullptr ) {
		std::map< std::string, core::Real > result;
		for ( auto e_map_pair : energy_totals ) {
			result[ e_map_pair.first] = e_map_pair.second - ref_totals[ e_map_pair.first];
		}
		return result;
	} else {
		return energy_totals;
	}
}

core::simple_metrics::SimpleMetricOP
CompositeEnergyMetric::clone() const {
	return core::simple_metrics::SimpleMetricOP(new CompositeEnergyMetric( *this ) );

}


void
CompositeEnergyMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	CompositeEnergyMetric::provide_xml_schema( xsd );
}

std::string
CompositeEnergyMetricCreator::keyname() const {
	return CompositeEnergyMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
CompositeEnergyMetricCreator::create_simple_metric() const {
	return core::simple_metrics::SimpleMetricOP( new CompositeEnergyMetric );

}

} //core
} //simple_metrics
} //metrics






