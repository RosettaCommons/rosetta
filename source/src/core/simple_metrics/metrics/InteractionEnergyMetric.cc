// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/metrics/InteractionEnergyMetric.cc
/// @brief Calculate the interaction energy between residues in two ResidueSelectors in a single pose, including long-range interactions.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
/// @author Original Logic from InterfaceDeltaEnergetics: John Karanicolas && Roland A. Pache

// Unit headers
#include <core/simple_metrics/metrics/InteractionEnergyMetric.hh>
#include <core/simple_metrics/simple_metric_creators.hh>

// Core headers
#include <core/simple_metrics/RealMetric.hh>
#include <core/simple_metrics/util.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/util.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/methods/Methods.hh> //for long range energies
#include <core/scoring/LREnergyContainer.hh> //long range energies
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/xml_util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <utility/vector1.srlz.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

// C++ headers
#include <utility>

static basic::Tracer TR( "core.simple_metrics.metrics.InteractionEnergyMetric" );

using namespace core::scoring;
using namespace core::select;
using namespace core::select::residue_selector;

namespace core {
namespace simple_metrics {
namespace metrics {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
InteractionEnergyMetric::InteractionEnergyMetric():
	core::simple_metrics::RealMetric()
{
	scorefxn_ = get_score_function();
}

InteractionEnergyMetric::InteractionEnergyMetric( ResidueSelectorCOP selector1, ResidueSelectorCOP selector2):
	core::simple_metrics::RealMetric()
{
	scorefxn_ = get_score_function();
	set_residue_selectors(selector1, selector2);
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
InteractionEnergyMetric::~InteractionEnergyMetric(){}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
InteractionEnergyMetric::InteractionEnergyMetric( InteractionEnergyMetric const &  ) = default;

core::simple_metrics::SimpleMetricOP
InteractionEnergyMetric::clone() const {
	return core::simple_metrics::SimpleMetricOP(new InteractionEnergyMetric( *this ) );

}

std::string
InteractionEnergyMetric::name() const {
	return name_static();
}

std::string
InteractionEnergyMetric::name_static() {
	return "InteractionEnergyMetric";

}
std::string
InteractionEnergyMetric::metric() const {
	return "interaction_energy";
}

///@brief Set the residue selectors we will use to calculate the interaction energies.
void
InteractionEnergyMetric::set_residue_selectors(
	select::residue_selector::ResidueSelectorCOP selector1,
	select::residue_selector::ResidueSelectorCOP selector2)
{
	selector1_ = selector1;
	selector2_ = selector2;
}

///@brief Set a scorefunction.  Only used if the pose is not scored.  Not recommended, this is our failsafe!
void
InteractionEnergyMetric::set_scorefunction( scoring::ScoreFunctionCOP scorefxn ){
	scorefxn_ = scorefxn;
}


void
InteractionEnergyMetric::set_include_only_scoretypes(utility::vector1<scoring::ScoreType> const & include_only_scoretypes){
	score_types_to_use_only_ = include_only_scoretypes;
}

void
InteractionEnergyMetric::set_ignore_scoretypes(utility::vector1<scoring::ScoreType> const & ignore_scoretypes){
	score_types_to_ignore_ = ignore_scoretypes;
}

void
InteractionEnergyMetric::set_include_rama_prepro_and_proclose( bool include ){
	include_rama_prepro_and_proclose_ = include;
}

void
InteractionEnergyMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap )
{
	SimpleMetric::parse_base_tag( tag );

	if ( tag->hasOption("scorefxn") ) {
		set_scorefunction(parse_score_function( tag, datamap ));
	}

	if ( tag->hasOption("residue_selector") && tag->hasOption("residue_selector2") ) {
		selector1_ = select::residue_selector::parse_residue_selector( tag, datamap );
	} else {
		utility_exit_with_message("InteractionEnergyMetric: Must set two residue_selectors");
	}

	if ( tag->hasOption("residue_selector2") ) {
		selector2_ = select::residue_selector::parse_residue_selector( tag, datamap, "residue_selector2" );
	}

	//////

	if ( tag->hasOption("scoretypes_only") ) {
		score_types_to_use_only_.clear();
		for ( std::string score_str : utility::string_split( tag->getOption<std::string>("scoretypes_only"), ',') ) {
			TR <<"Using Scoretype: "<< score_str << std::endl;
			scoring::ScoreType score_type = core::scoring::score_type_from_name( score_str );
			score_types_to_use_only_.push_back(score_type);
		}

	}

	if ( tag->hasOption("scoretypes_skip") ) {
		score_types_to_use_only_.clear();
		for ( std::string score_str : utility::string_split( tag->getOption<std::string>("scoretypes_skip"), ',') ) {
			TR <<"Skipping Scoretype: "<< score_str << std::endl;
			scoring::ScoreType score_type = core::scoring::score_type_from_name( score_str );
			score_types_to_ignore_.push_back(score_type);
		}

	}

	force_rescore_pose_ = tag->getOption< bool >("force_rescore", force_rescore_pose_ );
	set_include_rama_prepro_and_proclose(tag->getOption< bool >("include_rama_prepro_and_proclose", include_rama_prepro_and_proclose_));

}

void
InteractionEnergyMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;


	AttributeList attlist;
	attlist + XMLSchemaAttribute("scoretypes_only", xs_string, "Include only this list of ScoreTypes");
	attlist + XMLSchemaAttribute("scoretypes_skip", xs_string, "Always skip these score types");
	attlist + XMLSchemaAttribute::attribute_w_default("include_rama_prepro_and_proclose", xsct_rosetta_bool, "Include rama_prepro energy term AND pro_close? (Which are two-body backbone score terms) ", "false");

	attlist + XMLSchemaAttribute::attribute_w_default("force_rescore", xsct_rosetta_bool, "Should we rescore the pose, even if it has an energies object?  This will force a pose copy and rescore and will not be as efficent as having a scored pose", "false");

	attributes_for_parse_residue_selector( attlist, "residue_selector", "Selector specifying first set of residues." );
	attributes_for_parse_residue_selector( attlist, "residue_selector2", "Selector specifying second set of residues. Default is TrueResidueSelector" );

	attributes_for_get_score_function_name( attlist );

	core::simple_metrics::xsd_simple_metric_type_definition_w_attributes(xsd, name_static(),
		"Author: Jared Adolf-Bryfogle (jadolfbr@gmail.com)\n"
		"A metric for measuring the short and long range interaction energy between residues using two sets of residue selectors. Will use the current energies of the pose, unless force_rescore option is set.", attlist);
}

core::Real
InteractionEnergyMetric::calculate(const core::pose::Pose & pose ) const {

	if ( selector1_ == nullptr || selector2_ == nullptr ) {
		utility_exit_with_message("InteractionEnergyMetric requires two Residue Selectors to be set!");
	}

	utility::vector1< ScoreType > local_scoretypes_to_ignore = score_types_to_ignore_;
	if ( !include_rama_prepro_and_proclose_ ) {
		local_scoretypes_to_ignore.push_back(rama_prepro);
		local_scoretypes_to_ignore.push_back(pro_close);
	}

	// The pose must be scored for this metric to work. Here, we use a modified version of what Brian Conventry added to
	//  the Neighborhood Residue Selector.
	bool using_clone = ! pose.energies().energies_updated();
	pose::Pose pose_clone;
	if ( using_clone ) {
		pose_clone = pose;
		scorefxn_->score(pose_clone);
		if ( TR.Warning.visible() ) {
			TR.Warning << "################ Cloning pose and Scoring! ##############################" << std::endl;
			TR.Warning << "Ensure that pose is scored " << std::endl;
			TR.Warning << "before using InteractionEnergyMetric for maximum performance!" << std::endl;
			TR.Warning << "##########################################################################" << std::endl;
		}
	}
	const pose::Pose & pose_to_use = using_clone ? pose_clone : pose;

	EnergyGraph const & energy_graph( pose_to_use.energies().energy_graph() );

	// Loop over interactions across the selections

	utility::vector1< core::Size > set1 = get_residues_from_subset(selector1_->apply(pose));
	utility::vector1< core::Size > set2 = get_residues_from_subset(selector2_->apply(pose));
	utility::vector1< std::pair< core::Size, core::Size > > res_pairs;

	//Make our list, accounting for reflections
	for ( Size res_sel1 : set1 ) {
		for ( Size res_sel2 : set2 ) {
			if ( res_sel1 == res_sel2 ) continue;
			std::pair< core::Size, core::Size > res_pair     = std::make_pair(res_sel1, res_sel2);
			std::pair< core::Size, core::Size > res_pair_ref = std::make_pair(res_sel2, res_sel1);
			if ( res_pairs.contains( res_pair )  || res_pairs.contains(res_pair_ref) ) {
				continue;
			}

			res_pairs.push_back(res_pair);

		}
	}

	core::scoring::EnergyMap delta_energies_unweighted;
	core::scoring::EnergyMap final_energies_unweighted;

	core::scoring::EnergyMap weights = pose_to_use.energies().weights();
	core::Real weighted_total;

	//Short Range
	for ( auto const & res_pair : res_pairs ) {
		EnergyEdge const * edge = energy_graph.find_energy_edge( res_pair.first, res_pair.second );
		if ( edge == nullptr ) continue;
		delta_energies_unweighted += edge->fill_energy_map();
	}

	//Long Range
	for ( Size lr_enum = 1; lr_enum <= scoring::methods::n_long_range_types; ++lr_enum ) {
		auto lr_type = scoring::methods::LongRangeEnergyType( lr_enum );

		scoring::LREnergyContainerCOP lrec = pose_to_use.energies().long_range_container( lr_type );
		if ( !lrec ) continue;
		if ( lrec->empty() ) continue;

		//Have to double check that we are NOT double counting here!
		for ( auto const & pair : res_pairs ) {
			for ( scoring::ResidueNeighborConstIteratorOP
					rni = lrec->const_upper_neighbor_iterator_begin( pair.first );
					*rni != *( lrec->const_upper_neighbor_iterator_end( pair.first ) );
					++(*rni) ) {
				Size j = rni->upper_neighbor_id();
				if ( j == pair.second ) {
					scoring::EnergyMap emap;
					rni->retrieve_energy( emap );
					delta_energies_unweighted += emap;
				}
			}
			for ( scoring::ResidueNeighborConstIteratorOP
					rni = lrec->const_upper_neighbor_iterator_begin( pair.second );
					*rni != *( lrec->const_upper_neighbor_iterator_end( pair.second ) );
					++(*rni) ) {
				Size j = rni->upper_neighbor_id();
				if ( j == pair.first ) {
					scoring::EnergyMap emap;
					rni->retrieve_energy( emap );
					delta_energies_unweighted += emap;
				}
			}
		}
	}


	//Set to Ignore or Include only specific scoretypes.
	for ( core::Size ii = 1; ii <= n_score_types; ++ii ) {
		ScoreType score_type = static_cast<ScoreType>(ii);
		std::string s_name = name_from_score_type( score_type );
		if ( local_scoretypes_to_ignore.has_value( score_type ) &&  score_types_to_use_only_.has_value( score_type ) ) {
			utility_exit_with_message("InteractionEnergyMetric: A scoretype is both in the only include list and set to be ignored " + s_name);
		} else if ( score_types_to_use_only_.has_value( score_type) ) {
			final_energies_unweighted[ score_type ] = delta_energies_unweighted[ score_type ];
		} else if ( local_scoretypes_to_ignore.has_value( score_type ) ) {
			final_energies_unweighted[ score_type ] = 0.0;
		} else {
			final_energies_unweighted[ score_type ] = delta_energies_unweighted[ score_type ];
		}

		///Print some info on our piecemeal data
		if ( (weights[score_type] > 0.0) && (final_energies_unweighted[score_type] != 0) ) {
			TR<< std::endl << name_from_score_type(score_type) << " " << weights[score_type] << std::endl;
			TR<<"Unweighted: " << delta_energies_unweighted[score_type] << std::endl;
			TR<<"Weighted: " <<delta_energies_unweighted[score_type] * weights[score_type] << std::endl;
			TR.Debug<< "Dot: " << delta_energies_unweighted.dot(weights) << std::endl;
		}

	}

	//Save the total (weighted) delta score
	weighted_total = final_energies_unweighted.dot(weights);
	return weighted_total;
}

void
InteractionEnergyMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	InteractionEnergyMetric::provide_xml_schema( xsd );
}

std::string
InteractionEnergyMetricCreator::keyname() const {
	return InteractionEnergyMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
InteractionEnergyMetricCreator::create_simple_metric() const {
	return core::simple_metrics::SimpleMetricOP( new InteractionEnergyMetric );

}

} //core
} //simple_metrics
} //metrics


#ifdef    SERIALIZATION



template< class Archive >
void
core::simple_metrics::metrics::InteractionEnergyMetric::save( Archive & arc ) const {
	arc( cereal::base_class< core::simple_metrics::RealMetric>( this ) );
	arc( CEREAL_NVP( selector1_ ) );
	arc( CEREAL_NVP( selector2_ ) );
	arc( CEREAL_NVP( scorefxn_ ) );
	arc( CEREAL_NVP( score_types_to_ignore_));
	arc( CEREAL_NVP( score_types_to_use_only_));
	arc( CEREAL_NVP( include_rama_prepro_and_proclose_ ));
	arc( CEREAL_NVP( force_rescore_pose_ ));

}

template< class Archive >
void
core::simple_metrics::metrics::InteractionEnergyMetric::load( Archive & arc ) {
	arc( cereal::base_class< core::simple_metrics::RealMetric >( this ) );

	std::shared_ptr< core::select::residue_selector::ResidueSelector > local_selector;
	arc( local_selector ); // ResidueSelectorCOP
	selector1_ = local_selector; // copy the non-const pointer(s) into the const pointer(s)

	std::shared_ptr< core::select::residue_selector::ResidueSelector > local_selector2;
	arc( local_selector2 ); // ResidueSelectorCOP
	selector2_ = local_selector2; // copy the non-const pointer(s) into the const pointer(s)

	std::shared_ptr< core::scoring::ScoreFunction > local_scorefxn;
	arc( local_scorefxn ); // ResidueSelectorCOP
	scorefxn_ = local_scorefxn; // copy the non-const pointer(s) into the const pointer(s)

	arc( score_types_to_ignore_ );
	arc( score_types_to_use_only_ );
	arc( include_rama_prepro_and_proclose_ );
	arc( force_rescore_pose_);

}

SAVE_AND_LOAD_SERIALIZABLE( core::simple_metrics::metrics::InteractionEnergyMetric );
CEREAL_REGISTER_TYPE( core::simple_metrics::metrics::InteractionEnergyMetric )

CEREAL_REGISTER_DYNAMIC_INIT( core_simple_metrics_metrics_InteractionEnergyMetric )
#endif // SERIALIZATION




