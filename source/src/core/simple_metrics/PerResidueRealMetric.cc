// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/simple_metrics/PerResidueRealMetric.cc
///
/// @brief Main class for simple metrics.
/// @author Jared Adolf-Bryfogle ( jadolfbr@gmail.com )

// Unit Headers
#include <core/simple_metrics/PerResidueRealMetric.hh>
#include <core/simple_metrics/util.hh>

// Protocol Headers

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/simple_metrics/SimpleMetricData.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>

#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tag/util.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace core {
namespace simple_metrics {
using namespace core::select::residue_selector;

PerResidueRealMetric::PerResidueRealMetric():
	SimpleMetric("PerResidueRealMetric")
{
	selector_ = TrueResidueSelectorOP( new TrueResidueSelector());
}


PerResidueRealMetric::~PerResidueRealMetric() = default;

PerResidueRealMetric::PerResidueRealMetric( PerResidueRealMetric const & src ):
	SimpleMetric( src )
{


}

utility::vector1< std::string >
PerResidueRealMetric::get_metric_names() const {
	utility::vector1< std::string > names;
	names.push_back( metric() );
	return names;
}

void
PerResidueRealMetric::set_residue_selector(select::residue_selector::ResidueSelectorCOP selector){
	selector_ = selector;
}

void
PerResidueRealMetric::set_output_as_pdb_nums(bool output_as_pdb_nums){
	output_as_pdb_nums_ = output_as_pdb_nums;
}

select::residue_selector::ResidueSelectorCOP
PerResidueRealMetric::get_selector() const {
	return selector_;
}

void
PerResidueRealMetric::parse_per_residue_tag(utility::tag::TagCOP tag, basic::datacache::DataMap & datamap){
	set_output_as_pdb_nums(tag->getOption< bool >("output_as_pdb_nums", false));
	if ( tag->hasOption("residue_selector") ) {
		set_residue_selector(select::residue_selector::parse_residue_selector( tag, datamap ));
	}
}

///@brief Add options to the schema from this base class.
void
PerResidueRealMetric::add_schema( utility::tag::XMLSchemaComplexTypeGeneratorOP ct_gen) {
	add_per_residue_simple_metric_schema( ct_gen );
}

std::map< core::Size, core::Real>
PerResidueRealMetric::cached_calculate(
	pose::Pose const & pose,
	bool use_cache,
	std::string prefix,
	std::string suffix,
	bool fail_on_missing_cache,
	bool use_ref_pose_for_cache) const {

	std::string name = prefix + get_final_sm_type() + suffix;

	if ( use_cache && has_sm_data( pose ) ) {
		std::map< core::Size, core::Real > value;
		bool data_found = get_sm_data(pose)->get_value(name, value, pose, use_ref_pose_for_cache);
		if ( data_found ) {
			return value;
		} else if ( fail_on_missing_cache ) {
			utility_exit_with_message("Could not find PerResidueRealMetric: "+name+" in pose");
		} else {
			return calculate(pose);
		}
	} else {
		return calculate(pose);
	}
}

void
PerResidueRealMetric::apply( pose::Pose & pose, std::string prefix, std::string suffix, bool override_existing ) const {
	std::map< core::Size, core::Real > const values = calculate( pose ); //Index to value map

	std::string out_tag = prefix  + get_final_sm_type() + suffix;
	MetricKey mk;

	std::map< core::Size, core::Real > stored_value;

	if ( ( ! override_existing ) && get_sm_data(pose)->get_value(out_tag, stored_value) ) {
		throw_sm_override_error(out_tag, name());
	}

	get_sm_data(pose)->set_value(mk, pose, out_tag, values, output_as_pdb_nums_);

	/*
	//Output sum as well.
	core::Real sum = 0;
	std::string strip=" ";
	for ( auto value_pair : values ) {

	std::string out_number = utility::to_string(value_pair.first);

	if ( output_as_pdb_nums_ ) {

	out_number = utility::remove_from_string(pose.pdb_info()->pose2pdb(value_pair.first), strip);
	}

	sum += value_pair.second;
	std::string const out_tag = prefix  + get_final_sm_type() + "_"+ out_number + suffix;
	core::pose::setPoseExtraScore( pose, out_tag, value_pair.second);
	}

	std::string const out_tag = prefix + "SUM" + "_" + get_final_sm_type() + suffix;
	core::pose::setPoseExtraScore(pose, out_tag, sum);
	*/
}





} //namespace simple_metrics
} //namespace core


#ifdef    SERIALIZATION



template< class Archive >
void
core::simple_metrics::PerResidueRealMetric::save( Archive & arc ) const {
	arc( cereal::base_class< core::simple_metrics::SimpleMetric>( this ) );
	arc( CEREAL_NVP( selector_));
	arc( CEREAL_NVP( output_as_pdb_nums_ ) );

}

template< class Archive >
void
core::simple_metrics::PerResidueRealMetric::load( Archive & arc ) {
	arc( cereal::base_class< core::simple_metrics::SimpleMetric >( this ) );
	std::shared_ptr< core::select::residue_selector::ResidueSelector > local_selector;
	arc( local_selector ); // ResidueSelectorCOP
	selector_ = local_selector; // copy the non-const pointer(s) into the const pointer(s)
	arc( output_as_pdb_nums_ );


}

SAVE_AND_LOAD_SERIALIZABLE( core::simple_metrics::PerResidueRealMetric );
CEREAL_REGISTER_TYPE( core::simple_metrics::PerResidueRealMetric )

CEREAL_REGISTER_DYNAMIC_INIT( core_simple_metrics_PerResidueRealMetric )
#endif // SERIALIZATION


