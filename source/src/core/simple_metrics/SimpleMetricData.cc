// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/simple_metrics/SimpleMetricData.cc
/// @brief A container class for all Simple Metrics stored in the pose.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// unit headers
#include <core/simple_metrics/SimpleMetricData.hh>

#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/string_util.hh>
#include <basic/datacache/CacheableData.hh>


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/utility.hpp>



#endif // SERIALIZATION


static basic::Tracer TR( "core.simple_metrics.SimpleMetricData" );


namespace core {
namespace simple_metrics {

SimpleMetricData::SimpleMetricData() : basic::datacache::CacheableData() {}

SimpleMetricData::~SimpleMetricData() = default;

basic::datacache::CacheableDataOP
SimpleMetricData::clone() const {
	return utility::pointer::make_shared< SimpleMetricData >(*this);
}

SimpleMetricDataOP
SimpleMetricData::shared_from_this() {
	return utility::pointer::static_pointer_cast<SimpleMetricData>( CacheableData::shared_from_this() );
}

/////////////////////////////////////////////
/////           Data Access             /////
/////////////////////////////////////////////

void
SimpleMetricData::clear(){
	//Simple Metrics
	data_.string_data_.clear();
	data_.real_data_.clear();

	//CompositeMetrics
	data_.composite_string_data_.clear();
	data_.composite_real_data_.clear();

	//PerResidueMetrics
	data_.per_residue_string_data_.clear();
	data_.per_residue_real_data_.clear();

	data_.per_residue_string_output_.clear();
	data_.per_residue_real_output_.clear();
}

/////////////////////////////////////////////
/////          Simple Metrics           /////
/////////////////////////////////////////////

///@brief Get RealMetric data.  Return success status.
bool
SimpleMetricData::get_value(std::string const & name, Real & value) const{
	if ( data_.real_data_.count(name) ) {
		value = data_.real_data_.at(name);
		return true;
	} else {
		return false;
	}
}

///@brief Get StringMetric data.  Return success status.
bool
SimpleMetricData::get_value(std::string const & name, std::string & value) const{
	if ( data_.string_data_.count(name) ) {
		value = data_.string_data_.at(name);
		return true;
	} else {
		return false;
	}
}


/////////////////////////////////////////////
/////        Composite Metrics          /////
/////////////////////////////////////////////

///@brief Get CompositeRealMetric data.  Return success status.
bool
SimpleMetricData::get_value(std::string const & name, std::map< std::string, Real > & value) const {
	if ( data_.composite_real_data_.count(name) ) {
		value = data_.composite_real_data_.at(name);
		return true;
	} else {
		return false;
	}
}

///@brief Get CompositeStringMetric data.  Return success status.
bool
SimpleMetricData::get_value(std::string const & name, std::map< std::string, std::string > & value) const {
	if ( data_.composite_string_data_.count(name) ) {
		value = data_.composite_string_data_.at(name);
		return true;
	} else {
		return false;
	}
}


/////////////////////////////////////////////
/////       Per-Residue Metrics         /////
/////////////////////////////////////////////

///@brief Get PerResidueRealMetric data.  Return success status.
bool
SimpleMetricData::get_value( std::string const & name, std::map< Size, Real > & value) const {
	if ( data_.per_residue_real_data_.count(name) ) {
		value = data_.per_residue_real_data_.at(name);
		return true;
	} else {
		return false;
	}
}

///@brief Get PerResidueStringMetric data.  Return success status.
bool
SimpleMetricData::get_value( std::string const & name, std::map< Size, std::string > & value) const {
	if ( data_.per_residue_string_data_.count(name) ) {
		value = data_.per_residue_string_data_.at(name);
		return true;
	} else {
		return false;
	}
}


bool
SimpleMetricData::get_value( std::string const & name, std::map< Size, Real > & value, pose::Pose const & pose, bool use_ref_pose) const {

	if ( ! use_ref_pose ) {
		return get_value(name, value);
	}


	if ( data_.per_residue_real_data_.count(name) ) {
		std::map< Size, Real > ref_values;
		for ( auto const & res_pair: data_.per_residue_real_data_.at(name) ) {
			core::Size new_res = pose.corresponding_residue_in_current( res_pair.first, name);
			if ( new_res != 0 ) {
				ref_values[new_res] = res_pair.second;
			}
		}
		value = ref_values;
		return true;
	} else {
		return false;
	}
}


bool
SimpleMetricData::get_value( std::string const & name, std::map< Size, std::string > & value, pose::Pose const & pose, bool use_ref_pose) const {

	if ( ! use_ref_pose ) {
		return get_value(name, value);
	}

	if ( data_.per_residue_string_data_.count(name) ) {
		std::map< Size, std::string > ref_values;
		for ( auto const & res_pair: data_.per_residue_string_data_.at(name) ) {
			core::Size new_res = pose.corresponding_residue_in_current( res_pair.first, name);
			if ( new_res != 0 ) {
				ref_values[new_res] = res_pair.second;
			}
		}
		value = ref_values;
		return true;
	} else {
		return false;
	}
}


/////////////////////////////////////////////
/////         All Data Access           /////
/////////////////////////////////////////////

SimpleMetricStruct const &
SimpleMetricData::get_all_sm_data() const {
	return data_;
}

///@brief Get all RealMetric data
std::map< std::string, Real > const &
SimpleMetricData::get_real_metric_data() const {
	return data_.real_data_;
}

///@brief Get all StringMetric data
std::map< std::string, std::string > const &
SimpleMetricData::get_string_metric_data() const {
	return data_.string_data_;
}

///@brief Get all CompositeRealMetric data
std::map< std::string, std::map< std::string, Real >> const &
SimpleMetricData::get_composite_real_metric_data() const {
	return data_.composite_real_data_;
}

///@brief Get all CompositeStringMetric data
std::map< std::string, std::map< std::string, std::string >> const &
SimpleMetricData::get_composite_string_metric_data() const {
	return data_.composite_string_data_;
}

///@brief Get all PerResidueRealMetric data
std::map< std::string, std::map< core::Size, Real >> const &
SimpleMetricData::get_per_residue_real_metric_data() const {
	return data_.per_residue_real_data_;
}

///@brief Get all PerResidueStringMetric data
std::map< std::string, std::map< core::Size, std::string >> const &
SimpleMetricData::get_per_residue_string_metric_data() const {
	return data_.per_residue_string_data_;
}

///@brief Get all PerResidueRealMetric data
std::map< std::string, std::map< std::string, Real >> const &
SimpleMetricData::get_per_residue_real_metric_output() const {
	return data_.per_residue_real_output_;
}

///@brief Get all PerResidueStringMetric data
std::map< std::string, std::map< std::string, std::string >> const &
SimpleMetricData::get_per_residue_string_metric_output() const {
	return data_.per_residue_string_output_;
}

/////////////////////////////////////////////
/////          Simple Metrics           /////
/////////////////////////////////////////////

///@brief Set RealMetric data
void
SimpleMetricData::set_value( MetricKey , std::string const & name, Real value){
	data_.real_data_[name] = value;
}

///@brief Set StringMetric data
void
SimpleMetricData::set_value( MetricKey, std::string const & name, std::string const & value){
	data_.string_data_[name] = value;
}


/////////////////////////////////////////////
/////        Composite Metrics          /////
/////////////////////////////////////////////

///@brief Set CompositeRealMetric data
void
SimpleMetricData::set_value( MetricKey, std::string const & name, std::map< std::string, Real > const & value){
	data_.composite_real_data_[name] = value;
}

///@brief Set CompositeStringMetric data
void
SimpleMetricData::set_value( MetricKey, std::string const & name, std::map< std::string, std::string > const & value){
	data_.composite_string_data_[name] = value;
}


/////////////////////////////////////////////
/////       Per-Residue Metrics         /////
/////////////////////////////////////////////

///@brief Set PerResidueRealMetric data
/// Creates a ReferencePose with the given name for the pose to maintain data integrity
void
SimpleMetricData::set_value(
	MetricKey,
	pose::Pose & pose,
	std::string const & name,
	std::map< Size, Real > const & value,
	bool output_as_pdb_nums)
{

	std::string strip=" ";
	std::map<std::string, Real > output;
	for ( auto const & pair: value ) {
		std::string out_res;
		if ( output_as_pdb_nums ) {
			out_res = utility::remove_from_string(pose.pdb_info()->pose2pdb(pair.first), strip);
		} else {
			out_res = utility::to_string(pair.first);
		}
		output[out_res] = pair.second;
	}

	data_.per_residue_real_data_[name] = value;
	data_.per_residue_real_output_[name] = output;
	pose.reference_pose_from_current(name, true /*override*/);
	TR.Debug << "Created reference pose with name: " << name << std::endl;
}

///@brief Set PerResidueStringMetric data
/// Creates a ReferencePose with the given name for the pose to maintain data integrity
void
SimpleMetricData::set_value(
	MetricKey,
	pose::Pose & pose,
	std::string const & name,
	std::map< Size, std::string> const & value,
	bool output_as_pdb_nums)
{
	std::string strip=" ";
	std::map<std::string, std::string > output;
	for ( auto const & pair: value ) {
		std::string out_res;
		if ( output_as_pdb_nums ) {
			out_res = utility::remove_from_string(pose.pdb_info()->pose2pdb(pair.first), strip);
		} else {
			out_res = utility::to_string(pair.first);
		}
		output[out_res] = pair.second;
	}

	data_.per_residue_string_data_[name] = value;
	data_.per_residue_string_output_[name] = output;
	pose.reference_pose_from_current(name, true /*override*/);
	TR << "Created reference pose with name: " << name << std::endl;
}

void
SimpleMetricData::set_all_data(core::simple_metrics::MetricKey , SimpleMetricStruct const & data){
	data_ = data;
}

void
SimpleMetricData::show() const {

	TR << std::endl << "String Data" << std::endl;
	for ( auto const & kv : data_.real_data_ ) {
		TR << kv.first<<" " << kv.second << std::endl;
	}

	TR << std::endl << "Real Data" << std::endl;
	for ( auto const & kv : data_.string_data_ ) {
		TR << kv.first<< " " << kv.second << std::endl;
	}

	TR << std::endl << "CompositeString Data" << std::endl;
	for ( auto const & kv : data_.composite_string_data_ ) {
		TR << kv.first<< " " << utility::to_string(kv.second) << std::endl;
	}

	TR << std::endl << "CompositeReal Data" << std::endl;
	for ( auto const & kv : data_.composite_real_data_ ) {
		TR << kv.first<< " " << utility::to_string(kv.second) << std::endl;
	}

	TR << "PerResidueString Data" << std::endl;
	for ( auto const & kv : data_.per_residue_string_data_ ) {
		TR << kv.first<< " " << utility::to_string(kv.second) << std::endl;
	}

	TR << "PerResidueReal Data" << std::endl;
	for ( auto const & kv : data_.per_residue_real_data_ ) {
		TR << kv.first<< " " << utility::to_string(kv.second) << std::endl;
	}

}

} // namespace datacache
} // namespace basic


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::simple_metrics::SimpleMetricData::save( Archive & arc ) const {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc(data_);

}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::simple_metrics::SimpleMetricData::load( Archive & arc ) {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc(data_);

}

SAVE_AND_LOAD_SERIALIZABLE(core::simple_metrics::SimpleMetricData );
CEREAL_REGISTER_TYPE( core::simple_metrics::SimpleMetricData )

CEREAL_REGISTER_DYNAMIC_INIT( core_simple_metrics_SimpleMetricData )
#endif // SERIALIZATION
