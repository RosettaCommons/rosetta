// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/io/SimpleMetricStruct.hh
/// @brief A struct for passing around simple metric data for IO.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_simple_metrics_SimpleMetricStruct_hh
#define INCLUDED_simple_metrics_SimpleMetricStruct_hh

#include <platform/types.hh>
#include <core/types.hh>

#include <map>

#ifdef    SERIALIZATION
// Cereal headers
#include <utility/serialization/serialization.hh>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION

namespace core {
namespace simple_metrics {

///@brief Structure for all SimpleMetric Data.  Used for mmTF IO and multistage protocols.
struct SimpleMetricStruct {

	//Simple Metrics
	std::map< std::string, std::string > string_data_;
	std::map< std::string, core::Real > real_data_;

	//CompositeMetrics
	std::map< std::string, std::map< std::string, std::string >> composite_string_data_;
	std::map< std::string, std::map< std::string, core::Real >> composite_real_data_;

	//PerResidueMetrics
	std::map< std::string, std::map< core::Size, std::string >> per_residue_string_data_;
	std::map< std::string, std::map< core::Size, core::Real >> per_residue_real_data_;

	std::map< std::string, std::map< std::string, std::string >> per_residue_string_output_;
	std::map< std::string, std::map< std::string, core::Real >> per_residue_real_output_;

	///@brief Does this struct contain any data?
	bool has_data() const {
		if ( ! string_data_.empty() )                    return true;
		else if ( ! real_data_.empty() )                 return true;
		else if ( ! composite_string_data_.empty() )     return true;
		else if ( ! composite_real_data_.empty() )       return true;
		else if ( ! per_residue_string_data_.empty() )   return true;
		else if ( ! per_residue_real_data_.empty() )     return true;
		else return false;
	}


#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} //simple_metrics
} //core


#ifdef SERIALIZATION

template< class Archive >
void
core::simple_metrics::SimpleMetricStruct::save( Archive & arc ) const {
	arc( CEREAL_NVP( string_data_ ) );
	arc( CEREAL_NVP( real_data_ ) );
	arc( CEREAL_NVP( composite_string_data_ ) );
	arc( CEREAL_NVP( composite_real_data_ ) );
	arc( CEREAL_NVP( per_residue_string_data_ ) );
	arc( CEREAL_NVP( per_residue_real_data_ ) );
	arc( CEREAL_NVP( per_residue_string_output_ ) );
	arc( CEREAL_NVP( per_residue_real_output_ ) );
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::simple_metrics::SimpleMetricStruct::load( Archive & arc ) {
	arc( string_data_ );
	arc( real_data_ );
	arc( composite_string_data_ );
	arc( composite_real_data_ );
	arc( per_residue_string_data_);
	arc( per_residue_real_data_);
	arc( per_residue_string_output_);
	arc( per_residue_real_output_);
}

SAVE_AND_LOAD_SERIALIZABLE( core::simple_metrics::SimpleMetricStruct );

#endif // SERIALIZATION

#endif //INCLUDED_core_io_SimpleMetricStruct_hh





