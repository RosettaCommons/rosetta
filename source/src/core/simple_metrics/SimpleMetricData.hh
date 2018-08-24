// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/simple_metrics/SimpleMetricData.hh
/// @brief A container class for all SMs stored in the pose.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_basic_datacache_SimpleMetricData_hh
#define INCLUDED_basic_datacache_SimpleMetricData_hh

// unit headers
#include <core/simple_metrics/SimpleMetricData.fwd.hh>

// package headers
#include <basic/datacache/CacheableData.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

#include <core/simple_metrics/RealMetric.fwd.hh>
#include <core/simple_metrics/StringMetric.fwd.hh>
#include <core/simple_metrics/CompositeRealMetric.fwd.hh>
#include <core/simple_metrics/CompositeStringMetric.fwd.hh>
#include <core/simple_metrics/PerResidueRealMetric.fwd.hh>
#include <core/simple_metrics/PerResidueStringMetric.fwd.hh>

// C++ headers
#include <map>
#include <string>

#include <platform/types.hh>
#include <utility/down_cast.hh>
#include <cassert>
#include <cstddef>
#include <iosfwd>
#include <basic/datacache/CacheableData.fwd.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace simple_metrics {


///@brief The key needed to set data in our SimpleMetricData.
///  This allows only SimpleMetrics to set data
///
class MetricKey{

private:

	MetricKey(){}


	//Key Holders:

	friend RealMetric;
	friend StringMetric;
	friend CompositeRealMetric;
	friend CompositeStringMetric;
	friend PerResidueRealMetric;
	friend PerResidueStringMetric;

};


/// @brief A container class for all Simple Metrics stored in the pose.
class SimpleMetricData : public basic::datacache::CacheableData
{
public:
	SimpleMetricData();

	~SimpleMetricData() override;

	basic::datacache::CacheableDataOP
	clone() const override;

	SimpleMetricDataOP
	shared_from_this();


	/////////////////////////////////////////////
	/////           Data Access             /////
	/////////////////////////////////////////////

public:


	/////////////////////////////////////////////
	/////          Simple Metrics           /////
	/////////////////////////////////////////////

	///@brief Get RealMetric data. Return success status.
	///
	///@return
	///  Returns true if data was present.
	///
	bool
	get_value(std::string const & name, Real & value) const;

	///@brief Get StringMetric data.  Return success status.
	///
	///@return
	///  Returns true if data was present.
	///
	bool
	get_value(std::string const & name, std::string & value) const;


	/////////////////////////////////////////////
	/////        Composite Metrics          /////
	/////////////////////////////////////////////

	///@brief Get CompositeRealMetric data.  Return success status.
	///
	///@return
	///  Returns true if data was present.
	///
	bool
	get_value(std::string const & name, std::map< std::string, Real > & value) const;

	///@brief Get CompositeStringMetric data.  Return success status.
	///
	///@return
	///  Returns true if data was present.
	///
	bool
	get_value(std::string const & name, std::map< std::string, std::string > & value) const;


	/////////////////////////////////////////////
	/////       Per-Residue Metrics         /////
	/////////////////////////////////////////////

	///@brief Get PerResidueRealMetric data.  Return success status.
	///
	///@return
	///  Returns true if data was present.
	///
	bool
	get_value( std::string const & name, std::map< Size, Real > & value) const;

	///@brief Get PerResidueStringMetric data.  Return success status.
	///
	///@return
	///  Returns true if data was present.
	///
	bool
	get_value( std::string const & name, std::map< Size, std::string > & value) const;


	///@brief Get PerResidueRealMetric data, optionally convert using refpose.
	/// Any resnum not in current is removed.
	///
	///@return
	///  Returns true if data was present.
	///
	bool
	get_value( std::string const & name, std::map< Size, Real > & value, pose::Pose const & pose, bool use_ref_pose ) const;

	///@brief Get PerResidueStringMetric data, optionally convert using refpose.
	/// Any resnum not in current is removed.
	///
	///@return
	///  Returns true if data was present.
	///
	bool
	get_value( std::string const & name, std::map< Size, std::string > & value, pose::Pose const & pose, bool use_ref_pose ) const;

public:

	/////////////////////////////////////////////
	/////         All Data Access           /////
	/////////////////////////////////////////////

	///@brief Get all RealMetric data
	std::map< std::string, Real > const &
	get_real_metric_data() const;

	///@brief Get all StringMetric data
	std::map< std::string, std::string > const &
	get_string_metric_data() const;

	///@brief Get all CompositeRealMetric data
	std::map< std::string, std::map< std::string, Real >> const &
	get_composite_real_metric_data() const;

	///@brief Get all CompositeStringMetric data
	std::map< std::string, std::map< std::string, std::string >> const &
	get_composite_string_metric_data() const;

	///@brief Get all PerResidueRealMetric data
	/// Raw data - no ref-pose conversions.
	std::map< std::string, std::map< core::Size, Real >> const &
	get_per_residue_real_metric_data() const;

	///@brief Get all PerResidueStringMetric data
	/// Raw data - no ref-pose conversions.
	std::map< std::string, std::map< core::Size, std::string >> const &
	get_per_residue_string_metric_data() const;

	///@brief Get all PerResidueRealMetric output
	std::map< std::string, std::map< std::string, Real >> const &
	get_per_residue_real_metric_output() const;

	///@brief Get all PerResidueStringMetric output
	std::map< std::string, std::map< std::string, std::string >> const &
	get_per_residue_string_metric_output() const;


public:

	/////////////////////////////////////////////
	/////              SETTERS              /////
	/////////////////////////////////////////////

	/////////////////////////////////////////////
	/////          Simple Metrics           /////
	/////////////////////////////////////////////

	///@brief Set RealMetric data
	void
	set_value( MetricKey mk, std::string const & name, Real value);

	///@brief Set StringMetric data
	void
	set_value( MetricKey mk, std::string const & name, std::string const & value);


	/////////////////////////////////////////////
	/////        Composite Metrics          /////
	/////////////////////////////////////////////

	///@brief Set CompositeRealMetric data
	void
	set_value( MetricKey mk, std::string const & name, std::map< std::string, Real > const & value);

	///@brief Set CompositeStringMetric data
	void
	set_value( MetricKey mk, std::string const & name, std::map< std::string, std::string > const & value);


	/////////////////////////////////////////////
	/////       Per-Residue Metrics         /////
	/////////////////////////////////////////////

	///@brief Set PerResidueRealMetric data
	/// Creates a ReferencePose with the given name for the pose to maintain data integrity
	//  Creates an output map of the given data
	void
	set_value(
		MetricKey mk,
		pose::Pose & pose,
		std::string const & name,
		std::map< Size, Real > const & value,
		bool outut_as_pdb_nums=true);

	///@brief Set PerResidueStringMetric data
	/// Creates a ReferencePose with the given name for the pose to maintain data integrity
	/// Creates an ouput map of the given data
	void
	set_value(
		MetricKey mk,
		pose::Pose & pose,
		std::string const & name,
		std::map< Size, std::string> const & value,
		bool output_as_pdb_nums=true);


private:

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

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} // namespace simple_metrics
} // namespace core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_simple_metrics_SimpleMetricData )
#endif // SERIALIZATION


#endif /* INCLUDED_core_simple_metrics_SimpleMetricData_HH */
