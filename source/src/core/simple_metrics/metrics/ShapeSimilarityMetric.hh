// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file     core/simple_metrics/metrics/ShapeSimilarityMetric.hh
/// @brief    Header file for ShapeSimilarityMetric class
/// @details  The code modifies Luki Goldschmidt's
///           implementation of Lawrence & Coleman shape complementarity calculator
///           to allow for the comparison of similar surface shapes.
/// @author   Andreas Scheck (andreas.scheck@epfl.ch)

#ifndef INCLUDED_core_simple_metrics_metrics_ShapeSimilarityMetric_hh
#define INCLUDED_core_simple_metrics_metrics_ShapeSimilarityMetric_hh

// Unit Headers
#include <core/simple_metrics/metrics/ShapeSimilarityMetric.fwd.hh>

// Package Headers
#include <core/simple_metrics/RealMetric.hh>


// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/sc/ShapeSimilarityCalculator.hh>
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/vector1.hh>

// Parser headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/excn/Exceptions.hh>


namespace core {
namespace simple_metrics {
namespace metrics {


class ShapeSimilarityMetric : public RealMetric {

public:
	typedef core::scoring::sc::RESULTS ShapeSimilarityCalculatorResults;

public:// constructor/destructor
	// @brief default constructor
	ShapeSimilarityMetric();

	~ShapeSimilarityMetric() override= default;

	// @brief make clone
	SimpleMetricOP clone() const override { return SimpleMetricOP( new ShapeSimilarityMetric( *this ) ); }

public:
	///@brief Name of the class
	std::string
	name() const override;

	///@brief Name of the class for creator.
	static
	std::string
	name_static();

	///@brief Name of the metric
	std::string
	metric() const override;


public:// mutator
	void quick( core::Size const & quick ) { quick_ = quick; }
	void verbose( core::Size const & verbose ) { verbose_ = verbose; }
	void write_int_area( bool const & write_int_area ) { write_int_area_ = write_int_area; }


public:
	void parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap & data_map
	) override;

	/// @brief calc shape similarity, returns results of the ShapeSimilarityCalculator
	core::Real
	calculate( pose::Pose const & pose ) const override;

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	/// @brief Uses residue selectors to set up the ShapeSimilarityCalculator
	/// @param[in]  pose Pose to be analyzed
	/// @param[out] ssc Initialized, empty ShapeSimilarityCalculator, to which pose residues are added
	void
	setup_from_selectors( pose::Pose const & pose, scoring::sc::ShapeSimilarityCalculator & ssc ) const;

	/// @brief prints results to given tracer in a human-readable format
	/// @param[out] tr std::ostream object to write to
	/// @param[in]  r  ShapeSimilarityCalculatorResults object containing results
	void
	print_ss_result(
		std::ostream & tr,
		ShapeSimilarityCalculatorResults const & r,
		core::Real const nsubs_scalefactor ) const;


private:
	core::Size quick_;
	core::Size verbose_;
	core::pose::PoseCOP reference_pose_;
	core::select::residue_selector::ResidueSelectorCOP selector1_;
	core::select::residue_selector::ResidueSelectorCOP selector2_;
	bool write_int_area_;
	Real dist_weight_;
	bool median_;
};

/// @brief Super-simple exception to be thrown when we can't initialize the SS calculator
class EXCN_InitFailed : public utility::excn::Exception {
public:
	using utility::excn::Exception::Exception;
};

/// @brief Super-simple exception to be thrown when the SS calculator fails to compute
class EXCN_CalcFailed : public utility::excn::Exception {
public:
	using utility::excn::Exception::Exception;
};

} // core
} // simple_metrics
} // metrics

#endif //INCLUDED_core_simple_metrics_ShapeSimilarityMetric_hh
