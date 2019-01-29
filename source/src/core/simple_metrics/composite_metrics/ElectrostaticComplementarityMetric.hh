// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/simple_metrics/composite_metrics/ElectrostaticComplementarityMetric.hh
/// @brief  header file for ElectrostaticComplementarityMetric class
/// @author Brian Coventry (bcov@uw.edu)
/// @details Based on the method from
///          McCoy, A. J., Epa, V. C., & Colman, P. M. (1997).
///              Electrostatic complementarity at protein/protein interfaces1.
///              Journal of molecular biology, 268(2), 570-584.


#ifndef INCLUDED_core_simple_metrics_ElectrostaticComplementarityMetric_hh
#define INCLUDED_core_simple_metrics_ElectrostaticComplementarityMetric_hh

// Unit Headers
#include <core/simple_metrics/composite_metrics/ElectrostaticComplementarityMetric.fwd.hh>

// Package Headers
#include <core/simple_metrics/CompositeRealMetric.hh>

// Project Headers
#include <core/scoring/sc/ElectrostaticComplementarityCalculator.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Utility headers
#include <utility/vector1.fwd.hh>

// Parser headers

#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>


//// C++ headers

namespace core {
namespace simple_metrics {
namespace composite_metrics {

class ElectrostaticComplementarityMetric : public CompositeRealMetric {

public:// constructor/destructor
	// @brief default constructor
	ElectrostaticComplementarityMetric();

	~ElectrostaticComplementarityMetric() override = default;

	// @brief make clone
	SimpleMetricOP clone() const override { return SimpleMetricOP( new ElectrostaticComplementarityMetric( *this ) ); }

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

	///@brief Get the metric name(s) that this Metric will calculate
	utility::vector1< std::string >
	get_metric_names() const override;

public:
	void ignore_radius( Real ignore_radius ) { ignore_radius_ = ignore_radius; }
	void interface_trim_radius( Real interface_trim_radius ) { interface_trim_radius_ = interface_trim_radius; }
	void partially_solvated( bool partially_solvated ) { partially_solvated_ = partially_solvated; }
	void residue_selector1( select::residue_selector::ResidueSelectorCOP const & sel ) { selector1_ = sel; }
	void residue_selector2( select::residue_selector::ResidueSelectorCOP const & sel ) { selector2_ = sel; }
	void jump_id( Size jump ) { jump_id_ = jump;}
	void report_all_ec( bool report ) { report_all_ec_ = report; }


public:
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data ) override;

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	std::map< std::string, core::Real >
	calculate( pose::Pose const & pose ) const override;

private:
	void
	setup_from_selectors( pose::Pose const & pose, scoring::sc::ElectrostaticComplementarityCalculator & ecc ) const;


private:
	Real ignore_radius_;
	Real interface_trim_radius_;
	bool partially_solvated_;
	Size jump_id_;
	bool report_all_ec_;
	select::residue_selector::ResidueSelectorCOP selector1_;
	select::residue_selector::ResidueSelectorCOP selector2_;
};

/// @brief Super-simple exception to be thrown when we can't initialize the EC calculator
class EXCN_InitFailed : public utility::excn::Exception {
public:
	using utility::excn::Exception::Exception;
};

/// @brief Super-simple exception to be thrown when the EC calculator fails to compute
class EXCN_CalcFailed : public utility::excn::Exception {
public:
	using utility::excn::Exception::Exception;
};


} // namespace core
} // namespace simple_metrics
} // namespace composite_metrics

#endif
