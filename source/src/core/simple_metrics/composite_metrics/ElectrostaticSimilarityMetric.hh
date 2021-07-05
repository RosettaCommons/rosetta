// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/simple_metrics/composite_metrics/ElectrostaticSimilarityMetric.hh
/// @brief   Header file for ElectrostaticSimilarityMetric class
/// @details   The code closely follows Brian Coventry's implementation
///    of electrostatic complementarity with a added features
///    that in turn is based on the method:
///    McCoy, A. J., Epa, V. C., & Colman, P. M. (1997).
///    Electrostatic complementarity at protein/protein interfaces.
///    Journal of molecular biology, 268(2), 570-584.
/// @author  Andreas Scheck (andreas.scheck@epfl.ch)


#ifndef INCLUDED_core_simple_metrics_ElectrostaticSimilarityMetric_hh
#define INCLUDED_core_simple_metrics_ElectrostaticSimilarityMetric_hh

// Unit Headers
#include <core/simple_metrics/composite_metrics/ElectrostaticSimilarityMetric.fwd.hh>

// Package Headers
#include <core/simple_metrics/CompositeRealMetric.hh>

// Project Headers
#include <core/scoring/sc/ElectrostaticSimilarityCalculator.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Utility headers

// Parser headers
#include <utility/vector1.hh>


namespace core {
namespace simple_metrics {
namespace composite_metrics {

class ElectrostaticSimilarityMetric : public CompositeRealMetric {

public:// constructor/destructor
	// @brief default constructor
	ElectrostaticSimilarityMetric();

	~ElectrostaticSimilarityMetric() override = default;

	// @brief make clone
	SimpleMetricOP clone() const override { return SimpleMetricOP( new ElectrostaticSimilarityMetric( *this ) ); }

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
	void partially_solvated( bool partially_solvated ) { partially_solvated_ = partially_solvated; }
	void residue_selector1( select::residue_selector::ResidueSelectorCOP const & sel ) { selector1_ = sel; }
	void residue_selector2( select::residue_selector::ResidueSelectorCOP const & sel ) { selector2_ = sel; }
	void report_all_es( bool report ) { report_all_es_ = report; }


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
	setup_from_selectors( pose::Pose const & pose, scoring::sc::ElectrostaticSimilarityCalculator & esc ) const;


private:
	bool partially_solvated_;
	bool report_all_es_;
	core::pose::PoseCOP reference_pose_;
	select::residue_selector::ResidueSelectorCOP selector1_;
	select::residue_selector::ResidueSelectorCOP selector2_;

};


} // namespace core
} // namespace simple_metrics
} // namespace composite_metrics

#endif
