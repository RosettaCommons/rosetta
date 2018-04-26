// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/metrics/SelectedResidueCountMetric.hh
/// @brief A SimpleMetric that counts the number of residues in a residue selection.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

#ifndef INCLUDED_core_simple_metrics_metrics_SelectedResidueCountMetric_HH
#define INCLUDED_core_simple_metrics_metrics_SelectedResidueCountMetric_HH

#include <core/simple_metrics/metrics/SelectedResidueCountMetric.fwd.hh>
#include <core/simple_metrics/RealMetric.hh>

// Core headers
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.hh>

// Utility headers
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace core {
namespace simple_metrics {
namespace metrics {

///@brief A SimpleMetric that counts the number of residues in a residue selection.
class SelectedResidueCountMetric : public core::simple_metrics::RealMetric{

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	SelectedResidueCountMetric();

	/// @brief Copy constructor (not needed unless you need deep copies)
	SelectedResidueCountMetric( SelectedResidueCountMetric const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~SelectedResidueCountMetric() override;

	/////////////////////
	/// Metric Methods ///
	/////////////////////

	///Defined in RealMetric:
	///
	/// @brief Calculate the metric and add it to the pose as a score.
	///           labeled as prefix+metric+suffix.
	///
	/// @details Score is added through setExtraScorePose and is output
	///            into the score tables/file at pose output.
	//void
	//apply( pose::Pose & pose, prefix="", suffix="" ) override;

	///@brief Calculate the metric.
	core::Real
	calculate( core::pose::Pose const & pose ) const override;

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

public:

	/// @brief called by parse_my_tag -- should not be used directly
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data ) override;

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	core::simple_metrics::SimpleMetricOP
	clone() const override;

	/// @brief Set the residue selector.
	/// @details  Copies the input pointer; doesn't clone the object.
	void set_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector_in );

	/// @brief Remove the residue selector.
	/// @details In the absence of a residue selector, this metric returns the number of residues in a pose.
	void remove_residue_selector();

private:

	/// @brief A residue selector used for counting.
	/// @details If left nullptr, the whole pose is counted.
	core::select::residue_selector::ResidueSelectorCOP residue_selector_;
};

} //protocols
} //analysis
} //simple_metrics



#endif //core_simple_metrics_metrics_SelectedResidueCountMetric_HH





