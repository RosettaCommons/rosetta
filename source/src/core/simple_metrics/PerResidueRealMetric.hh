// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/simple_metrics/PerResidueRealMetric.hh
///
/// @brief  Base class for core::Real PerResidueRealMetrics
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_simple_metrics_PerResidueRealMetric_hh
#define INCLUDED_core_simple_metrics_PerResidueRealMetric_hh

// Unit headers
#include <core/simple_metrics/SimpleMetric.hh>
#include <core/simple_metrics/PerResidueRealMetric.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Core headers

#include <core/pose/Pose.fwd.hh>

#include <core/types.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ includes
#include <map>

namespace core {
namespace simple_metrics {

/// @brief A class that is derived to calculate a set of core::Real values for each Residue.
///  Apply(pose) method calculates this metric and adds it to the pose score for output.
///
///
/// @details
///  Calculate(pose) method calculates core::Real values and returns them as a map<Size, Real>
///   Resnum:Value
///
class PerResidueRealMetric : public core::simple_metrics::SimpleMetric {

public: // constructors / destructors

	PerResidueRealMetric();

	~PerResidueRealMetric() override;

	PerResidueRealMetric( PerResidueRealMetric const & other );

	/// @brief Calculate the metric and add it to the pose as a score.
	///           labeled as prefix+metric+suffix.
	///
	/// @details Score is added through setPoseExtraScore and is output
	///            into the score table/ score file at pose output.
	void
	apply( pose::Pose & pose, std::string prefix="", std::string suffix="" ) const override;

	///@brief Set a ResidueSelector for which we will calculate values over.
	void
	set_residue_selector( select::residue_selector::ResidueSelectorCOP selector );

	///@brief Set to output in PDB numbering instead of Rosetta during the Apply function,
	/// which adds the data to pose as extra scores.
	void
	set_output_as_pdb_nums( bool output_as_pdb_nums );

	///@brief Calculate the metric.
	/// This map is Rosetta Resnum->value and includes only those residues selected.
	///
	///@details
	/// Return by value as this function can not STORE the result, it only calculates.
	/// Store the result in the pose by using the apply method, which calls this method and stores the result
	/// in the pose as ExtraScoreValues.
	virtual std::map< core::Size, core::Real >
	calculate( pose::Pose const & pose ) const = 0;

public:

	///@brief Name of the class
	virtual std::string
	name() const override = 0;

	///@brief Name of the metric
	virtual std::string
	metric() const override = 0;

	///@brief Get the submetric names that this Metric will calculate
	utility::vector1< std::string >
	get_metric_names() const override;

	///@brief Get the set residue selector of this class.
	select::residue_selector::ResidueSelectorCOP
	get_selector() const;

public:
	/// @brief called by parse_my_tag -- should not be used directly
	virtual void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data ) override = 0;

	virtual SimpleMetricOP
	clone() const override = 0;

	///@brief Add options to the schema from this base class.
	static
	void
	add_schema( utility::tag::XMLSchemaComplexTypeGeneratorOP complex_schema);

	///Parse the base class tag.  Keep required interface for parse_my_tag.
	virtual void
	parse_per_residue_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data );


private:

	select::residue_selector::ResidueSelectorCOP selector_; //Set as a TrueResidueSelector.
	bool output_as_pdb_nums_ = false;

}; //class PerResidueRealMetrics

} //namespace simple_metrics
} //namespace core


#endif
