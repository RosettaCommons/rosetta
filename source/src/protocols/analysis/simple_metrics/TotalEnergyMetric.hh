// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/metrics/TotalEnergyMetric.hh
/// @brief A metric to report the total energy of the system or the delta total energy between another input pose or the set native.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_analysis_simple_metrics_TotalEnergyMetric_HH
#define INCLUDED_protocols_analysis_simple_metrics_TotalEnergyMetric_HH

#include <protocols/analysis/simple_metrics/TotalEnergyMetric.fwd.hh>
#include <core/simple_metrics/RealMetric.hh>

// Core headers
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Utility headers
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace protocols {
namespace analysis {
namespace simple_metrics {

///@brief A metric to report the total energy of the system
/// or the delta total energy between another input pose or the set native.
///
class TotalEnergyMetric : public core::simple_metrics::RealMetric{

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	TotalEnergyMetric();

	TotalEnergyMetric(
		core::select::residue_selector::ResidueSelectorCOP selector );

	TotalEnergyMetric(
		core::select::residue_selector::ResidueSelectorCOP selector,
		core::scoring::ScoreFunctionCOP scorefxn );

	/// @brief Copy constructor (not needed unless you need deep copies)
	TotalEnergyMetric( TotalEnergyMetric const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~TotalEnergyMetric() override;

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
	/// Returns the total score from the scorefunction or the total score of each residue in the residue selector.
	///
	///@details
	/// If the native pose is given on the command line or a comparison pose is given
	///  will calculate the delta between them (pose - comparison_pose).
	core::Real
	calculate( core::pose::Pose const & pose ) const override;

public:

	///@brief Set a scorefunction.  Will use default Rosetta scorefunction if not set.
	void
	set_scorefunction( core::scoring::ScoreFunctionCOP scorefxn );

	///@brief Set a residue selector to calculate total energy of a subset of residues.
	void
	set_residue_selector( core::select::residue_selector::ResidueSelectorCOP residue_selector );

	///@brief Set a pose into to calculate/report delta of total energy.
	/// (apply_pose - comparison_pose)
	///
	void
	set_comparison_pose( core::pose::Pose const & pose );

	///@brief Load any set in:file:native as the reference pose.
	///  This is opt-in to save on loading time.
	void
	load_native_pose_as_reference();

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
		basic::datacache::DataMap  & data ) override;

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	core::simple_metrics::SimpleMetricOP
	clone() const override;

private:

	core::scoring::ScoreFunctionCOP scorefxn_ = nullptr;
	core::select::residue_selector::ResidueSelectorCOP residue_selector_ = nullptr;
	core::pose::PoseOP ref_pose_ = nullptr;

};

} // metrics
} // simple_metrics
} // core



#endif //protocols_analysis_simple_metrics_TotalEnergyMetric_HH





