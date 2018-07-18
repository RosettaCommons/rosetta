// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/per_residue_metrics/PerResidueEnergyMetric.hh
/// @brief A per-residue metric that will calculate/output per residue total energies or a specific score component.  Correctly decomposes energies to per-residue.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_simple_metrics_per_residue_metrics_PerResidueEnergyMetric_HH
#define INCLUDED_core_simple_metrics_per_residue_metrics_PerResidueEnergyMetric_HH

#include <core/simple_metrics/per_residue_metrics/PerResidueEnergyMetric.fwd.hh>
#include <core/simple_metrics/PerResidueRealMetric.hh>

// Core headers
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreType.hh>

// Utility headers
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers
#include <map>

namespace core {
namespace simple_metrics {
namespace per_residue_metrics {

///@brief A per-residue metric that will calculate/output per residue total energies or a specific score component. WEIGHTED.
///  Correctly decomposes energies to per-residue.
///
class PerResidueEnergyMetric : public core::simple_metrics::PerResidueRealMetric{

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	PerResidueEnergyMetric();

	/// @brief Copy constructor (not needed unless you need deep copies)
	PerResidueEnergyMetric( PerResidueEnergyMetric const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~PerResidueEnergyMetric() override;

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


	///Defined in PerResidueRealMetric
	//void
	//set_residue_selector( select::residue_selector::ResidueSelectorCOP selector );

	///@brief Calculate the metric.
	/// This map is Rosetta Resnum->value and includes only those residues selected.
	///
	///@details
	/// Return by value as this function can not STORE the result, it only calculates.
	/// Store the result in the pose by using the apply method, which calls this method and stores the result
	/// in the pose as ExtraScoreValues.
	std::map< core::Size, core::Real >
	calculate( core::pose::Pose const & pose ) const override;

public:

	///@brief Set a pose into to calculate/report delta of total energy or the energy term provided.
	/// (apply_pose - comparison_pose)
	///
	void
	set_comparison_pose( pose::PoseCOP pose );

	///@brief Set a scorefunction.  Will use default Rosetta scorefunction if not set.
	void
	set_scorefunction( scoring::ScoreFunctionCOP scorefxn );

	///@brief Set the scoretype data that we return.  Default is total_score.
	void
	set_scoretype( scoring::ScoreType scoretype );

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

private:

	scoring::ScoreFunctionCOP scorefxn_ = nullptr;
	pose::PoseCOP ref_pose_ = nullptr;
	scoring::ScoreType scoretype_ = scoring::total_score;
};

} //core
} //simple_metrics
} //per_residue_metrics



#endif //core_simple_metrics_per_residue_metrics_PerResidueEnergyMetric_HH





