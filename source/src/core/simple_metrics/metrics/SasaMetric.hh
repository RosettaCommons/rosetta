// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/metrics/SasaMetric.hh
/// @brief A Metric to cacluate overall sasa of a pose.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
/// @modified Vikram K. Mulligan (vmulligan@flatironinstitute.org) -- Added support for polar or hydrophobic SASA.


#ifndef INCLUDED_core_simple_metrics_metrics_SasaMetric_hh
#define INCLUDED_core_simple_metrics_metrics_SasaMetric_hh

#include <core/simple_metrics/metrics/SasaMetric.fwd.hh>
#include <core/simple_metrics/RealMetric.hh>

#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/scoring/sasa/SasaCalc.fwd.hh>
#include <core/scoring/sasa/SasaMethod.hh>

// Utility headers
#include <core/types.hh>

#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace core {
namespace simple_metrics {
namespace metrics {

/// @brief A Metric to cacluate overall sasa of a pose.
/// @details Use the Apply method to calculate and add the metric to the pose score to be written
class SasaMetric : public core::simple_metrics::RealMetric{

public:

	/// @brief Default constructor.
	SasaMetric() = default;

	SasaMetric( select::residue_selector::ResidueSelectorCOP selector, core::scoring::sasa::SasaMethodHPMode const mode = core::scoring::sasa::SasaMethodHPMode::ALL_SASA );

	SasaMetric( SasaMetric const & ) = default;

	~SasaMetric() override;

	///Base Class Interface:
	///
	/// @brief Calculate the metric and add it to the pose as a score.
	///           labeled as prefix+metric+suffix.
	///
	/// @details Score is added through setExtraScorePose and is output
	///            into the score tables/file at pose output.
	//void
	//apply( pose::Pose & pose, prefix="", suffix="" ) override;

	///@brief Calculate The total sasa of the pose or the SASA of the residues given in a residue selector.
	core::Real
	calculate( pose::Pose const & pose ) const override;

public:

	///@brief Set a residue selector to calculate total sasa of residues in the selector.
	void
	set_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector );

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

	/// @brief Set the behaviour of this metric (count all SASA, count polar SASA, count hydrophobic SASA, etc.).
	/// @details Default is to compute all SASA.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	void set_sasa_metric_mode( core::scoring::sasa::SasaMethodHPMode const mode_in );

	/// @brief Set the behaviour of this metric (count all SASA, count polar SASA, count hydrophobic SASA, etc.).
	/// @details Default is to compute all SASA.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	void set_sasa_metric_mode( std::string const & mode_in );

	/// @brief Get the behaviour of this metric (count all SASA, count polar SASA, count hydrophobic SASA, etc.).
	/// @details Default is to compute all SASA.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	inline core::scoring::sasa::SasaMethodHPMode sasa_metric_mode() const { return sasa_metric_mode_; }

public:
	/// @brief called by parse_my_tag -- should not be used directly
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data ) override;

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


	SimpleMetricOP
	clone() const override;

private:

	core::select::residue_selector::ResidueSelectorCOP selector_;

	/// @brief The mode for this metric.  Default is to count all solvent-accessible
	/// surface area.
	core::scoring::sasa::SasaMethodHPMode sasa_metric_mode_ = core::scoring::sasa::SasaMethodHPMode::ALL_SASA;

};

} //core
} //simple_metrics
} //metrics



#endif //INCLUDED_core_metrics_SasaMetric_hh





