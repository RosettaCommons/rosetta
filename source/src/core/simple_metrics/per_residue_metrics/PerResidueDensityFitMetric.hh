// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/per_residue_metrics/PerResidueDensityFitMetric.hh
/// @brief A per-residue metric that will calculate the density fit for each residue using either a correlation or a zscore.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_simple_metrics_per_residue_metrics_PerResidueDensityFitMetric_HH
#define INCLUDED_core_simple_metrics_per_residue_metrics_PerResidueDensityFitMetric_HH

#include <core/simple_metrics/per_residue_metrics/PerResidueDensityFitMetric.fwd.hh>
#include <core/simple_metrics/PerResidueRealMetric.hh>

// Core headers
#include <core/types.hh>

// Utility headers
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers
#include <map>

namespace core {
namespace simple_metrics {
namespace per_residue_metrics {

/// @brief A per-residue metric that will calculate the density fit for each residue using either a correlation or a zscore.
///
/// @details Zscore Uses weighted sum of density, density-compared-to-neighbors, rama (where applicable) and cart_bonded to compute
///
class PerResidueDensityFitMetric : public core::simple_metrics::PerResidueRealMetric{

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	PerResidueDensityFitMetric();

	/// @brief Copy constructor (not needed unless you need deep copies)
	PerResidueDensityFitMetric( PerResidueDensityFitMetric const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~PerResidueDensityFitMetric() override;

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
	std::map< core::Size, core::Real >
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

public:

	///@brief Use density correlation instead of a zscore to fit to density
	/// Default False.  Useful for comparison to other densities of different resolution.
	void
	set_match_mode( bool match_mode );

	///@brief Use the selector as true mask to calculate the Zscore.  Otherwise, use it just as a selection for result.  Default true.
	void
	set_use_selector_as_zscore_mask( bool selector_as_mask );


	///@brief Set the sliding window size.  Default is 3.  If you have anything other than protein, you probably want 1.
	/// Only used for Z-score mode.
	void
	set_sliding_window_size( core::Size window_size );

	///@brief Set the code to use a mixed sliding window depending on the residue.  3 for protein, 1 for everything else.
	/// Only used for Z-score mode.
	void
	set_mixed_sliding_window( bool mixed_sliding_window );

private:

	///@brief Compute All per-residue scores
	void
	compute_scores(
		pose::Pose & pose,

		std::map< Size, Real > & per_rsd_dens,
		std::map< Size, Real > & per_rsd_nbrdens,
		std::map< Size, Real > & per_rsd_rama,
		std::map< Size, Real > & per_rsd_geometry
	) const;


private:

	Size sliding_window_size_ = 1;
	bool mixed_sliding_window_ = false;
	bool match_res_ = false;
	bool use_selector_as_zscore_mask_ = true;
	pose::PoseCOP rs_native_ = nullptr; //Only used in RS.
};

} //core
} //simple_metrics
} //per_residue_metrics



#endif //core_simple_metrics_per_residue_metrics_PerResidueDensityFitMetric_HH





