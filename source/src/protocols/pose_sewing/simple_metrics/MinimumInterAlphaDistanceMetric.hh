// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pose_sewing/simple_metrics/MinimumInterAlphaDistanceMetric.hh
/// @brief min dist
/// @author Frank Teets (frankdt@email.unc.edu)

#ifndef INCLUDED_protocols_pose_sewing_simple_metrics_MinimumInterAlphaDistanceMetric_HH
#define INCLUDED_protocols_pose_sewing_simple_metrics_MinimumInterAlphaDistanceMetric_HH

#include <protocols/pose_sewing/simple_metrics/MinimumInterAlphaDistanceMetric.fwd.hh>
#include <core/simple_metrics/RealMetric.hh>

// Core headers
#include <core/types.hh>

// Utility headers
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION
#include <core/select/residue_selector/util.hh>
#include <core/select/residue_selector/ResidueSelector.hh>

namespace protocols {
namespace pose_sewing {
namespace simple_metrics {

/// @brief min dist
class MinimumInterAlphaDistanceMetric : public core::simple_metrics::RealMetric{

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	MinimumInterAlphaDistanceMetric();

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~MinimumInterAlphaDistanceMetric() override;


public:

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
	//apply( core::pose::Pose & pose, prefix="", suffix="" ) override;

	/// @brief Calculate the metric.
	core::Real
	calculate( core::pose::Pose const & pose ) const override;

	void
	set_selectors( core::select::residue_selector::ResidueSelectorCOP A_selector, core::select::residue_selector::ResidueSelectorCOP B_selector, core::select::residue_selector::ResidueSelectorCOP all_atom_selector);

public:

	/// @brief Name of the class
	std::string
	name() const override;

	/// @brief Name of the class for creator.
	static
	std::string
	name_static();

	/// @brief Name of the metric
	std::string
	metric() const override;

	/// @brief Provide authorship information for an unpublished Rosetta module.
	void
	provide_citation_info(basic::citation_manager::CitationCollectionList & citations) const override;


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

	core::Real sequence_gap_ = 6;
	core::select::residue_selector::ResidueSelectorCOP A_selector_;
	core::select::residue_selector::ResidueSelectorCOP B_selector_;
	core::select::residue_selector::ResidueSelectorCOP all_atom_selector_;
};

} //simple_metrics
} //pose_sewing
} //protocols


#endif //protocols_pose_sewing_simple_metrics_MinimumInterAlphaDistanceMetric_HH





