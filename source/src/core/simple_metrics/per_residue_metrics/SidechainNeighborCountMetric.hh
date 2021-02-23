// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/per_residue_metrics/SidechainNeighborCountMetric.hh
/// @brief A metric for calculating each sidechains neighbors based on cones.  This metric uses the same core code as the LayerSelector.  It can be combined with the SimpleMetricSelector to select any set of residues based on burial.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_simple_metrics_per_residue_metrics_SidechainNeighborCountMetric_HH
#define INCLUDED_core_simple_metrics_per_residue_metrics_SidechainNeighborCountMetric_HH

#include <core/simple_metrics/per_residue_metrics/SidechainNeighborCountMetric.fwd.hh>
#include <core/simple_metrics/PerResidueRealMetric.hh>

// Core headers
#include <core/types.hh>

// Utility headers
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers
#include <map>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace simple_metrics {
namespace per_residue_metrics {

///@brief A metric for calculating each sidechains neighbors based on cones.  This metric uses the same core code as the LayerSelector.  It can be combined with the SimpleMetricSelector to select any set of residues based on burial.
class SidechainNeighborCountMetric : public core::simple_metrics::PerResidueRealMetric{

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	SidechainNeighborCountMetric();

	/// @brief Copy constructor (not needed unless you need deep copies)
	SidechainNeighborCountMetric( SidechainNeighborCountMetric const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~SidechainNeighborCountMetric() override;


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

	/// @brief Set the exponent for the angle term, which affects how close other atoms have to be to the CA-CB line to be counted fully.
	/// @details Residues in the cone are counted and the count is multiplied by the cosine of the angle between the CA-CB vector and the CA-other atom
	/// vector.  The shift factor is then added, and the resulting value is raised to the angle_exponent (default 2.0) and multiplied by the distance factor.
	void
	set_angle_exponent( core::Real angle_exponent );

	/// @brief Set the shift factor for the angle term.
	/// @details Residues in the cone are counted and the count is multiplied by the cosine of the angle between the CA-CB vector and the CA-other atom
	/// vector.  The shift factor (default 0.5) is then added, and the resulting value is raised to the angle_exponent and multiplied by the distance factor.
	void
	set_angle_shift_factor( core::Real angle_shift_factor );

	/// @brief Set the exponent for the distance term, which affects how sharp the falloff is with distance.
	/// @details The distance term is: 1/(1+exp(n*(d - m))), where d is the distance, n is the exponent set by this term,
	/// and m is the midpoint of the falloff.  The n value sets the sharpness.  Defaults to 1.0.
	void
	set_dist_exponent( core::Real dist_exponent);

	/// @brief Midpoint of the distance falloff sigmoid.
	/// @details Defaults to 9.0.  Only used by the sidchain_neighbors code.
	void
	set_dist_midpoint( core::Real dist_midpoint );

	/// @brief Factor by which number of residue neighbors is divided.
	/// @details Defaults to 1.0.
	void
	set_res_denominator( core::Real res_denom );

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

	core::Real angle_exponent_ = 2.0;
	core::Real angle_shift_factor_ = 0.5;
	core::Real dist_exponent_ = 1.0;
	core::Real dist_midpoint_ = 9.0;
	core::Real res_denominator_ = 1.0;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} //per_residue_metrics
} //simple_metrics
} //core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_simple_metrics_per_residue_metrics_SidechainNeighborCountMetric )
#endif // SERIALIZATION

#endif //core_simple_metrics_per_residue_metrics_SidechainNeighborCountMetric_HH





