// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/per_residue_metrics/HbondMetric.hh
/// @brief A metric to report the total h-bonds of a residue, or from a set of residues to another set of residues.  Use the SummaryMetric to get total hbonds of a selection or between selections. See the WaterMediatedBridgedHBondMetric for water-mediated h-bonds.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_simple_metrics_per_residue_metrics_HbondMetric_HH
#define INCLUDED_core_simple_metrics_per_residue_metrics_HbondMetric_HH

#include <core/simple_metrics/per_residue_metrics/HbondMetric.fwd.hh>
#include <core/simple_metrics/PerResidueRealMetric.hh>

// Core headers
#include <core/types.hh>
#include <core/scoring/hbonds/HBondSet.fwd.hh>

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

///@brief A metric to report the total h-bonds residues from a selection to all [OTHER] residues, or from a set of residues to another set of residues.  If No selection is given, will report ALL vs ALL.
///
///Tips:
///  Use the SummaryMetric to get total hbonds of a selection or total number of residues having some number of hbonds. . See the WaterMediatedBridgedHBondMetric for water-mediated h-bonds.
///
///  By default does not report self-self hbonds
///
class HbondMetric : public core::simple_metrics::PerResidueRealMetric{

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	HbondMetric();

	/// @brief Copy constructor (not needed unless you need deep copies)
	HbondMetric( HbondMetric const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~HbondMetric() override;


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

	///@brief Set to include hbonds as so:
	///  resi - resi
	/// Default false.
	///
	/// Will only find these if selector1 and selector2 contain the same residue.
	void
	set_include_self( bool include_self);

	///@brief Optionally set a second residue selector to get bridged hbonds between
	///  bothselections.
	void
	set_residue_selector2( core::select::residue_selector::ResidueSelectorCOP selector);



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

	bool include_self_ = false;
	core::select::residue_selector::ResidueSelectorCOP selector_two_ = nullptr;


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
CEREAL_FORCE_DYNAMIC_INIT( core_simple_metrics_per_residue_metrics_HbondMetric )
#endif // SERIALIZATION

#endif //core_simple_metrics_per_residue_metrics_HbondMetric_HH





