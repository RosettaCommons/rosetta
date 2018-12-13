// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/per_residue_metrics/PerResidueClashMetric.hh
/// @brief A SimpleMetric that calculates the number of atomic clashes per residue using the LJ radius (at 0).
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_simple_metrics_per_residue_metrics_PerResidueClashMetric_HH
#define INCLUDED_core_simple_metrics_per_residue_metrics_PerResidueClashMetric_HH

#include <core/simple_metrics/per_residue_metrics/PerResidueClashMetric.fwd.hh>
#include <core/simple_metrics/PerResidueRealMetric.hh>

// Core headers
#include <core/types.hh>

// Utility headers
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers
#include <map>
#include <cmath>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace simple_metrics {
namespace per_residue_metrics {

///@brief A SimpleMetric that calculates the total number of atom-atom clashes from a residue in a residue selector to all other
/// residues defined in a second residue selector using the LJ radius of each atom.
///
/// Can use a soft radius, which reduces it by 33%.
///
///@details
/// Does NOT calculate INTRA-RESIDUE clashes!
/// Default is to calculate clashes of each residue to ALL OTHER residues (including those in the primary selector)
///  Set a secondary selector if you want other behavior
///
/// We use a SCORED pose to only calculate energies if fa_rep between a pair of residues under consideration is NOT 0
///   This is in order to speed up the calculation.
///
/// INTRA-residue clashes can be done through the PUBLIC is_clashing method to calculate atomic-clashes.
///
class PerResidueClashMetric : public core::simple_metrics::PerResidueRealMetric{

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	PerResidueClashMetric();

	///@brief Constructor with other residue selector set.  The other selector is what we use to calculate clashes TO
	/// Default is to calculate clashes of each residue to ALL OTHER residues (including those in the primary selector)
	PerResidueClashMetric(
		select::residue_selector::ResidueSelectorCOP inner_selector,
		select::residue_selector::ResidueSelectorCOP outer_selector );

	/// @brief Copy constructor (not needed unless you need deep copies)
	PerResidueClashMetric( PerResidueClashMetric const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~PerResidueClashMetric() override;


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

	///@brief Set the Residue Selector to use for TO other residues.
	///
	void
	set_secondary_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector );

	///@brief Should we calculate only heavy-heavy atom clashes
	/// Default True
	void
	set_use_hydrogens( bool use_hydrogens );

	///@brief When we calculate atom-atom distances using LJ distances,
	/// Default True
	///
	///@details
	///  clash if distance < (atomI_LJ + atomJ_LJ)*(1 - soft_clash)
	///
	void
	set_use_soft_clash( bool soft_clash_check );

	///@brief Set the dampening of the LJ to use for soft-clash.
	/// Default=.33
	void
	set_soft_dampening( core::Real dampening );

public:

	///@brief Are these atoms clashing using soft or hard definition?
	bool
	is_clashing( core::pose::Pose const & pose, core::Size resA, core::Size atomA, core::Size resB, core::Size atomB ) const ;

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

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

private:

	///@brief Return the radius at which the energy is 0 (IE where energy goes from attraction to repulsion).
	core::Real
	lj_radius_to_zero_e_radius( core::Real num) const;

private:

	bool use_soft_clash_ = true;
	bool use_hydrogens_ = false;
	core::Real soft_clash_dampening_ = .33;
	core::select::residue_selector::ResidueSelectorCOP selector2_ = nullptr;
	core::Real lj_n_ = pow( 2, 1.0/6 );


};

} //core
} //simple_metrics
} //per_residue_metrics

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_simple_metrics_per_residue_metrics_PerResidueClashMetric )
#endif // SERIALIZATION

#endif //core_simple_metrics_per_residue_metrics_PerResidueClashMetric_HH





