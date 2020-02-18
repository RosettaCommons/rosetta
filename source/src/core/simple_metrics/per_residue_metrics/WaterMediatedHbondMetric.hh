// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/per_residue_metrics/WaterMediatedHbondMetric.hh
/// @brief A metric to measure hydrogen bonds between a set of residues that are water-mediated.  Depth of 1 is default where one water mediates the interface between these residues.  Depth can be set.  Make sure to use the -include_waters flag to have Rosetta not ignore HOH

/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_simple_metrics_per_residue_metrics_WaterMediatedHbondMetric_HH
#define INCLUDED_core_simple_metrics_per_residue_metrics_WaterMediatedHbondMetric_HH

#include <core/simple_metrics/per_residue_metrics/WaterMediatedHbondMetric.fwd.hh>
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


///@brief A metric to measure hydrogen bonds between a set of residues that are water-mediated.
///
///@details
/// DEPTH:
///    We only traverse a depth of 1 by default. Make sure to set the -ignore_waters flag to false in order to have Rosetta include the HOH residues.
///
/// SELECTION:
///   If one residue selector is given, will calculate bridged waters between residues of the selection and all [OTHER] residues, otherwise it will calculate bridges between one selection and another.
///
///   If NO SELECTION is give, will report ALL bridged hbonds in a pose [by default without self-water-self hbonds].
///
/// HBOND:
///  Since these are bridged hbonds, and h-bond networks can be rather complex, the numbers
///   reported here are the unique h-bond paths from sele1 to sele2.  If you give only a single residue selector, these are bridged hbonds from sele1 to OTHER residues in sele1.
///   By default we do not include mediated hbonds back on itself, but this is an option.
///
/// TIPS:
///  It is generally recommended to repack the waters using the OptH TaskOperation before input into this metric, especially if only Oxygens were present.  Using the option -include_vrt false will keep all waters present in the resulting structure.  Use the option -corrections::water::wat_rot_sampling 10 to decrease the angle of sampling from a default of 30 to 10.  This will result in many more rotamers, but will improve networks.
///
///  Hydration shells can be calculated by passing selection1 as all waters in the pose, selection2 as the protein or chain, and then reporting only a single depth (with max_depth at 0 and 1 for the first and second shell waters."
///
///  It is also recommended to use the Beta scorefunction as this can improve hydrogen bond detection.
class WaterMediatedHbondMetric : public core::simple_metrics::PerResidueRealMetric{

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	WaterMediatedHbondMetric();

	/// @brief Copy constructor (not needed unless you need deep copies)
	WaterMediatedHbondMetric( WaterMediatedHbondMetric const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~WaterMediatedHbondMetric() override;


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


	///@brief After packing, some waters become HOH_v if they don't match certain ref energy.
	/// The option -include_vrt false will have them all be HOH during packing.
	/// Set this option to true to incude them in our calculation if they are still virts.
	void
	set_include_virt_waters(bool include_virt_waters);

	///@brief Optionally set a second residue selector to get bridged hbonds between
	///  both selections.
	void
	set_residue_selector2( core::select::residue_selector::ResidueSelectorCOP selector);

	///@brief Set the depth we will go to find hbonds.
	/// Default is 1. IE resi - water -resj
	///
	void
	set_depth( core::Size const depth);

	///@brief Set to include only the depth that is set instead of Up to the max depth.
	/// IE - if depth is set to 2, we only include hbond paths mediated by two waters.
	///  Default false
	void
	set_include_only_set_depth( bool include_only_depth);

	///@brief Set to include hbonds as so:
	///  resi - water - resi
	/// Default false.
	///
	/// Will only find these if selector1 and selector2 contain the same residue.
	void
	set_include_self( bool include_self);


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

	///@brief Recursive function to find all unique HBond paths involving water and ending
	///  with a residue that is not in our water list.
	///
	///@details Populates paths as the reference.
	/// Path is a splittalbe string of the connections.
	///  ex: 191-221-333 where residue 221 is a water in our list.
	///
	void
	find_hb_paths(
		core::scoring::hbonds::HBondSet const & local_hb_set,
		utility::vector1< core::Size > const & local_waters,
		std::map<std::string, core::Size> & paths,
		core::Size const current_res,
		core::Size const max_depth=1,
		std::string const & current_path = "",
		core::Size const current_depth=0,
		core::scoring::hbonds::HBondCOP prev_hb= nullptr) const;

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

	bool include_virt_waters_ = false;
	bool include_self_ = false;
	bool include_only_set_depth_ = false;

	core::Size depth_ = 1;

	core::select::residue_selector::ResidueSelectorCOP selector_two_ = nullptr;

	//bool at_depth = false;

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
CEREAL_FORCE_DYNAMIC_INIT( core_simple_metrics_per_residue_metrics_WaterMediatedHbondMetric )
#endif // SERIALIZATION

#endif //core_simple_metrics_per_residue_metrics_WaterMediatedHbondMetric_HH





