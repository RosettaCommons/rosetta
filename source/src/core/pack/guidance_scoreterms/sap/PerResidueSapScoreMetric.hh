// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/guidance_scoreterms/sap/PerResidueSapScoreMetric.hh
/// @brief A per-residue metric that will calculate the SapScore for each residue.
/// @author Brian Coventry (bcov@uw.edu)

#ifndef INCLUDED_core_pack_guidance_scoreterms_sap_PerResidueSapScoreMetric_HH
#define INCLUDED_core_pack_guidance_scoreterms_sap_PerResidueSapScoreMetric_HH

#include <core/pack/guidance_scoreterms/sap/PerResidueSapScoreMetric.fwd.hh>
#include <core/simple_metrics/PerResidueRealMetric.hh>

// Core headers
#include <core/types.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Utility headers
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/pointer/memory.hh>

// C++ headers
#include <map>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace sap {

///@brief A per-residue metric that will calculate SASA for each residue given in a selector.
class PerResidueSapScoreMetric : public core::simple_metrics::PerResidueRealMetric{

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	PerResidueSapScoreMetric(
		core::select::residue_selector::ResidueSelectorCOP score_selector = utility::pointer::make_shared<core::select::residue_selector::TrueResidueSelector>(),
		core::select::residue_selector::ResidueSelectorCOP sap_calculate_selector = nullptr,
		core::select::residue_selector::ResidueSelectorCOP sasa_selector = nullptr
	);

	/// @brief Copy constructor (not needed unless you need deep copies)
	PerResidueSapScoreMetric( PerResidueSapScoreMetric const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~PerResidueSapScoreMetric() override;

	PerResidueSapScoreMetric & operator=( PerResidueSapScoreMetric const & ot );

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

	void
	set_score_selector( core::select::residue_selector::ResidueSelectorCOP const & selector );

	void
	set_sap_calculate_selector( core::select::residue_selector::ResidueSelectorCOP const & selector );

	void
	set_sasa_selector( core::select::residue_selector::ResidueSelectorCOP const & selector );


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

	core::select::residue_selector::ResidueSelectorCOP score_selector_;
	core::select::residue_selector::ResidueSelectorCOP sap_calculate_selector_;
	core::select::residue_selector::ResidueSelectorCOP sasa_selector_;


#ifdef    SERIALIZATION
protected:
	friend class cereal::access;

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} //sap
} //guidance_scoreterms
} //pack
} //core


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_pack_guidance_scoreterms_sap_PerResidueSapScoreMetric )
#endif // SERIALIZATION


#endif //core_pack_guidance_scoreterms_sap_PerResidueSapScoreMetric_HH





