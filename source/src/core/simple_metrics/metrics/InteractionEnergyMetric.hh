// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/metrics/InteractionEnergyMetric.hh
/// @brief Calculate the interaction energy between residues in two ResidueSelectors in a single pose, including long-range interactions.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
/// @author Original Logic from InterfaceDeltaEnergetics: John Karanicolas && Roland A. Pache


#ifndef INCLUDED_core_simple_metrics_metrics_InteractionEnergyMetric_HH
#define INCLUDED_core_simple_metrics_metrics_InteractionEnergyMetric_HH

#include <core/simple_metrics/metrics/InteractionEnergyMetric.fwd.hh>
#include <core/simple_metrics/RealMetric.hh>

// Core headers
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Utility headers
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace simple_metrics {
namespace metrics {

///@brief Calculate the interaction energy between residues in two ResidueSelectors in a single pose, including long-range interactions.
/// Pose should be scored!
///
/// By default, we skip rama_prepro and pro_close
///
///@details
///
/// Scoretypes from REF:
///
///   fa_atr
///   fa_rep
///   fa_sol
///   fa_elec
///   lk_ball_wtd
///   pro_close #ignored by default
///   rama_prepro #ignored by default
///   hbond_sr_bb
///   hbond_lr_bb
///   hbond_bb_sc
///   hbond_sc
///
class InteractionEnergyMetric : public core::simple_metrics::RealMetric{

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	InteractionEnergyMetric();

	InteractionEnergyMetric(
		select::residue_selector::ResidueSelectorCOP selector1,
		select::residue_selector::ResidueSelectorCOP selector2
	);

	/// @brief Copy constructor (not needed unless you need deep copies)
	InteractionEnergyMetric( InteractionEnergyMetric const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~InteractionEnergyMetric() override;

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

	///@brief Calculate the metric.
	core::Real
	calculate( core::pose::Pose const & pose ) const override;

public:

	///@brief Set the residue selectors we will use to calculate the interaction energies.
	void
	set_residue_selectors(
		select::residue_selector::ResidueSelectorCOP selector1,
		select::residue_selector::ResidueSelectorCOP selector2);

	///@brief Set a scorefunction.  Only used if the pose is not scored.  Not recommended, this is our failsafe!
	void
	set_scorefunction( scoring::ScoreFunctionCOP scorefxn );

	///@brief Set this to include only these score types.
	void
	set_include_only_scoretypes( utility::vector1< scoring::ScoreType > const & include_only_scoretypes );

	///@brief Always ignore these scoretypes.
	void
	set_ignore_scoretypes( utility::vector1< scoring::ScoreType > const & ignore_scoretypes );

	///@brief Include Rama PrePro and ProClose?
	/// Default False
	///
	void
	set_include_rama_prepro_and_proclose( bool include );


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
	select::residue_selector::ResidueSelectorCOP selector1_ = nullptr;
	select::residue_selector::ResidueSelectorCOP selector2_ = nullptr;
	scoring::ScoreFunctionCOP scorefxn_ = nullptr;

	utility::vector1<scoring::ScoreType> score_types_to_ignore_;
	utility::vector1<scoring::ScoreType> score_types_to_use_only_;
	bool include_rama_prepro_and_proclose_ = false;
	bool force_rescore_pose_ = false;
};

} //core
} //simple_metrics
} //metrics

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_simple_metrics_metrics_InteractionEnergyMetric )
#endif // SERIALIZATION

#endif //core_simple_metrics_metrics_InteractionEnergyMetric_HH





