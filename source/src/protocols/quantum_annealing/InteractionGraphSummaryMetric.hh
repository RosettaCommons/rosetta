// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/quantum_annealing/InteractionGraphSummaryMetric.hh
/// @brief A simple metric that allows an InteractionGraph to be written out in a format that external annealers (including quantum annealers) can read.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_protocols_quantum_annealing_InteractionGraphSummaryMetric_HH
#define INCLUDED_protocols_quantum_annealing_InteractionGraphSummaryMetric_HH

#include <protocols/quantum_annealing/InteractionGraphSummaryMetric.fwd.hh>
#include <core/simple_metrics/StringMetric.hh>

// Core headers
#include <core/types.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Utility headers
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace quantum_annealing {

///@brief A simple metric that allows an InteractionGraph to be written out in a format that external annealers (including quantum annealers) can read.
class InteractionGraphSummaryMetric : public core::simple_metrics::StringMetric{

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	InteractionGraphSummaryMetric();

	/// @brief Copy constructor (not needed unless you need deep copies)
	InteractionGraphSummaryMetric( InteractionGraphSummaryMetric const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~InteractionGraphSummaryMetric() override;

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

	///@brief Calculate the metric.
	std::string
	calculate( core::pose::Pose const & pose ) const override;

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

	/// @brief Set the task factory to be used for interaxction graph setup.
	/// @details Doesn't clone input.
	void set_task_factory( core::pack::task::TaskFactoryCOP task_factory_in );

	/// @brief Set the scorefunction used for scoring.
	/// @details Copied, not cloned.
	void set_scorefunction( core::scoring::ScoreFunctionCOP sfxn_in );

	/// @brief Set whether we're storing a short summary, or a long one.
	/// @details If true, then only the interaction graph summary is printed.  If false, then the
	/// interaction graph summary plus the full information needed to reconstruct the pose is included.
	/// @details False by default.
	inline void set_short_version( bool const setting ) { short_version_ = setting; }

	/// @brief Get whether we're storing a short summary, or a long one.
	/// @details If true, then only the interaction graph summary is printed.  If false, then the
	/// interaction graph summary plus the full information needed to reconstruct the pose is included.
	/// @details False by default.
	inline bool short_version() const { return short_version_; }

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

private: //Variables

	/// @brief The task factory that we'll use to set up the packing.  We won't
	/// actually DO the packing, mind you.  We'll just use this to set up the
	/// interaction graph.
	core::pack::task::TaskFactoryCOP task_factory_;

	/// @brief The scorefunction that we'll use.
	core::scoring::ScoreFunctionCOP scorefxn_;

	/// @brief If true, then only the interaction graph summary is printed.  If false, then the
	/// interaction graph summary plus the full information needed to reconstruct the pose is included.
	/// @details False by default.
	bool short_version_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION


};

} //protocols
} //quantum_annealing

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_quantum_annealing_InteractionGraphSummaryMetric )
#endif // SERIALIZATION


#endif //protocols_quantum_annealing_InteractionGraphSummaryMetric_HH

