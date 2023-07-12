// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/simple_metrics/PerResidueProbabilitiesMetric.hh
///
/// @brief  Base class for core::Real PerResidueProbabilitiesMetrics
/// @author Moritz Ertelt (moritz.ertelt@gmail.com)

#ifndef INCLUDED_core_simple_metrics_PerResidueProbabilitiesMetric_hh
#define INCLUDED_core_simple_metrics_PerResidueProbabilitiesMetric_hh

// Unit headers
#include <core/simple_metrics/SimpleMetric.hh>
#include <core/simple_metrics/PerResidueProbabilitiesMetric.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Core headers
#include <core/chemical/AA.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ includes
#include <map>
#include "core/sequence/SequenceProfile.hh"

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace simple_metrics {

/// @brief A class that is derived to calculate probabilities for different residue types for each Residue.
///  Apply(pose) method calculates this metric and adds it to the pose score for output.
///
///
/// @details
///  Calculate(pose) method calculates core::Real values and returns them as a map<Size, map< AA, Real>>
///   Resnum:(Restype:Value)
///
class PerResidueProbabilitiesMetric : public core::simple_metrics::SimpleMetric {

public: // constructors / destructors

	PerResidueProbabilitiesMetric();

	~PerResidueProbabilitiesMetric() override;

	PerResidueProbabilitiesMetric( PerResidueProbabilitiesMetric const & other );

	PerResidueProbabilitiesMetric & operator=( PerResidueProbabilitiesMetric const & );

	using SimpleMetric::apply; // Add back the base class implementation

	/// @brief Calculate the metric and add it to the pose as a score.
	///           labeled as out_label
	///
	/// @details Score is added to the SimpleMetricData cache in the pose
	///            A ReferencePose is created with out_label as a name for further access.
	///            Data is output to the final scorefile.
	void
	apply(
		std::string const & out_label,
		pose::Pose & pose,
		bool override_existing_data = false) const override;

	///@brief Set a ResidueSelector for which we will calculate values over.
	void
	set_residue_selector( select::residue_selector::ResidueSelectorCOP selector );

	///@brief Set to output in PDB numbering instead of Rosetta during the Apply function,
	/// which adds the data to pose as extra scores.
	void
	set_output_as_pdb_nums( bool output_as_pdb_nums );

	///@brief Calculate the metric.
	/// This map is Rosetta Resnum->(ResType->value) and includes only those residues selected.
	///
	///@details
	/// Return by value as this function can not STORE the result, it only calculates.
	/// Store the result in the pose by using the apply method, which calls this method and stores the result
	/// in the pose as ExtraScoreValues.
	virtual std::map< core::Size, std::map< core::chemical::AA, core::Real >>
	calculate( pose::Pose const & pose ) const = 0;

	///@brief Grab the data from the pose if it exists or calculate the metric
	///
	///@details If use_cache is true, we will attempt to pull the data from the pose.
	/// If fail_on_missing_cache is true, we will fail, otherwise, we will calculate the metric.
	///
	/// This function is meant to support caching metrics, so values do not need to be calculated twice,
	///  for example in SimpleMetricFilter/Features
	///  or code-wise where data takes a while to calculate and can be reused.
	///
	/// If we cached the data, we have created a ref-pose and can match the current resnums with our refpose resnums
	///  using the use_ref_pose_for_cache option.
	///  This allows us to delete residues and still retain the given data to match.
	///
	std::map< core::Size, std::map< core::chemical::AA, core::Real >>
	cached_calculate(
		pose::Pose const & pose,
		bool use_cache,
		std::string const & prefix="",
		std::string const & suffix="",
		bool fail_on_missing_cache=true,
		bool use_ref_pose_for_cache=true) const;

public:

	///@brief Name of the class
	std::string
	name() const override = 0;

	///@brief Name of the metric
	std::string
	metric() const override = 0;

	///@brief Get the submetric names that this Metric will calculate
	utility::vector1< std::string >
	get_metric_names() const override;

	///@brief Get the set residue selector of this class.
	select::residue_selector::ResidueSelectorCOP
	get_selector() const;

public:
	/// @brief called by parse_my_tag -- should not be used directly
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data ) override = 0;

	SimpleMetricOP
	clone() const override = 0;

	///@brief Add options to the schema from this base class.
	static
	void
	add_schema( utility::tag::XMLSchemaComplexTypeGeneratorOP complex_schema);

	///Parse the base class tag.  Keep required interface for parse_my_tag.
	virtual void
	parse_per_residue_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data );


private:

	select::residue_selector::ResidueSelectorCOP selector_; //Set as a TrueResidueSelector.
	bool output_as_pdb_nums_ = false;

	/// @brief Output the sequence_profile as file
	/// @param[in] profile A SequenceProfile filled with the logits
	/// @param[in] output_filename A string defining the name of the output file
	static void write_profile(core::sequence::SequenceProfile &profile, std::string const &output_filename) ;

protected:

	/// @brief Format the probabilities in psi-blast position-specific-scoring-matrix (PSSM) format and write to file
	/// @param[in] sequence The sequence of the pose
	/// @param[in] logit_map A map containing the predicted logits for each position
	/// @param[in] output_filename A string defining the name of the output file
	static void output_sequence_profile(std::string const & sequence, std::map<core::Size, std::map<core::chemical::AA, core::Real>> const &logit_map,
		std::string const &output_filename) ;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION


}; //class PerResidueProbabilitiesMetrics

} //namespace simple_metrics
} //namespace core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_simple_metrics_PerResidueProbabilitiesMetric )
#endif // SERIALIZATION

#endif

