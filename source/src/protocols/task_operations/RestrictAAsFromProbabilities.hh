// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/task_operations/RestrictAAsFromProbabilities
/// @brief A class to restrict designable amino acids depending on the probabilities from a PerResidueProbabilitiesMetric
/// @author Moritz Ertelt (moritz.ertelt@gmail.com)


#ifndef INCLUDED_protocols_task_operations_RestrictAAsFromProbabilities_hh
#define INCLUDED_protocols_task_operations_RestrictAAsFromProbabilities_hh

#include <protocols/task_operations/RestrictAAsFromProbabilities.fwd.hh>

// core headers
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/simple_metrics/PerResidueProbabilitiesMetric.hh>


// Utility headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace protocols {
namespace task_operations {

///@brief A class to restrict designable amino acids depending on the probabilities from a PerResidueProbabilitiesMetric
class RestrictAAsFromProbabilities: public core::pack::task::operation::TaskOperation {
public:

	RestrictAAsFromProbabilities();

	~RestrictAAsFromProbabilities() override;

	core::pack::task::operation::TaskOperationOP
	clone() const override;

	/// @brief Configure from a RosettaScripts XML tag.
	void
	parse_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & ) override;

	//////////////////////

	void
	apply(core::pose::Pose const & pose, core::pack::task::PackerTask & task) const override;

	/// @brief Return the name used to construct this TaskOperation from an XML file
	static std::string keyname();

	/// @brief Describe the format of XML file used to initialize this TaskOperation
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	///@brief Set the PerResidueProbabilitiesMetric that will be used to calculate the pseudo-perplexity.
	void
	set_metric( core::simple_metrics::PerResidueProbabilitiesMetricCOP metric );

	///@brief Set a boolean to attempt to find cached data matching the name/custom_type of the passed in simple_metric.
	/// Optionally pass any set prefix/suffix.
	///
	/// This will allow the filter to re-use previously calculated data.
	///
	void
	set_use_cached_data( bool use_cache, std::string const & prefix="", std::string const & suffix="");

	///@brief If use_cache is set to false, do we fail if no data is found in the pose?
	/// Default True
	void
	set_fail_on_missing_cache(bool fail);

	///@brief Set the probability cutoff option
	void set_prob_cutoff( core::Real prob_cutoff );

	///@brief Set the delta probability cutoff option
	void set_delta_prob_cutoff( core::Real delta_prob_cutoff );

	///@brief disable all amino acids that do not meet user provided probability cutoffs
	void disable_AAs(
		std::map<core::Size, std::map<core::chemical::AA, core::Real>>& probabilities,
		core::pose::Pose const& pose,
		core::pack::task::PackerTask& task
	) const;

private:

	// metric to be used
	core::simple_metrics::PerResidueProbabilitiesMetricCOP metric_ = nullptr;

	// probability cutoffs set by the user
	core::Real prob_cutoff_ = 0.0001;
	core::Real delta_prob_cutoff_ = 0.0;

	//Accessing cached data from the pose
	bool use_cache_ = false;
	std::string cache_prefix_;
	std::string cache_suffix_;
	bool fail_on_missing_cache_ = true;
};

} //task_operations
} //protocols

#endif //INCLUDED_protocols/task_operations_RestrictAAsFromProbabilities_fwd_hh
