// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/denovo_design/filters/SSPredictionFilter.hh
/// @brief header file for filter to determine agreement with psipred for secondary structure prediction
/// @details
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_denovo_design_filters_sspredictionfilter_hh
#define INCLUDED_protocols_denovo_design_filters_sspredictionfilter_hh

// Unit Headers
#include <protocols/denovo_design/filters/SSPredictionFilter.fwd.hh>

// Package headers

// Project headers
#include <core/io/external/PsiPredInterface.fwd.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/parser/BluePrint.fwd.hh>
#include <protocols/ss_prediction/SS_predictor.fwd.hh>

namespace protocols {
namespace denovo_design {
namespace filters {

// main SSPredictionFilter class
class SSPredictionFilter : public protocols::filters::Filter {
public:
	SSPredictionFilter();
	SSPredictionFilter( core::Real const threshold,
		std::string const & cmd,
		std::string const & blueprint_filename,
		bool const use_probability,
		bool const mismatch_probability,
		bool const use_scratch_dir = false );
	~SSPredictionFilter() override;
	bool apply( core::pose::Pose const & pose ) const override;
	protocols::filters::FilterOP clone() const override;
	protocols::filters::FilterOP fresh_instance() const override;
	void report( std::ostream & out, core::pose::Pose const & pose ) const override;
	core::Real report_sm( core::pose::Pose const & pose ) const override;
	core::Real compute( core::pose::Pose const & pose ) const;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data_map
	) override;

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	/// @brief Sets the secondary structure to be compared with the psipred result. If neither this nor the blueprint is set, DSSP will be used to determine the pose secondary structure
	void
	set_secstruct( std::string const secstruct );

	/// @brief Sets the psipred executable (has no effect if "use_svm" is true)
	void
	set_cmd( std::string const & cmd );

	/// @brief Sets a blueprint object which can be used to determine the desired secondary structure
	void
	set_blueprint( protocols::parser::BluePrintOP blueprint );

	/// @brief Sets a blueprint file which can be used to determine the desired secondary structure
	void
	set_blueprint_file( std::string const & blueprint_file );

	/// @brief Sets the threshold. If use_probablility is true, the filter will pass if the calculated value is <= threshold.  If use_probability is false, the filter will pass if the calcualted value is >= threshold.
	void
	set_threshold( core::Real const threshold );

	/// @brief Sets whether or not to use the predicted SS probabilities to compute a value, or to use the predicted secondary structure (match or not matching) for each residue. Values are combined into a boltzmann sum-like-value using the provided temperature, so that poorly predicted residues have a larger impact on the score
	void
	set_use_probability( bool const use_probability );

	/// @brief assumes use_probability, if set this will compute the cumulative probability of having correct secondary structure at all residues
	void
	set_mismatch_probability( bool const mismatch_probability );

	/// @brief tells whether to use the psipred pass2 confidence values -- overrrides use_probability.
	void
	set_use_confidence( bool const use_confidence );

	/// @brief Sets the temperature which is used to compute probability-based values
	void
	set_temperature( core::Real const temp );

	/// @brief If set, a trained SVM will be used instead of psipred to predict the secondary structure
	void
	set_use_svm( bool const use_svm );

	/// @brief Sets the psipred interface object that will be used to call psipred
	void
	set_psipred_interface( core::io::external::PsiPredInterfaceOP psipred ) { psipred_interface_ = psipred; }

private:
	/// @brief computes the weighted boltzmann sum of the passed vector
	core::Real compute_boltz_sum( utility::vector1< core::Real > const & probabilities ) const;

private:
	/// @brief computes one minus the geometric mean of the passed vector
	core::Real compute_mismatch_prob( utility::vector1< core::Real > const & probabilities ) const;

	/// @brief set the scratch dir to out::path::scratch or error if that's not set
	void set_scratch_dir();


private:
	core::Real threshold_;
	std::string cmd_;
	protocols::parser::BluePrintOP blueprint_;
	// tells whether to just count residues that don't match, or use the probability generated by psipred
	bool use_probability_;
	// assumes use_probability, tells whether to use probabilities directly or in a (weird) weighted way
	bool mismatch_probability_;
	/// @brief tells whether to use the psipred pass2 confidence values -- overrrides use_probability.
	bool use_confidence_;
	/// @brief should we use svm to estimate secondary structure?
	bool use_svm_;
	/// @brief what temperature should we be using doing a boltzmann sum?
	core::Real temp_;
	/// @brief the object which predicts the secondary structure
	protocols::ss_prediction::SS_predictorOP ss_predictor_;
	/// @brief the object which communicates with psipred and interprets its output
	core::io::external::PsiPredInterfaceOP psipred_interface_;
	std::string secstruct_;
	/// @brief the scratch folder to use with psipred (or "")
	std::string scratch_dir_;
};  //SSPredictionFilter

//namespaces
}
}
}

#endif
