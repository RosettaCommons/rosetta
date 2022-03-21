// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rosetta_scripts/ParsedProtocol.hh
/// @author Sarel Fleishman (sarelf@u.washington.edu)
/// @author Vikram K. Mulligan (vmullig@uw.edu) -- Modified this to facilitate use of ParsedProtocols to combine movers and filters in code outside of a RosettaScripts context.

#ifndef INCLUDED_protocols_rosetta_scripts_ParsedProtocol_HH
#define INCLUDED_protocols_rosetta_scripts_ParsedProtocol_HH

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/jd2/JobOutputterObserver.hh>

#include <core/types.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <core/simple_metrics/SimpleMetric.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <utility/exit.hh>
#include <protocols/moves/ResId.hh>

// C++ headers
#include <string>

// Unit headers
#include <protocols/rosetta_scripts/ParsedProtocol.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace rosetta_scripts {

class ParsedProtocol :
	public protocols::moves::Mover,
	public protocols::moves::ResId,
	public protocols::jd2::JobOutputterObserver
{
public:
	typedef core::Real Real;
	typedef core::pose::Pose Pose;

	//@brief This enum is used by ParsedProtocolStep to label filter behavior for rerunning checks
	enum class FilterReportTime
	{
		AT_END,      //rerun the filter at the end of ParsedProtocol
		AFTER_APPLY, //rerun the filter immediately after applying the filter
		NONE         //never rerun the filter after filter::apply
	};

	/// @brief Represents a step in the ParsedProtocol
	/// Note that one or more of mover/filter/metrics may be null/empty.
	struct ParsedProtocolStep {
		ParsedProtocolStep();
		ParsedProtocolStep( ParsedProtocolStep const & ) = default;

		/// @brief Assignment operator.
		inline ParsedProtocolStep &
		operator=( ParsedProtocolStep const & ) = default;

		ParsedProtocolStep (
			moves::MoverOP mover_in,
			std::string const & mover_name,
			filters::FilterOP filter_in = nullptr,
			FilterReportTime frt = FilterReportTime::AT_END
		);

		protocols::moves::MoverOP mover;
		std::string mover_user_name;
		protocols::filters::FilterOP filter;
		utility::vector1< core::simple_metrics::SimpleMetricCOP > metrics;
		utility::vector1< std::string > metric_labels;

		//Interpret report_time_ to decide when to report filters
		bool report_at_end() const;
		bool report_after_apply() const;

	private:
		//Making report_time private so that commandline options can intervene
		FilterReportTime report_time_;
		bool never_rerun_filters_; //local cache of basic::options::option[ basic::options::OptionKeys::parser::never_rerun_filters ]()
	};

	typedef utility::vector1< ParsedProtocolStep > ParsedProtocolStepVector;
	typedef ParsedProtocolStepVector::iterator iterator;
	typedef ParsedProtocolStepVector::const_iterator const_iterator;

public:
	ParsedProtocol();
	~ParsedProtocol() override;
	void apply( Pose & pose ) override;
	core::pose::PoseOP get_additional_output( ) override;
	void final_scorefxn( core::scoring::ScoreFunctionCOP scorefxn );
	core::scoring::ScoreFunctionCOP final_scorefxn() const;
	void final_score(core::pose::Pose & pose) const;
	void report_all( Pose const & pose ) const; // cycles over all filter->report methods to output their values to a common stream.
	void report_filters_to_pose( Pose & pose ) const; // as above but reports to pose DataCache

	// Called directly from JobOutputter via Observer pattern
	void add_values_to_job( Pose const & pose, protocols::jd2::Job & ) const override;


	// void report_all_sm( std::map< std::string, core::Real > & score_map, Pose const & pose ) const; // ditto, but outputs filter values into score_map object
	protocols::moves::MoverCOP get_mover( core::Size const mover_number ) const {
		runtime_assert( steps_.size() >= mover_number && mover_number > 0 );
		return( steps_[ mover_number ].mover );
	}
	ParsedProtocolStep get_step( core::Size const step_number ) const {
		runtime_assert( steps_.size() >= step_number && step_number > 0 );
		return( steps_[ step_number ] );
	}

	/// @brief Add a mover-filter pair.
	/// @details Indended for use OUTSIDE of a RosettaScripts context.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void add_step(
		protocols::moves::MoverOP mover,
		std::string const &mover_name,
		protocols::filters::FilterOP filter,
		bool const report_filter_at_end=false
	);

	/// @brief Add a step to the protocol
	/// @details Indended for use OUTSIDE of a RosettaScripts context.
	void add_step(
		ParsedProtocolStep const & step
	);

	void set_resid( core::Size const resid ) override;
	void set_resid( core::pose::ResidueIndexDescriptionCOP r ) override;


	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override { return utility::pointer::make_shared< ParsedProtocol >(); }
	void parse_my_tag( utility::tag::TagCOP, basic::datacache::DataMap & ) override; // this is defined as public here, b/c I need to circumvent the name-check, since this is called both by the Movers section (as ParsedProtocol) and the PROTOCOLS section.

	void clear() { steps_.clear(); }
	std::string mode() const{ return mode_; }
	iterator begin();
	const_iterator begin() const;
	iterator end();
	const_iterator end() const;
	void apply_probability( utility::vector1< core::Real > const & a );
	utility::vector1< core::Real > apply_probability();
	core::Size size() { return steps_.size(); }
	core::Size last_attempted_mover_idx() { return last_attempted_mover_idx_; }
	void last_attempted_mover_idx( core::Size const s ){ last_attempted_mover_idx_ = s;}
	bool report_call_order() const { return report_call_order_; }
	void report_call_order( bool const c ) { report_call_order_ = c; }
	std::string call_order() const{ return call_order_; }

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	/// @brief Provide the citation.
	void provide_citation_info(basic::citation_manager::CitationCollectionList & ) const override;

	core::Size n_steps_passed_in_previous_run() const {
		return n_steps_passed_in_previous_run_;
	}

	///@brief Get the status of any filters after running.
	bool
	get_filter_status() const { return filter_status_;}

private:
	void finish_protocol(Pose & pose);

	/// @brief apply the components of the step
	/// @details Returns false on failure
	bool apply_step(Pose & pose, ParsedProtocolStep const & step, bool skip_mover = false);

	/// @brief apply the mover of the step
	bool apply_mover(Pose & pose, ParsedProtocolStep const & step);

	/// @brief apply the filter of the step
	/// @details Returns false on failure
	bool apply_filter(Pose & pose, ParsedProtocolStep const & step);

	/// @brief apply the metric of the step
	void apply_metrics(Pose & pose, ParsedProtocolStep const & step);

	void sequence_protocol(Pose & pose, ParsedProtocolStepVector::const_iterator mover_it_in);

	void random_order_protocol(Pose & pose);
	void random_single_protocol(Pose & pose);

private:

	ParsedProtocolStepVector steps_;
	core::scoring::ScoreFunctionCOP final_scorefxn_ = nullptr;
	std::string mode_;
	utility::vector1< core::Real > apply_probability_; // if mode_="single_random", assigns a probability of execution to each mover/filter pair. Defaults to equal probabilities to all.
	core::Size last_attempted_mover_idx_; //index to last attempted mover; useful for adaptive monte carlo
	bool report_call_order_; //dflt false; At the end of the run, write to out the sequence of mover/filter calls (good for stochastic application
	std::string call_order_; // saved call order, not writeable
	protocols::moves::MoverOP last_mover_;
	bool resume_support_;

	core::Size n_steps_passed_in_previous_run_ = 0;
	bool filter_status_ = true; //JAB - used by MultiStageRosettaScripts to properly fail if filter failed.
};

} // rosetta_scripts
} // protocols

#endif //INCLUDED_protocols_rosetta_scripts_ParsedProtocol_HH
