// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/helical_bundle/BundleReporterFilter.hh
/// @brief Header files for a filter that is mainly intended for reporting helical bundle parameters and an associated energy.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_protocols_helical_bundle_BundleReporterFilter_hh
#define INCLUDED_protocols_helical_bundle_BundleReporterFilter_hh

//unit headers
#include <protocols/helical_bundle/BundleReporterFilter.fwd.hh>

// Project Headers
#include <core/scoring/ScoreFunction.hh>
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/scoring/ScoreType.hh>

//Auto Headers


namespace protocols {
namespace helical_bundle {

/// @brief The behaviour of this filter.
/// @details The filter can return true always, false always, or can filter by energy.
/// If you update this enum, add new behaviours to the two set_filter_behaviour()
/// functions (setting by string and by enum).
enum BundleReporterFilterBehaviour {
	brf_always_true = 1,
	brf_always_false,
	brf_filter_by_energy,
	brf_undefined //keep this last
};

class BundleReporterFilter : public filters::Filter
{
public:

	/// @brief Constructor
	///
	BundleReporterFilter();
	
	/// @brief Copy constructor
	///
	BundleReporterFilter( BundleReporterFilter const &src );

	/// @brief Destructor.
	///
	virtual ~BundleReporterFilter();

	/// @brief Returns an owning pointer to a new instance of this filter, with copied
	/// variables (a copy of this filter).
	filters::FilterOP clone() const {
		return filters::FilterOP( new BundleReporterFilter( *this ) );
	}
	
	/// @brief Returns an owning pointer to a new instance of this filter, with default
	/// initialization (NOT a copy).
	filters::FilterOP fresh_instance() const{
		return filters::FilterOP( new BundleReporterFilter() );
	}
	
	/// @brief Set the score threshold.
	///
	void set_score_type_threshold( core::Real const &val ) { score_type_threshold_ = val; return; }

	/// @brief Get the score threshold.
	///
	inline core::Real score_type_threshold( ) const { return score_type_threshold_; }
	
	/// @brief Set the score type.
	///
	void set_score_type( core::scoring::ScoreType const &type ) { score_type_=type; return; }

	/// @brief Get the score type.
	///
	inline core::scoring::ScoreType score_type( ) const { return score_type_; }

	/// @brief Set the scorefunction.
	///
	void set_scorefxn( core::scoring::ScoreFunctionOP scorefxn_in )
	{
		runtime_assert_string_msg( scorefxn_in, "Error in protocols::helical_bundle::BundleReporterFilter::set_scorefxn(): NULL pointer passed to this function.  Consult a developer -- this is a programming error, most likely." );
		scorefxn_ = scorefxn_in;
		return;
	}
	
	/// @brief Get the scorefunction (by owning pointer).
	/// @details Non-const access.
	inline core::scoring::ScoreFunctionOP scorefxn() { return scorefxn_; }
	
	/// @brief Set the behaviour by string.
	/// @details Options are "ALWAYS_TRUE", "ALWAYS_FALSE", or "FILTER".
	void set_filter_behaviour( std::string const &behaviour_string );

	/// @brief Set the behaviour by enum.
	/// @details Options are brf_always_true, brf_always_false, or brf_filter_by_energy.
	void set_filter_behaviour( BundleReporterFilterBehaviour const behaviour_in ) {
		runtime_assert_string_msg( behaviour_in > 0 && behaviour_in < brf_undefined,
			"Error in protocols::helical_bundle::BundleReporterFilter::set_filter_behaviour(): The filter can only be set to always true (\"ALWAYS_TRUE\"), always false (\"ALWAYS_FALSE\"), or filter by energy  (\"FILTER\")." );
		behaviour_ = behaviour_in;
		return;
	}
	
	/// @brief Get the filter behaviour.
	///
	inline BundleReporterFilterBehaviour filter_behaviour() const { return behaviour_; }

	/// @brief Actually apply the filter to a pose.
	/// @details This scores the pose with the scorefunction, writes out the energy and the
	/// bundle parameters (if any) to the REPORT tracer, then applies the selected behaviour
	/// of the filter (always true, always false, or actually filtering by score).
	bool apply( core::pose::Pose const & pose ) const;

  /// @brief Allows reporting of filter values to a stream.
  ///
	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	
	/// @brief Allows reporting of the filter value to a float.
	///
	core::Real report_sm( core::pose::Pose const & pose ) const;

	/// @brief Generates the text of the full report that the filter writes out to the REPORT tracer.
	/// @details Called by the APPLY function.  Jobno is the RosettaScripts job number; ignored if set to 0.
	std::string generate_full_tracer_report(
		core::Size const jobno,
		core::Real const &score,
		core::pose::Pose const &pose
	) const; 

	/// @brief Computes the energy of the pose.
	/// @details the energy function must be suitable for residue type set of the pose,
	/// and must be symmetric if this is a symmetric pose.
	core::Real compute( core::pose::Pose const &pose ) const;
	
	/// @brief Parse XML (RosettaScripts) setup.
	///
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );

	/// @brief Set whether we report sequences.
	///
	void set_report_sequence( bool const val=true ) { report_sequence_=val; return; }

	/// @brief Set whether we report three-letter codes.
	///
	void set_use_threeletter( bool const val=true ) { report_three_letter_codes_=val; return; }
	
	/// @brief Get whether we report sequences.
	///
	inline bool report_sequence( ) const { return report_sequence_; }

	/// @brief Get whether we report three-letter codes.
	///
	inline bool use_threeletter( ) const { return report_three_letter_codes_; }


private:

	/// @brief If the score is used to filter, what's the cutoff threshold above which the filter returns false?
	/// @details 0.0 by default.
	core::Real score_type_threshold_;
	
	/// @brief What score term should be reported/used to filter?
	/// @details Set to total_score by default.
	core::scoring::ScoreType score_type_;
	
	/// @brief What scorefunction should be reported/used to filter?
	/// @details Must be specified by the user or in the code prior to invoking the filter; defaults to a null pointer.
	core::scoring::ScoreFunctionOP scorefxn_;
	
	/// @brief The behaviour of this filter.
	/// @details Defaults to always returning true, but can be set to always return false or to filter by energy.
	BundleReporterFilterBehaviour behaviour_;
	
	/// @brief Report sequence?  False by default.
	bool report_sequence_;
	
	/// @brief Report one- or three-letter codes?  False (one-letter codes) by default.
	bool report_three_letter_codes_;

};

} //helical_bundle
} //protocols

#endif //INCLUDED_protocols_helical_bundle_BundleReporterFilter_hh
