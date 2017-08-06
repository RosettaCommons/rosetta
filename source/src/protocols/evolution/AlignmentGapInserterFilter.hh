// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/evolution/AlignmentGapInserterFilter.hh

#ifndef INCLUDED_protocols_evolution_AlignmentGapInserterFilter_hh
#define INCLUDED_protocols_evolution_AlignmentGapInserterFilter_hh

//unit headers
#include <protocols/evolution/AlignmentGapInserterFilter.fwd.hh>

// Project Headers
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector0.hh>


typedef utility::vector1< bool > bools;


namespace protocols {
namespace evolution {

/// @brief test whether a pose contains a comment that evaluates to a predefined value. This is useful in controlling execution flow in RosettaScripts.
class AlignmentGapInserter : public filters::Filter
{
public:
	AlignmentGapInserter();
	~AlignmentGapInserter() override;
	filters::FilterOP clone() const override {
		return filters::FilterOP( new AlignmentGapInserter( *this ) );
	}
	filters::FilterOP fresh_instance() const override{
		return filters::FilterOP( new AlignmentGapInserter() );
	}

	utility::vector0< core::Size > find_char_location_in_string(std::string const & string, char const findIt) const;

	bool apply( core::pose::Pose const & p ) const override;

	void report( std::ostream & out, core::pose::Pose const & pose ) const override;

	core::Real report_sm( core::pose::Pose const & pose ) const override;

	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &data, filters::Filters_map const &filters, moves::Movers_map const &, core::pose::Pose const & ) override;

	std::string alignment_file() const { return alignment_file_; }
	void alignment_file( std::string const & a ) { alignment_file_ = a; }

	std::string available_AAs_file() const { return available_AAs_file_; }
	void available_AAs_file( std::string const & a ) { available_AAs_file_ = a; }

	std::string cleaned_alignment_file() const { return cleaned_alignment_file_; }
	void cleaned_alignment_file( std::string const & a ) { cleaned_alignment_file_ = a; }

	core::Real nbr_e_threshold() const { return nbr_e_threshold_; }
	void nbr_e_threshold( core::Real const & t ) { nbr_e_threshold_ = t; }

	utility::vector1< core::Real > loop_seqid_thresholds() const { return loop_seqid_thresholds_; };
	void loop_seqid_thresholds( utility::vector1< core::Real > const & v ) { loop_seqid_thresholds_ = v;};

	core::Size indel_motif_radius() const { return indel_motif_radius_; }
	void indel_motif_radius( core::Size const & r ) { indel_motif_radius_ = r; }

	core::Size only_clean_seq_num() const { return only_clean_seq_num_; }
	void only_clean_seq_num( core::Size const & t ) { only_clean_seq_num_ = t; }

	core::scoring::ScoreFunctionOP scorefxn() const{ return scorefxn_; }
	void scorefxn( core::scoring::ScoreFunctionOP scorefxn ) { scorefxn_ = scorefxn; }

	utility::vector1< core::Real > max_score_diffs() const { return max_score_diffs_; };
	void max_score_diffs( utility::vector1< core::Real > const & v ) { max_score_diffs_ = v;};

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	std::string alignment_file_;
	std::string available_AAs_file_;
	std::string cleaned_alignment_file_;
	core::Real nbr_e_threshold_;
	utility::vector1< core::Real > loop_seqid_thresholds_;
	core::Size indel_motif_radius_;
	core::Real only_clean_seq_num_;
	core::scoring::ScoreFunctionOP scorefxn_;
	utility::vector1< core::Real > max_score_diffs_;

};
}
}

#endif
