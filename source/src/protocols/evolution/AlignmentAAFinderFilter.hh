// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/evolution/AlignmentAAFinderFilter.hh

#ifndef INCLUDED_protocols_evolution_AlignmentAAFinderFilter_hh
#define INCLUDED_protocols_evolution_AlignmentAAFinderFilter_hh

//unit headers
#include <protocols/evolution/AlignmentAAFinderFilter.fwd.hh>

// Project Headers
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/vector1.hh>


typedef utility::vector1< bool > bools;


namespace protocols {
namespace evolution {

/// @brief test whether a pose contains a comment that evaluates to a predefined value. This is useful in controlling execution flow in RosettaScripts.
class AlignmentAAFinder : public filters::Filter
{
public:
	AlignmentAAFinder();
	~AlignmentAAFinder() override;
	filters::FilterOP clone() const override {
		return filters::FilterOP( new AlignmentAAFinder( *this ) );
	}
	filters::FilterOP fresh_instance() const override{
		return filters::FilterOP( new AlignmentAAFinder() );
	}

	bool apply( core::pose::Pose const & p ) const override;

	void report( std::ostream & out, core::pose::Pose const & pose ) const override;

	core::Real report_sm( core::pose::Pose const & pose ) const override;

	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &data, filters::Filters_map const &filters, moves::Movers_map const &movers, core::pose::Pose const & ) override;

	core::Real exclude_AA_threshold() const { return exclude_AA_threshold_; }
	void exclude_AA_threshold( core::Real const & t ) { exclude_AA_threshold_ = t; }

	std::string alignment_file() const { return alignment_file_; }
	void alignment_file( std::string const & a ) { alignment_file_ = a; }

	std::string available_AAs_file() const { return available_AAs_file_; }
	void available_AAs_file( std::string const & a ) { available_AAs_file_ = a; }

	core::Size indel_motif_radius() const { return indel_motif_radius_; }
	void indel_motif_radius( core::Size const & r ) { indel_motif_radius_ = r; }

	core::Real loop_seqid_threshold() const { return loop_seqid_threshold_; }
	void loop_seqid_threshold( core::Real const & t ) { loop_seqid_threshold_ = t; }

	core::scoring::ScoreFunctionOP scorefxn() const{ return scorefxn_; }
	void scorefxn( core::scoring::ScoreFunctionOP scorefxn ) { scorefxn_ = scorefxn; }

	protocols::moves::MoverOP relax_mover() const { return relax_mover_; };
	void relax_mover( protocols::moves::MoverOP mover ) { relax_mover_ = mover; };

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	core::Real exclude_AA_threshold_;
	std::string alignment_file_;
	std::string available_AAs_file_;
	core::Size indel_motif_radius_;
	core::Real loop_seqid_threshold_;
	core::scoring::ScoreFunctionOP scorefxn_;
	protocols::moves::MoverOP relax_mover_;

};
}
}

#endif
