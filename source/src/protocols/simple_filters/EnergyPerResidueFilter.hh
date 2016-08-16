// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/EnergyPerResidueFilter.hh
/// @brief definition of filter class EnergyPerResidueFilter.
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

#ifndef INCLUDED_protocols_simple_filters_EnergyPerResidueFilter_hh
#define INCLUDED_protocols_simple_filters_EnergyPerResidueFilter_hh

#include <protocols/simple_filters/EnergyPerResidueFilter.fwd.hh>
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>

namespace protocols {
namespace simple_filters {

class EnergyPerResidueFilter : public filters::Filter
{
public:
	EnergyPerResidueFilter() : filters::Filter( "EnergyPerResidue" ) {}

	EnergyPerResidueFilter( core::Size const resnum, core::scoring::ScoreFunctionCOP scorefxn, core::scoring::ScoreType const score_type, core::Real const threshold, bool const whole_interface = false, bool const whole_protein = false, bool const select_resnums = false, bool const select_around_resnums = false, core::Real aournd_shell=8.0, std::string string_resnums="", std::string string_around_resnums="", core::Size const rb_jump = 1, core::Real const interface_distance_cutoff =  8.0 , bool const bb_bb = false );

	EnergyPerResidueFilter( EnergyPerResidueFilter const &init );
	void apply_helper( std::string name, core::pose::Pose const & pose, utility::vector1< bool > & use_all_residues ) const;
	bool apply( core::pose::Pose const & pose ) const;
	filters::FilterOP clone() const {
		return filters::FilterOP( new EnergyPerResidueFilter( *this ) );
	}
	filters::FilterOP fresh_instance() const{
		return filters::FilterOP( new EnergyPerResidueFilter() );
	}

	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	core::Real report_sm( core::pose::Pose const & pose ) const;
	core::Real compute( core::pose::Pose const & pose, core::Size const resid ) const;
	core::Real compute( core::pose::Pose const & pose ) const;
	virtual ~EnergyPerResidueFilter();
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
	core::Size resnum() const;
	core::scoring::ScoreFunctionOP scorefxn() const;
	core::scoring::ScoreType score_type() const;
	core::Real threshold() const;
	bool bb_bb() const;
	void resnum( core::Size const rn );
	void scorefxn( core::scoring::ScoreFunctionOP scorefxn );
	void score_type( core::scoring::ScoreType score_type );
	void threshold( core::Real const th );
	void bb_bb( bool const b_b );
private:
	core::Size resnum_;
	core::scoring::ScoreFunctionOP scorefxn_; //?
	core::scoring::ScoreType score_type_;
	core::Real threshold_;
	bool whole_interface_;
	bool whole_protein_;
	bool select_resnums_;
	bool select_around_resnums_;
	std::string string_resnums_;
	std::string string_around_resnums_;
	core::Real around_shell_;
	core::Size rb_jump_;
	core::Real interface_distance_cutoff_;
	bool bb_bb_;
};

}
}

#endif
