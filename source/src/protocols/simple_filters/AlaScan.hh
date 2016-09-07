// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/AlaScan.hh
/// @brief definition of filter class AlsScan.
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

#ifndef INCLUDED_protocols_simple_filters_AlaScan_hh
#define INCLUDED_protocols_simple_filters_AlaScan_hh

#include <protocols/simple_filters/AlaScan.fwd.hh>


// Project Headers
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <utility/exit.hh>

#include <utility/vector1.hh>

namespace protocols {
namespace simple_filters {

class AlaScan : public filters::Filter
{
public :
	AlaScan();
	AlaScan( bool const chain1, bool const chain2, core::Size const repeats, core::Real const dist, core::scoring::ScoreFunctionCOP scorefxn, core::Size const jump, bool const symmetry );
	bool apply( core::pose::Pose const & ) const override{ return true; }
	filters::FilterOP clone() const override {
		return filters::FilterOP( new AlaScan( *this ) );
	}
	filters::FilterOP fresh_instance() const override{
		return filters::FilterOP( new AlaScan() );
	}

	void report( std::ostream & out, core::pose::Pose const & pose ) const override;
	void report_symmetry( std::ostream & out, core::pose::Pose const & pose ) const;
	core::Real report_sm( core::pose::Pose const & ) const override { return (0.0); };
	void chain1( bool const c1 ){ chain1_ = c1; }
	void chain2( bool const c2 ){ chain2_ = c2; }
	void dist( core::Real const d ){ distance_threshold_ = d; }
	void scorefxn( core::scoring::ScoreFunctionOP scorefxn );
	core::Real ddG_for_single_residue( core::pose::Pose const & pose, core::Size const resi ) const;
	bool chain1() const{ return chain1_; }
	bool chain2() const{ return chain2_; }
	core::Size repeats() const { return repeats_; }
	void repeats( core::Size const r ) { repeats_ = r; }
	core::Size jump() const { return jump_; }
	void jump( core::Size const j ) { jump_ = j; }
	core::Real dist() const{ return distance_threshold_; }
	~AlaScan() override;
	void repack( bool const repack );
	bool repack() const;
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;
private:
	bool chain1_, chain2_;
	core::Size repeats_;
	core::Real distance_threshold_;
	core::Size jump_;
	bool symmetry_;
	core::scoring::ScoreFunctionOP scorefxn_;
	bool repack_; //dflt true; do you want to repack the partners in the bound and unbound states?
};


}
}
#endif
