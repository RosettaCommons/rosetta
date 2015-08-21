// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/AngleToVectorFilter.hh

#ifndef INCLUDED_protocols_simple_filters_AngleToVectorFilter_hh
#define INCLUDED_protocols_simple_filters_AngleToVectorFilter_hh

//unit headers
#include <protocols/simple_filters/AngleToVectorFilter.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

namespace protocols {
namespace simple_filters {

class AngleToVector : public filters::Filter
{
public:
	AngleToVector();
	virtual ~AngleToVector();
	filters::FilterOP clone() const {
		return filters::FilterOP( new AngleToVector( *this ) );
	}
	filters::FilterOP fresh_instance() const{
		return filters::FilterOP( new AngleToVector() );
	}

	virtual bool apply( core::pose::Pose const & pose ) const;
	virtual core::Real compute( core::pose::Pose const & pose ) const;
	virtual void report( std::ostream & out, core::pose::Pose const & pose ) const;
	virtual core::Real report_sm( core::pose::Pose const & pose ) const;
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &filters, moves::Movers_map const &, core::pose::Pose const & );

	core::Size chain() const{ return chain_; }
	void chain( core::Size const r ){ chain_ = r; }
	core::Real min_angle() const{ return min_angle_; }
	void min_angle( core::Real const r ){ min_angle_ = r; }
	core::Real max_angle() const{ return max_angle_; }
	void max_angle( core::Real const r ){ max_angle_ = r; }

	core::Real refx() const{ return refx_; }
	void refx( core::Real const r ){ refx_ = r; }
	core::Real refy() const{ return refy_; }
	void refy( core::Real const r ){ refy_ = r; }
	core::Real refz() const{ return refz_; }
	void refz( core::Real const r ){ refz_ = r; }

	std::string atm1() const{ return atm1_; }
	void atm1( std::string const a ){ atm1_ = a; }
	std::string atm2() const{ return atm2_; }
	void atm2( std::string const a ){ atm2_ = a; }

private:
	core::Real min_angle_, max_angle_; //dflt 0, +90
	core::Real refx_,refy_,refz_;// a reference vector for calculation
	core::Size chain_;//dflt 2; which chain to compute
	std::string atm1_, atm2_; //atoms for which to compute the vectors
};

}
}

#endif
