// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/ResiduesInInterfaceFilter.hh
/// @brief Reports to Tracer which residues are designable in a taskfactory
/// @author Sarel Fleishman (sarelf@uw.edu)

#ifndef INCLUDED_protocols_simple_filters_ResiduesInInterfaceFilter_hh
#define INCLUDED_protocols_simple_filters_ResiduesInInterfaceFilter_hh

#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

namespace protocols {
namespace simple_filters{

class ResiduesInInterfaceFilter : public filters::Filter
{
public:
	ResiduesInInterfaceFilter() : Filter( "ResInInterface" ) {}
	ResiduesInInterfaceFilter( core::Size const residues_in_interface_threshold, core::Size const rb_jump ) : Filter( "ResInInterface" ) {
		residues_in_interface_threshold_ = residues_in_interface_threshold;
		rb_jump_ = rb_jump;
	}
	bool apply( core::pose::Pose const & pose ) const;
	filters::FilterOP clone() const {
		return filters::FilterOP( new ResiduesInInterfaceFilter( *this ) );
	}
	filters::FilterOP fresh_instance() const{
		return filters::FilterOP( new ResiduesInInterfaceFilter() );
	}
  	ResiduesInInterfaceFilter( ResiduesInInterfaceFilter const & init ) :
	//utility::pointer::ReferenceCount(),
	Filter( init ) {
    residues_in_interface_threshold_ = init.residues_in_interface_threshold_;
		rb_jump_ = init.rb_jump_;
  	}
	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	core::Real report_sm( core::pose::Pose const & pose ) const;
	core::Size compute( core::pose::Pose const & pose ) const;
	virtual ~ResiduesInInterfaceFilter();
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
private:
	core::Size residues_in_interface_threshold_, rb_jump_;
};

}
}

#endif
