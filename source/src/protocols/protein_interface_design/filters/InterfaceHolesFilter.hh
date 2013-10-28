// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/filters/InterfaceHolesFilter.hh
/// @brief Find packing defects at an interface using packstat score terms
/// @author Jacob Corn (jecorn@u.washington.edu)

#ifndef INCLUDED_protocols_protein_interface_design_filters_InterfaceHolesFilter_hh
#define INCLUDED_protocols_protein_interface_design_filters_InterfaceHolesFilter_hh


#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace protein_interface_design {
namespace filters {

class InterfaceHolesFilter : public protocols::filters::Filter
{
public:
	InterfaceHolesFilter() : protocols::filters::Filter( "InterfaceHoles" ) {}
	InterfaceHolesFilter( core::Size const rb_jump, core::Real const threshold ) : protocols::filters::Filter( "InterfaceHoles" ) {
		rb_jump_ = rb_jump;
		threshold_ = threshold;
	}
	bool apply( core::pose::Pose const & pose ) const;
	protocols::filters::FilterOP clone() const {
		return new InterfaceHolesFilter( *this );
	}
	protocols::filters::FilterOP fresh_instance() const{
		return new InterfaceHolesFilter();
	}
  InterfaceHolesFilter( InterfaceHolesFilter const & init ) : 
	//utility::pointer::ReferenceCount(), 
	protocols::filters::Filter( init ) 
{
		rb_jump_ = init.rb_jump_;
		threshold_ = init.threshold_;
  };
	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	core::Real report_sm( core::pose::Pose const & pose ) const;
	core::Real compute( core::pose::Pose const & pose ) const;
	virtual ~InterfaceHolesFilter();
	void parse_my_tag( utility::tag::TagCOP const tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
private:
	core::Size rb_jump_;
	core::Real threshold_;
};
} // filters
} // protein_interface_design
} // protocols

#endif //INCLUDED_protocols_protein_interface_design_filters_InterfaceHolesFilter_HH_

