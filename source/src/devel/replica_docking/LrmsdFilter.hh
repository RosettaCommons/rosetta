// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Zhe Zhang

#ifndef INCLUDED_devel_replica_docking_LrmsdFilter_hh
#define INCLUDED_devel_replica_docking_LrmsdFilter_hh

#include <devel/replica_docking/LrmsdFilter.fwd.hh>

#include <protocols/filters/Filter.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/docking/types.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

#include <utility/tag/Tag.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
//#include <utility/vector1.hh>


namespace devel {
namespace replica_docking {

class LrmsdFilter : public protocols::filters::Filter
{
public:
	LrmsdFilter( core::Size const rb_jump=1, core::Real const lower_threshold=0, core::Real const upper_threshold=9999 );
	bool apply( core::pose::Pose const & pose ) const;
	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	core::Real report_sm( core::pose::Pose const & pose ) const;
	core::Real compute( core::pose::Pose const & pose ) const; // which residue numbers are neighbors
	protocols::filters::FilterOP clone() const;
	protocols::filters::FilterOP fresh_instance() const;

	void register_options();

	virtual ~LrmsdFilter();
	void jump( core::Size const jump_id );

	void parse_my_tag( utility::tag::TagCOP const tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
private:
	core::Real lower_threshold_;
	core::Real upper_threshold_;
	//  utility::vector1< core::Size > movable_jumps_;
 	protocols::docking::DockJumps movable_jumps_;
	core::pose::PoseCOP native_pose_;
	bool options_registered_;
};

}
}

#endif
