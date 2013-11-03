// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/simple_moves/sidechain_moves/SetChiMover.hh
/// @brief  A mover to change one chi angle
/// @author Noah Ollikanen

#ifndef INCLUDED_protocols_simple_moves_sidechain_moves_SetChiMover_hh
#define INCLUDED_protocols_simple_moves_sidechain_moves_SetChiMover_hh

#include <protocols/simple_moves/sidechain_moves/SetChiMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

//parsing
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh> //Movers_map
#include <protocols/filters/Filter.fwd.hh> //Filters_map

namespace protocols {
namespace simple_moves {
namespace sidechain_moves {

/// @brief A mover to change one chi angle
class SetChiMover : public protocols::moves::Mover
{
private:
	typedef protocols::moves::Mover parent;
public:
	///@brief default ctor
	SetChiMover();
	virtual ~SetChiMover();

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
	virtual protocols::moves::MoverOP clone() const {
		return (protocols::moves::MoverOP( new protocols::simple_moves::sidechain_moves::SetChiMover( *this ) ) );
	}
	virtual protocols::moves::MoverOP fresh_instance() const {
		return protocols::moves::MoverOP( new SetChiMover );
	}

	void parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & );
	core::Real angle() const { return angle_; }
	core::Size resnum() const { return resnum_; }
	core::Size chinum() const { return chinum_; }
	void angle( core::Real const a ){ angle_ = a; }
	void resnum( core::Size const r ){ resnum_ = r; }
	void chinum( core::Size const c ){ chinum_ = c; }

private:
	core::Real angle_;
	core::Size resnum_;
	core::Size chinum_;
};

} // sidechain_moves
} // simple_moves
} // protocols

#endif //INCLUDED_protocols_simple_moves_sidechain_moves_SetChiMover_HH_
