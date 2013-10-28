// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file
/// @brief

#ifndef INCLUDED_protocols_simple_moves_SetTorsion_hh
#define INCLUDED_protocols_simple_moves_SetTorsion_hh

#include <protocols/simple_moves/SetTorsion.fwd.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
// AUTO-REMOVED #include <core/chemical/AA.hh>

//parsing
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh> //Movers_map
#include <protocols/filters/Filter.fwd.hh> //Filters_map

#include <utility/vector1.hh>


// Utility headers

// C++ headers

// Unit headers

namespace protocols {
namespace simple_moves {

/// @brief A mover to change one torsion angle
class SetTorsion : public protocols::moves::Mover
{
private:
	typedef protocols::moves::Mover parent;
public:
	///@brief default ctor
	SetTorsion();
	virtual ~SetTorsion();

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
	virtual protocols::moves::MoverOP clone() const {
		return (protocols::moves::MoverOP( new protocols::simple_moves::SetTorsion( *this ) ) );
	}
	virtual protocols::moves::MoverOP fresh_instance() const {
		return protocols::moves::MoverOP( new SetTorsion );
	}

	void parse_my_tag( utility::tag::TagCOP const tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & );
	core::Real angle() const { return angle_; }
	core::Size resnum() const { return resnum_; }
	std::string torsion_name() const { return torsion_name_; }
	void angle( core::Real const a ){ angle_ = a; }
	void resnum( core::Size const r ){ resnum_ = r; }
  void torsion_name( std::string const s ){ torsion_name_ = s; }

private:
	core::Real angle_;
	core::Size resnum_;
	std::string torsion_name_; // phi/psi etc.
};

} // moves
} // protocols

#endif //INCLUDED_protocols_simple_moves_SetTorsion_HH_

