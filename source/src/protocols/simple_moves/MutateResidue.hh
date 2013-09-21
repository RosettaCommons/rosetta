// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file
/// @brief
/// @author Spencer Bliven <blivens@u.washington.edu>
/// @date 6/26/2009

#ifndef INCLUDED_protocols_simple_moves_MutateResidue_hh
#define INCLUDED_protocols_simple_moves_MutateResidue_hh

#include <protocols/simple_moves/MutateResidue.fwd.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
// AUTO-REMOVED #include <core/chemical/AA.hh>

//parsing
#include <utility/tag/Tag.fwd.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh> //Movers_map
#include <protocols/filters/Filter.fwd.hh> //Filters_map

#include <utility/vector1.hh>


// Utility headers

// C++ headers

// Unit headers

namespace protocols {
namespace simple_moves {

/// @brief A mover to mutate a single residue
class MutateResidue : public protocols::moves::Mover
{
private:
	typedef protocols::moves::Mover parent;
public:
	///@brief default ctor
	MutateResidue();
	///@brief copy ctor
	MutateResidue(MutateResidue const& dm);
	///@brief Mutate a single residue to a new amino acid
	MutateResidue( core::Size const target, std::string const new_res );
	MutateResidue( core::Size const target, int const new_res/*one letter code*/);  // Changing char --> int so PyRosetta could use overloaded function
	virtual ~MutateResidue() {};

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
	virtual protocols::moves::MoverOP clone() const {
		return (protocols::moves::MoverOP( new protocols::simple_moves::MutateResidue( *this ) ) );
	}
	virtual protocols::moves::MoverOP fresh_instance() const {
		return protocols::moves::MoverOP( new MutateResidue );
	}

	void parse_my_tag( utility::tag::TagPtr const tag,
		protocols::moves::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & );

  virtual void parse_def( utility::lua::LuaObject const & def,
		utility::lua::LuaObject const & score_fxns,
		utility::lua::LuaObject const & tasks,
		protocols::moves::MoverCacheSP cache );
private:
	core::Size rb_jump_;
	core::Size target_;
	std::string res_name_;
};

} // moves
} // protocols

#endif //INCLUDED_protocols_simple_moves_MutateResidue_HH_
