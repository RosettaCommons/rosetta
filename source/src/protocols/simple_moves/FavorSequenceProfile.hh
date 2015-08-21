// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/FavorSequenceProfile.hh
/// @brief Add a SequenceProfileConstraint to a pose.
/// @author Rocco Moretti (rmoretti@u.washington.edu)

#ifndef INCLUDED_protocols_simple_moves_FavorSequenceProfile_hh
#define INCLUDED_protocols_simple_moves_FavorSequenceProfile_hh


#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/sequence/Sequence.fwd.hh>
#include <core/sequence/SequenceProfile.hh>

#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>

#include <utility/vector1.hh>

// #include <protocols/protein_interface_design/design_utils.hh>

namespace protocols {
namespace simple_moves {

class FavorSequenceProfile : public protocols::moves::Mover
{
public:
	FavorSequenceProfile();
	virtual std::string get_name() const { return "FavorSequenceProfile"; }
	protocols::moves::MoverOP clone() const { return( protocols::moves::MoverOP( new protocols::simple_moves::FavorSequenceProfile( *this ) ) ); }
	protocols::moves::MoverOP fresh_instance() const { return protocols::moves::MoverOP( new FavorSequenceProfile ); }
	void set_sequence( core::sequence::Sequence & seq, std::string matrix);
	/// @brief Set the profile object to use.
	/// Remember to set set_scaling() appropriately for the profile matrix you pass in.
	void set_profile( core::sequence::SequenceProfile & profile);
	/// @brief What type of manipulation to apply to the profile prior to using it for constraints.
	void set_scaling( std::string const & scaling );
	/// @brief What multiplication factor to apply to the profile prior to using it for constraints.
	void set_weight( core::Real weight );
	void apply( core::pose::Pose & pose );
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
	virtual void parse_def( utility::lua::LuaObject const & def,
		utility::lua::LuaObject const & score_fxns,
		utility::lua::LuaObject const & tasks,
		protocols::moves::MoverCacheSP cache );
	virtual ~FavorSequenceProfile() {};
private:
	core::Real weight_;
	bool use_current_;
	std::string matrix_;
	std::string scaling_;
	core::Size chain_;
	std::string string_exclude_resnums_;
	core::sequence::SequenceProfileOP ref_profile_;
};

} // moves
} // protocols


#endif /*INCLUDED_protocols_simple_moves_FavorSequenceProfile_hh*/
