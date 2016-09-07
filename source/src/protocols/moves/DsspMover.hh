// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ./protocols/moves/DsspMover.hh
/// @brief  header file of DsspMover.cc
/// @author Nobuyasu Koga ( nobuyasau@uw.edu )

#ifndef INCLUDED_protocols_moves_DsspMover_hh
#define INCLUDED_protocols_moves_DsspMover_hh

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

#include <utility/vector1.hh>

#ifdef WIN32
#include <utility/tag/Tag.hh>
#endif


namespace protocols {
namespace moves {


class DsspMover : public Mover {
public:

	typedef protocols::moves::MoverOP MoverOP;
	typedef utility::tag::TagCOP TagCOP;
	typedef protocols::filters::Filters_map Filters_map;
	typedef basic::datacache::DataMap DataMap;
	typedef protocols::moves::Movers_map Movers_map;

public:

	// default constructor
	DsspMover();

	// @brief destructor
	~DsspMover() override;

	/// @brief clone this object
	MoverOP clone() const override;

	/// @brief create this type of object
	MoverOP fresh_instance() const override;

	// @brief virtual main operation
	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;

	void parse_my_tag( TagCOP, basic::datacache::DataMap &, Filters_map const &, Movers_map const &, Pose const & ) override;

private:
	bool reduced_IG_as_L_;


};


} // moves
} // protocols


#endif
