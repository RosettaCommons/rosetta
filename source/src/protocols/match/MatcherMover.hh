// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/MatherMover.hh
/// @brief  mover wrapper for the matcher
/// @author Florian Richter, floric@u.washington.edu, june 2010

#ifndef INCLUDED_protocols_match_MatcherMover_hh
#define INCLUDED_protocols_match_MatcherMover_hh

// Unit headers
#include <protocols/match/MatcherMover.fwd.hh>

// Package headers

// Project headers
#include <protocols/moves/Mover.hh>

#include <core/types.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>

// Utility headers
#include <utility/vector1.hh>

namespace protocols {
namespace match {

class MatcherMover : public protocols::moves::Mover {
public:

	typedef core::Real Real;
	typedef core::Size Size;

	typedef protocols::moves::MoverOP MoverOP;

	typedef utility::tag::TagCOP TagCOP;
	typedef protocols::filters::Filters_map Filters_map;
	typedef basic::datacache::DataMap DataMap;
	typedef protocols::moves::Movers_map Movers_map;

public:


	/// Construction and Destruction
	MatcherMover( bool incorporate_matches_into_pose = true );
	virtual ~MatcherMover();

	/// @brief copy constructor
	MatcherMover( MatcherMover const & rval );


public: // virtual constructors


	/// @brief clone this object
	virtual
	MoverOP clone() const;


	/// @brief create this type of object
	virtual
	MoverOP fresh_instance() const;

	virtual
	void parse_my_tag( TagCOP tag,
										 basic::datacache::DataMap &,
									   Filters_map const &,
										 Movers_map const &,
										 Pose const & );


public:

	void
	apply(
		core::pose::Pose & pose );

	std::string
	get_name() const;

	void
	set_ligres(
		core::conformation::ResidueCOP ligres );

	void
	set_match_positions(
		utility::vector1< core::Size > const & match_positions );


private:

	//dictates whether matches will be output to disk
	//or one of them will be incorporated into the pose
	bool incorporate_matches_into_pose_;

	core::conformation::ResidueCOP ligres_;
	utility::vector1< core::Size > match_positions_;

};

}
}

#endif
