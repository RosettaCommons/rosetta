// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
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
#include <protocols/rosetta_scripts/MultiplePoseMover.hh>

#include <core/types.hh>

#include <core/pack/task/residue_selector/ResidueSelector.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>

// Utility headers
#include <utility/vector1.hh>

namespace protocols {
namespace match {

class MatcherMover : public protocols::rosetta_scripts::MultiplePoseMover {
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

	void
	generate_match_pos( core::pose::Pose const & pose ) const;

	void
	clean_match_pos() const;

	/// @brief if set to true, a single random match will be returned.
	/// @details The default behavior is to use the MultiplePoseMover framework
	/// to return all matches
	void
	set_return_single_random_match( bool const single_random );

protected:
	virtual bool process_pose( core::pose::Pose &, utility::vector1 < core::pose::PoseOP > & );

private:

	//dictates whether matches will be output to disk
	//or one of them will be incorporated into the pose
	bool incorporate_matches_into_pose_;
	bool return_single_random_match_;

	core::conformation::ResidueCOP ligres_;
	utility::vector1< core::Size > match_positions_;
	utility::vector1< core::pack::task::residue_selector::ResidueSelectorCOP > selectors_;

	static std::string const MATCH_POS_FILE;
};

void
set_ligpose_rotamer( core::pose::Pose & ligpose );

}
}

#endif
