// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/movers/FoldTreeFromFoldGraphMover.hh
/// @brief Creates and sets a new fold tree for the pose by traversing a FoldGraph
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_denovo_design_movers_FoldTreeFromFoldGraphMover_hh
#define INCLUDED_protocols_denovo_design_movers_FoldTreeFromFoldGraphMover_hh

// Unit headers
#include <protocols/denovo_design/movers/FoldTreeFromFoldGraphMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/denovo_design/components/FoldGraph.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/loops/Loops.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>

namespace protocols {
namespace denovo_design {
namespace movers {

struct CutAndJump {
public:
	CutAndJump( core::Size const cut_resid, int const jump_id ):
		cutpoint( cut_resid ), jump( jump_id ) {}

	core::Size cutpoint;
	int jump;

private:
	CutAndJump();
};
typedef utility::vector1< CutAndJump > JumpInfo;

///@brief Creates and sets a new fold tree for the pose by traversing a FoldGraph
class FoldTreeFromFoldGraphMover : public protocols::moves::Mover {
public:

	FoldTreeFromFoldGraphMover();

	/// @brief if loops is not NULL, this will create jumps/cuts as appropriate from the
	///        loops object
	FoldTreeFromFoldGraphMover(
		SegmentNames const & roots,
		protocols::loops::Loops const & loops );

	// destructor (important for properly forward-declaring smart-pointer members)
	virtual ~FoldTreeFromFoldGraphMover();

	static std::string
	class_name();

public:
	// mover virtual API
	virtual void
	apply( core::pose::Pose & pose );

	virtual void
	show( std::ostream & output = std::cout ) const;

	virtual std::string
	get_name() const;

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	virtual void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose );

	//FoldTreeFromFoldGraphMover & operator=( FoldTreeFromFoldGraphMover const & src );

	/// @brief required in the context of the parser/scripting scheme
	virtual protocols::moves::MoverOP
	fresh_instance() const;

	/// @brief required in the context of the parser/scripting scheme
	virtual protocols::moves::MoverOP
	clone() const;

public:
	void
	set_roots( SegmentNames const & roots );

private:
	/// @brief sets terminal variants for broken-chain folding using remodel
	///        this method also sets the last_jump_info_
	void
	prepare_termini_for_remodel( core::pose::Pose & pose );

private:
	SegmentNames roots_;
	protocols::loops::Loops loops_;
	JumpInfo last_jump_info_;
};

std::ostream &
operator<<( std::ostream & os, FoldTreeFromFoldGraphMover const & mover );

} //protocols
} //denovo_design
} //movers

#endif //protocols/denovo_design/movers_FoldTreeFromFoldGraphMover_hh
