// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/denovo_design/movers/SealFoldTreeMover.hh
/// @brief Creates a sealed foldtree, and removes all cutpoints from the pose
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_denovo_design_movers_SealFoldTreeMover_hh
#define INCLUDED_protocols_denovo_design_movers_SealFoldTreeMover_hh

// Unit headers
#include <protocols/denovo_design/movers/SealFoldTreeMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/denovo_design/components/FoldGraph.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>

namespace protocols {
namespace denovo_design {
namespace movers {

///@brief Creates a sealed foldtree, and removes all cutpoints from the pose
class SealFoldTreeMover : public protocols::moves::Mover {
public:
	typedef components::Cutpoints Cutpoints;

public:

	SealFoldTreeMover();

	SealFoldTreeMover( components::StructureData const & sd, protocols::loops::Loops const & loops );

	// destructor (important for properly forward-declaring smart-pointer members)
	virtual ~SealFoldTreeMover();

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

	//SealFoldTreeMover & operator=( SealFoldTreeMover const & src );

	/// @brief required in the context of the parser/scripting scheme
	virtual protocols::moves::MoverOP
	fresh_instance() const;

	/// @brief required in the context of the parser/scripting scheme
	virtual protocols::moves::MoverOP
	clone() const;

private:
	void
	remove_cutpoints( core::pose::Pose & pose ) const;

private:
	Cutpoints cutpoints_;
	components::FoldGraphCOP fg_;
	SegmentNames roots_;

};

std::ostream &
operator<<( std::ostream & os, SealFoldTreeMover const & mover );

} //protocols
} //denovo_design
} //movers

#endif //protocols/denovo_design/movers_SealFoldTreeMover_hh
