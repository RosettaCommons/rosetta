// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /protocols/simple_moves/LoadUnboundRotMover.hh
/// @brief
/// @author Steven Lewis smlewi@gmail.com, Brian Weitzner brian.weitzner@gmail.com

#ifndef INCLUDED_protocols_simple_moves_LoadUnboundRotMover_hh
#define INCLUDED_protocols_simple_moves_LoadUnboundRotMover_hh

// Unit Headers
#include <protocols/simple_moves/LoadUnboundRotMover.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>

// Utility Headers

namespace protocols {
namespace simple_moves {

/// @details This BS class exists to route around a hard-coded pseudo-Mover in the APPLY_TO_POSE section of the Parser.  At XRW2 we needed to move junk out of Parser to use the protocols::moves::MoverFactory instead.
class LoadUnboundRotMover : public protocols::moves::Mover {

public:

	LoadUnboundRotMover();

	virtual ~LoadUnboundRotMover();

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

	virtual protocols::moves::MoverOP fresh_instance() const;
	virtual protocols::moves::MoverOP clone() const;

	/// @brief parse XML (specifically in the context of the parser/scripting scheme)
	virtual void parse_my_tag(
		utility::tag::TagCOP,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		Pose const & );

private:
	//woo, no data

};//end LoadUnboundRotMover

}//namespace simple_moves
}//namespace protocols

#endif // INCLUDED_protocols_simple_moves_LoadUnboundRotMover_HH
