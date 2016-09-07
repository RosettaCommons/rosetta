// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/simple_moves/VirtualRootMover.hh
/// @brief  Manipulate virtual roots on poses.
/// @author Rocco Moretti (rmoretti@u.washington.edu)

#ifndef INCLUDED_protocols_simple_moves_VirtualRootMover_hh
#define INCLUDED_protocols_simple_moves_VirtualRootMover_hh

// Unit headers
#include <protocols/simple_moves/VirtualRootMover.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <core/types.hh>

namespace protocols {
namespace simple_moves {

class VirtualRootMover : public moves::Mover {

public:
	VirtualRootMover();
	~VirtualRootMover() override;
	moves::MoverOP clone() const override;
	moves::MoverOP fresh_instance() const override;

	void apply( Pose & pose ) override;
	std::string get_name() const override;


	void parse_my_tag(
		TagCOP,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		moves::Movers_map const &,
		Pose const & ) override;

	void set_remove( bool const remove );
	void set_removable( bool const removable );
private:
	bool remove_;
	bool removable_;

};

} // simple_moves
} // protocols

#endif //INCLUDED_protocols_simple_moves_VirtualRootMover_HH
