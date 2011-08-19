// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/NullDesignMover.hh
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

#ifndef INCLUDED_protocols_protein_interface_design_movers_NullDesignMover_hh
#define INCLUDED_protocols_protein_interface_design_movers_NullDesignMover_hh

#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>

namespace protocols {
namespace protein_interface_design {
namespace movers {

/// @brief does absolutely nothing. useful as if you want to havea final filter without making a move.
/// this can now be assimilated with NullMover
class NullDesignMover : public protocols::moves::Mover
{
public:
	NullDesignMover() : protocols::moves::Mover( "Null" ) {}
	protocols::moves::MoverOP clone() const {
		return( protocols::moves::MoverOP( new NullDesignMover() ));
	}
	protocols::moves::MoverOP fresh_instance() const { return protocols::moves::MoverOP( new NullDesignMover); }
	void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
	virtual ~NullDesignMover();
};

} // movers
} // protein_interface_design
} // protocols


#endif /*INCLUDED_protocols_protein_interface_design_movers_NullDesignMover_HH*/
