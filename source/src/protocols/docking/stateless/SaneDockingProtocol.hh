// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   SaneDockingProtocol.hh
///
/// @brief DockingProtocol class that doesn't use the more esoteric features of
/// the Mover base class (e.g get_input_pose()).
/// @author James Thompson

#ifndef INCLUDED_protocols_docking_stateless_SaneDockingProtocol_hh
#define INCLUDED_protocols_docking_stateless_SaneDockingProtocol_hh

#include <protocols/docking/DockingProtocol.hh>
#include <protocols/docking/stateless/SaneDockingProtocol.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace docking {
namespace stateless {

class SaneDockingProtocol : public protocols::docking::DockingProtocol {
public:
	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
};

} // stateless
} // docking
} // protocols

#endif
