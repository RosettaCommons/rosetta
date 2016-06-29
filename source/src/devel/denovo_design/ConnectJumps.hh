// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/devel/denovo_design/ConnectJumps.hh
/// @brief The ConnectJumps Protocol
/// @details
/// @author Tom Linsky


#ifndef INCLUDED_devel_denovo_design_ConnectJumps_hh
#define INCLUDED_devel_denovo_design_ConnectJumps_hh

// Unit headers
#include <devel/denovo_design/ConnectJumps.fwd.hh>

// Protocol headers
#include <protocols/denovo_design/movers/BridgeChainsMover.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

namespace devel {
namespace denovo_design {

class ConnectJumps : public protocols::denovo_design::movers::BridgeChainsMover {
public:
	/// @brief default constructor
	ConnectJumps();

	/// @brief virtual constructor to allow derivation
	virtual ~ConnectJumps();

	/// @brief Return the name of this mover.
	virtual std::string get_name() const;

	/// @brief return a fresh instance of this class in an owning pointer
	virtual protocols::moves::MoverOP clone() const;

	/// @brief Apply the ConnectJumps. Overloaded apply function from mover base class.
	virtual void apply( core::pose::Pose & pose );
};


} // denovo_design
} // devel

#endif
