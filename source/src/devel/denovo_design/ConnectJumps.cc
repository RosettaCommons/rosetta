// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/devel/denovo_design/ConnectJumps.cc
/// @brief The ConnectJumps
/// @details
/// @author Tom Linsky


// Unit Headers
#include <devel/denovo_design/ConnectJumps.hh>
#include <devel/denovo_design/ConnectJumpsCreator.hh>

// Basic Headers
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "devel.denovo_design.ConnectJumps" );

////////////////////////////////////////////////////////////////////////////////////////////////////

namespace devel {
namespace denovo_design {
////////////////////////////////////////////////////////////////////////////////////////////////////

std::string
ConnectJumpsCreator::keyname() const
{
	return ConnectJumpsCreator::mover_name();
}

protocols::moves::MoverOP
ConnectJumpsCreator::create_mover() const {
	return protocols::moves::MoverOP( new ConnectJumps() );
}

std::string
ConnectJumpsCreator::mover_name()
{
	return "ConnectJumps";
}

///  ---------------------------------------------------------------------------------
///  ConnectJumps main code:
///  ---------------------------------------------------------------------------------

/// @brief default constructor
ConnectJumps::ConnectJumps() :
	protocols::denovo_design::movers::BridgeChainsMover()
{
}

/// @brief destructor - this class has no dynamic allocation, so
ConnectJumps::~ConnectJumps() = default;


/// Return a copy of ourselves
protocols::moves::MoverOP
ConnectJumps::clone() const
{
	return protocols::moves::MoverOP( new ConnectJumps(*this) );
}

std::string
ConnectJumps::get_name() const
{
	return ConnectJumpsCreator::mover_name();
}

void ConnectJumps::apply( core::pose::Pose & pose )
{
	TR << "*****************************************************************************" << std::endl;
	TR << "WARNING: ConnectJumps is deprecated. Please use BridgeChains instead." << std::endl;
	TR << "*****************************************************************************" << std::endl;
	protocols::denovo_design::movers::BridgeChainsMover::apply( pose );
}

} // namespace denovo_design
} // namespace devel
