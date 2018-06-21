// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/network/ui_mover.hh
/// @brief: UIMover, apply to send in-progress results to UI
///
/// @author Sergey Lyskov

#pragma once

#include <protocols/network/ui_mover.fwd.hh>

#include <protocols/moves/Mover.hh>

#ifdef ZEROMQ
#include <libzmq/include/zmq.hpp>
#endif // ZEROMQ

namespace protocols {
namespace network {


/// Mover to send in-progress data to UI client
/// Note: you will need to build Rosetta with extra=zmq for this Mover to do any real work
class UIMover : public protocols::moves::Mover
{
public:
	/// @brief ctor
	explicit UIMover();

	/// @brief cctor
	explicit UIMover(UIMover const &);

	/// @brief dtor
	~UIMover() override;

	void apply( Pose & ) override;

	std::string get_name() const override { return "UIMover"; };


#ifdef ZEROMQ

private:
	zmq::socket_t hal;

#endif // ZEROMQ

};

} // namespace network
} // namespace protocols
