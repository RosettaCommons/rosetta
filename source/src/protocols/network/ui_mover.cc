// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/network/ui_mover.cc
/// @brief: UIMover, apply to send in-progress results to UI
///
/// @author Sergey Lyskov

#include <protocols/network/ui_mover.hh>

#include <protocols/network/util.hh>

#include <utility/json_utilities.hh>

namespace protocols {
namespace network {

using std::string;

#if defined(ZEROMQ)

UIMover::UIMover() : hal(zmq_context(), ZMQ_DEALER)
{
	hal.connect(_hal_address_);
}

UIMover::UIMover(UIMover const &) : UIMover()
{
}

UIMover::~UIMover() = default;

void UIMover::apply(Pose & pose)
{
	auto pose_binary = protocols::network::pose_to_bytes(pose);

	nlohmann::json result;

	result["pose"] = pose_binary;

	string binary_result;
	nlohmann::json::basic_json::to_msgpack(result, binary_result);

	send_message(hal, _m_progress_, ZMQ_SNDMORE);
	send_message(hal, binary_result);
}


#else // defined(ZEROMQ)

UIMover::UIMover() = default;

UIMover::UIMover(UIMover const &) = default;

UIMover::~UIMover() = default;

void UIMover::apply( Pose & ) {}

#endif // defined(ZEROMQ)

} // namespace network
} // namespace protocols
