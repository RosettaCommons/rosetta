// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
/// @file protocols/network/cloud.hh
/// @brief: RosettaCloud integration
///
/// @author Sergey Lyskov

///
/// For example of using RosettaCloud API please see main/source/src/apps/pilot/sergey/cloud_demo.cc
///
/// https://ui.graylab.jhu.edu/execution/summaries - list of all summaries already created on your account
///
/// https://ui.graylab.jhu.edu/execution/es/0      - `magic` link to _the latest_ summary (re-bound latest ES on page refresh)
///
///
/// Relevant command line options:
///
/// -cloud:auth  - specify user name and password (visit https://ui.graylab.jhu.edu/settings to get credentials for your account)
///
/// -cloud:key   - specify `key` to use when process query Cloud server for new ExecutionSummary ID.
///                Default is empty string (no key) which mean always create a new ES instance
///                Specifying `key` is a way to have somewhat permanent ES instance. The downside of
///                using it is that each new run override the previous one (but this could be desirable
///                when debugging)
///
/// -cloud:clean - delete all previously posted files upon ES instance creation. If you using `key`
///                option above for debugging then you most likely want to set `-clean` to true as well.
///
/// -cloud:block - specify what to do in `post_*` when network queue is full: blocking will pause main thread execution until
///                some network operation is finished while setting block to `false` will instruct `post_*` function to drop
///                new request without delaying main thread.
///

#pragma once

#include <core/pose/Pose.fwd.hh>

#include <string>

namespace protocols {
namespace network {

void post_file(std::string const &file_name, std::string const &data, bool append = false);

void post_decoy(std::string const &file_name, core::pose::Pose const &pose);

void post_decoy(core::pose::Pose const &pose);

} // namespace network
} // namespace protocols
