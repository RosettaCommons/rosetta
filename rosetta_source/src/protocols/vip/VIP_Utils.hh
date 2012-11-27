// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/vip/VIP_Utils.hh
/// @brief

#ifndef INCLUDED_protocols_vip_VIP_Utils_HH
#define INCLUDED_protocols_vip_VIP_Utils_HH

#include "core/kinematics/MoveMap.fwd.hh"
#include "core/types.hh"
#include "core/pose/Pose.hh"
#include "core/scoring/packstat/types.hh"
#include <string>

namespace protocols {
namespace vip {

using core::Real;

using namespace core::scoring::packstat;

std::string base_name(const std::string& str);
std::string get_out_tag(std::string fname);
core::Real output_packstat( core::pose::Pose & );
void set_local_movemap( core::pose::Pose & pose, core::Size position, core::kinematics::MoveMapOP mmap );

}}
#endif
