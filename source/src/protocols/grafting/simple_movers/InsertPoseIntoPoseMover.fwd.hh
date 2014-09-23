// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/grafting/InsertPoseIntoPoseMover.fwd.hh
/// @brief   Forward declarations for InsertPoseIntoPoseMover class
/// @author  Jared Adolf-Bryfogle jadolfbr@gmail.com

#ifndef INCLUDED_my_namespace_name_InsertPoseIntoPoseMover_FWD_HH
#define INCLUDED_my_namespace_name_InsertPoseIntoPoseMover_FWD_HH

// Utility header
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace grafting {

/// @brief  
class InsertPoseIntoPoseMover;

typedef utility::pointer::shared_ptr<InsertPoseIntoPoseMover> InsertPoseIntoPoseMoverOP;
typedef utility::pointer::shared_ptr<InsertPoseIntoPoseMover const> InsertPoseIntoPoseMoverCOP;

}}  // namespace grafting_protocls

#endif  // INCLUDED_grafting_protocols_InsertPoseIntoPoseMover.FWD_HH
