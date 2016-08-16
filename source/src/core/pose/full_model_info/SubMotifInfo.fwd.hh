// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/full_model_info/SubMotifInfo.fwd.hh
/// @brief  Stores information about submotifs in a pose.
/// @author Caleb Geniesse

#ifndef INCLUDED_core_pose_full_model_info_SubMotifInfo_fwd_hh
#define INCLUDED_core_pose_full_model_info_SubMotifInfo_fwd_hh


#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>


namespace core {
namespace pose {
namespace full_model_info {

class SubMotifInfo;
typedef utility::pointer::shared_ptr< SubMotifInfo > SubMotifInfoOP;
typedef utility::pointer::shared_ptr< SubMotifInfo const > SubMotifInfoCOP;


} //full_model_info
} //pose
} //core

#endif
