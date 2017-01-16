// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/recces/setup_util.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_recces_setup_util_HH
#define INCLUDED_protocols_recces_setup_util_HH

#include <protocols/recces/options/RECCES_Options.fwd.hh>
#include <core/pose/Pose.fwd.hh>

namespace protocols {
namespace recces {

//////////////////////////////////////////////////////////////////////////////
core::pose::PoseOP
recces_pose_setup( options::RECCES_Options const & options );

//////////////////////////////////////////////////////////////////////////////
core::pose::PoseOP
pose_setup_turner(
	std::string const & seq1,
	std::string const & seq2);


//////////////////////////////////////////////////////////////////////////////
core::pose::PoseOP
pose_setup_from_file ( options::RECCES_Options const & options );


} //recces
} //protocols

#endif
