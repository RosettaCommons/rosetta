// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/kinematics/Edge.fwd.hh
/// @brief  Kinematics Edge forward declarations header
/// @author Phil Bradley


#ifndef INCLUDED_core_kinematics_Edge_fwd_hh
#define INCLUDED_core_kinematics_Edge_fwd_hh


// Utility headers
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>


namespace core {
namespace kinematics {


// Forward
class Edge;


// Types
typedef  utility::pointer::weak_ptr< Edge >  EdgeAP;
typedef  utility::pointer::shared_ptr< Edge >  EdgeOP;


} // namespace kinematics
} // namespace core

#endif // INCLUDED_core_kinematics_Edge_FWD_HH
