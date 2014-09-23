// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file		core/conformation/membrane/MembraneInfo.fwd.hh
///
/// @brief		MembraneInfo - Membrane Pose Definition Object
/// @details	The membrane info object is responsible for storing non-coordinate
///				derived information, extending the traditional definiton of a pose
///				to describe a membrane protein in Rosetta. This information includes:
///				 - the resiue number of the membrane virtual residue containing
///					the posiiton of the membrane
///				 - membrane spanning topology
///				 - membrane lipophilicity info (user-specified)
///
///				This object belongs to the conformation and should be accessed
///				via pose.conformation().membrane_info().
///
///	     		Last Modified: 6/21/14
///
/// @author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_conformation_membrane_MembraneInfo_fwd_hh
#define INCLUDED_core_conformation_membrane_MembraneInfo_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace conformation {
namespace membrane {
		
/// @brief Class: Membrane Info
/// @details Handles memrbane conformation, foldtree, and maintains membrane info
class MembraneInfo;
typedef utility::pointer::shared_ptr< MembraneInfo > MembraneInfoOP;
typedef utility::pointer::shared_ptr< MembraneInfo const > MembraneInfoCOP;
		
} // membrane
} // conformation
} // core

#endif // INCLUDED_core_conformation_membrane_MembraneInfo_fwd_hh

