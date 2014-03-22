// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 	 core/conformation/membrane/MembraneInfo.fwd.hh
///
/// @brief 	 Membrane Conformation Info
/// @details The Membrane Conformation Info object is responsible for:
///             - maintaining a correct membrane foldtree
///             - maintaining references to the membrane and embedding residues
///             - providing access to membrane related data
///
/// @note    Last Modified 3/12/14
/// @author  Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_conformation_membrane_MembraneInfo_fwd_hh
#define INCLUDED_core_conformation_membrane_MembraneInfo_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace conformation {
namespace membrane {
		
	/// @brief Class: Membrane Conformation Info
	/// @details Handles memrbane conformation, foldtree, and maintains memrbane info
	class MembraneInfo;
	typedef utility::pointer::owning_ptr< MembraneInfo > MembraneInfoOP;
	typedef utility::pointer::owning_ptr< MembraneInfo const > MembraneInfoCOP;
		
} // membrane
} // conformation
} // core

#endif // INCLUDED_core_conformation_membrane_MembraneInfo_fwd_hh

