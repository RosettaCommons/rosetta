// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/forge/components/VarLengthBuild.fwd.hh
/// @brief  forward declaration for protocols::forge::components::VarLengthBuild
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_protocols_forge_components_VarLengthBuild_fwd_hh
#define INCLUDED_protocols_forge_components_VarLengthBuild_fwd_hh

// utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>


namespace protocols {
namespace forge {
namespace components {


/// @brief forward declaration for protocols::forge::components::VarLengthBuild
class VarLengthBuild;


/// @brief access pointer for VarLengthBuild
typedef utility::pointer::weak_ptr< VarLengthBuild > VarLengthBuildAP;


/// @brief const access pointer for VarLengthBuild
typedef utility::pointer::weak_ptr< VarLengthBuild const > VarLengthBuildCAP;


/// @brief owning pointer for VarLengthBuild
typedef utility::pointer::shared_ptr< VarLengthBuild > VarLengthBuildOP;


/// @brief const owning pointer for VarLengthBuild
typedef utility::pointer::shared_ptr< VarLengthBuild const > VarLengthBuildCOP;


namespace VLB_VallMemoryUsage {

/// @brief enum dictating memory usage wrt VallLibrary in VarLengthBuild
/// @details Dictates whether to keep the VallLibrary in memory, or whether
///  to clear it under particular circumstances after picking fragments.
///  Default for VarLengthBuild is currently KEEP_IN_MEMORY.
enum Enum {
	KEEP_IN_MEMORY,
	CLEAR_IF_CACHING_FRAGMENTS,
	ALWAYS_CLEAR
};

}


} // components
} // forge
} // protocols


#endif /* INCLUDED_protocols_forge_components_VarLengthBuild_FWD_HH */
