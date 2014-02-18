// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		core/membrane/properties/LipidAccInfo.fwd.hh
///
/// @brief      Membrane Lipid Accessibility Data
/// @details    Stores lipid accessibility data derived from OCTOPUS spanning file
///             and psiblast search using run_lips.pl script
///
/// @note       Last Modified: 1/1/14
/// @author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_properties_LipidAccInfo_fwd_hh
#define INCLUDED_core_membrane_properties_LipidAccInfo_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace membrane {
namespace properties {
            
    /// @brief      Class: Membrane Lipid Accessibility Data
    /// @details    Stores lipid accessibility data derived from OCTOPUS spanning file
    ///             and psiblast search using run_lips.pl script
    class LipidAccInfo;
    typedef utility::pointer::owning_ptr< LipidAccInfo > LipidAccInfoOP;
    typedef utility::pointer::owning_ptr< LipidAccInfo const > LipidAccInfoCOP;
    
} // properties
} // membrane
} // core

#endif // INCLUDED_core_membrane_properties_LipidAccInfo_fwd_hh

