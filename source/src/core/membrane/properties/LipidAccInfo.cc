// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		core/membrane/properties/LipidAccInfo.hh
///
/// @brief      Membrane Lipid Accessibility Data
/// @details    Stores lipid accessibility data derived from OCTOPUS spanning file
///             and psiblast search using run_lips.pl script
///
/// @note       Last Modified: 1/1/14
/// @author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_properties_LipidAccInfo_cc
#define INCLUDED_core_membrane_properties_LipidAccInfo_cc

// Unit headers
#include <core/membrane/properties/LipidAccInfo.hh>

// Package Headers
#include <core/types.hh>

#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <numeric/xyzVector.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/FArray4D.hh>

#include <utility/vector1.hh>
#include <utility/vector0.hh>

// Platform headers
#include <platform/types.hh>

// C++ Headers
#include <cstddef>
#include <cstdlib>
#include <string>

/// @brief      Membrane Lipid Accessibility Data
/// @details    Stores lipid accessibility data derived from OCTOPUS spanning file
///             and psiblast search using run_lips.pl script

namespace core {
namespace membrane {
namespace properties {
    
    /// @brief Constructor
    LipidAccInfo::LipidAccInfo() :
        utility::pointer::ReferenceCount()
    {}

    /// @brief Conpy Constructor
    LipidAccInfo::LipidAccInfo( LipidAccInfo const & src ) :
        utility::pointer::ReferenceCount()
    {
        copy_data(src, *this);
    }

    /// @brief Destructor
    LipidAccInfo::~LipidAccInfo() {}

    /// @brief Getters
    // Lipid burial and exposure
    utility::vector1< core::Real > LipidAccInfo::lipid_exposure() { return lipid_exposure_; }
    utility::vector1< core::Real > LipidAccInfo::lipid_burial() { return lipid_burial_; }

    /// @brief Setters
    // Lipid burial and exposure
    void LipidAccInfo::set_lipid_exposure( utility::vector1< core::Real > exp ) { lipid_exposure_ = exp; }
    void LipidAccInfo::set_lipid_burial( utility::vector1< core::Real > buried ) { lipid_burial_ = buried; }

    /// @brief Copy Data
    void
    LipidAccInfo::copy_data( LipidAccInfo src, LipidAccInfo copy ) {
        
        src.lipid_exposure_ = copy.lipid_exposure_;
        src.lipid_burial_ = copy.lipid_burial_;
    }
    
    
} // properties
} // membrane
} // core

#endif // INCLUDED_core_membrane_properties_LipidAccInfo_cc



