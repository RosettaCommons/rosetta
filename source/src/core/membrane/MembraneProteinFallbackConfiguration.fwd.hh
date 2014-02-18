// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/membrane/MembraneProteinFallbackConfiguration.fwd.hh
/// @brief  Membrane Protein Fallback Configuration
/// @author Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_MembraneProteinFallbackConfiguration_fwd_hh
#define INCLUDED_core_membrane_MembraneProteinFallbackConfiguration_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace membrane {

    /// @brief Membrane Protein Fallback Configuraiton Class
    class MembraneProteinFallbackConfiguration;
    typedef utility::pointer::owning_ptr< MembraneProteinFallbackConfiguration > MembraneProteinFallbackConfigurationOP;
    typedef utility::pointer::owning_ptr< MembraneProteinFallbackConfiguration const > MembraneProteinFallbackConfigurationCOP;
        
} // namespace membrane
} // namespace core

#endif // INCLUDED_core_membrane_MembraneProteinFallbackConfiguration_fwd_hh

