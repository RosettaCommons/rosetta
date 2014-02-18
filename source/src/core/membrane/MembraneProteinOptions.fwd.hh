// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/membrane/MembraneProteinOptions.fwd.hh
/// @brief  loader creator for membrane protein options
/// @author Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_MembraneProteinOptions_fwd_hh
#define INCLUDED_core_membrane_MembraneProteinOptions_fwd_hh

//utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace membrane {
    
    /// @brief Membrane Protein Options Class
    class MembraneProteinOptions;
    typedef utility::pointer::owning_ptr< MembraneProteinOptions > MembraneProteinOptionsOP;
    typedef utility::pointer::owning_ptr< MembraneProteinOptions const > MembraneProteinOptionsCOP;
        
} // namespace membrane
} // namespace core

#endif // INCLUDED_core_membrane_MembraneProteinOptions_fwd_hh

