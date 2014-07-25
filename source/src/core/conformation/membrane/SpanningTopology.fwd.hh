// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		core/conformation/membrane/SpanningTopology.fwd.hh
///
/// @brief      Membrane Spanning Topology Data
/// @details    Stores information describing the membrane spanning
///             topology of a pose. This definition is a dependency for embedding definitions
///             and requires a spanningfile from OCTOPUS for initialization
///
/// @note       Last Modified: 1/1/14
/// @author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_conformation_membrane_SpanningTopology_fwd_hh
#define INCLUDED_core_conformation_membrane_SpanningTopology_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
    namespace conformation {
        namespace membrane {
            
            /// @brief      Class: Membrane Spanning Topology
            /// @details    Stores information describing the membrane spanning
            ///             topology of a pose. This definition is a dependency for embedding definitions
            ///             and requires a spanningfile from OCTOPUS for initialization
            class SpanningTopology;
            typedef utility::pointer::owning_ptr< SpanningTopology > SpanningTopologyOP;
            typedef utility::pointer::owning_ptr< SpanningTopology const > SpanningTopologyCOP;
            
            
        } // membrane
    } // conformation
} // core

#endif // INCLUDED_core_membrane_properties_SpanningTopology_fwd_hh
