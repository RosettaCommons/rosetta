// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		core/membrane/geometry/EmbeddingFactory.fwd.hh
///
/// @brief 		Factory Method for creating a membrane embedding residue from user options and tags
/// @details 	Creates reference center, normal, depth/spanning definition and topology
///				definition from the user's specifications
///
/// @author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_geometry_EmbeddingFactory_fwd_hh
#define INCLUDED_core_membrane_geometry_EmbeddingFactory_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace membrane {
namespace geometry {

    /// @brief Create an Embedding Definition
    /// @details Generate the correct protein emebdding using a user-specified
    /// method
    class EmbeddingFactory;
    typedef utility::pointer::owning_ptr< EmbeddingFactory > EmbeddingFactoryOP;
    typedef utility::pointer::owning_ptr< EmbeddingFactory const > EmbeddingFactoryCOP;
    
} // geometry
} // membrane
} // core

#endif // INCLUDED_core_membrane_geometry_EmbeddingFactory_fwd_hh

