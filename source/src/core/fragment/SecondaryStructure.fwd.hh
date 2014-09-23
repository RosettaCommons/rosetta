// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Oliver Lange


#ifndef INCLUDED_core_fragment_SecondaryStructure_fwd_hh
#define INCLUDED_core_fragment_SecondaryStructure_fwd_hh


// Utility headers
// AUTO-REMOVED #include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>

namespace core {
namespace fragment {

// Forward
class SecondaryStructure;

// Types
typedef  utility::pointer::shared_ptr< SecondaryStructure >  SecondaryStructureOP;
typedef  utility::pointer::shared_ptr< SecondaryStructure const >  SecondaryStructureCOP;

} // namespace fragment
} // namespace core

#endif
