// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragments/FragSet.fwd.hh
/// @brief  set of fragments for a certain alignment frame
/// @author Oliver Lange (olange@u.washington.edu)
/// @date   Wed Oct 20 12:08:31 2007
///


#ifndef INCLUDED_core_fragment_FragSet_fwd_hh
#define INCLUDED_core_fragment_FragSet_fwd_hh

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.fwd.hh>

// C++ headers
#include <utility>

namespace core {
namespace fragment {

// Forward
class FragSet;

typedef utility::pointer::shared_ptr< FragSet > FragSetOP;
typedef utility::pointer::shared_ptr< FragSet const > FragSetCOP;

typedef std::pair< Size, Size > Range; //does such a type already exist ?

//silly set that contains residues that can be sampled
// in exact terms: the InsertMap contains the start() numbers of all applicable Frames
typedef utility::vector1< Size > InsertMap;
typedef utility::vector1< Size > InsertSize;


} // namespace fragment
} // namespace core

#endif
