// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/graph/ArrayPool.fwd.hh
/// @brief  Forward declaration for ArrayPool classes
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_graph_ArrayPool_fwd_hh
#define INCLUDED_core_graph_ArrayPool_fwd_hh

/// Project headers
#include <platform/types.hh>

/// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace graph {

template< class T >
class Array0;

template < class T >
class NegSpaceElement;

template < class T >
class ArrayPoolElement;

template < class T >
class ArrayPool;


typedef utility::pointer::shared_ptr< ArrayPool< platform::Real > > RealArrayPoolOP;
typedef utility::pointer::shared_ptr< ArrayPool< platform::Real > const > RealArrayPoolCOP;

}
}

#endif
