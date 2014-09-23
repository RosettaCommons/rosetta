// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/rotamer_set/RotamerSetOperation.fwd.hh
/// @brief  rotamer set operation forward declaration
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)


#ifndef INCLUDED_core_pack_rotamer_set_RotamerSetOperation_fwd_hh
#define INCLUDED_core_pack_rotamer_set_RotamerSetOperation_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

// C++ Headers
#include <list>

namespace core {
namespace pack {
namespace rotamer_set {

class RotamerOperation;
class RotamerSetOperation;

typedef utility::pointer::shared_ptr< RotamerOperation > RotamerOperationOP;
typedef utility::pointer::shared_ptr< RotamerOperation const > RotamerOperationCOP;

typedef std::list< RotamerOperationOP > RotamerOperations;

typedef utility::pointer::shared_ptr< RotamerSetOperation > RotamerSetOperationOP;
typedef utility::pointer::shared_ptr< RotamerSetOperation const > RotamerSetOperationCOP;

typedef std::list< RotamerSetOperationOP > RotSetOperationList;
typedef std::list< RotamerSetOperationOP >::const_iterator RotSetOperationListIterator;

}
}
}

#ifdef USEBOOSTSERIALIZE
#include <boost/serialization/access.hpp>
#endif


#endif
