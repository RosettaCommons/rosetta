// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/chemical/IdealBondLengthSet.fwd.hh
/// @author P. Douglas Renfrew (renfrew@nyu.edu)


#ifndef INCLUDED_core_chemical_IdealBondLengthSet_fwd_hh
#define INCLUDED_core_chemical_IdealBondLengthSet_fwd_hh

#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>

namespace core {
namespace chemical {

class IdealBondLengthSet;

typedef  utility::pointer::weak_ptr< IdealBondLengthSet > IdealBondLengthSetAP;
typedef  utility::pointer::weak_ptr< IdealBondLengthSet const > IdealBondLengthSetCAP;
typedef  utility::pointer::shared_ptr< IdealBondLengthSet > IdealBondLengthSetOP;
typedef  utility::pointer::shared_ptr< IdealBondLengthSet const > IdealBondLengthSetCOP;

}
}

#endif
