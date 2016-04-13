// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/SASAPotential.fwd.hh
/// @brief
/// @author


#ifndef INCLUDED_core_scoring_SASAPotential_fwd_hh
#define INCLUDED_core_scoring_SASAPotential_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {

///

class PDatom;
typedef utility::pointer::shared_ptr< PDatom > PDatomOP;
typedef utility::pointer::shared_ptr< PDatom const > PDatomCOP;

class PDvertex;
typedef utility::pointer::shared_ptr< PDvertex > PDvertexOP;
typedef utility::pointer::shared_ptr< PDvertex const > PDvertexCOP;

class PDinter;
typedef utility::pointer::shared_ptr< PDinter > PDinterOP;
typedef utility::pointer::shared_ptr< PDinter const > PDinterCOP;

class SASAPotential;
typedef  utility::pointer::shared_ptr< SASAPotential > SASAPotentialOP;
typedef  utility::pointer::shared_ptr< SASAPotential const > SASAPotentialCOP;

} // scoring
} // core

#endif
