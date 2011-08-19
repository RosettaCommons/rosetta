// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @file   core/scoring/methods/EnvPairPotential.fwd.hh
/// @brief   Membrane Potential
/// @author Patrick Barth


#ifndef INCLUDED_core_scoring_Membrane_FAPotential_fwd_hh
#define INCLUDED_core_scoring_Membrane_FAPotential_fwd_hh


// Unit headers
#include <utility/pointer/owning_ptr.hh>

// Package headers

// Project headers

// Utility headers

// C++


namespace core {
namespace scoring {

class Membrane_FAEmbed;
typedef utility::pointer::owning_ptr< Membrane_FAEmbed > Membrane_FAEmbedOP;

class Membrane_FAPotential;

} // ns scoring
} // ns core

#endif
