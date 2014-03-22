// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file		core/scoring/Membrane_FAPotential.cc
///
/// @brief		Membrane FA Potential - Class for Fullatom Membrane Scoring Methods
/// @details	Compute High resolution energy terms and high resolution embedding corrections
///				for penalties. Also contains pass-through methods for accessing and updating
///				mp framework supported data in a membrane conformation.
///				Last Modified: 3/11/14
///
///	@author		Rebecca Faye Alford (rfalford12@gmail.com)
/// @author		Julia Koehler Leman (julia.koehler1982@gmail.com)
/// @author		Patrick Barth (original)


#ifndef INCLUDED_core_scoring_Membrane_FAPotential_fwd_hh
#define INCLUDED_core_scoring_Membrane_FAPotential_fwd_hh

// Unit headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
	namespace scoring {
		
		class Membrane_FAEmbed;
		typedef utility::pointer::owning_ptr< Membrane_FAEmbed > Membrane_FAEmbedOP;
		typedef utility::pointer::owning_ptr< Membrane_FAEmbed const > Membrane_FAEmbedCOP;
		
		class Membrane_FAPotential;
		typedef  utility::pointer::owning_ptr< Membrane_FAPotential > Membrane_FAPotentialOP;
		typedef  utility::pointer::owning_ptr< Membrane_FAPotential const > Membrane_FAPotentialCOP;
		
	} // scoring
} // core

#endif // INCLUDED_core_scoring_Membrane_FAPotential_fwd_hh