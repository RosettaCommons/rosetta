// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file		core/scoring/MembranePotential.fwd.hh
///
/// @brief		Membrane Potential - Base Scoring Methods for LowRes Energy Function
/// @details	Compute Low Res membrane energy terms: Menv, MPair, MCBeta and Membrane
///				penalties. Also contains pass-through methods for accessing and updating
///				mp framework supported data in a membrane conformation.
///				Last Modified: 3/11/14
///
///	@author		Rebecca Faye Alford (rfalford12@gmail.com)
/// @author		Julia Koehler Leman (julia.koehler1982@gmail.com)
/// @author		Bjorn Wallner (Original)

#ifndef INCLUDED_core_scoring_MembranePotential_fwd_hh
#define INCLUDED_core_scoring_MembranePotential_fwd_hh

// Unit headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
	namespace scoring {
		
		class MembraneEmbed;
		typedef utility::pointer::owning_ptr< MembraneEmbed > MembraneEmbedOP;
		typedef utility::pointer::owning_ptr< MembraneEmbed const > MembraneEmbedCOP;
		
		
		class MembranePotential;
		typedef utility::pointer::owning_ptr< MembranePotential > MembranePotentialOP;
		typedef utility::pointer::owning_ptr< MembranePotential const > MembranePotentialCOP;
		
	} // scoring
} // core

#endif // INCLUDED_core_scoring_MembranePotential_fwd_hh

