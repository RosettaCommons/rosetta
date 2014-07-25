// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file		core/scoring/membrane/FaMPSolvEnergy.cc
///
/// @brief		LK-Type Membrane Solvation Energy
/// @details	Last Modified: 5/13/14
///
/// @author		Patrick Barth (Original)
/// @author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_scoring_membrane_FaMPSolvEnergy_cc
#define INCLUDED_core_scoring_membrane_FaMPSolvEnergy_cc

// Unit Headers
#include <core/scoring/membrane/FaMPSolvEnergy.hh>
#include <core/scoring/methods/ContextDependentTwoBodyEnergy.hh>

// Project Headers
#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/SpanningTopology.hh>

// Package headers
#include <core/conformation/Atom.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/memb_etable/MembEtable.hh>
#include <core/scoring/ScoreFunction.hh>

#include <core/pose/Pose.hh>

// Utility headers
#include <utility/vector1.hh>

#include <ObjexxFCL/FArray3D.hh>

// C++ Headers
#include <cstdlib>

namespace core {
namespace scoring {
namespace membrane {
	
class FaMPSolvEnergy;
typedef utility::pointer::owning_ptr< FaMPSolvEnergy > FaMPSolvEnergyOP;
typedef utility::pointer::owning_ptr< FaMPSolvEnergy const > FaMPSolvEnergyCOP;
	
} // membrane
} // scoring
} // core

#endif // INCLUDED_core_scoring_methods_FaMPSolvEnergy_cc
