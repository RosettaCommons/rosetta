// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   	core/membrane/scoring/MembraneSearch.fwd.hh
///
/// @brief  	Membrane Search
/// @detail		Performs a search method for optimal membrane normal and center from scoring
///				this is used in both conformation and scoring
///
/// @author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_scoring_MembraneSearch_fwd_hh
#define INCLUDED_core_membrane_scoring_MembraneSearch_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

/// @brief 	Class: MembraneSearch
/// @detail Search for membrane embedding normal and center - can be called on both
///			a membrane definition and an embedding definition's parameters
class MembraneSearch;
typedef utility::pointer::owning_ptr< MembraneSearch > MembraneSearchOP;
typedef utility::pointer::owning_ptr< MembraneSearch const > MembraneSearchCOP;


#endif // INCLUDED_core_membrane_scoring_MembraneSearch_fwd_hh