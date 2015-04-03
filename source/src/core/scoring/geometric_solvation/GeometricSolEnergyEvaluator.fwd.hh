// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/scoring/geometric_solvation/GeometricSolEnergyEvaluator.fwd.hh
/// @brief 
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_scoring_geometric_solvation_GeometricSolEnergyEvaluator_FWD_HH
#define INCLUDED_core_scoring_geometric_solvation_GeometricSolEnergyEvaluator_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace geometric_solvation {
	
	class GeometricSolEnergyEvaluator;
	typedef utility::pointer::shared_ptr< GeometricSolEnergyEvaluator > GeometricSolEnergyEvaluatorOP;
	typedef utility::pointer::shared_ptr< GeometricSolEnergyEvaluator const > GeometricSolEnergyEvaluatorCOP;
	
} //geometric_solvation 
} //scoring 
} //core 

#endif
