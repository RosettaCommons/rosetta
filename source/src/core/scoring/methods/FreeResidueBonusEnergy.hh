// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/FreeResidueBonusEnergy.hh
/// @brief  Score bonus for residues that form 0 interactions (i.e. bulges) - rough approximation to entropic bonus for flexible residues
/// @author Arvind Kannan


#ifndef INCLUDED_core_scoring_methods_FreeResidueBonusEnergy_hh
#define INCLUDED_core_scoring_methods_FreeResidueBonusEnergy_hh


// Package headers
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>



// Utility headers


namespace core {
namespace scoring {
namespace methods {


class FreeResidueBonusEnergy : public methods::ContextIndependentOneBodyEnergy  {
public:
	typedef methods::ContextIndependentOneBodyEnergy  parent;

public:

	/// @brief Created in order to stabilize bulges during SWM calculations.
	FreeResidueBonusEnergy();

	/// clone
	virtual
	methods::EnergyMethodOP
	clone() const;
	
	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////
	virtual
	void
	residue_energy(
				   conformation::Residue const & rsd,
				   pose::Pose const &,
				   EnergyMap & emap
				   ) const;
	
	void indicate_required_context_graphs(utility::vector1< bool > &) const {}
	
private:
	///@brief Return the version of the energy method
	virtual
	core::Size version() const;

	
	
/////////////////////////////////////////////////////////////////////////////
// data
/////////////////////////////////////////////////////////////////////////////
private:


};


}
}
}

#endif // INCLUDED_core_scoring_methods_FreeResidueBonusEnergy_HH
