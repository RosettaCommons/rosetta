// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/MembraneLipo.hh
/// @brief  RMS Energy function. Used to optimize the RMSD between two structures.
/// @author James Thompson


#ifndef INCLUDED_core_scoring_methods_MembraneLipo_hh
#define INCLUDED_core_scoring_methods_MembraneLipo_hh


// Package headers
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/MembranePotential.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <utility/vector1.hh>


// Project headers

// Utility headers


namespace core {
namespace scoring {
namespace methods {


class MembraneLipo : public WholeStructureEnergy  {
public:
	typedef WholeStructureEnergy  parent;

public:


	MembraneLipo();

	/// clone
	virtual
	EnergyMethodOP
	clone() const;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	void
	finalize_total_energy(
		pose::Pose & pose,
		ScoreFunction const &,
		EnergyMap & totals
	) const;
	virtual
	void
	setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const;
	void
	indicate_required_context_graphs(
		utility::vector1< bool > & /*context_graphs_required*/
	) const {}

private:
	// const-ref to scoring database
	MembranePotential const & potential_;
virtual
core::Size version() const;
};


} // methods
} // scoring
} // core

#endif
