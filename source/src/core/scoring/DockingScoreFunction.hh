// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/symmetry/DockingScoreFunction.hh
/// @brief  Symmetric Score function class
/// @author Ingemar Andre


#ifndef INCLUDED_core_scoring_DockingScoreFunction_hh
#define INCLUDED_core_scoring_DockingScoreFunction_hh

// Unit headers
#include <core/scoring/DockingScoreFunction.fwd.hh>

// Package headers
#include <core/scoring/ScoreFunction.hh>

// Project headers
#include <core/conformation/symmetry/SymmetricConformation.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {

class DockingScoreFunction : public ScoreFunction
{
public:
	typedef ScoreFunction parent;

public:

	/// ctor
	DockingScoreFunction();

private:
	/// @brief Copy constructor and assignment operators private for *ScoreFunctions as they discard subtype info.
	DockingScoreFunction &
	operator=( DockingScoreFunction const & );

	DockingScoreFunction( DockingScoreFunction const & );

public:

	/// @brief INTERNAL USE ONLY
	virtual void assign( ScoreFunction const & src);

	/// @brief INTERNAL USE ONLY
	virtual void assign( DockingScoreFunction const & src);

	virtual ScoreFunctionOP clone() const;

  /////////////////////////////////////////////////////////////////////////////
  // score
  /////////////////////////////////////////////////////////////////////////////

	virtual Real
	operator ()( pose::Pose & pose ) const;


};


} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_DockingScoreFunction_HH
