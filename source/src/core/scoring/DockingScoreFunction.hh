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

/// @file   core/scoring/symmetry/DockingScoreFunction.hh
/// @brief  Symmetric Score function class
/// @author Ingemar Andre


#ifndef INCLUDED_core_scoring_DockingScoreFunction_hh
#define INCLUDED_core_scoring_DockingScoreFunction_hh

// Unit headers
#include <core/scoring/DockingScoreFunction.fwd.hh>

// Package headers
#include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <core/scoring/hbonds/HBondSet.fwd.hh>

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
	///@brief Copy constructor and assignment operators private for *ScoreFunctions as they discard subtype info.
	DockingScoreFunction &
	operator=( DockingScoreFunction const & );

	DockingScoreFunction( DockingScoreFunction const & );

public:

	///@brief INTERNAL USE ONLY
	virtual void assign( ScoreFunction const & src);

	///@brief INTERNAL USE ONLY
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
