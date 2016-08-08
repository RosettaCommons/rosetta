// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file --path--/--class--.hh
/// @brief --brief--
/// @author --name-- (--email--)

#ifndef INCLUDED_--path_underscore--_--class--_hh
#define INCLUDED_--path_underscore--_--class--_hh

// Unit headers
#include <--path--/--class--.fwd.hh>
#include <core/scoring/methods/WholeStructureEnergy.hh>

// Project Headers
#include <core/scoring/methods/EnergyMethod.hh> 

#include <core/scoring/ScoreFunction.fwd.hh> 
#include <core/scoring/ScoreType.hh>
#include <core/pose/Pose.fwd.hh>

// Basic/Utility headers
#include <utility/vector1.hh> 
#include <core/types.hh> 

--namespace--

///@brief --brief--
class --class-- : public core::scoring::methods::WholeStructureEnergy {

	typedef core::scoring::methods::WholeStructureEnergy parent;

public:

	--class--();

	// copy constructor (not needed unless you need deep copies)
	//--class--( --class-- const & src );

	// destructor (important for properly forward-declaring smart-pointer members)
	virtual ~--class--();

	/// @brief Clone: create a copy of this object, and return an owning pointer
	/// to the copy.
	virtual core::scoring::methods::EnergyMethodOP clone() const;

	/// @brief Indicate required setup steps for scoring
	virtual
	void setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const;

	/// @brief Is the score context dependent or context independent? 
	virtual void indicate_required_context_graphs( utility::vector1< bool > &context_graphs_required ) const;

	/// @brief Indicates the current version of this score term
	virtual core::Size version() const;

	/// @brief Actually calculate the total energy
	/// @details Called by the scoring machinery.
	virtual void finalize_total_energy( core::pose::Pose & pose, ScoreFunction const &, EnergyMap & totals ) const;

	virtual
	Distance atomic_interaction_cutoff() const; 

private:

};

--end_namespace--

#endif //--path--_--class--_hh
