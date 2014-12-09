// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/downstream/SecMatchResiduePairEvaluator.hh
/// @brief
/// @author Florian Richter, floric@u.washington.edu, june 09


#ifndef INCLUDED_protocols_match_downstream_SecMatchResiduePairEvaluator_hh
#define INCLUDED_protocols_match_downstream_SecMatchResiduePairEvaluator_hh

// Unit headers
#include <protocols/match/downstream/SecMatchResiduePairEvaluator.fwd.hh>


// Project headers
#include <core/types.hh>
// AUTO-REMOVED #include <core/chemical/ResidueType.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
// AUTO-REMOVED #include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>


// C++ headers

namespace protocols {
namespace match {
namespace downstream {


/// @brief base for class used by secondary matcher
/// to determine whether a given residue to be matched
/// interacts satisfactorily with a target residue
class SecMatchResiduePairEvaluator : public utility::pointer::ReferenceCount
{
public:
	typedef core::Size Size;
	typedef core::Real Real;

public:
	SecMatchResiduePairEvaluator();

	virtual ~SecMatchResiduePairEvaluator();

	/// @brief Returns true if the interaction between the two residues satisifies the
	/// secondary match requirement.
	/// candidate_res: the rotamer of the residue trying to be placed
	/// target_res: the previously placed residue
	virtual
	bool
	evaluate_residues(
		core::conformation::Residue const & candidate_res,
		core::conformation::Residue const & target_res
	) const = 0;

	/// @brief Returns true if all coordinates of the target residue are required
	/// in order to evaluate the interaction between the candidate and the target residues.
	virtual
	bool
	require_all_target_residue_atom_coordinates() const = 0;

	/// @brief If require_all_target_residue_atom_coordinates() returns false, then
	/// this method should return true for the atoms on the target residue that the
	/// evaluator requires.
	virtual
	bool
	require_target_atom_coordinate( Size target_atom_id ) const = 0;

	/// @brief Are there atoms of the candidate residue that must be within
	/// some cutoff distance of a given atom on the target residue?  Base
	/// class implementation returns false.
	virtual
	bool
	require_candidate_residue_atoms_to_lie_near_target_atom( Size target_atom_id ) const;

	/// @brief Return a list of atom indices on the candidate residue; if all atoms
	/// in this list are further than max_separation_dist_to_target_atom() away
	/// from the target_atom_id atom for a given pair of conformations of the
	/// target_residue and the candidate_residue, then this evaluator will return false
	/// in the call to evaluate( candidate_residue, target_residue ).
	/// This list will allow the SecondaryMatcher to (conservatively!) prune
	/// conformations of the candidate_residue from consideration.  The base class implements
	/// a noop -- it returns an empty list.
	virtual
	utility::vector1< Size >
	candidate_res_atoms_reqd_near_target_atom(
		Size target_atom_id
	) const;

	/// @brief Return the maximum separation distance that any of the match-residue
	/// atoms identified by the function match_atoms_reqd_near_target_atom
	/// may be from a particular atom on the target residue.  Returns a negative value
	/// if there is no requirement that any atom be within a certain radius of the
	/// target atom.  The base class implementation returns -1.0.
	virtual
	Real
	max_separation_dist_to_target_atom( Size target_atom_id ) const;

};


}
}
}

#endif
