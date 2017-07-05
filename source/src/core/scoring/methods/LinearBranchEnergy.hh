// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/EtableEnergy.hh
/// @brief  Etable energy method class declaration
/// @author Phil Bradley
/// @author Andrew Leaver-Fay
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDED_core_scoring_methods_LinearBranchEnergy_hh
#define INCLUDED_core_scoring_methods_LinearBranchEnergy_hh

// Unit headers
#include <core/scoring/methods/LinearBranchEnergy.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/kinematics/ShortestPathInFoldTree.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/methods/WholeStructureEnergy.hh>

// Project headers

// Third-party headers
#include <boost/scoped_ptr.hpp>

// C++ headers

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {

/// @brief LinearBranchEnergy class iterates across all residues in finalize()
/// and determines the penalty between branch-connected residues by how much their
/// psueduo atoms do not align (if they have them).
///
/// @details Calculates linear_branch_conn.
///  linear_branch_conn measures 3 distances (branch variants with cutpoint-like atoms
///  must be added to pose)
///
///  For ideal poses, this score should be very close to 0.  Real PDBs, however have bond length and angle
///   deviations that will cause this score to be fairly high.
///
class LinearBranchEnergy : public WholeStructureEnergy {
public:
	typedef WholeStructureEnergy parent;

	// @brief Creates a new LinearBranchEnergy with the default allowable sequence
	// separation (+inf)
	LinearBranchEnergy();

	/// @brief The auto-generated copy constructor does not properly handle smart
	/// pointer types, so we must explicitly define our own.
	LinearBranchEnergy(const LinearBranchEnergy&);

	/// @brief The auto-generated operator=() method does not properly handle pointer types.
	LinearBranchEnergy& operator=(const LinearBranchEnergy&);

	// @brief Releases resources associated with an instance.
	~LinearBranchEnergy();

	/// clone
	virtual EnergyMethodOP clone() const {
		return EnergyMethodOP( new LinearBranchEnergy(*this) );
	}

	/// called at the end of energy evaluation
	virtual void finalize_total_energy(
		pose::Pose & pose,
		ScoreFunction const &,
		EnergyMap & totals) const;

	/// called during gradient-based minimization inside dfunc
	/**
	F1 and F2 are not zeroed -- contributions from this atom are
	just summed in
	**/
	virtual void eval_atom_derivative(id::AtomID const & id,
		pose::Pose const & pose,
		kinematics::DomainMap const & domain_map,
		ScoreFunction const & sfxn,
		EnergyMap const & weights,
		Vector & F1,
		Vector & F2) const;

	virtual void indicate_required_context_graphs( utility::vector1< bool > & ) const;

private:

	// Initialization routine common to both constructor
	void initialize(Size allowable_sequence_sep);

	// Maximum allowable sequence separation permitted for scoring
	Size allowable_sequence_sep_;

	virtual
	core::Size version() const;
};

} // methods
} // scoring
} // core
#endif
