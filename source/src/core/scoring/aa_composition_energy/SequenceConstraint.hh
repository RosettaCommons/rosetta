// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/aa_composition_energy/SequenceConstraint.hh
/// @brief Headers for a constraint for constraining sequences to have a desired amino acid composition, analogous to a geometric constraint.
/// @details The corresponding energy term for this constraint is the AACompositionEnergy (aa_composition in wts files).
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_core_scoring_aa_composition_energy_SequenceConstraint_hh
#define INCLUDED_core_scoring_aa_composition_energy_SequenceConstraint_hh

#include <core/scoring/aa_composition_energy/SequenceConstraint.fwd.hh>

#include <core/scoring/constraints/Constraint.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/scoring/aa_composition_energy/AACompositionEnergySetup.fwd.hh>
#include <core/scoring/func/Func.fwd.hh>
#include <core/scoring/func/XYZ_Func.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>

#include <core/scoring/ScoreType.hh>
#include <core/id/AtomID.hh>

#include <string>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace aa_composition_energy {

class SequenceConstraint : public core::scoring::constraints::Constraint {

public: //Constructor, destructor, copy, clone:

	/// @brief Constructor
	SequenceConstraint( core::scoring::ScoreType const & t );

	/// @brief Copy constructor
	SequenceConstraint( SequenceConstraint const &src );

	/// @brief Destructor
	~SequenceConstraint();


public: //Functions that all Constraints must implement, but which aren't really applicable for a seqeunce constraint:

	/// @brief Implemented because the base class requires it; not used by a sequence constraint.
	virtual void score( core::scoring::func::XYZ_Func const &, EnergyMap const &, EnergyMap & ) const { return; }

	/// @brief Implemented because the base class requires it; not used by a sequence constraint.
	virtual void fill_f1_f2( core::id::AtomID const &, core::scoring::func::XYZ_Func const &, Vector &, Vector &, EnergyMap const & ) const { return; }

	/// @brief Implemented because the base class requires it; not used by a sequence constraint.
	/// @details Always returns zero.
	virtual core::Size natoms() const { return 0; }

	/// @brief Implemented because the base class requires it; not used by a sequence constraint.
	/// @details Always returns AtomID(0,0).
	virtual core::id::AtomID const &atom( core::Size const ) const { return dummy_atomid_; }


private:
	// Member variables

	/// @brief Silly fake atomID necessary for the atom() function that all Constraints have to implement.
	/// @details Not used for anything.
	core::id::AtomID dummy_atomid_;

};

} // aa_composition_energy
} // scoring
} // core

#endif
