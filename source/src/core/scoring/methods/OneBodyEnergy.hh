// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author


#ifndef INCLUDED_core_scoring_methods_OneBodyEnergy_hh
#define INCLUDED_core_scoring_methods_OneBodyEnergy_hh

// Unit headers
#include <core/scoring/methods/OneBodyEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/MinimizationData.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Project headers
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>

// Utility headers
#include <utility/vector1.fwd.hh>

#include <core/scoring/DerivVectorPair.fwd.hh>
#include <utility/vector1.hh>


#ifdef PYROSETTA
	#include <core/scoring/DerivVectorPair.hh>
#endif


namespace core {
namespace scoring {
namespace methods {

class OneBodyEnergy : public EnergyMethod {
public:
	typedef EnergyMethod parent;

public:
	/// @brief Constructor with an EnergyMethodCreator to inform the EnergyMethod
	/// parent which ScoreTypes this EnergyMethod is responsible for computing.
	OneBodyEnergy( EnergyMethodCreatorOP );

	// @brief dstor;
	virtual ~OneBodyEnergy();

	/// @brief Evaluate the one-body energies for a particular residue, in the context of a
	/// given Pose, and increment those energies in the input Emap (do not overwrite them).
	virtual
	void
	residue_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		EnergyMap & emap
	) const = 0;

	/// @brief During minimization, energy methods are allowed to decide that they say nothing
	/// about a particular residue (e.g. no non-zero energy) and as a result they will not be queried for
	/// a derivative or an energy.  The default behavior is to return "true" for all residues.
	virtual
	bool
	defines_score_for_residue(
		conformation::Residue const &
	) const;

	/// @brief Rely on the extended version of the residue_energy function during score-function
	/// evaluation in minimization? The extended version (below) takes a ResSingleMinimizationData.
	/// Return 'true' for the extended version.  The default method implemented in this class returns 'false'
	virtual
	bool
	use_extended_residue_energy_interface() const;

	/// @brief Evaluate the one-body energies for a particular residue, in the context of a
	/// given Pose, and with the help of a piece of cached data for minimization, increment those
	/// one body energies into the input EnergyMap.  The calling function must guarantee that this
	/// EnergyMethod has had the opportunity to update the input ResSingleMinimizationData object
	/// for the given residue in a call to setup_for_minimizing_for_residue before this function is
	/// invoked. This function should not be called unless the use_extended_residue_energy_interface()
	/// method returns "true".  Default implementation provided by this base class calls
	/// utility::exit(). The Pose merely serves as context, and the input residue is not required
	/// to be a member of the Pose.
	virtual
	void
	residue_energy_ext(
		conformation::Residue const & rsd,
		ResSingleMinimizationData const & min_data,
		pose::Pose const & pose,
		EnergyMap & emap
	) const;

	/// @brief Called at the beginning of minimization, allowing this energy method to cache data
	/// pertinent for a single residue in the the ResSingleMinimizationData that is used for a
	/// particular residue in the context of a particular Pose.  This base class provides a noop
	/// implementation for this function if there is nothing that the derived class needs to perform
	/// in this setup phase.   The Pose merely serves as context, and the input residue is not
	/// required to be a member of the Pose.
	virtual
	void
	setup_for_minimizing_for_residue(
		conformation::Residue const & rsd,
		pose::Pose const & ,
		ScoreFunction const & ,
		kinematics::MinimizerMapBase const & ,
		ResSingleMinimizationData &
	) const;

	/// @brief Does this EnergyMethod require the opportunity to examine the residue before scoring begins?  Not
	/// all energy methods would.  The ScoreFunction will not ask energy methods to examine residues that are uninterested
	/// in doing so.
	virtual
	bool
	requires_a_setup_for_scoring_for_residue_opportunity( pose::Pose const & pose ) const;

	/// @brief Do any setup work should the coordinates of this residue, who is still guaranteed to be
	/// of the same residue type as when setup_for_minimizing_for_residue was called, have changed so dramatically
	/// as to possibly require some amount of setup work before scoring should proceed
	virtual
	void
	setup_for_scoring_for_residue(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		ResSingleMinimizationData & min_data
	) const;


	/// @brief Does this EnergyMethod require the opportunity to examine the residue before derivative evaluation begins?  Not
	/// all energy methods would.  The ScoreFunction will not ask energy methods to examine residues that are uninterested
	/// in doing so.
	virtual
	bool
	requires_a_setup_for_derivatives_for_residue_opportunity( pose::Pose const & pose ) const;

	/// @brief Do any setup work necessary before evaluating the derivatives for this residue
	virtual
	void
	setup_for_derivatives_for_residue(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		ResSingleMinimizationData & min_data
	) const;

	/// @brief Evaluate the derivative for an atom in a residue in the context of a particular pose,
	/// and increment the F1 and F2 vectors.  This base class provides a default noop implementation
	/// of this function. The calling function must guarantee that this EnergyMethod has had the
	/// opportunity to update the input ResSingleMinimizationData object for the given residue
	/// in a call to prepare_for_minimization before this function is invoked.   The Pose merely
	/// serves as context, and the input residue is not required to be a member of the Pose.
	/// DEPRECATED -- too slow.
	/*virtual
	void
	eval_atom_derivative_for_residue(
	Size const atom_index,
	conformation::Residue const & rsd,
	ResSingleMinimizationData const & min_data,
	pose::Pose const & pose, // provides context
	kinematics::DomainMap const & domain_map,
	ScoreFunction const & sfxn,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
	) const;*/

	/// @brief Evaluate the derivatives for all atoms on this residue and increment them
	/// into the input atom_derivs vector1.  The calling function must guarantee that
	/// setup for derivatives is called before this function is, and that the atom_derivs
	/// vector contains at least as many entries as there are atoms in the input Residue.
	/// This base class provides a default noop implementation of this function.
	virtual
	void
	eval_residue_derivatives(
		conformation::Residue const & rsd,
		ResSingleMinimizationData const & min_data,
		pose::Pose const & pose,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & atom_derivs
	) const;

	/// @brief Use the dof_derivative interface for this energy method when
	/// calculating derivatives?  It is possible to define both dof_derivatives and
	/// atom-derivatives; they are not mutually exclusive.
	virtual
	bool
	defines_dof_derivatives( pose::Pose const & p ) const;

	/// @brief Evaluate the DOF derivative for a particular residue.  The Pose merely serves as context,
	/// and the input residue is not required to be a member of the Pose.
	virtual
	Real
	eval_residue_dof_derivative(
		conformation::Residue const & rsd,
		ResSingleMinimizationData const & min_data,
		id::DOF_ID const & dof_id,
		id::TorsionID const & torsion_id,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap const & weights
	) const;


};

} // methods
} // scoring
} // core


#endif // INCLUDED_core_scoring_EtableEnergy_HH
