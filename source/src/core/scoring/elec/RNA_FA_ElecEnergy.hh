// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/RNA_FA_ElecEnergy.hh
/// @brief  Electrostatics for RNA
/// @author Rhiju Das


#ifndef INCLUDED_core_scoring_elec_RNA_FA_ElecEnergy_hh
#define INCLUDED_core_scoring_elec_RNA_FA_ElecEnergy_hh

/// Unit Headers
#include <core/scoring/elec/RNA_FA_ElecEnergy.fwd.hh>

/// Package Headers
#include <core/scoring/elec/FA_ElecEnergy.hh>

#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/scoring/methods/EnergyMethod.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>

#include <ObjexxFCL/FArray2D.fwd.hh>

#include <utility/vector1.hh>



namespace core {
namespace scoring {
namespace elec {

///
class RNA_FA_ElecEnergy : public FA_ElecEnergy  {
public:
	typedef FA_ElecEnergy parent;
	typedef ContextIndependentTwoBodyEnergy grandparent;

public:

	///
	RNA_FA_ElecEnergy( methods::EnergyMethodOptions const & options );

	///
	RNA_FA_ElecEnergy( RNA_FA_ElecEnergy const & src );


	/// clone
	virtual
	methods::EnergyMethodOP
	clone() const;

	virtual
	void
	setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const;

	virtual
	void
	setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const;

	///
	virtual
	void
	setup_for_packing( pose::Pose & pose, utility::vector1< bool > const &, utility::vector1< bool > const & ) const;


	/// @brief overrides parent class implementation which would have
	/// created several tries
	virtual
	void
	prepare_rotamers_for_packing(
		pose::Pose const & pose,
		conformation::RotamerSetBase & set ) const;

	/// @brief overrides parent class implementation which would have
	/// updated a trie
	virtual
	void
	update_residue_for_packing( pose::Pose & pose, Size resid ) const;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	virtual
	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const &,
		EnergyMap & emap
	) const;

	/// @brief Returns "true" because this energy method has not been updated to
	/// use the new derivative evaluation machinery.  Note that this class requires
	/// the definition of this method because it's parent class, FA_ElecEnergy,
	/// HAS been updated to use the new derivative evaluation machinery, and,
	/// if this class did not return "true", it would be asked to evaluate derivatives
	/// in ways it cannot yet evaluate them in.
	bool
	minimize_in_whole_structure_context( pose::Pose const & ) const { return true; }


	/// @brief Jan 10, 2012. Parin Sripakdeevon (sripakpa@stanford.edu)
	/// Returns "false" to overwrite the behavior in the parent class (FA_ElecEnergy)! 
	virtual
	bool
	use_extended_residue_pair_energy_interface() const { return false; }

	virtual
	void
	eval_intrares_energy(
		conformation::Residue const &,
		pose::Pose const &,
		ScoreFunction const &,
		EnergyMap &
	) const {}

	//@brief overrides default rotamer/rotamer energy calculation
	// and overrides the parent class trie implementatoin
	virtual
	void
	evaluate_rotamer_pair_energies(
		conformation::RotamerSetBase const & set1,
		conformation::RotamerSetBase const & set2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap const & weights,
		ObjexxFCL::FArray2D< core::PackerEnergy > & energy_table
	) const;


	//@brief overrides default rotamer/background energy calculation
	// and overrides the parent class trie implementatoin
	virtual
	void
	evaluate_rotamer_background_energies(
		conformation::RotamerSetBase const & set,
		conformation::Residue const & residue,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap const & weights,
		utility::vector1< core::PackerEnergy > & energy_vector
	) const;


 	virtual
 	void
 	eval_atom_derivative(
 		id::AtomID const & atom_id,
 		pose::Pose const & pose,
		kinematics::DomainMap const & domain_map,
 		ScoreFunction const &,
 		EnergyMap const & weights,
 		Vector & F1,
 		Vector & F2
 	) const;

	virtual
	bool
	defines_intrares_energy( EnergyMap const & /*weights*/ ) const { return false; }

	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & context_graphs_required ) const;

public:


	Real
	residue_pair_energy_RNA(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		EnergyMap & emap
	) const;

	void
	eval_atom_derivative_RNA(
	   conformation::Residue const & rsd1,
		 Size const & i,
		 conformation::Residue const & rsd2,
		 EnergyMap const & weights,
		 Vector & F1,
		 Vector & F2
	) const;
virtual
core::Size version() const;

};


}
}
}

#endif
