// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/vdwaals/VDW_Energy.hh
/// @brief  Low-resolution (centroid) repulsive energy
/// @author Phil Bradley


#ifndef INCLUDED_core_scoring_vdwaals_VDW_Energy_hh
#define INCLUDED_core_scoring_vdwaals_VDW_Energy_hh

// Unit Headers
#include <core/scoring/vdwaals/VDW_Energy.fwd.hh>
#include <core/scoring/vdwaals/VDWTrie.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/AtomVDW.hh>
#include <core/scoring/ScoreFunction.fwd.hh>


// Project headers
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


// Utility headers


namespace core {
namespace scoring {
namespace vdwaals {


class VDW_Energy : public methods::ContextIndependentTwoBodyEnergy  {
public:
	typedef methods::ContextIndependentTwoBodyEnergy  parent;
public:

	/// @brief  C-tor, requires options to tell us the atom_type_set_name for the AtomVDW data
	VDW_Energy( methods::EnergyMethodOptions const & options );

private:

	/// @brief Private unimplemented default constructor - must be initialized with options object.
	VDW_Energy();

public:

	/// clone
	virtual
	methods::EnergyMethodOP
	clone() const;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	virtual
	void
	setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const;

	virtual
	void
	setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const;

	// Creates a rotamer trie for the input set of rotamers and stores the trie
	// in the rotamer set.
	virtual
	void
	prepare_rotamers_for_packing(
		pose::Pose const & pose,
		conformation::RotamerSetBase & set
	) const;

	virtual
	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const &,
		EnergyMap & emap
	) const;

	// requires that residue-pair derivative evaluation be implemented, first.
	//bool
	//minimize_in_whole_structure_context( pose::Pose const & ) const { return false; }

	virtual
	void
	eval_intrares_energy(
		conformation::Residue const &,
		pose::Pose const &,
		ScoreFunction const &,
		EnergyMap &
	) const {}


	virtual
	void
	eval_atom_derivative(
		id::AtomID const & atom_id,
		pose::Pose const & pose,
		kinematics::DomainMap const & domain_map,
		ScoreFunction const & scorefxn,
		EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	) const;

	//  virtual
	//  void
	//  eval_atom_derivative(
	//   id::AtomID const & atom_id,
	//   pose::Pose const & pose,
	//   ScoreFunction const &,
	//   EnergyMap const & weights,
	//   Vector & F1,
	//   Vector & F2
	//  ) const;

	virtual
	bool
	defines_intrares_energy( EnergyMap const & /*weights*/ ) const { return false; }

	virtual
	Distance
	atomic_interaction_cutoff() const;

	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & context_graphs_required ) const;

	/// Trie related functions


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


	virtual
	core::Size version() const;


	/// @brief inlined version of the atom-pair energy evaluation for use in the trie
	inline
	Energy
	atom_pair_energy( VDWAtom const & at1, VDWAtom const & at2, Real & d2 ) const {
		Real const vdw12 = atom_vdw_( at1.atom_type() )[ at2.atom_type() ];
		d2 = at1.xyz().distance_squared( at2.xyz() );
		Real const clash = vdw12 - d2;
		//std::cout << "VDW atom_pair_energy: " << at1 << " " << at2 << " d2: " << d2 << " vdw12: " << vdw12 << " clash: " << clash << std::endl;
		if ( clash > 0.0 ) {
			return vdw_scale_factor_ * ( clash*clash ) / vdw12;
		} else {
			return 0.0;
		}
	}


private:

	void
	calculate_hydrogen_interaction_cutoff();

	VDWRotamerTrieOP
	create_rotamer_trie(
		conformation::RotamerSetBase const & rotset
	) const;

	trie::TrieCountPairBaseOP
	get_count_pair_function_trie(
		conformation::RotamerSetBase const & set1,
		conformation::RotamerSetBase const & set2,
		pose::Pose const & pose
	) const;

	trie::TrieCountPairBaseOP
	get_count_pair_function_trie(
		conformation::Residue const & res1,
		conformation::Residue const & res2,
		trie::RotamerTrieBaseCOP trie1,
		trie::RotamerTrieBaseCOP trie2
	) const;

	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////

private:

	// const-ref to scoring database
	AtomVDW const & atom_vdw_;

	std::string const atom_type_set_name_;

	Real const vdw_scale_factor_;
	Real hydrogen_interaction_cutoff2_; // square distance at which the hydrogen-hydrogen interaction energy goes to zero

};

class VDWTrieEvaluator
{
public:
	VDWTrieEvaluator(
		VDW_Energy const & vdw_,
		Real const atomic_interaction_cutoff_,
		Real const hydrogen_interaction_cutoff2_,
		Real const vdw_weight_ );

	Distance
	atomic_interaction_cutoff() const;

	Real
	hydrogen_interaction_cutoff2() const;

	Real
	vdw_weight() const;

	inline
	Energy heavyatom_heavyatom_energy(
		VDWAtom const & at1,
		VDWAtom const & at2,
		DistanceSquared & d2,
		Size & /*path_dist*/
	) const
	{
		return weighted_atom_pair_energy( at1,at2,d2 );
	}

	inline
	Energy heavyatom_hydrogenatom_energy(
		VDWAtom const & at1,
		VDWAtom const & at2,
		Size & /*path_dist*/
	) const
	{
		Real d2;
		return weighted_atom_pair_energy( at1,at2,d2 );
	}

	inline
	Energy hydrogenatom_heavyatom_energy(
		VDWAtom const & at1,
		VDWAtom const & at2,
		Size & /*path_dist*/
	) const
	{
		Real d2;
		return weighted_atom_pair_energy( at1,at2,d2 );
	}

	inline
	Energy hydrogenatom_hydrogenatom_energy(
		VDWAtom const & at1,
		VDWAtom const & at2,
		Size & /*path_dist*/
	) const
	{
		Real d2;
		return weighted_atom_pair_energy( at1,at2,d2 );
	}

	inline
	Energy
	weighted_atom_pair_energy(
		VDWAtom const & at1,
		VDWAtom const & at2,
		Real & d2
	) const
	{
		return vdw_weight_ * vdw_.atom_pair_energy( at1, at2, d2 );
	}

private:
	VDW_Energy const & vdw_;
	Real const atomic_interaction_cutoff_;
	Real const hydrogen_interaction_cutoff2_;
	Real const vdw_weight_;

};


}
}
}

#endif
