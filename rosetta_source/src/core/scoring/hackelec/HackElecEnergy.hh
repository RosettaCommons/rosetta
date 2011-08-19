// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/HackElecEnergy.hh
/// @brief  Statistically derived rotamer pair potential class declaration
/// @author Phil Bradley, modifed by James Gleixner


#ifndef INCLUDED_core_scoring_hackelec_HackElecEnergy_hh
#define INCLUDED_core_scoring_hackelec_HackElecEnergy_hh

// Unit Headers
#include <core/scoring/hackelec/HackElecEnergy.fwd.hh>

// Package headers
#include <core/scoring/hackelec/ElecAtom.hh>

// Project headers
#include <core/chemical/AtomType.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/trie/RotamerTrieBase.fwd.hh>
// AUTO-REMOVED #include <core/scoring/trie/RotamerTrie.fwd.hh>
// AUTO-REMOVED #include <core/scoring/trie/TrieCountPairBase.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.fwd.hh>

//Auto Headers
#include <core/scoring/trie/TrieCountPairBase.fwd.hh>


// Utility headers


namespace core {
namespace scoring {
namespace hackelec {

struct weight_triple
{
	Real wbb_bb_;
	Real wbb_sc_;
	Real wsc_sc_;
};

///
class HackElecEnergy : public methods::ContextIndependentTwoBodyEnergy  {
public:
	typedef methods::ContextIndependentTwoBodyEnergy  parent;
public:

	///
	HackElecEnergy( methods::EnergyMethodOptions const & options );

	///
	HackElecEnergy( HackElecEnergy const & src );

	void
	initialize(); //james added in option to change max-dis

	/// clone
	virtual
	methods::EnergyMethodOP
	clone() const;

	/// stashes nblist if use_nblist is true
	virtual
	void
	setup_for_minimizing(
		pose::Pose & pose,
		ScoreFunction const & sfxn,
		kinematics::MinimizerMapBase const & min_map
	) const;

	virtual
	void
	setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const;

	virtual
	void
	setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const;

	///
	virtual
	void
	setup_for_packing(
		pose::Pose & pose,
		utility::vector1< bool > const &,
		utility::vector1< bool > const &
	) const;

	// Creates a rotamer trie for the input set of rotamers and stores the trie
	// in the rotamer set.
	virtual
	void
	prepare_rotamers_for_packing(
		pose::Pose const & pose,
		conformation::RotamerSetBase & set ) const;

	// Updates the cached rotamer trie for a residue if it has changed during the course of
	// a repacking
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

	/// @brief Returns true if we're using neighborlist-autoupdate
	virtual
	bool
	minimize_in_whole_structure_context( pose::Pose const & pose ) const;

	virtual
	bool
	defines_score_for_residue_pair(
		conformation::Residue const & res1,
		conformation::Residue const & res2,
		bool res_moving_wrt_eachother
	) const;

	virtual
	bool
	use_extended_residue_pair_energy_interface() const;

	virtual
	void
	residue_pair_energy_ext(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		ResPairMinimizationData const & min_data,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;

	void
	setup_for_minimizing_for_residue_pair(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		kinematics::MinimizerMapBase const & minmap,
		ResSingleMinimizationData const & res1_data_cache,
		ResSingleMinimizationData const & res2_data_cache,
		ResPairMinimizationData & data_cache
	) const;

	virtual
	void
	eval_residue_pair_derivatives(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		ResSingleMinimizationData const &,
		ResSingleMinimizationData const &,
		ResPairMinimizationData const & min_data,
		pose::Pose const & pose, // provides context
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & r1_atom_derivs,
		utility::vector1< DerivVectorPair > & r2_atom_derivs
	) const;

	virtual
	void
	backbone_backbone_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;

	virtual
	void
	backbone_sidechain_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;

	virtual
	void
	sidechain_sidechain_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;

	void
	finalize_total_energy(
		pose::Pose & pose,
		ScoreFunction const &,
		EnergyMap & totals
	) const;

	// was private, but outside world wants to call this one:
	inline
	Real
	eval_atom_atom_hack_elecE(
		Vector const & i_xyz,
		Real const i_charge,
		Vector const & j_xyz,
		Real const j_charge
	) const;


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
	evaluate_rotamer_pair_energies(
		conformation::RotamerSetBase const & set1,
		conformation::RotamerSetBase const & set2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap const & weights,
		ObjexxFCL::FArray2D< core::PackerEnergy > & energy_table
	) const;


	//@brief overrides default rotamer/background energy calculation and uses
	// the trie-vs-trie algorithm instead
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


	/*virtual
	void
	eval_atom_derivative(
		id::AtomID const & atom_id,
		pose::Pose const & pose,
		kinematics::DomainMap const & domain_map,
		ScoreFunction const &,
		EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	) const;*/

	virtual
	bool
	defines_intrares_energy( EnergyMap const & /*weights*/ ) const { return false; }

	/// @brief Interface function for class NeighborList.
	etable::count_pair::CountPairFunctionCOP
	get_intrares_countpair(
		conformation::Residue const &,
		pose::Pose const &,
		ScoreFunction const &
	) const;

	/// @brief Interface function for class NeighborList.
	etable::count_pair::CountPairFunctionCOP
	get_count_pair_function(
		Size const,
		Size const,
		pose::Pose const &,
		ScoreFunction const &
	) const;

	etable::count_pair::CountPairFunctionCOP
	get_count_pair_function(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2
	) const;


	virtual
	Distance
	atomic_interaction_cutoff() const;

	virtual
	bool
	divides_backbone_and_sidechain_energetics() const
	{
		return true;
	}

	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & context_graphs_required ) const;

public:
	 // Hooks for trie algorithms

	///  How close two heavy atoms have to be such that their hydrogen atoms might interact, squared.
	Real
	hydrogen_interaction_cutoff2() const
	{
		return ( hydrogen_interaction_cutoff() )*( hydrogen_interaction_cutoff() );
	}

	/// How close two heavy atoms have to be such that their hydrogen atoms might interact
	Real
	hydrogen_interaction_cutoff() const
	{
		return ( max_dis + 2*core::chemical::MAX_CHEMICAL_BOND_TO_HYDROGEN_LENGTH );
	}

	inline
	Energy heavyatom_heavyatom_energy(
		ElecAtom const & at1,
		ElecAtom const & at2,
		DistanceSquared & d2,
		Size & /*path_dist*/
	) const
	{
		return hackelec_weight(at1.isbb(),at2.isbb()) * eval_atom_atom_hack_elecE( at1.xyz(), at1.charge(), at2.xyz(), at2.charge(), d2  );
	}

	inline
	Energy heavyatom_hydrogenatom_energy(
		ElecAtom const & at1,
		ElecAtom const & at2,
		Size & /*path_dist*/
	) const
	{
		return hackelec_weight(at1.isbb(),at2.isbb()) * eval_atom_atom_hack_elecE( at1.xyz(), at1.charge(), at2.xyz(), at2.charge()  );
	}

	inline
	Energy hydrogenatom_heavyatom_energy(
		ElecAtom const & at1,
		ElecAtom const & at2,
		Size & /*path_dist*/
	) const
	{
		return hackelec_weight(at1.isbb(),at2.isbb()) * eval_atom_atom_hack_elecE( at1.xyz(), at1.charge(), at2.xyz(), at2.charge()  );
	}

	inline
	Energy hydrogenatom_hydrogenatom_energy(
		ElecAtom const & at1,
		ElecAtom const & at2,
		Size & /*path_dist*/
	) const
	{
		return hackelec_weight(at1.isbb(),at2.isbb()) * eval_atom_atom_hack_elecE( at1.xyz(), at1.charge(), at2.xyz(), at2.charge()  );
	}

	/// @brief This has to go
	inline
	Real
	hackelec_weight( bool at1isbb, bool at2isbb ) const {
		return ( at1isbb ? ( at2isbb ? wbb_bb_ : wbb_sc_ ) : ( at2isbb ? wbb_sc_ : wsc_sc_ ) );
	}


	inline
	Real
	hackelec_weight( bool at1isbb, bool at2isbb, weight_triple const & wts ) const {
		return ( at1isbb ? ( at2isbb ? wts.wbb_bb_ : wts.wbb_sc_ ) : ( at2isbb ? wts.wbb_sc_ : wts.wsc_sc_ ) );
	}


/// Private methods
private:

	void
	setup_weight_triple(
		EnergyMap const & weights,
		weight_triple & wttrip
	) const;

	inline
	Real
	eval_atom_atom_hack_elecE(
		Vector const & i_xyz,
		Real const i_charge,
		Vector const & j_xyz,
		Real const j_charge,
		DistanceSquared & d2
	) const;

protected:
	inline
	Real
	eval_dhack_elecE_dr_over_r(
		Real const dis2,
		Real const q1,
		Real const q2
	) const;

private:
	/* ALL RNA specific electrostatics have been moved to the RNAHackElec class
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
	*/

	trie::RotamerTrieBaseOP
	create_rotamer_trie(
		conformation::RotamerSetBase const & rotset,
		pose::Pose const & pose
	) const;

	trie::RotamerTrieBaseOP
	create_rotamer_trie(
		conformation::Residue const & res,
		pose::Pose const & // will be need to create tries for disulfides
	) const;

	trie::TrieCountPairBaseOP
	get_count_pair_function_trie(
		conformation::RotamerSetBase const & set1,
		conformation::RotamerSetBase const & set2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn
	) const;


	trie::TrieCountPairBaseOP
	get_count_pair_function_trie(
		conformation::Residue const & res1,
		conformation::Residue const & res2,
		trie::RotamerTrieBaseCOP trie1,
		trie::RotamerTrieBaseCOP trie2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn
	) const;

	inline
	Real
	score_atom_pair(
		conformation::Residue const & res1,
		conformation::Residue const & res2,
		Size const at1,
		Size const at2,
		EnergyMap & emap,
		Real const cpweight,
		Real & d2
	) const;

	void
	set_nres_mono(
		core::pose::Pose const & pose
	) const;
	
	bool
	monomer_test(
		Size irsd,
		Size jrsd
	) const;
	

private:

	/////////////////////////////////////////////////////////////////////////////
	// data
  // These variables were moved here from HackElecEnergy.cc and made public so
	// that distance cuttoff could be changed and the r dependence of the dielectric
	// could be changed by option flags.  They are all defined in HackElecEnergy.cc
	//		James Gleixner and Liz Kellogg - April 2011
	/////////////////////////////////////////////////////////////////////////////

	/**
	static Real const max_dis;
	static Real const max_dis2;
	static Real const min_dis;
	static Real const min_dis2;
	**/

	Real max_dis;
	Real max_dis2;
	Real min_dis;
	Real min_dis2;
	bool r_option;

	bool exclude_protein_protein_;
	bool exclude_monomer_;
	bool exclude_DNA_DNA_;

	//mutable Real hackelec_weight_; // used during trie-vs-trie algorithm
	mutable Real wbb_bb_;
	mutable Real wbb_sc_;
	mutable Real wsc_sc_;

	/**
	static Real const C0_;
	static Real const die_;
	static Real const C1_;
	static Real const C2_;
	static Real const min_dis_score_;
	static Real const dEfac_;
	**/
	Real C0_;
	Real die_;
	Real C1_;
	Real C2_;
	Real min_dis_score_;
	Real dEfac_;
	
	mutable Size nres_monomer_;

	virtual
	core::Size version() const;

};

inline
Real
HackElecEnergy::eval_atom_atom_hack_elecE(
	Vector const & i_xyz,
	Real const i_charge,
	Vector const & j_xyz,
	Real const j_charge
) const
{

	Real d2;
	if( r_option ) {
		d2 = i_xyz.distance( j_xyz );
		//debug_output
		//std::cout << "Using distance " << std::endl;
	} else {
		d2 = i_xyz.distance_squared( j_xyz );
		//debug_output
		//std::cout << "Using distance squared " << std::endl;
	}
	Real energy(0.0);
	//debug_output
	//std::cout << "d2: " << d2 << "max_dis2: " << max_dis2 << std::endl;
	if ( d2 <= max_dis2 ) {   //max_dis2 is defined in HackElecEnergy.cc
		if ( d2 < min_dis2 ) energy = i_charge * j_charge * min_dis_score_;
		else energy = i_charge * j_charge * ( C1_ / d2 - C2_ );

	}
//debug_output
	//std::cout << "Energy: " << energy << std::endl;
	return energy;
}

inline
Real
HackElecEnergy::eval_atom_atom_hack_elecE(
	Vector const & i_xyz,
	Real const i_charge,
	Vector const & j_xyz,
	Real const j_charge,
	DistanceSquared & d2
) const
{
	if (r_option) {
		d2 = i_xyz.distance( j_xyz );
	}
	else {
		d2 = i_xyz.distance_squared( j_xyz );
	}
	Real energy(0.0);
	if ( d2 <= max_dis2 ) {
		if ( d2 < min_dis2 ) energy = i_charge * j_charge * min_dis_score_;
		else energy = i_charge * j_charge * ( C1_ / d2 - C2_ );
	}
	return energy;
}

}
}
}

#endif
