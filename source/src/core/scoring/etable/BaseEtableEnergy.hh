// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/etable/EtableEnergy.hh
/// @brief  Etable energy method class declaration
/// @author Phil Bradley
/// @author Andrew Leaver-Fay
/// @author Oliver Lange

///***********************************************************************
//++++++++++++++++++++++++++++ WARNING +++++++++++++++++++++++++++++++++++
//
//the CoarseEtable is currently not threadsafe, since it has the
//"mutable" element seq_dist...
//
//could move seq_dist into CoarseEtableEnergy --> then at least the
//CoarseEtable can be shared between threads
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//************************************************************************/


#ifndef INCLUDED_core_scoring_etable_BaseEtableEnergy_hh
#define INCLUDED_core_scoring_etable_BaseEtableEnergy_hh

// Unit headers
// #include <core/scoring/etable/BaseEtableEnergy.fwd.hh> -- file doesn't exist?

// Package headers
#include <core/scoring/etable/count_pair/types.hh>
#include <core/scoring/trie/TrieCountPairBase.fwd.hh>
#include <core/scoring/etable/etrie/EtableTrie.fwd.hh>


#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/methods/EnergyMethodCreator.fwd.hh>
#include <core/scoring/MinimizationData.hh>
#include <core/scoring/rna/RNA_FullAtomVDW_BasePhosphateCreator.fwd.hh>

// Project headers
#include <core/conformation/Atom.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>

#include <core/kinematics/DomainMap.fwd.hh>

#include <ObjexxFCL/FArray3D.hh>

#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.fwd.hh>
#include <core/scoring/etable/etrie/EtableAtom.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <utility/vector1.hh>

namespace core {
namespace scoring {
namespace etable {

// PyRosetta: Python can't pass reference to value, so we have to create data holding structure for atom_pair_energy* functions.
struct AtomPairEnergy {
	Energy attractive;
	Energy repulsive;
	Energy solvation;
	Energy bead_bead_interaction;
	Real distance_squared;
};


template < class Derived >
class BaseEtableEnergy : public methods::ContextIndependentTwoBodyEnergy
{
public:
	typedef methods::ContextIndependentTwoBodyEnergy parent;
	typedef methods::EnergyMethodOP EnergyMethodOP;

public:

	/// construction with an etable
	BaseEtableEnergy(
		methods::EnergyMethodCreatorOP creator,
		Etable const & etable_in,
		methods::EnergyMethodOptions const & options,
		bool const do_classic_intrares_ = false
	);

	/// @brief explicit copy constructor, since the BaseEtableEnergy now contains OP data
	BaseEtableEnergy( BaseEtableEnergy< Derived > const & );


	/////////////////////////////////////////////////////////////////////////////
	// methods for ContextIndependentTwoBodyEnergies
	/////////////////////////////////////////////////////////////////////////////

	/// @brief The neighborlist-autoupdate algorithm requires that the EtableEnergy
	/// be able to control the definition and update for its atom-neighbors.  This
	/// will bypass the standard neighborlist evaluation inside the ScoreFunction,
	/// avoiding the use the MinimizationGraph.
	virtual
	bool
	minimize_in_whole_structure_context( pose::Pose const & ) const;

	/// @brief stashes nblist if pose.energies().use_nblist_auto_update() is true
	/// This is only invoked, now, if the neighborlist-autoupdate flag is on.
	virtual
	void
	setup_for_minimizing(
		pose::Pose & pose,
		ScoreFunction const & sfxn,
		kinematics::MinimizerMapBase const & min_map
	) const;

	/// @brief check compatibility with atomtypeset
	virtual
	void
	setup_for_scoring( pose::Pose &pose, ScoreFunction const &scfxn ) const;

	virtual
	void
	setup_for_derivatives(
		pose::Pose &pose,
		ScoreFunction const &scfxn
	) const;

	// The EtableEnergy method stores a vector of rotamer trie objects in the Energies
	// object for use in rapid rotamer/background energy calculations.  Overrides default
	// do-nothing behavior.
	virtual
	void
	setup_for_packing( pose::Pose &, utility::vector1< bool > const &, utility::vector1< bool > const & ) const;

	// Creates a rotamer trie for the input set of rotamers and stores the trie
	// in the rotamer set.
	virtual
	void
	prepare_rotamers_for_packing(
		pose::Pose const & pose,
		conformation::RotamerSetBase & set
	) const;

	// Updates the cached rotamer trie for a residue if it has changed during the course of
	// a repacking
	virtual
	void
	update_residue_for_packing( pose::Pose & pose, Size resid ) const;


	count_pair::CountPairFunctionCOP
	get_count_pair_function(
		Size res1,
		Size res2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn
	) const;

	count_pair::CountPairFunctionCOP
	get_count_pair_function(
		conformation::Residue const & res1,
		conformation::Residue const & res2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn
	) const;

	count_pair::CountPairFunctionOP
	get_intrares_countpair(
		conformation::Residue const & res,
		pose::Pose const & pose,
		ScoreFunction const & sfxn
	) const;


	virtual
	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;

	/// APL -- note, new
	virtual
	bool
	use_extended_residue_pair_energy_interface() const;

	/// APL -- note, new
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


	/// APL -- note, new
	virtual
	void
	setup_for_minimizing_for_residue_pair(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		kinematics::MinimizerMapBase const & minmap,
		ResSingleMinimizationData const & res1_data_cache,
		ResSingleMinimizationData const & res2_data_cache,
		ResPairMinimizationData & min_data
	) const;

	/// @brief Does this EnergyMethod require the opportunity to examine the residue before scoring begins?  Not
	/// all energy methods would.  The ScoreFunction will not ask energy methods to examine residues that are uninterested
	/// in doing so.
	virtual
	bool
	requires_a_setup_for_scoring_for_residue_opportunity( pose::Pose const & pose ) const;

	/// @brief Do any setup work should the coordinates of this residue (who is still guaranteed to be
	/// of the same residue type as when setup_for_minimizing_for_residue was called) have changed so dramatically
	/// as to possibly require some amount of setup work before scoring should proceed.
	/// This function is used for both intra-residue setup and pre-inter-residue setup
	virtual
	void
	setup_for_scoring_for_residue(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		ScoreFunction const &,
		ResSingleMinimizationData & min_data
	) const;

	/// @brief Does this EnergyMethod require the opportunity to examine each residue before derivative evaluation begins?  Not
	/// all energy methods would.  The ScoreFunction will not ask energy methods to examine residue pairs that are uninterested
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
		ScoreFunction const &,
		ResSingleMinimizationData & min_data
	) const;

	/// @brief Does this EnergyMethod require the opportunity to examine each residue pair before scoring begins?  Not
	/// all energy methods would.  The ScoreFunction will not ask energy methods to examine residue pairs that are uninterested
	/// in doing so.
	virtual
	bool
	requires_a_setup_for_scoring_for_residue_pair_opportunity( pose::Pose const & pose ) const;

	/// @brief Do any setup work should the coordinates of a pair of residues, who are still guaranteed to be
	/// of the same residue type as when setup_for_minimizing_for_residue was called, have changed so dramatically
	/// as to possibly require some amount of setup work before scoring should proceed
	virtual
	void
	setup_for_scoring_for_residue_pair(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		ResSingleMinimizationData const & minsingle_data1,
		ResSingleMinimizationData const & minsingle_data2,
		pose::Pose const & pose,
		ScoreFunction const &,
		ResPairMinimizationData & data_cache
	) const;

	/// @brief Does this EnergyMethod require the opportunity to examine each residue pair before derivative evaluation begins?  Not
	/// all energy methods would.  The ScoreFunction will not ask energy methods to examine residue pairs that are uninterested
	/// in doing so.
	virtual
	bool
	requires_a_setup_for_derivatives_for_residue_pair_opportunity( pose::Pose const & pose ) const;

	/// @brief Do any setup work necessary before evaluating the derivatives for this residue pair
	virtual
	void
	setup_for_derivatives_for_residue_pair(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		ResSingleMinimizationData const & minsingle_data1,
		ResSingleMinimizationData const & minsingle_data2,
		pose::Pose const & pose,
		ScoreFunction const &,
		ResPairMinimizationData & data_cache
	) const;

	/// APL -- note, new
	/*virtual
	void
	eval_atom_derivative_for_residue_pair(
	Size const atom_index,
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResSingleMinimizationData const & minsingle_data1,
	ResSingleMinimizationData const & minsingle_data2,
	ResPairMinimizationData const & min_data,
	pose::Pose const & pose, // provides context
	kinematics::DomainMap const & domain_map,
	ScoreFunction const & sfxn,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
	) const;*/

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
		utility::vector1< DerivVectorPair > & r1_at_derivs,
		utility::vector1< DerivVectorPair > & r2_at_derivs
	) const;

	/// @brief Evaluates the interaction between the backbone of rsd1 and the
	/// backbone of rsd2 and accumulates the unweighted energy.
	virtual
	void
	backbone_backbone_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;


	/// @brief Evaluates the interaction between the backbone of rsd1 and the
	/// sidechain of rsd2 and accumulates the unweighted energy.
	virtual
	void
	backbone_sidechain_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;

	/// @brief Evaluates the interaction between the sidechain of rsd1 and the
	/// sidechain of rsd2 and accumulates the unweighted energy.
	virtual
	void
	sidechain_sidechain_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;

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

	/// APL -- note, new
	virtual
	bool
	use_extended_intrares_energy_interface() const;

	/// APL -- note, new
	virtual
	void
	eval_intrares_energy_ext(
		conformation::Residue const & rsd,
		ResSingleMinimizationData const & min_data,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;

	/// APL -- note, new
	virtual
	void
	setup_for_minimizing_for_residue(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		kinematics::MinimizerMapBase const & minmap,
		ResSingleMinimizationData & min_data
	) const;

	void
	eval_intrares_derivatives(
		conformation::Residue const & rsd,
		ResSingleMinimizationData const & min_data,
		pose::Pose const & pose,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & atom_derivs
	) const;

	virtual
	void
	bump_energy_full(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;

	virtual
	void
	bump_energy_backbone(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;


	/// called at the end of energy evaluation
	virtual
	void
	finalize_total_energy(
		pose::Pose & pose,
		ScoreFunction const &,
		EnergyMap & totals
	) const;


	/// called during gradient-based minimization inside dfunc
	/**
	F1 and F2 are not zeroed -- contributions from this atom are
	just summed in
	**/
	virtual
	void
	eval_atom_derivative(
		id::AtomID const & id,
		pose::Pose const & pose,
		kinematics::DomainMap const &, // domain_map,
		ScoreFunction const & sfxn,
		EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	) const;


	/////////////////////////////////////////////////////////////////////////////
	// methods specific to EtableEnergy:
	/////////////////////////////////////////////////////////////////////////////
	/// Should/Could these inlined" in a separate .inline.hh?

	// for some reason the virtual-function call to residue_pair_energy does not get passed thru to
	// CoarseEnergyEtable... Use template based call instead...
	///
	inline
	void
	prepare_for_residue_pair(
		Size const res1,
		Size const res2,
		pose::Pose const & pose
	) const {
		static_cast< Derived const* > (this) -> derived_prepare_for_residue_pair(res1,res2,pose);
	}


	inline
	void
	atom_pair_energy(
		conformation::Atom const & atom1,
		conformation::Atom const & atom2,
		Real const weight,
		EnergyMap & emap,
		Real & dsq
	) const {
		static_cast< Derived const & > (*this).interres_evaluator().atom_pair_energy(atom1,atom2,weight,emap,dsq);
	}

	/// @brief for the trie-vs-trie algorithm; could test if the other
	/// atom pair energy function could inline this function to avoid
	/// the table reading code duplication.
	inline
	void
	atom_pair_energy(
		conformation::Atom const & atom1,
		conformation::Atom const & atom2,
		Real const weight,
		Energy & atr,
		Energy & rep,
		Energy & solv,
		Energy & bb,
		Real & dsq
	) const  {
		//  std::cerr << __FILE__<< ' ' << __LINE__ << std::endl;
		///static_cast< Derived const* > (this) -> atom_pair_energy_(atom1,atom2,weight,atr,rep,solv,bb,dsq);
		bb = 0; /// <-- need to kill this parameter
		static_cast< Derived const & > (*this).interres_evaluator().atom_pair_energy( atom1, atom2, weight, atr, rep, solv, dsq );
	}

	/// APL -- consider reinstating this function!
	// PyRosetta friendly version
	inline void atom_pair_energy(
		conformation::Atom const & atom1,
		conformation::Atom const & atom2,
		Real const weight,
		AtomPairEnergy & ape)
	const {
		atom_pair_energy(atom1, atom2, weight, ape.attractive, ape.repulsive, ape.solvation, ape.bead_bead_interaction, ape.distance_squared);
	}


	//inline
	//void
	//pair_energy_H(
	// conformation::Atom const & atom1,
	// conformation::Atom const & atom2,
	// Real weight,
	// Energy & atr,
	// Energy & rep,
	// Energy & solv,
	// Energy & bb
	//) const  {
	// return static_cast< Derived const* > (this) -> pair_energy_H_(atom1,atom2,weight,atr,rep,solv,bb);
	//};

	/// APL -- Consider reinstating this function!
	// PyRosetta friendly version
	// inline
	// void
	// pair_energy_H(
	//  conformation::Atom const & atom1,
	//  conformation::Atom const & atom2,
	//  Real weight,
	//  AtomPairEnergy & ape
	// ) const {
	//  ape.distance_squared = 0.0;
	//  pair_energy_H(atom1, atom2, weight, ape.attractive, ape.repulsive, ape.solvation, ape.bead_bead_interaction);
	// }


	inline
	void
	pair_energy_H(
		conformation::Atom const & atom1,
		conformation::Atom const & atom2,
		Real weight,
		EnergyMap &emap
	) const
	{
		Real d2;
		(static_cast< Derived const * > (this))->interres_evaluator().atom_pair_energy( atom1, atom2, weight, emap, d2 );
	}


	Real
	hydrogen_interaction_cutoff2() const
	{
		return hydrogen_interaction_cutoff2_;
	}

	/// @brief Etable atomic distance cutoff is 5.5 A
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
	void indicate_required_context_graphs( utility::vector1< bool > & /*context_graphs_required*/ ) const;

	inline
	EnergiesCacheableDataType::Enum
	nblist_type() const {
		return ( do_classic_intrares_ ?
			EnergiesCacheableDataType::ETABLE_CLASSIC_INTRARES_NBLIST :
			EnergiesCacheableDataType::ETABLE_NBLIST );
	}

	inline
	min_single_data
	min_single_data_type() const {
		return ( do_classic_intrares_ ? etab_classic_intrares_single_nblist : etab_single_nblist );
	}

	inline
	min_pair_data
	min_pair_data_type() const {
		return ( do_classic_intrares_ ? etab_classic_intrares_pair_nblist : etab_pair_nblist );
	}


	inline
	bool
	calculate_interres( core::conformation::Residue const & rsd1,
		core::conformation::Residue const & rsd2 ) const {
		if ( exclude_DNA_DNA_ && rsd1.is_DNA() && rsd2.is_DNA() ) return false;
		if ( do_classic_intrares_ ) return false;
		return true;
	}

	inline
	bool
	calculate_intrares( core::conformation::Residue const & rsd ) const {
		if ( rsd.is_DNA() && exclude_DNA_DNA_ ) return false;
		if ( !do_classic_intrares_ && rsd.is_protein() && exclude_intra_res_protein_ /* default true */ ) return false;
		return true;
	}

protected: //protected methods that may be used by derived classes

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

	count_pair::CPCrossoverBehavior
	determine_crossover_behavior(
		conformation::Residue const & res1,
		conformation::Residue const & res2,
		pose::Pose const &,
		ScoreFunction const & sfxn
	) const;

	etrie::EtableRotamerTrieOP
	create_rotamer_trie(
		conformation::RotamerSetBase const & rotset,
		pose::Pose const & pose // will be need to create tries for disulfides
	) const;

	etrie::EtableRotamerTrieOP
	create_rotamer_trie(
		conformation::Residue const & residue,
		pose::Pose const & // will be need to create tries for disulfides
	) const;

	etrie::EtableRotamerTrieOP
	create_rotamer_trie_2(
		conformation::RotamerSetBase const & rotset
	) const;

	etrie::EtableRotamerTrieOP
	create_rotamer_trie_1(
		conformation::RotamerSetBase const & rotset,
		Size connection_type // HACK: 1 for lower connect, 2 for upper connect.  Replace this with an enum(?)
	) const;

protected:


	inline
	Real
	eval_dE_dR_over_r_(
		conformation::Atom const & atom1,
		conformation::Atom const & atom2,
		EnergyMap const & weights,
		Vector & f1,
		Vector & f2
	) const;

	void
	derived_prepare_for_residue_pair(
		Size const,
		Size const,
		pose::Pose const &
	) const
	{ //do nothing if the derived class does not override this method
	}

	// this is the only function that uses Residue -- should this be somewhere else?
	inline
	bool
	calculate_intra_xover4( core::conformation::Residue const & rsd ) const {
		if ( rsd.is_protein() /*&& options_.exclude_intra_res_protein()*/ ) return false;
		if ( rsd.is_DNA() && exclude_DNA_DNA_ /* default true */ ) return false;
		return true;
	}

	//little helper methods for interpolation:
	bool interpolate_bins(
		conformation::Atom const & atom1,
		conformation::Atom const & atom2,
		Real &d2,
		int &disbin,
		Real &frac
	) const;

	ScoreType
	rep_scoretype() const;

	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////

protected:

	Etable const & etable() const { return etable_; }

	inline bool do_classic_intrares() const { return do_classic_intrares_; }
	inline bool put_intra_into_total() const { return put_intra_into_total_; }

private:

	Etable const & etable_; // held as a const reference instead of as a pointer for fast access.

	Real safe_max_dis2;
	Real hydrogen_interaction_cutoff2_;

	// temporary hack -- make this configurable/cleaner, Phil
	bool exclude_DNA_DNA_;
	bool do_classic_intrares_;
	bool put_intra_into_total_;
	bool exclude_intra_res_protein_;
};


///////////////////////////////////////////////////////////////////////////////
// inline methods
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
// inline methods
///////////////////////////////////////////////////////////////////////////////
//template < class Derived >
//inline
//bool
//BaseEtableEnergy< Derived >::interpolate_bins(
// conformation::Atom const & atom1,
// conformation::Atom const & atom2,
// Real &d2,
// int &disbin,
// Real &frac
//) const
//{
// d2 = atom1.xyz().distance_squared( atom2.xyz() );
//
// if ( ( d2 >= safe_max_dis2 ) || ( d2 == Real(0.0) ) ) {
//  return false;
// }
//
// // bin by distance:
// Real const d2_bin = d2 * etable_bins_per_A2;
// disbin = static_cast< int >( d2_bin ) + 1;
// // int const disbin2 = disbin + 1;
// frac = d2_bin - ( disbin - 1 );
// return true;
// //ctsa
// //ctsa  tables have been hacked so that if disbin2 = lastbin, all values = 0.
// //ctsa
//}

//template <class Derived>
//inline
//void
//BaseEtableEnergy< Derived>::atom_pair_energy_(
// conformation::Atom const & atom1,
// conformation::Atom const & atom2,
// Real const weight,
// Real &atr,
// Real &rep,
// Real &solv,
// Real &bb,
// Real & d2
//) const
//{
// bb = 0;
// etable_evaluator_->atom_pair_energy( atom1, atom2, weight, atr, rep, solv, d2 );
//
//debug_assert( ljatr_.active() );
// bb = 0.0; //bead-bead interaction energy only in CoarseTable
// int disbin; Real frac;
// atr = rep = solv = bb = 0.0;
//
// //etable_.interpolated_analytic_etable_evaluation( atom1, atom2, atr, rep, solv, d2 );
// etable_.analytic_etable_evaluation( atom1, atom2, atr, rep, solv, d2 );
// atr  *= weight;
// rep  *= weight;
// solv *= weight;
// return;
//
// if (interpolate_bins(atom1,atom2,d2,disbin,frac)) {
//
//
//   //  std::cerr << "atom_pair_energy... " << disbin << ' ' << d2 << ' ' << frac << ' ' << ljatr.size() << std::endl;
//   // l1 and l2 are FArray LINEAR INDICES for fast lookup:
//   // [ l1 ] == (disbin,attype2,attype1)
//   // [ l2 ] == (disbin2,attype2,attype1)
//
//   int const l1 = ljatr_.index( disbin, atom1.type(), atom2.type()),
//     l2 = l1 + 1;
//
//   Real e1 = ljatr_[ l1 ];
//   atr = weight * ( e1 + frac * ( ljatr_[ l2 ] - e1 ) );
//
//   e1 = ljrep_[ l1 ];
//   rep = weight * ( e1 + frac * ( ljrep_[ l2 ] - e1 ) );
//
//   e1 = solv1_[ l1 ] + solv2_[ l1 ];
//   solv = weight * ( e1 + frac * ( solv1_[ l2 ] + solv2_[l2] - e1 ) );
//   //  std::cout << "solv " << solv << std::endl;
//   //  std::cerr << "finished evaluating atom _pair energy " << std::endl;
//
// } //if within cutoff
//}

/////////////////////////////////////////////////////////////////////////////////
//template <class Derived>
//inline
//Real
//BaseEtableEnergy< Derived >::eval_dE_dR_over_r_(
// conformation::Atom const & atom1,
// conformation::Atom const & atom2,
// EnergyMap const & weights,
// Vector & f1,
// Vector & f2
//) const
//{
// return etable_evaluator_->eval_dE_dR_over_r( atom1, atom2, weights, f1, f2 );
//
// Real d2,frac;
// int disbin;
//
// if ( atom1.xyz().distance_squared( atom2.xyz() ) > safe_max_dis2 ) return 0.;
//
// f1 = atom1.xyz().cross( atom2.xyz() );
// f2 = atom1.xyz() - atom2.xyz();
//
// Real datr, drep, dsol, invd;
// etable_.analytic_etable_derivatives( atom1, atom2, datr, drep, dsol, invd );
// return ( weights[ fa_atr ] * datr + weights[ fa_rep ] * drep + weights[ fa_sol ] * dsol ) * invd;
//
//
// if ( interpolate_bins(atom1,atom2,d2,disbin,frac) ) {
//
//
//   f1 = atom1.xyz().cross( atom2.xyz() );
//   f2 = atom1.xyz() - atom2.xyz();
//
//   // l1 and l2 are FArray LINEAR INDICES for fast lookup:
//   // [ l1 ] == (disbin  ,attype2,attype1)
//   // [ l2 ] == (disbin+1,attype2,attype1)
//
//   /// BEGIN DERIVATIVE INTERPOLATION
//   Real deriv = 0.0;
//
//   int const l1 = dljatr_.index( disbin, atom1.type(), atom2.type()),
//     l2 = l1 + 1;
//
//   Real e1 = dljatr_[ l1 ];
//   deriv = weights[ st_atr_] * ( e1 + frac * ( dljatr_[ l2 ] - e1 ) );
//
//   e1 = dljrep_[ l1 ];
//   deriv += weights[ st_rep_ ] * ( e1 + frac * ( dljrep_[ l2 ] - e1 ) );
//
//   e1 = dsolv_[ l1 ];
//   deriv += weights[ st_sol_ ] * ( e1 + frac * ( dsolv_[ l2 ] - e1 ) );
//
//   return deriv / std::sqrt( d2 );
//
//   /// BEGIN EXACT DERIVATIVE CALCULATION
//
//   /*Real deriv( 0.0 );
//  int const l1 = ljatr_.index( disbin, atom1.type(), atom2.type()),
//   l2 = l1 + 1;
//
//  // d g(x) / dx with g(x) = x^2 ---> 2x;  x is the distance
//  // we want to avoid the sqrt. Since at1-at2 (f2) already is the right length, consider it pre-multiplied by sqrt(d2).
//  // so what we have below is 2x / ( d2step * x ) = 2 / d2step = 2 * inv(d2step).  This will be multiplied by the un-normalized f1 and f2 vectors.
//  Real dxsquared_dx_times_x2step_over_x =  2 * etable_bins_per_A2;
//
//  Real const atr1 = ljatr_[ l1 ];
//  Real const atr2 = ljatr_[ l2 ];
//  Real const rep1 = ljrep_[ l1 ];
//  Real const rep2 = ljrep_[ l2 ];
//  Real const sol1 = solv1_[ l1 ] + solv2_[ l1 ];
//  Real const sol2 = solv1_[ l2 ] + solv2_[ l2 ];
//
//  deriv = weights[ st_atr_ ] * ( atr2 - atr1 );
//  deriv += weights[ st_rep_ ] * ( rep2 - rep1 );
//  deriv += weights[ st_sol_ ] * ( sol2 - sol1 );
//
//  return deriv * dxsquared_dx_times_x2step_over_x;*/
//
//   /// TEMP
//   /*
//  std::cout << "Testing numerically: ";
//  Real f11( 0.0 );
//  Real step = 0.00001;
//  {// scope
//   Real altd = std::sqrt( d2 ) - step;
//   Real altd2 = altd*altd;
//   Real altfrac = ( altd2 * etable_bins_per_A2 - ( disbin - 1 ) );
//   int const altl1 = ljatr_.index( disbin, atom1.type(), atom2.type()),
//    altl2 = altl1 + 1;
//
//
//   Real e1 = ljatr_[ altl1 ];
//   f11 = weights[ st_atr_ ] * ( e1 + altfrac * ( ljatr_[ altl2 ] - e1 ) );
//
//   e1 = ljrep_[ altl1 ];
//   f11 += weights[ st_rep_ ] * ( e1 + altfrac * ( ljrep_[ altl2 ] - e1 ) );
//
//   e1 = solv1_[ altl1 ] + solv2_[ altl1 ];
//   f11 += weights[ st_sol_ ] * ( e1 + altfrac * ( solv1_[ altl2 ] + solv2_[altl2] - e1 ) );
//  }
//
//  Real f22(0.0);
//  {// scope
//   Real altd = std::sqrt( d2 ) + step;
//   Real altd2 = altd*altd;
//   Real altfrac = ( altd2 * etable_bins_per_A2 - ( disbin - 1 ) );
//   int const altl1 = ljatr_.index( disbin, atom1.type(), atom2.type()),
//    altl2 = altl1 + 1;
//
//
//   Real e1 = ljatr_[ altl1 ];
//   f22 = weights[ st_atr_ ] * ( e1 + altfrac * ( ljatr_[ altl2 ] - e1 ) );
//
//   e1 = ljrep_[ altl1 ];
//   f22 += weights[ st_rep_ ] * ( e1 + altfrac * ( ljrep_[ altl2 ] - e1 ) );
//
//   e1 = solv1_[ altl1 ] + solv2_[ altl1 ];
//   f22 += weights[ st_sol_ ] * ( e1 + altfrac * ( solv1_[ altl2 ] + solv2_[altl2] - e1 ) );
//  }
//  std::cout << "Deriv discrep: " << ( f22 - f11 ) / (2 * step ) << " vs " << deriv * dxsquared_dx_times_x2step_over_x * std::sqrt( d2 ) << std::endl;
//    */
//
// } else {
//  return 0.0;
// }
//}


//template < class Derived >
//inline
//void
//BaseEtableEnergy< Derived >::pair_energy_H(
// conformation::Atom const & atom1,
// conformation::Atom const & atom2,
// Real weight,
// EnergyMap &emap
//) const {
// Energy atr(0.0);
// Energy rep(0.0);
// Energy solv(0.0);
// Energy bb(0.0);
// pair_energy_H(atom1,atom2,weight,atr,rep,solv,bb);
// emap[st_atr_]+=atr;
// emap[st_rep_]+=rep;
// emap[st_sol_]+=solv;
// emap[ coarse_beadlj ]+=bb;
//}


///////////////////////////////////////////////////////////////////////////////
//template < class Derived >
//inline
//void
//BaseEtableEnergy< Derived >::pair_energy_H_(
// conformation::Atom const & atom1,
// conformation::Atom const & atom2,
// Real const weight,
// Real &atr,
// Real &rep,
// Real &solv,
// Real &bb
//) const
//{
// bb = 0;
// Real dis2;
// return etable_evaluator_->atom_pair_energy( atom1, atom2, weight, atr, rep, solv, dis2 );
//debug_assert( ljrep_.active() );
// Real d2,frac;
// int disbin;
// atr = rep = solv = bb = 0.0;
// //etable_.interpolated_analytic_etable_evaluation( atom1, atom2, atr, rep, solv, d2 );
// etable_.analytic_etable_evaluation( atom1, atom2, atr, rep, solv, d2 );
// atr  *= weight;
// rep  *= weight;
// solv *= weight;
// return;
//
// if (interpolate_bins(atom1,atom2,d2,disbin,frac)) {
//
//
//
//   //ctsa
//   //ctsa  tables have been hacked so that if disbin2 = lastbin, all values = 0.
//   //ctsa
//
//   // l is an FArray LINEAR INDICES for fast lookup:
//   // [ ll ] == (disbin,attype2,attype1)
//
//   int l1 = ljrep_.index( disbin, atom1.type(), atom2.type() );
//   int l2 = l1+1;
//
//   Real const rep_e1( ljrep_[ l1 ] );
//   rep =  weight * ( rep_e1 + frac * ( ljrep_[ l2 ] - rep_e1 ) );
//
//   Real const atr_e1 = ljatr_[ l1 ];
//   atr = weight * ( atr_e1 + frac * ( ljatr_[ l2 ] - atr_e1 ) );
//   //  std::cerr << __FILE__<< ' ' << __LINE__ << std::endl;
//
// }
//}


} // etable
} // scoring
} // core


#endif // INCLUDED_core_scoring_EtableEnergy_HH
