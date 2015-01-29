// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file		core/scoring/etable/BaseMembEtableEnergy.hh
///
/// @brief		Base class for membrane fullatom energy method classs
/// @details	Warning: This class is not threadsafe because it has mutable elements
///
/// @author		(Original) Phil Bradley, Andrew Leaver-Fay, Oliver Lange
/// @author		Patrick Barth
/// @author		(updates) Rebecca Alford (rfalford12@gmail.com)


#ifndef INCLUDED_core_scoring_memb_etable_BaseMembEtableEnergy_hh
#define INCLUDED_core_scoring_memb_etable_BaseMembEtableEnergy_hh

// Unit headers
#include <core/scoring/etable/EtableEnergy.fwd.hh>

// Package headers
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/etable/etrie/EtableAtom.hh>
#include <core/scoring/etable/count_pair/types.hh>
#include <core/scoring/trie/TrieCountPairBase.fwd.hh>
#include <core/scoring/etable/etrie/EtableTrie.fwd.hh>

#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/ContextGraphTypes.hh>

// Project headers
#include <core/conformation/Atom.hh>
#include <core/pose/Pose.fwd.hh>

#include <ObjexxFCL/FArray3D.hh>

namespace core {
namespace scoring {
namespace etable {

template < class Derived >
class BaseEtableEnergy : public methods::ContextIndependentTwoBodyEnergy
{

public:
	typedef methods::EnergyMethodOP EnergyMethodOP;

public:

	/// @brief Construct given an Etable
	BaseEtableEnergy(
		Etable const & etable_in,
		methods::EnergyMethodOptions const & options,
		ScoreType st_atr,
		ScoreType st_rep,
		ScoreType st_sol
	);


	/////////////////////////////////////////////////////////////////////////////
	// methods for ContextIndependentTwoBodyEnergies
	/////////////////////////////////////////////////////////////////////////////

	/// @brief Setup TwoBody Energy for Minimization
	/// @details Stash nblist if use_nblist is true
	virtual
	void
	setup_for_minimizing(
		pose::Pose & pose,
		ScoreFunction const & sfxn,
		kinematics::MinimizerMapBase const & min_map
	) const;

	/// @brief Setup for Scoring Two Body Term
	/// @details Check compatibility with atom typeset
 	virtual
 	void
 	setup_for_scoring( pose::Pose &pose, ScoreFunction const &scfxn ) const;

	/// @brief Setup two-body terms for derivatives
 	virtual
 	void
	setup_for_derivatives(
		pose::Pose &pose,
		ScoreFunction const &scfxn
	) const;

	/// @brief Setup Two body Term for Packing
	/// @details Etable Energy method scores a vector of rotamer trie objects in the
	/// energies object for use in rapid rotamer/background energy calculations. Overrides
	/// default do-nothing behavior.
	virtual
	void
	setup_for_packing( pose::Pose &, utility::vector1< bool > const &, utility::vector1< bool > const & ) const;

	/// @brief Create Rotamer Trie for the input set of rotamers and scores the trie in the rotamer set.
	virtual
	void
	prepare_rotamers_for_packing(
		pose::Pose const & pose,
		conformation::RotamerSetBase & set
	) const;

	/// @brief Updates the cached rotamer trie for a residue if it has changed during the course of a repacking
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

	/// @brief Compute Two-Body residue Pair energy
	virtual
	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
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


	/// @brief Overrides default rotamer/background energy calculation and uses
	/// the trie-vs-trie algorithm instead
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

	/// @brief Called at the end of an energy evaluation
	virtual
	void
	finalize_total_energy(
		pose::Pose & pose,
		ScoreFunction const &,
		EnergyMap & totals
	) const;

	/// @brief Called during graident-based minimization inside dfunc
	/// @details F1 and F2 are not zeroed -- contributions from this
	/// atom are just summed in.
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


	/// @brief Standard Atom Pair Energy
	inline
	void
	atom_pair_energy(
		conformation::Atom const & atom1,
		conformation::Atom const & atom2,
		Real const weight,
		EnergyMap & emap,
		Real & dsq
	) const {
		Energy atr, rep, solv, bb;
		atr = rep = solv = bb = 0;
		atom_pair_energy(atom1,atom2,weight,atr,rep,solv,bb,dsq);
		emap[st_atr_]+=atr;
		emap[st_rep_]+=rep;
		emap[st_sol_]+=solv;
		emap[coarse_beadlj]+=bb;
	}

	/// @brief Membrane Specific Atom Pair Energy
	inline
	void
	memb_atom_pair_energy(
		conformation::Atom const & atom1,
		conformation::Atom const & atom2,
		Real const weight,
		EnergyMap & emap,
		Real & dsq
	) const {
		Energy atr, rep, solv, bb;
		atr = rep = solv = bb = 0;
		memb_atom_pair_energy(atom1,atom2,weight,atr,rep,solv,bb,dsq);
		emap[st_atr_]+=atr;
		emap[st_rep_]+=rep;
		emap[st_sol_]+=solv;
		emap[coarse_beadlj]+=bb;
	}

	/// @brief Speed Optimized Membrane Environment Energy
	inline
	void
	fast_memb_env_energy(
		conformation::Atom const & atom1,
		int const attype1,
		Real & mbenvE
	) const {
		fast_memb_env_energy(atom1,attype1,mbenvE);
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
		static_cast< Derived const* > (this) -> atom_pair_energy_(atom1,atom2,weight,atr,rep,solv,bb,dsq);
	}

	/// @brief for the trie-vs-trie algorithm; could test if the other
	/// atom pair energy function could inline this function to avoid
	/// the table reading code duplication.
	inline
	void
	memb_atom_pair_energy(
		conformation::Atom const & atom1,
		conformation::Atom const & atom2,
		Real const weight,
		Energy & atr,
		Energy & rep,
		Energy & solv,
		Energy & bb,
		Real & dsq
	) const  {
		static_cast< Derived const* > (this) -> memb_atom_pair_energy_(atom1,atom2,weight,atr,rep,solv,bb,dsq);
	}

	/// @brief Compute fast membrane environment energy
	inline
	void
	fast_memb_env_energy(
		conformation::Atom const & atom1,
		int const attype1,
		Real & mbenvE
	) const {
		static_cast< Derived const* > (this) -> fast_memb_env_energy_(atom1,attype1,mbenvE);
	}

	///
	inline
	void
	pair_energy_H(
		conformation::Atom const & atom1,
		conformation::Atom const & atom2,
		Real weight,
		Energy & atr,
		Energy & rep,
		Energy & solv,
		Energy & bb
	) const  {
		return static_cast< Derived const* > (this) -> pair_energy_H_(atom1,atom2,weight,atr,rep,solv,bb);
	};

	///
	inline
	void
	pair_energy_H(
		conformation::Atom const & atom1,
		conformation::Atom const & atom2,
		Real weight,
		EnergyMap &emap
	) const;



	///
	inline
	Real
	eval_dE_dR_over_r(
		conformation::Atom const & atom1,
		conformation::Atom const & atom2,
		EnergyMap const & weights,
		Vector & f1,
		Vector & f2
	) const  {
		return static_cast< Derived const* > (this) -> eval_dE_dR_over_r_(atom1,atom2,weights,f1,f2);
	};

	///pba
	inline
	Real
	memb_eval_dE_dR_over_r(
		conformation::Atom const & atom1,
		conformation::Atom const & atom2,
		EnergyMap const & weights,
		Vector & f1,
		Vector & f2
	) const  {
		return static_cast< Derived const* > (this) -> memb_eval_dE_dR_over_r_(atom1,atom2,weights,f1,f2);
	};

	///
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
	void indicate_required_context_graphs( utility::vector1< bool > & /*context_graphs_required*/ ) const;

	/// Inline Methods For Trie-vs-Trie Algorithm
	inline
	Energy sum_energies ( Real atr, Real rep, Real solv, Real bb ) const {
		return weights_[ st_atr_ ] * atr
			+ weights_[ st_rep_ ] * rep
			+ weights_[ st_sol_ ] * solv
			+ weights_ [ coarse_beadlj ] * bb;
	}

	inline
	Energy heavyatom_heavyatom_energy(
		etrie::EtableAtom const & at1,
		etrie::EtableAtom const & at2,
		DistanceSquared & d2
	) const
	{
		Energy atr(0.0), rep(0.0), solv(0.0), bb(0.0);
		atom_pair_energy( at1, at2, 1.0, atr, rep, solv, bb, d2 );
		return sum_energies( atr, rep, solv, bb );
	}

	inline
	Energy heavyatom_hydrogenatom_energy(
		etrie::EtableAtom const & at1,
		etrie::EtableAtom const & at2
	) const
	{
		Energy atr(0.0), rep(0.0), solv(0.0), bb(0.0);
		pair_energy_H( at1, at2, 1.0, atr, rep, solv, bb );
		return sum_energies( atr, rep, solv, bb );
	}

	inline
	Energy hydrogenatom_heavyatom_energy(
		etrie::EtableAtom const & at1,
		etrie::EtableAtom const & at2
	) const
	{
		Energy atr(0.0), rep(0.0), solv(0.0), bb(0.0);
		pair_energy_H( at1, at2, 1.0, atr, rep, solv, bb );
		return sum_energies( atr, rep, solv, bb );
	}

	inline
	Energy hydrogenatom_hydrogenatom_energy(
		etrie::EtableAtom const & at1,
		etrie::EtableAtom const & at2
	) const
	{
		Energy atr(0.0), rep(0.0), solv(0.0), bb(0.0);
		pair_energy_H( at1, at2, 1.0, atr, rep, solv, bb );
		return sum_energies( atr, rep, solv, bb );
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

	/*count_pair::CPResidueConnectionType
	determine_residue_connection(
		conformation::Residue const & res1,
		conformation::Residue const & res2,
		pose::Pose const &
	) const;*/

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
	create_rotamer_trie_2(
		conformation::RotamerSetBase const & rotset
	) const;

	etrie::EtableRotamerTrieOP
	create_rotamer_trie_1(
		conformation::RotamerSetBase const & rotset,
		Size connection_type // HACK: 1 for lower connect, 2 for upper connect.  Replace this with an enum(?)
	) const;


protected:
	// implementation for quasi-virtual functions
	// rename this derived_atom_pair_energy
	inline
	void
	atom_pair_energy_(
		conformation::Atom const & atom1,
		conformation::Atom const & atom2,
		Real const weight,
		Real &atr,
		Real &rep,
		Real &solv,
		Real &bb,
		Real & dsq
	) const;

  //pba implementation for quasi-virtual functions
  // rename this derived_atom_pair_energy
  inline
  void
  memb_atom_pair_energy_(
    conformation::Atom const & atom1,
    conformation::Atom const & atom2,
    Real const weight,
    Real &atr,
    Real &rep,
    Real &solv,
    Real &bb,
    Real & dsq
  ) const;

  //pba
  inline
  void
  fast_memb_env_energy(
  conformation::Atom const & atom1,
  int const attype1,
  Real & mbenvE
  ) const;

	///
	inline
	void
	pair_energy_H_(
		conformation::Atom const & atom1,
		conformation::Atom const & atom2,
		Real const weight,
		Real &atr,
		Real &rep,
		Real &solv,
		Real &bb
	) const;

	///
	inline
	Real
	eval_dE_dR_over_r_(
		conformation::Atom const & atom1,
		conformation::Atom const & atom2,
		EnergyMap const & weights,
		Vector & f1,
		Vector & f2
	) const;

  ///pba
  inline
  Real
  memb_eval_dE_dR_over_r_(
    conformation::Atom const & atom1,
    conformation::Atom const & atom2,
    EnergyMap const & weights,
    Vector & f1,
    Vector & f2
  ) const;

	void
	derived_prepare_for_residue_pair(
		Size const res1,
		Size const res2,
		pose::Pose const &
	) const
	{ //do nothing if the derived class does not override this method
	}


	//little helper methods for interpolation:
	bool interpolate_bins(
		conformation::Atom const & atom1,
		conformation::Atom const & atom2,
		Real &d2,
		int &disbin,
		Real &frac
	) const;

	void
	set_scoretypes(
		ScoreType atr_type,
		ScoreType rep_type,
		ScoreType sol_type
	) const;

	ScoreType
	rep_scoretype() const;

	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////
protected:
	Etable const & etable_; // shouldn't this be a pointer? Reference count information is (dangerously) lost when
	//a reference is taken, instead of a smart pointer.  There's the potential for a dangling reference with this.

private:
	/// these guys are taken from the etable
	ObjexxFCL::FArray3D< Real > const & ljatr_;
	ObjexxFCL::FArray3D< Real > const & ljrep_;
	ObjexxFCL::FArray3D< Real > const & solv1_;
	ObjexxFCL::FArray3D< Real > const & solv2_;
	ObjexxFCL::FArray3D< Real > const & dljatr_;
	ObjexxFCL::FArray3D< Real > const & dljrep_;
	ObjexxFCL::FArray3D< Real > const & dsolv_;
  ObjexxFCL::FArray3D< Real > const & memb_solv1_; //pba
  ObjexxFCL::FArray3D< Real > const & memb_solv2_; //pba
  ObjexxFCL::FArray3D< Real > const & memb_dsolv_; //pba
	Real safe_max_dis2;

	int etable_bins_per_A2;

	Real hydrogen_interaction_cutoff2_;

	mutable ScoreType st_rep_;
	mutable ScoreType st_atr_;
	mutable ScoreType st_sol_;

	// For use in the trie-vs-trie algorithm
	// written to at the beginning of the evaluate_rotamer_pair_energies call (which is const)
	mutable EnergyMap weights_;

	// temporary hack -- make this configurable/cleaner, Phil
	bool exclude_DNA_DNA;

	// not clear that we need these:
	//chemical::AtomTypeSet const * atom_set_ptr_; // for atomtype info
};


///////////////////////////////////////////////////////////////////////////////
// inline methods
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
// inline methods
///////////////////////////////////////////////////////////////////////////////
template < class Derived >
inline
bool
BaseEtableEnergy< Derived >::interpolate_bins(
	conformation::Atom const & atom1,
	conformation::Atom const & atom2,
	Real &d2,
	int &disbin,
	Real &frac
) const
{
	d2 = atom1.xyz().distance_squared( atom2.xyz() );

	if ( ( d2 >= safe_max_dis2 ) || ( d2 == Real(0.0) ) ) {
		return false;
	}

	// bin by distance:
	Real const d2_bin = d2 * etable_bins_per_A2;
	disbin = static_cast< int >( d2_bin ) + 1;
//	int const disbin2 = disbin + 1;
	frac = d2_bin - ( disbin - 1 );
	return true;
//ctsa
//ctsa  tables have been hacked so that if disbin2 = lastbin, all values = 0.
//ctsa
}

template <class Derived>
inline
void
BaseEtableEnergy< Derived>::atom_pair_energy_(
  conformation::Atom const & atom1,
  conformation::Atom const & atom2,
  Real const weight,
  Real &atr,
  Real &rep,
  Real &solv,
  Real &bb,
  Real & d2
) const
{
 debug_assert( ljatr_.active() );
  bb = 0.0; //bead-bead interaction energy only in CoarseTable
  int disbin; Real frac;
  atr = rep = solv = bb = 0.0;
  if (interpolate_bins(atom1,atom2,d2,disbin,frac)) {
    //    std::cerr << "atom_pair_energy... " << disbin << ' ' << d2 << ' ' << frac << ' ' << ljatr.size() << std::endl;
    // l1 and l2 are FArray LINEAR INDICES for fast lookup:
    // [ l1 ] == (disbin,attype2,attype1)
    // [ l2 ] == (disbin2,attype2,attype1)

    int const l1 = ljatr_.index( disbin, atom1.type(), atom2.type()),
      l2 = l1 + 1;


    Real e1 = ljatr_[ l1 ];
    atr = weight * ( e1 + frac * ( ljatr_[ l2 ] - e1 ) );

    e1 = ljrep_[ l1 ];
    rep = weight * ( e1 + frac * ( ljrep_[ l2 ] - e1 ) );

    e1 = solv1_[ l1 ] + solv2_[ l1 ];
    solv = weight * ( e1 + frac * ( solv1_[ l2 ] + solv2_[l2] - e1 ) );
    //    std::cout << "solv " << solv << std::endl;
    //    std::cerr << "finished evaluating atom _pair energy " << std::endl;
  } //if within cutoff
}

///////////////////////////////////////////////////////////////////////////////

//pba
template <class Derived>
inline
void
BaseEtableEnergy< Derived>::fast_memb_env_energy_(
  conformation::Atom const & atom1,
  int const attype1,
  Real & mbenvE
) const
{

    Vector const normal(MembraneEmbed_from_pose( pose ).normal());
    Vector const center(MembraneEmbed_from_pose( pose ).center());
    Real const center_FA(MembraneEmbed_from_pose( pose ).center_FA);
    Real const thickness(MembraneEmbed_from_pose( pose ).thickness);
    Real const steepness(MembraneEmbed_from_pose( pose ).steepness);

    Real z = dot(atom1.xyz(),normal)+30;
    z = fabs(z-center_FA);
    z /= thickness;
    Real zn = std::pow( z, steepness );
    Real f = zn/(1 + zn);
    mbenvE = (1 - f) * (memb_lk_dgrefce(attype1) - lk_dgrefce(attype1));
}
///////////////////////////////////////////////////////////////////////////////

//pba
template <class Derived>
inline
void
BaseEtableEnergy< Derived>::memb_atom_pair_energy_(
	conformation::Atom const & atom1,
	conformation::Atom const & atom2,
	Real const weight,
	Real &atr,
	Real &rep,
	Real &solv,
	Real &bb,
	Real & d2
) const
{
debug_assert( ljatr_.active() );
	bb = 0.0; //bead-bead interaction energy only in CoarseTable
	int disbin; Real frac;
	atr = rep = solv = bb = 0.0;
	if (interpolate_bins(atom1,atom2,d2,disbin,frac)) {
		//		std::cerr << "atom_pair_energy... " << disbin << ' ' << d2 << ' ' << frac << ' ' << ljatr.size() << std::endl;
		// l1 and l2 are FArray LINEAR INDICES for fast lookup:
		// [ l1 ] == (disbin,attype2,attype1)
		// [ l2 ] == (disbin2,attype2,attype1)

		int const l1 = ljatr_.index( disbin, atom1.type(), atom2.type()),
			l2 = l1 + 1;


		Real e1 = ljatr_[ l1 ];
		atr = weight * ( e1 + frac * ( ljatr_[ l2 ] - e1 ) );

		e1 = ljrep_[ l1 ];
		rep = weight * ( e1 + frac * ( ljrep_[ l2 ] - e1 ) );

    //pba Membrane specific solvation

    Vector const normal(MembraneEmbed_from_pose( pose ).normal());
    Vector const center(MembraneEmbed_from_pose( pose ).center());
    Real const center_FA(MembraneEmbed_from_pose( pose ).center_FA);
    Real const thickness(MembraneEmbed_from_pose( pose ).thickness);
    Real const steepness(MembraneEmbed_from_pose( pose ).steepness);

    //pba solvation of atom1 based on its distance from the membrane center on the membrane normal

    Real z = dot(atom1.xyz(),normal)+30;
    z = fabs(z-center_FA);
    z /= thickness;
    Real zn = std::pow( z, steepness );
    Real f = zn/(1 + zn);

    e11 = f * solv1_[ l1 ] + (1 - f) * memb_solv1_[ l1 ];
    e12 = f * solv1_[ l2 ] + (1 - f) * memb_solv1_[ l2 ];

    //pba solvation of atom2 based on its distance from the membrane center on the membrane normal

    z = dot(atom2.xyz(),normal)+30;
    z = fabs(z-center_FA);
    z /= thickness;
    zn = std::pow( z, steepness );
    f = zn/(1 + zn);

    e21 = f * solv2_[ l1 ] + (1 - f) * memb_solv2_[ l1 ];
    e22 = f * solv2_[ l2 ] + (1 - f) * memb_solv2_[ l2 ];

    e1 = e11 + e21;
    Real e2 = e12 + e22;

    solv = weight * ( e1 + frac * ( e2 - e1 ) );

		//		std::cout << "solv " << solv << std::endl;
		//		std::cerr << "finished evaluating atom _pair energy " << std::endl;
	} //if within cutoff
}

///////////////////////////////////////////////////////////////////////////////
template <class Derived>
inline
Real
BaseEtableEnergy< Derived >::memb_eval_dE_dR_over_r_(
  conformation::Atom const & atom1,
  conformation::Atom const & atom2,
  EnergyMap const & weights,
  Vector & f1,
  Vector & f2
) const
{
 debug_assert( dljatr_.active() );

  f1 = atom1.xyz().cross( atom2.xyz() );
  f2 = atom1.xyz() - atom2.xyz();
  Real d2,frac;
  int disbin;

  if ( interpolate_bins(atom1,atom2,d2,disbin,frac) ) {
    // l1 and l2 are FArray LINEAR INDICES for fast lookup:
    // [ l1 ] == (disbin  ,attype2,attype1)
    // [ l2 ] == (disbin+1,attype2,attype1)

    Real deriv = 0.0;

    int const l1 = dljatr_.index( disbin, atom1.type(), atom2.type()),
      l2 = l1 + 1;

    Real e1 = dljatr_[ l1 ];
    deriv = weights[ st_atr_] * ( e1 + frac * ( dljatr_[ l2 ] - e1 ) );

    e1 = dljrep_[ l1 ];
    deriv += weights[ st_rep_ ] * ( e1 + frac * ( dljrep_[ l2 ] - e1 ) );

    //pba Membrane specific solvation

    Vector const normal(MembraneEmbed_from_pose( pose ).normal());
    Vector const center(MembraneEmbed_from_pose( pose ).center());
    Real const center_FA(MembraneEmbed_from_pose( pose ).center_FA);
    Real const thickness(MembraneEmbed_from_pose( pose ).thickness);
    Real const steepness(MembraneEmbed_from_pose( pose ).steepness);

    //pba solvation of atom1 based on its distance from the membrane center on the membrane normal

    Real z = dot(atom1.xyz(),normal)+30;
    z = fabs(z-center_FA);
    z /= thickness;
    Real zn = std::pow( z, steepness );
    Real f = zn/(1 + zn);

    e1 = f * dsolv1_[ l1 ] + (1 - f) * memb_dsolv1_[ l1 ];

    //pba solvation of atom2 based on its distance from the membrane center on the membrane normal

    z = dot(atom2.xyz(),normal)+30;
    z = fabs(z-center_FA);
    z /= thickness;
    zn = std::pow( z, steepness );
    f = zn/(1 + zn);

    Real e2 = f * dsolv2_[ l2 ] + (1 - f) * memb_dsolv2_[ l2 ];

    deriv += weights[ st_sol_ ] * ( e1 + frac * ( e2 - e1 ) );

    return deriv / std::sqrt( d2 );
  } else {
    return 0.0;
  }
}


///////////////////////////////////////////////////////////////////////////////
template <class Derived>
inline
Real
BaseEtableEnergy< Derived >::eval_dE_dR_over_r_(
	conformation::Atom const & atom1,
	conformation::Atom const & atom2,
	EnergyMap const & weights,
	Vector & f1,
	Vector & f2
) const
{
debug_assert( dljatr_.active() );

	f1 = atom1.xyz().cross( atom2.xyz() );
	f2 = atom1.xyz() - atom2.xyz();
	Real d2,frac;
	int disbin;

	if ( interpolate_bins(atom1,atom2,d2,disbin,frac) ) {
		// l1 and l2 are FArray LINEAR INDICES for fast lookup:
		// [ l1 ] == (disbin  ,attype2,attype1)
		// [ l2 ] == (disbin+1,attype2,attype1)

		Real deriv = 0.0;

		int const l1 = dljatr_.index( disbin, atom1.type(), atom2.type()),
			l2 = l1 + 1;

		Real e1 = dljatr_[ l1 ];
		deriv = weights[ st_atr_] * ( e1 + frac * ( dljatr_[ l2 ] - e1 ) );

		e1 = dljrep_[ l1 ];
		deriv += weights[ st_rep_ ] * ( e1 + frac * ( dljrep_[ l2 ] - e1 ) );

		e1 = dsolv_[ l1 ];
		deriv += weights[ st_sol_ ] * ( e1 + frac * ( dsolv_[ l2 ] - e1 ) );


		return deriv / std::sqrt( d2 );
	} else {
		return 0.0;
	}
}


template < class Derived >
inline
void
BaseEtableEnergy< Derived >::pair_energy_H(
	conformation::Atom const & atom1,
	conformation::Atom const & atom2,
	Real weight,
	EnergyMap &emap
) const {
	Energy atr(0.0);
	Energy rep(0.0);
	Energy solv(0.0);
	Energy bb(0.0);
	pair_energy_H(atom1,atom2,weight,atr,rep,solv,bb);
	emap[st_atr_]+=atr;
	emap[st_rep_]+=rep;
	emap[st_sol_]+=solv;
	emap[ coarse_beadlj ]+=bb;
}



///////////////////////////////////////////////////////////////////////////////
template < class Derived >
inline
void
BaseEtableEnergy< Derived >::pair_energy_H_(
	conformation::Atom const & atom1,
	conformation::Atom const & atom2,
	Real const weight,
	Real &atr,
	Real &rep,
	Real &solv,
	Real &bb
) const
{
debug_assert( ljrep_.active() );
	Real d2,frac;
	int disbin;
	atr = rep = solv = bb = 0.0;
	if (interpolate_bins(atom1,atom2,d2,disbin,frac)) {
//ctsa
//ctsa  tables have been hacked so that if disbin2 = lastbin, all values = 0.
//ctsa

	// l is an FArray LINEAR INDICES for fast lookup:
	// [ ll ] == (disbin,attype2,attype1)

		int l1 = ljrep_.index( disbin, atom1.type(), atom2.type() );
		int l2 = l1+1;

		Real const rep_e1( ljrep_[ l1 ] );
		rep =  weight * ( rep_e1 + frac * ( ljrep_[ l2 ] - rep_e1 ) );

		Real const atr_e1 = ljatr_[ l1 ];
		atr = weight * ( atr_e1 + frac * ( ljatr_[ l2 ] - atr_e1 ) );
		//		std::cerr << __FILE__<< ' ' << __LINE__ << std::endl;
	}
}

} // etable
} // scoring
} // core


#endif // INCLUDED_core_scoring_EtableEnergy_HH
