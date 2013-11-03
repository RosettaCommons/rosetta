// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/ScoreFunction.hh
/// @brief  Score function class
/// @author Phil Bradley


#ifndef INCLUDED_core_scoring_etable_EtableEnergy_hh
#define INCLUDED_core_scoring_etable_EtableEnergy_hh

// Unit headers
#include <core/scoring/etable/EtableEnergy.fwd.hh>

#include <core/scoring/methods/EnergyMethod.fwd.hh>


// Package headers
#include <core/scoring/etable/BaseEtableEnergy.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/scoring/EnergyMap.hh>

#include <utility/vector1.hh>



namespace core {
namespace scoring {
namespace etable {


//////////////////// Evaluators ///////////////////////////////


class EtableEvaluator : public utility::pointer::ReferenceCount {

public:
	EtableEvaluator( Etable const & etable );
	~EtableEvaluator();

	void
	set_weights(
		EnergyMap const & weights
	) {
		atr_weight_ = weights[ st_atr_ ];
		rep_weight_ = weights[ st_rep_ ];
		sol_weight_ = weights[ st_sol_ ];
	}

	void
	set_scoretypes(
		ScoreType st_atr_in,
		ScoreType st_rep_in,
		ScoreType st_sol_in
	)
	{
		st_atr_ = st_atr_in;
		st_rep_ = st_rep_in;
		st_sol_ = st_sol_in;
	}


	inline Real atr_weight() const { return atr_weight_; }
	inline Real rep_weight() const { return rep_weight_; }
	inline Real sol_weight() const { return sol_weight_; }

	inline ScoreType st_atr() const { return st_atr_; }
	inline ScoreType st_rep() const { return st_rep_; }
	inline ScoreType st_sol() const { return st_sol_; }

	inline
	Energy
	sum_energies( Real atr, Real rep, Real solv ) const {
		return
			atr_weight_ * atr +
			rep_weight_ * rep +
			sol_weight_ * solv;
	}

	inline
	Real
	hydrogen_interaction_cutoff2() const
	{
		return hydrogen_interaction_cutoff2_;
	}

	virtual
	void
	atom_pair_energy_v(
    conformation::Atom const & atom1,
    conformation::Atom const & atom2,
    Real const weight,
		Real & atrE,
		Real & repE,
		Real & solE,
		Real & d2
	) const = 0;

	virtual
	void
	atom_pair_lk_energy_and_deriv_v(
    conformation::Atom const & atom1,
    conformation::Atom const & atom2,
		Real & solE1,
		Real & dsolE1,
		bool const eval_deriv = false
		) const = 0;

    virtual
    void
    atom_pair_lk_energy_and_deriv_v_efficient(
        conformation::Atom const & atom1,
        conformation::Atom const & atom2,
        Real & solE1,
        Real & solE2,
        Real & dsolE1,
        bool const eval_deriv
    ) const;

private:

	Real atr_weight_;
	Real rep_weight_;
	Real sol_weight_;

	/// score types: could be either fa_atr/fa_atr_intra, etc.
	ScoreType st_atr_;
	ScoreType st_rep_;
	ScoreType st_sol_;

	Real hydrogen_interaction_cutoff2_;
};



class AnalyticEtableEvaluator : public EtableEvaluator
{
public:
	AnalyticEtableEvaluator( Etable const & etable );
	~AnalyticEtableEvaluator();

	/// Atom pair energy inline type resolution functions
	void
	residue_atom_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		count_pair::CountPairFunction const & cp,
		EnergyMap & emap
	) const;

	void
	residue_atom_pair_energy_sidechain_backbone(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		count_pair::CountPairFunction const & cp,
		EnergyMap & emap
	) const;

	void
	residue_atom_pair_energy_sidechain_whole(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		count_pair::CountPairFunction const & cp,
		EnergyMap & emap
	) const;

	void
	residue_atom_pair_energy_backbone_backbone(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		count_pair::CountPairFunction const & cp,
		EnergyMap & emap
	) const;

	void
	residue_atom_pair_energy_sidechain_sidechain(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		count_pair::CountPairFunction const & cp,
		EnergyMap & emap
	) const;

	/// Trie vs trie / trie vs path type resolution functions

	void
	trie_vs_trie(
		trie::RotamerTrieBase const & trie1,
		trie::RotamerTrieBase const & trie2,
		trie::TrieCountPairBase & cp,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const;

	void
	trie_vs_path(
		trie::RotamerTrieBase const & trie1,
		trie::RotamerTrieBase const & trie2,
		trie::TrieCountPairBase & cp,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const;

	inline
	void
	atom_pair_energy(
    conformation::Atom const & atom1,
    conformation::Atom const & atom2,
    Real const weight,
		EnergyMap & emap,
		Real & d2
	) const
	{
		Real atr(0),rep(0),solv(0);
		atom_pair_energy( atom1, atom2, weight, atr, rep, solv, d2 );
		emap[st_atr()]+=atr;
    emap[st_rep()]+=rep;
    emap[st_sol()]+=solv;
	}

	virtual
	void
	atom_pair_energy_v(
    conformation::Atom const & atom1,
    conformation::Atom const & atom2,
    Real const weight,
		Real & atrE,
		Real & repE,
		Real & solE,
		Real & d2
	) const {
		atom_pair_energy( atom1, atom2, weight, atrE, repE, solE, d2 );
	}

	virtual
	void
	atom_pair_lk_energy_and_deriv_v(
    conformation::Atom const & atom1,
    conformation::Atom const & atom2,
		Real & solE1,
		Real & dsolE1,
		bool const eval_deriv = false ) const;

    virtual
    void
    atom_pair_lk_energy_and_deriv_v_efficient(
        conformation::Atom const & atom1,
        conformation::Atom const & atom2,
        Real & solE1,
        Real & solE2,
        Real & dsolE1,
        bool const eval_deriv
    ) const;


	inline
	void
	atom_pair_energy(
		conformation::Atom const & atom1,
		conformation::Atom const & atom2,
		Real const weight,
		Real &atr,
		Real &rep,
		Real &solv,
		Real & d2
	) const;

  inline
  void
  pair_energy_H(
    conformation::Atom const & atom1,
    conformation::Atom const & atom2,
    Real const weight,
    EnergyMap & emap
  ) const
  {
    Real atr(0), rep(0), sol(0), d2(0);
    atom_pair_energy( atom1, atom2, weight, atr, rep, sol, d2 );
    emap[ st_atr() ] += atr;
    emap[ st_rep() ] += rep;
    emap[ st_sol() ] += sol;
  }

	inline
	Real
	eval_dE_dR_over_r(
		conformation::Atom const & atom1,
		conformation::Atom const & atom2,
		EnergyMap const & weights,
		Vector & f1,
		Vector & f2
	) const;

	/// Inline Methods For Trie-vs-Trie Algorithm
	inline
	Energy heavyatom_heavyatom_energy(
		etrie::EtableAtom const & at1,
		etrie::EtableAtom const & at2,
		DistanceSquared & d2,
		Size & /*path_dist*/
	) const
	{
		Energy atr(0.0), rep(0.0), solv(0.0);
		atom_pair_energy( at1, at2, 1.0, atr, rep, solv, d2 );
		return sum_energies( atr, rep, solv );
	}

	inline
	Energy heavyatom_hydrogenatom_energy(
		etrie::EtableAtom const & at1,
		etrie::EtableAtom const & at2,
		Size & /*path_dist*/
	) const
	{
		Energy atr(0.0), rep(0.0), solv(0.0), d2dummy(0.0);
		atom_pair_energy( at1, at2, 1.0, atr, rep, solv, d2dummy );
		return sum_energies( atr, rep, solv );
	}

	inline
	Energy hydrogenatom_heavyatom_energy(
		etrie::EtableAtom const & at1,
		etrie::EtableAtom const & at2,
		Size & /*path_dist*/
	) const
	{
		Energy atr(0.0), rep(0.0), solv(0.0), d2dummy(0.0);
		atom_pair_energy( at1, at2, 1.0, atr, rep, solv, d2dummy );
		return sum_energies( atr, rep, solv );
	}

	inline
	Energy hydrogenatom_hydrogenatom_energy(
		etrie::EtableAtom const & at1,
		etrie::EtableAtom const & at2,
		Size & /*path_dist*/
	) const
	{
		Energy atr(0.0), rep(0.0), solv(0.0), d2dummy(0.0);
		atom_pair_energy( at1, at2, 1.0, atr, rep, solv, d2dummy );
		return sum_energies( atr, rep, solv );
	}


private:

	Etable const & etable_;
	Real safe_max_dis2_;

};


class TableLookupEvaluator : public EtableEvaluator
{
public:
	TableLookupEvaluator( Etable const & etable_in );
	~TableLookupEvaluator();

	/// Atom pair energy inline type resolution functions
	void
	residue_atom_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		count_pair::CountPairFunction const & cp,
		EnergyMap & emap
	) const;

	void
	residue_atom_pair_energy_sidechain_backbone(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		count_pair::CountPairFunction const & cp,
		EnergyMap & emap
	) const;

	void
	residue_atom_pair_energy_sidechain_whole(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		count_pair::CountPairFunction const & cp,
		EnergyMap & emap
	) const;

	void
	residue_atom_pair_energy_backbone_backbone(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		count_pair::CountPairFunction const & cp,
		EnergyMap & emap
	) const;

	void
	residue_atom_pair_energy_sidechain_sidechain(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		count_pair::CountPairFunction const & cp,
		EnergyMap & emap
	) const;

	/// Trie vs trie / trie vs path type resolution functions
	void
	trie_vs_trie(
		trie::RotamerTrieBase const & trie1,
		trie::RotamerTrieBase const & trie2,
		trie::TrieCountPairBase & cp,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const;

	void
	trie_vs_path(
		trie::RotamerTrieBase const & trie1,
		trie::RotamerTrieBase const & trie2,
		trie::TrieCountPairBase & cp,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const;

	virtual
	void
	atom_pair_energy_v(
    conformation::Atom const & atom1,
    conformation::Atom const & atom2,
    Real const weight,
		EnergyMap & emap,
		Real & d2
	) const
	{
		atom_pair_energy( atom1, atom2, weight, emap, d2 );
	}

	virtual
	void
	atom_pair_lk_energy_and_deriv_v(
    conformation::Atom const & atom1,
    conformation::Atom const & atom2,
		Real & solE1,
		Real & dsolE1,
		bool const eval_deriv = false ) const;

    virtual
    void
    atom_pair_lk_energy_and_deriv_v_efficient(
                                              conformation::Atom const & atom1,
                                              conformation::Atom const & atom2,
                                              Real & solE1,
                                              Real & solE2,
                                              Real & dsolE1,
                                              bool const eval_deriv
                                              ) const;

	inline
	void
	atom_pair_energy(
    conformation::Atom const & atom1,
    conformation::Atom const & atom2,
    Real const weight,
		EnergyMap & emap,
		Real & d2
	) const
	{
		Real atr(0),rep(0),solv(0);
		atom_pair_energy( atom1, atom2, weight, atr, rep, solv, d2 );
		emap[st_atr()]+=atr;
		emap[st_rep()]+=rep;
		emap[st_sol()]+=solv;
	}

	virtual
	void
	atom_pair_energy_v(
    conformation::Atom const & atom1,
    conformation::Atom const & atom2,
    Real const weight,
		Real & atrE,
		Real & repE,
		Real & solE,
		Real & d2
	) const {
		atom_pair_energy( atom1, atom2, weight, atrE, repE, solE, d2 );
	}

	inline
	void
	atom_pair_energy(
		conformation::Atom const & atom1,
		conformation::Atom const & atom2,
		Real const weight,
		Real &atr,
		Real &rep,
		Real &solv,
		Real & d2
	) const;

	inline
	Real
	eval_dE_dR_over_r(
		conformation::Atom const & atom1,
		conformation::Atom const & atom2,
		EnergyMap const & weights,
		Vector & f1,
		Vector & f2
	) const;


	/// Inline Methods For Trie-vs-Trie Algorithm
	inline
	Energy heavyatom_heavyatom_energy(
		etrie::EtableAtom const & at1,
		etrie::EtableAtom const & at2,
		DistanceSquared & d2,
		Size & /*path_dist*/
	) const
	{
		Energy atr(0.0), rep(0.0), solv(0.0);
		atom_pair_energy( at1, at2, 1.0, atr, rep, solv, d2 );
		return sum_energies( atr, rep, solv );
	}

	inline
	Energy heavyatom_hydrogenatom_energy(
		etrie::EtableAtom const & at1,
		etrie::EtableAtom const & at2,
		Size & /*path_dist*/
	) const
	{
		Energy atr(0.0), rep(0.0), solv(0.0);
		pair_energy_H( at1, at2, 1.0, atr, rep, solv );
		return sum_energies( atr, rep, solv );
	}

	inline
	Energy hydrogenatom_heavyatom_energy(
		etrie::EtableAtom const & at1,
		etrie::EtableAtom const & at2,
		Size & /*path_dist*/
	) const
	{
		Energy atr(0.0), rep(0.0), solv(0.0);
		pair_energy_H( at1, at2, 1.0, atr, rep, solv );
		return sum_energies( atr, rep, solv );
	}

	inline
	Energy hydrogenatom_hydrogenatom_energy(
		etrie::EtableAtom const & at1,
		etrie::EtableAtom const & at2,
		Size & /*path_dist*/
	) const
	{
		Energy atr(0.0), rep(0.0), solv(0.0);
		pair_energy_H( at1, at2, 1.0, atr, rep, solv );
		return sum_energies( atr, rep, solv );
	}

	inline
	void
	pair_energy_H(
		conformation::Atom const & atom1,
		conformation::Atom const & atom2,
		Real const weight,
		EnergyMap & emap
	) const
	{
		Real atr(0), rep(0), sol(0);
		pair_energy_H( atom1, atom2, weight, atr, rep, sol );
		emap[ st_atr() ] += atr;
		emap[ st_rep() ] += rep;
		emap[ st_sol() ] += sol;
	}

	inline
	void
	pair_energy_H(
		conformation::Atom const & atom1,
		conformation::Atom const & atom2,
		Real const weight,
		Real &atr,
		Real &rep,
		Real &solv
	) const;

private:
	bool
	interpolate_bins(
		conformation::Atom const & atom1,
		conformation::Atom const & atom2,
		Real & d2,
		int &  disbin,
		Real & frac
	) const;

private:

	ObjexxFCL::FArray3D< Real > const & ljatr_;
	ObjexxFCL::FArray3D< Real > const & ljrep_;
	ObjexxFCL::FArray3D< Real > const & solv1_;
	ObjexxFCL::FArray3D< Real > const & solv2_;
	ObjexxFCL::FArray3D< Real > const & dljatr_;
	ObjexxFCL::FArray3D< Real > const & dljrep_;
	ObjexxFCL::FArray3D< Real > const & dsolv_;

	Real safe_max_dis2_;

	int etable_bins_per_A2_;
	Real dis2_step_;

};


///
class TableLookupEtableEnergy : public BaseEtableEnergy< TableLookupEtableEnergy >  {
public:
	typedef BaseEtableEnergy< TableLookupEtableEnergy >  parent;
	typedef TableLookupEvaluator Evaluator;

public:
	using parent::eval_residue_pair_derivatives;
	using parent::eval_intrares_derivatives;

public:

	/// @brief construction with an etable
	TableLookupEtableEnergy(
		Etable const & etable_in,
		methods::EnergyMethodOptions const & options
	);

	/// @brief explicit copy constructor
	TableLookupEtableEnergy( TableLookupEtableEnergy const & );

	/// clone
	EnergyMethodOP
	clone() const;

	void
	setup_for_scoring_(pose::Pose const& pose, scoring::ScoreFunction const&) const;


	//
	/////////////////////////////////////////////////////////////////////////////
	// methods for ContextIndependentTwoBodyEnergies
	/////////////////////////////////////////////////////////////////////////////


	virtual
	bool
	defines_intrares_energy( EnergyMap const & /*weights*/ ) const;

	virtual
	void
	eval_intrares_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;

public:

//	inline
//	void
//	derived_prepare_for_residue_pair(
//		Size const res1,
//		Size const res2,
//		pose::Pose const &
//	) const;
//
	TableLookupEvaluator const & intrares_evaluator() const { return intrares_evaluator_; }
	TableLookupEvaluator const & interres_evaluator() const { return interres_evaluator_; }

	TableLookupEvaluator & intrares_evaluator() { return intrares_evaluator_; }
	TableLookupEvaluator & interres_evaluator() { return interres_evaluator_; }

private:
	// bookkeeping variable -- are intra- or interresidue weights being using at the moment?
	// don't bother to reset the BaseEtable's st_atr_, st_rep_, st_sol_ if possible.
	//mutable bool using_interres_scoretypes_;

	virtual
	core::Size version() const;

	TableLookupEvaluator intrares_evaluator_;
	TableLookupEvaluator interres_evaluator_;

};

class AnalyticEtableEnergy : public BaseEtableEnergy< AnalyticEtableEnergy >  {
public:
	typedef BaseEtableEnergy< AnalyticEtableEnergy >  parent;
	typedef AnalyticEtableEvaluator Evaluator;

public:
	using parent::eval_residue_pair_derivatives;
	using parent::eval_intrares_derivatives;

public:

	/// @brief construction with an etable
	AnalyticEtableEnergy(
		Etable const & etable_in,
		methods::EnergyMethodOptions const & options
	);

	/// @brief explicit copy constructor
	AnalyticEtableEnergy( AnalyticEtableEnergy const & );

	/// clone
	EnergyMethodOP
	clone() const;

	void
	setup_for_scoring_(pose::Pose const& pose, scoring::ScoreFunction const&) const;


	//
	/////////////////////////////////////////////////////////////////////////////
	// methods for ContextIndependentTwoBodyEnergies
	/////////////////////////////////////////////////////////////////////////////


	virtual
	bool
	defines_intrares_energy( EnergyMap const & /*weights*/ ) const;

	virtual
	void
	eval_intrares_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;

public:

//	inline
//	void
//	derived_prepare_for_residue_pair(
//		Size const res1,
//		Size const res2,
//		pose::Pose const &
//	) const;

	AnalyticEtableEvaluator const & intrares_evaluator() const { return intrares_evaluator_; }
	AnalyticEtableEvaluator const & interres_evaluator() const { return interres_evaluator_; }

	AnalyticEtableEvaluator & intrares_evaluator() { return intrares_evaluator_; }
	AnalyticEtableEvaluator & interres_evaluator() { return interres_evaluator_; }

private:
	// bookkeeping variable -- are intra- or interresidue weights being using at the moment?
	// don't bother to reset the BaseEtable's st_atr_, st_rep_, st_sol_ if possible.
	//mutable bool using_interres_scoretypes_;

	virtual
	core::Size version() const;

	mutable AnalyticEtableEvaluator intrares_evaluator_;
	mutable AnalyticEtableEvaluator interres_evaluator_;


};


/*class EtableEnergyEvaluator
{
public:
	// default ctor -- set inter-residue score types
	EtableEnergyEvaluator( EtableEnergy const & etable_energy ) :
		etable_energy_( etable_energy ),
		st_atr_( fa_atr ),
		st_rep_( fa_rep ),
		st_sol_( fa_sol )
	{}

	// switch to intra-res mode
	void set_intrares() {
		st_atr_ = fa_intra_atr;
		st_rep_ = fa_intra_rep;
		st_sol_ = fa_intra_sol;
	}

public:

	Real
	hydrogen_interaction_cutoff2() const {
		return etable_energy_.hydrogen_interaction_cutoff2();
	}

	void
	pair_energy_H(
		Atom const & at1,
		Atom const & at2,
		Real const & cp_weight,
		EnergyMap & emap
	) const;

	void
	atom_pair_energy(
		Atom const & at1,
		Atom const & at2,
		Real const & cp_weight,
		EnergyMap & emap,
		Real & dsq
	) const;


private:
	EtableEnergy const & etable_energy_;
	ScoreType st_atr_;
	ScoreType st_rep_;
	ScoreType st_sol_;
};*/

///////////////////////////////////////////////////////////////////////////////
// inline methods
///////////////////////////////////////////////////////////////////////////////

/// @brief decide between intra and interresidue score terms.  Pseudo polymorphic call
/// from base class
//inline
//void
//TableLookupEtableEnergy::derived_prepare_for_residue_pair(
//	Size const res1,
//	Size const res2,
//	pose::Pose const &
//) const {
//
//	if ( res1 == res2 ) {
//		if ( using_interres_scoretypes_ ) {
//			//std::cout << "setting intaresidue scoretypes for residue " << res1 << std::endl;
//			//set_scoretypes( fa_intra_atr, fa_intra_rep, fa_intra_sol );
//			etable_evaluator()->set_scoretypes( fa_intra_atr, fa_intra_rep, fa_intra_sol );
//			using_interres_scoretypes_ = false;
//		}
//		assert( etable_evaluator()->st_rep() == fa_intra_rep );
//	} else {
//		if ( ! using_interres_scoretypes_ ) {
//			//std::cout << "setting interresidue scoretypes for residue pair " << res1 << " & " << res2 << std::endl;
//			//set_scoretypes( fa_atr, fa_rep, fa_sol );
//			etable_evaluator()->set_scoretypes( fa_atr, fa_rep, fa_sol );
//			using_interres_scoretypes_ = true;
//		}
//		assert( etable_evaluator()->st_rep() == fa_rep );
//	}
//
//}

// TableLookupEvaluator's inline methods
inline
bool
TableLookupEvaluator::interpolate_bins(
	conformation::Atom const & atom1,
	conformation::Atom const & atom2,
	Real & d2,
	int  & disbin,
	Real & frac
) const
{
	d2 = atom1.xyz().distance_squared( atom2.xyz() );

	if ( ( d2 >= safe_max_dis2_ ) || ( d2 == Real(0.0) ) ) {
		return false;
	}

	// bin by distance:
	Real const d2_bin = d2 * etable_bins_per_A2_;
	disbin = static_cast< int >( d2_bin ) + 1;
	//	int const disbin2 = disbin + 1;
	frac = d2_bin - ( disbin - 1 );
	return true;
	//ctsa
	//ctsa  tables have been hacked so that if disbin2 = lastbin, all values = 0.
	//ctsa
}

inline
void
TableLookupEvaluator::atom_pair_energy(
 conformation::Atom const & atom1,
 conformation::Atom const & atom2,
 Real const weight,
 Real & atr,
 Real & rep,
 Real & solv,
 Real & d2
) const
{
	assert( ljatr_.active() );
	int disbin; Real frac;
	atr = rep = solv = 0.0;

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

		e1 = solv1_[ l1 ] + solv2_[ l1 ];
		solv = weight * ( e1 + frac * ( solv1_[ l2 ] + solv2_[l2] - e1 ) );
		//		std::cout << "solv " << solv << std::endl;
		//		std::cerr << "finished evaluating atom _pair energy " << std::endl;

	} //if within cutoff

}



inline
void
TableLookupEvaluator::pair_energy_H(
	conformation::Atom const & atom1,
	conformation::Atom const & atom2,
	Real const weight,
	Real &atr,
	Real &rep,
	Real &solv
) const
{
	assert( ljrep_.active() );
	Real d2,frac;
	int disbin;
	atr = rep = solv = 0.0;

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

inline
void
TableLookupEvaluator::atom_pair_lk_energy_and_deriv_v(
                                                      conformation::Atom const & atom1,
                                                      conformation::Atom const & atom2,
                                                      Real & solv1,
                                                      Real & dsolv1,
                                                      bool const eval_deriv /* = false */
                                                      ) const
{

    int disbin;
    Real frac, d2;

    if (interpolate_bins(atom1,atom2,d2,disbin,frac)) {

        int const l1 = solv1_.index( disbin, atom2.type(), atom1.type()),
        l2 = l1 + 1;

        Real const e1 = solv1_[ l1 ];
        solv1 = ( e1 + frac * ( solv1_[ l2 ] - e1 ) );

        if ( eval_deriv ){
            // Following (commented out) is used in dE_dR_over_R below,
            //  but its a mistake, I think -- rhiju.
            //			Real e1 = dsolv1_[ l1 ];
            //			deriv = ( e1 + frac * ( dsolv1_[ l2 ] - e1 ) );
            dsolv1 = ( solv1_[ l2 ] - solv1_[ l1 ] ) * etable_bins_per_A2_ * std::sqrt( d2 ) * 2;
        }

    } //if within cutoff

}

inline
void
TableLookupEvaluator::atom_pair_lk_energy_and_deriv_v_efficient(
                                                          conformation::Atom const & atom1,
                                                          conformation::Atom const & atom2,
                                                          Real & solv1,
                                                          Real & solv2,
                                                          Real & dsolv1,
                                                          bool const eval_deriv /* = false */
                                                          ) const
{
    int disbin;
    Real frac, d2;

    if (interpolate_bins(atom1,atom2,d2,disbin,frac)) {

        int const l1 = solv1_.index( disbin, atom2.type(), atom1.type()),
        l2 = l1 + 1;

        Real const e1 = solv1_[ l1 ];
        solv1 = ( e1 + frac * ( solv1_[ l2 ] - e1 ) );

        if ( eval_deriv ){
            // Following (commented out) is used in dE_dR_over_R below,
            //  but its a mistake, I think -- rhiju.
            //			Real e1 = dsolv1_[ l1 ];
            //			deriv = ( e1 + frac * ( dsolv1_[ l2 ] - e1 ) );
            dsolv1 = ( solv1_[ l2 ] - solv1_[ l1 ] ) * etable_bins_per_A2_ * std::sqrt( d2 ) * 2;
        }

    } //if within cutoff

	disbin = 0;
	frac = 0.0;
	d2 = 0.0;

	if (interpolate_bins(atom2,atom1,d2,disbin,frac)) {

        int const l1 = solv1_.index( disbin, atom1.type(), atom2.type()),
        l2 = l1 + 1;

        Real const e1 = solv1_[ l1 ];
        solv2 = ( e1 + frac * ( solv1_[ l2 ] - e1 ) );

        if ( eval_deriv ){
            // Following (commented out) is used in dE_dR_over_R below,
            //  but its a mistake, I think -- rhiju.
            //			Real e1 = dsolv1_[ l1 ];
            //			deriv = ( e1 + frac * ( dsolv1_[ l2 ] - e1 ) );
            dsolv1 = ( solv1_[ l2 ] - solv1_[ l1 ] ) * etable_bins_per_A2_ * std::sqrt( d2 ) * 2;
        }

    } //if within cutoff
}

inline
void
EtableEvaluator::atom_pair_lk_energy_and_deriv_v_efficient(
	 conformation::Atom const &,
	 conformation::Atom const &,
	 Real &,
	 Real &,
	 Real &,
	 bool const /* = false */
) const {
	//AnalyticEtableEvaluator::atom_pair_lk_energy_and_deriv_v_efficient(atom1, atom2, solv1, solv2, dsolv1, eval_deriv);
	return;
}

Real
TableLookupEvaluator::eval_dE_dR_over_r(
 conformation::Atom const & atom1,
 conformation::Atom const & atom2,
 EnergyMap const & weights,
 Vector & f1,
 Vector & f2
) const
{
	Real d2,frac;
	int disbin;

	if ( atom1.xyz().distance_squared( atom2.xyz() ) > safe_max_dis2_ ) return 0.;

	if ( interpolate_bins(atom1,atom2,d2,disbin,frac) ) {


		f1 = atom1.xyz().cross( atom2.xyz() );
		f2 = atom1.xyz() - atom2.xyz();

		// l1 and l2 are FArray LINEAR INDICES for fast lookup:
		// [ l1 ] == (disbin  ,attype2,attype1)
		// [ l2 ] == (disbin+1,attype2,attype1)

		/// BEGIN DERIVATIVE INTERPOLATION
		Real deriv = 0.0;

		int const l1 = dljatr_.index( disbin, atom1.type(), atom2.type()),
				l2 = l1 + 1;

		Real e1 = dljatr_[ l1 ];
		deriv = weights[ st_atr() ] * ( e1 + frac * ( dljatr_[ l2 ] - e1 ) );

		e1 = dljrep_[ l1 ];
		deriv += weights[ st_rep() ] * ( e1 + frac * ( dljrep_[ l2 ] - e1 ) );

		e1 = dsolv_[ l1 ];
		deriv += weights[ st_sol() ] * ( e1 + frac * ( dsolv_[ l2 ] - e1 ) );

		return deriv / std::sqrt( d2 );
		/// END DERIVATIVE INTERPOLATION

		/// BEGIN EXACT DERIVATIVE CALCULATION

		/*Real deriv( 0.0 );
		int const l1 = ljatr_.index( disbin, atom1.type(), atom2.type()),
			l2 = l1 + 1;

		// d g(x) / dx with g(x) = x^2 ---> 2x;  x is the distance
		// we want to avoid the sqrt. Since at1-at2 (f2) already is the right length, consider it pre-multiplied by sqrt(d2).
		// so what we have below is 2x / ( d2step * x ) = 2 / d2step = 2 * inv(d2step).  This will be multiplied by the un-normalized f1 and f2 vectors.
		Real dxsquared_dx_times_x2step_over_x =  2 * etable_bins_per_A2;

		Real const atr1 = ljatr_[ l1 ];
		Real const atr2 = ljatr_[ l2 ];
		Real const rep1 = ljrep_[ l1 ];
		Real const rep2 = ljrep_[ l2 ];
		Real const sol1 = solv1_[ l1 ] + solv2_[ l1 ];
		Real const sol2 = solv1_[ l2 ] + solv2_[ l2 ];

		deriv = weights[ st_atr_ ] * ( atr2 - atr1 );
		deriv += weights[ st_rep_ ] * ( rep2 - rep1 );
		deriv += weights[ st_sol_ ] * ( sol2 - sol1 );

		return deriv * dxsquared_dx_times_x2step_over_x;*/

		/// TEMP
		/*
		std::cout << "Testing numerically: ";
		Real f11( 0.0 );
		Real step = 0.00001;
		{// scope
			Real altd = std::sqrt( d2 ) - step;
			Real altd2 = altd*altd;
			Real altfrac = ( altd2 * etable_bins_per_A2 - ( disbin - 1 ) );
			int const altl1 = ljatr_.index( disbin, atom1.type(), atom2.type()),
				altl2 = altl1 + 1;


			Real e1 = ljatr_[ altl1 ];
			f11 = weights[ st_atr_ ] * ( e1 + altfrac * ( ljatr_[ altl2 ] - e1 ) );

			e1 = ljrep_[ altl1 ];
			f11 += weights[ st_rep_ ] * ( e1 + altfrac * ( ljrep_[ altl2 ] - e1 ) );

			e1 = solv1_[ altl1 ] + solv2_[ altl1 ];
			f11 += weights[ st_sol_ ] * ( e1 + altfrac * ( solv1_[ altl2 ] + solv2_[altl2] - e1 ) );
		}

		Real f22(0.0);
		{// scope
			Real altd = std::sqrt( d2 ) + step;
			Real altd2 = altd*altd;
			Real altfrac = ( altd2 * etable_bins_per_A2 - ( disbin - 1 ) );
			int const altl1 = ljatr_.index( disbin, atom1.type(), atom2.type()),
				altl2 = altl1 + 1;


			Real e1 = ljatr_[ altl1 ];
			f22 = weights[ st_atr_ ] * ( e1 + altfrac * ( ljatr_[ altl2 ] - e1 ) );

			e1 = ljrep_[ altl1 ];
			f22 += weights[ st_rep_ ] * ( e1 + altfrac * ( ljrep_[ altl2 ] - e1 ) );

			e1 = solv1_[ altl1 ] + solv2_[ altl1 ];
			f22 += weights[ st_sol_ ] * ( e1 + altfrac * ( solv1_[ altl2 ] + solv2_[altl2] - e1 ) );
		}
		std::cout << "Deriv discrep: " << ( f22 - f11 ) / (2 * step ) << " vs " << deriv * dxsquared_dx_times_x2step_over_x * std::sqrt( d2 ) << std::endl;
		*/

	} else {
		return 0.0;
	}
}


/// AnalyticEtableEvaluator inline methods.

void
AnalyticEtableEvaluator::atom_pair_energy(
 conformation::Atom const & atom1,
 conformation::Atom const & atom2,
 Real const weight,
 Real & atr,
 Real & rep,
 Real & solv,
 Real & d2
) const
{
	atr = rep = solv = 0.0;

	//etable_.interpolated_analytic_etable_evaluation( atom1, atom2, atr, rep, solv, d2 );
	etable_.analytic_etable_evaluation( atom1, atom2, atr, rep, solv, d2 );
	atr  *= weight;
	rep  *= weight;
	solv *= weight;
	return;
}

Real
AnalyticEtableEvaluator::eval_dE_dR_over_r(
 conformation::Atom const & atom1,
 conformation::Atom const & atom2,
 EnergyMap const & weights,
 Vector & f1,
 Vector & f2
) const
{
	if ( atom1.xyz().distance_squared( atom2.xyz() ) > safe_max_dis2_ ) return 0.;

	f1 = atom1.xyz().cross( atom2.xyz() );
	f2 = atom1.xyz() - atom2.xyz();

	Real datr, drep, dsol, invd;
	etable_.analytic_etable_derivatives( atom1, atom2, datr, drep, dsol, invd );
	return ( weights[ st_atr() ] * datr + weights[ st_rep() ] * drep + weights[ st_sol() ] * dsol ) * invd;

}


} // methods
} // scoring
} // core


#endif // INCLUDED_core_scoring_EtableEnergy_HH
