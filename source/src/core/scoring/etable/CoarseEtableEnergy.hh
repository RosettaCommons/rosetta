// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/CoarseEtableEnergy.hh
/// @brief  Energy method for coarse grained residue representation class header
/// @author Oliver


#ifndef INCLUDED_core_scoring_etable_CoarseEtableEnergy_hh
#define INCLUDED_core_scoring_etable_CoarseEtableEnergy_hh

// Unit headers
#include <core/scoring/etable/CoarseEtableEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/etable/BaseEtableEnergy.hh>
//#include <core/scoring/etable/BaseEtableEnergy.tmpl.hh>
// AUTO-REMOVED #include <core/scoring/etable/Etable.hh>
#include <core/scoring/ContextGraphTypes.hh>
//#include <core/coarse/CoarseEtable.hh>

// Project headers
#include <core/conformation/Atom.hh>
#include <core/pose/Pose.fwd.hh>

#include <ObjexxFCL/FArray3D.hh>

namespace core {
namespace scoring {
namespace etable {


///

class CoarseEtableEnergy;


class CoarseEtableEnergy : public BaseEtableEnergy< CoarseEtableEnergy > {
public:
	typedef BaseEtableEnergy< CoarseEtableEnergy > parent;
public:
	/// @brief Constructor that belongs in a .cc file!
	CoarseEtableEnergy (
		Etable const & etable_in,
		methods::EnergyMethodOptions const & options,
		coarse::CoarseEtableCAP coarse_etable_in
	);

	//clone
	EnergyMethodOP
	clone() const {
		return new CoarseEtableEnergy( *this );
	}

	void
	derived_prepare_for_residue_pair(
			Size const res1,
		Size const res2,
		pose::Pose const & pose
	) const;

	void indicate_required_context_graphs( utility::vector1< bool > & context_graphs_required ) const
	{
		context_graphs_required[ ten_A_neighbor_graph ] = true;
	}

	//
	void
	setup_for_scoring_(pose::Pose const& pose, scoring::ScoreFunction const&) const;

public:
	// implementation for quasi-virtual functions
	inline
	void
	atom_pair_energy_(
		conformation::Atom const & atom1,
		conformation::Atom const & atom2,
		Real const weight,
		Energy &atr,
		Energy &rep,
		Energy &solv,
		Energy &bb,
		Real & dsq
	) const;

	//
	inline
	void
	pair_energy_H_(
		 conformation::Atom const & atom1,
		 conformation::Atom const & atom2,
		 Real weight,
		 Energy &atr,
		 Energy &rep,
		 Energy &solv,
		 Energy &bb
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

	virtual
	bool
	defines_intrares_energy( EnergyMap const & /*weights*/ ) const
	{
		return false;
	}

	virtual
	void
	eval_intrares_energy(
		conformation::Residue const &,
		pose::Pose const &,
		ScoreFunction const &,
		EnergyMap &
	) const {}

private:
	coarse::CoarseEtableCAP coarse_etable_;

};

	//inline methods

inline
void
CoarseEtableEnergy::atom_pair_energy_(
	conformation::Atom const & atom1,
	conformation::Atom const & atom2,
	Real const weight,
	Energy &atr,
	Energy &rep,
	Energy &solv,
	Energy &bb,
	Real & d2
) const {
	//	std::cerr << __FILE__<< ' ' << __LINE__ << std::endl;
	if ( coarse_etable_->handles( atom1, atom2 ) ) {
		int disbin;
		Real frac;
		if (interpolate_bins( atom1, atom2, d2, disbin, frac )) {
			coarse_etable_->atom_pair_energy( disbin, frac, atom1, atom2, bb);
			atr = rep = solv = 0.0;
		}
	} else {
		parent::atom_pair_energy_( atom1, atom2, weight, atr, rep, solv, bb, d2);
	}
}
///////////////////////////////////////////////////////////////////////////////

/// How does the coarse_etable_ know which residue pair this is?
/// when is coarse_etable_->set_residue_pair( rsd1, rsd2 ); called?
inline
Real
CoarseEtableEnergy::eval_dE_dR_over_r_(
	conformation::Atom const & atom1,
	conformation::Atom const & atom2,
	EnergyMap const & weights,
	Vector & f1,
	Vector & f2
) const
{
	if ( coarse_etable_->handles(atom1,atom2) ) {
		f1 = atom1.xyz().cross( atom2.xyz() );
		f2 = atom1.xyz() - atom2.xyz();
		Real d2,frac;
		int disbin;
		if (interpolate_bins(atom1,atom2,d2,disbin,frac)) {
			return coarse_etable_->eval_dE_dR(disbin,frac,atom1,atom2,weights);
		} else {
			return 0.0;
		}
	} else {
		return parent::eval_dE_dR_over_r_(atom1,atom2,weights,f1,f2);
	}
}


///////////////////////////////////////////////////////////////////////////////

inline
void
CoarseEtableEnergy::pair_energy_H_(
	conformation::Atom const & atom1,
	conformation::Atom const & atom2,
	Real weight,
	Energy &atr,
	Energy &rep,
	Energy &solv,
	Energy &bb
) const
{

	if ( coarse_etable_->handles( atom1, atom2) ) {
		Real d2,frac; int disbin;
		if ( interpolate_bins( atom1, atom2, d2, disbin, frac) ) {
			coarse_etable_->atom_pair_energy( disbin, frac, atom1, atom2, bb);
			atr = rep = solv = 0.0;
		}
	} else {
		parent::pair_energy_H_( atom1, atom2, weight, atr, rep, solv, bb );
	}
}



}
}
}

#endif
