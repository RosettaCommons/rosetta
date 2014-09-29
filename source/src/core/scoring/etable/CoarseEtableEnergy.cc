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


// Unit headers
#include <core/scoring/etable/CoarseEtableEnergy.hh>
#include <core/scoring/etable/CoarseEtableEnergyCreator.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/etable/BaseEtableEnergy.hh>
#include <core/scoring/etable/BaseEtableEnergy.tmpl.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/ContextGraphTypes.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
//#include <core/coarse/CoarseEtable.hh>

// Project headers
#include <core/conformation/Atom.hh>
#include <core/pose/Pose.hh>

#include <ObjexxFCL/FArray3D.hh>

//Auto Headers
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/scoring/etable/etrie/CountPairData_1_1.hh>
#include <core/scoring/etable/etrie/CountPairData_1_2.hh>
#include <core/scoring/etable/etrie/CountPairData_1_3.hh>
#include <core/scoring/etable/etrie/EtableAtom.hh>
#include <core/scoring/trie/trie.functions.hh>


namespace core {
namespace scoring {
namespace etable {


/// @details This must return a fresh instance of the CoarseEtableEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
CoarseEtableEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {
	return new CoarseEtableEnergy( *( ScoringManager::get_instance()->etable( options.etable_type() )),
      options,
      ScoringManager::get_instance()->coarse_etable( chemical::COARSE_TWO_BEAD ) );
}

ScoreTypes
CoarseEtableEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( coarse_fa_atr );
	sts.push_back( coarse_fa_rep );
	sts.push_back( coarse_fa_sol );
	return sts;
}


///

CoarseEtableEnergy::CoarseEtableEnergy(
	Etable const & etable_in,
	methods::EnergyMethodOptions const & options,
	coarse::CoarseEtableCAP coarse_etable_in
) :
 	BaseEtableEnergy< CoarseEtableEnergy >(
		methods::EnergyMethodCreatorOP( new CoarseEtableEnergyCreator ),
		etable_in, options,
		coarse_fa_atr, coarse_fa_rep, coarse_fa_sol ),
	coarse_etable_( coarse_etable_in )
{
	if ( coarse_etable_in ) {
		coarse_etable_in->check_atomset_compatibility( etable_in.atom_set() );
	} else {
		utility_exit_with_message("coarse etable not present ");
	}
}


void
CoarseEtableEnergy::derived_prepare_for_residue_pair(
		Size const res1,
	Size const res2,
	pose::Pose const & pose
) const {
	coarse_etable_->set_residue_pair( pose.residue( res1 ), pose.residue( res2 ) );
}



void
CoarseEtableEnergy::setup_for_scoring_(
	pose::Pose const & pose,
	scoring::ScoreFunction const &
) const {
	if (pose.total_residue()) {
		if ( pose.residue(1).type().atom_type_set_ptr() != coarse_etable_->atom_set() ) {
			utility_exit_with_message( "Illegal attempt to score with non-identical atom set between pose and etable " );
		}
	}
}


/*
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
	) const {
		coarse_etable_->set_residue_pair( pose.residue( res1 ), pose.residue( res2 ) );
	}

	void indicate_required_context_graphs( utility::vector1< bool > & context_graphs_required ) const
	{
		context_graphs_required[ ten_A_neighbor_graph ] = true;
	}

	//
	void
	setup_for_scoring_(pose::Pose const& pose, scoring::ScoreFunction const&) const {
		if (pose.total_residue()) {
			if ( pose.residue(1).type().atom_type_set_ptr() != coarse_etable_->atom_set() ) {
				utility_exit_with_message( "Illegal attempt to score with non-identical atom set between pose and etable " );
			}
		}
	}


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
	defines_intrares_energy( EnergyMap const & weights ) const
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

*/

}
}
}

