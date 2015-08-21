// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//////////////////////////////////////////////////////////////////////
///
/// @brief
/// A class for scoring fa_atr, fa_rep, fa_sol
///
/// @details
/// This class is invoked when scoring. The terms that it is responsible for are fa_atr,
/// fa_rep, fa_sol. The class is highly optimized for speed and will break if you are not
/// careful. It calls functions within BaseEtableEnergy, which actually does the scoring.
/// It passes an energy map, which contains the energy for that residue for the atr, rep, and sol
/// terms. This is modified once the scores are tallied in BaseEtableEnergy.
///
///
/// @author
/// Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// Kevin P. Hinshaw (KevinHinshaw@gmail.com)
/// Andrew Leaver-Fay (leaverfa@email.unc.edu)
/// Steven Combs - comments and skipping of virtual atoms
///
/////////////////////////////////////////////////////////////////////////
#ifndef INCLUDED_core_scoring_etable_atom_pair_energy_inline_hh
#define INCLUDED_core_scoring_etable_atom_pair_energy_inline_hh

#include <core/scoring/types.hh>

#include <core/conformation/Atom.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AtomType.hh>


#include <numeric/xyzVector.hh>

#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreType.hh>
#include <utility/vector1.hh>


namespace core {
namespace scoring {


///////////////////////////////////////////////////////////////////////////////

/// class T must define
// operator() ( int, int, Real & ) const;
// class T_Etable must define
// pair_energy_H( Atom const &, Atom const &, Real, EnergyMap & )
template< class T, class T_Etable >
inline
void
residue_fast_pair_energy_attached_H(
	conformation::Residue const & res1,
	int const atomno1,
	conformation::Residue const & res2,
	Size const atomno2,
	Size const at1hbegin, //at1hbegin and at1hend define a range of hydrogen atom indices -- those h's bound to at1
	Size const at1hend,
	Size const at2hbegin,
	Size const at2hend,
	T const & count_pair,
	T_Etable const & etable_energy,
	EnergyMap &emap
)
{
	using conformation::Atom;

	Weight weight( 1.0 );
	Size path_dist( 0 );
	Atom const & atom1( res1.atom( atomno1 ) );
	Atom const & atom2( res2.atom( atomno2 ) );


	// Heavy Atom in res1 to Hs in res2
	for ( Size i = at2hbegin; i<= at2hend; ++i ) {
		Atom const & H2( res2.atom( i ) );
		weight = 1.0;
		path_dist = 0;
		if ( count_pair( atomno1, i, weight, path_dist ) ) {
			etable_energy.pair_energy_H( atom1, H2, weight, emap );
		}
	}


	// Hs in res1 to heavy Atom and Hs in res2
	for ( Size i = at1hbegin; i<= at1hend; ++i ) {
		Atom const & H1( res1.atom(i) );
		weight = 1.0;
		path_dist = 0;
		// H in res1 to heavy Atom in res2
		if ( count_pair( i, atomno2, weight, path_dist ) ) {
			etable_energy.pair_energy_H( H1, atom2, weight, emap );
		}

		// H in res1 to Hs in res2
		for ( Size j = at2hbegin; j<= at2hend; ++j ) {
			Atom const & H2( res2.atom(j) );
			weight = 1.0f;
			path_dist = 0;
			if ( count_pair( i, j, weight, path_dist ) ) {
				etable_energy.pair_energy_H( H1, H2, weight, emap );
			}
		}
	}
}


///////////////////////////////////////////////////////////////////////////////


/// @brief templated atom pair energy calculations
///
/// loops over the heavy atoms of residue1 and the heavy atoms of residue2,
/// evaluates their energies, and if a pair of heavy atoms is close enough,
/// descendes into the attached hydrogen atoms for each.
///
/// Templates are for count_pair type resolution and etable type resolution: there are
/// no polymorphic lookups within these functions
///
/// class T must define
// operator() ( int, int, Real & ) const;
///
/// class T_Etable must define
/// atom_pair_energy( Atom const &, Atom const &, Real, EnergyMap &, Distance ) and
//. pair_energy_H( Atom const &, Atom const &, Real, EnergyMap & )
template < class T, class T_Etable >
inline
void
inline_residue_atom_pair_energy(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	T_Etable const & etable_energy,
	T const & count_pair,
	EnergyMap & emap,
	int res1_start,
	int res1_end,
	int res2_start,
	int res2_end
)
{
	using conformation::Atom;

	//std::cout << "inline_residue_atom_pair_energy( res" << res1.seqpos() << ", res";
	//std::cout << res2.seqpos() << ");" << std::endl;
	//debug_assert ( !( res1.is_DNA() && res2.is_DNA() ) ); // pb for testing, will remove

	DistanceSquared dsq;
	Weight weight;
	Size path_dist;
	// get hydrogen interaction cutoff
	Real const Hydrogen_interaction_cutoff2
		( etable_energy.hydrogen_interaction_cutoff2() );

	typedef utility::vector1< Size > const & vect;

	vect r1hbegin( res1.attached_H_begin() );
	vect r1hend(   res1.attached_H_end()   );
	vect r2hbegin( res2.attached_H_begin() );
	vect r2hend(   res2.attached_H_end()   );

	// Atom pairs
	for ( int i = res1_start, i_end = res1_end; i <= i_end; ++i ) {
		Atom const & atom1( res1.atom(i) );
		//get virtual information
		bool atom1_virtual(res1.atom_type(i).is_virtual());
		for ( int j=res2_start, j_end = res2_end; j <= j_end; ++j ) {
			//check if virtual
			bool atom2_virtual(res2.atom_type(j).is_virtual());
			Atom const & atom2( res2.atom(j) );
			weight = 1.0;
			path_dist = 0;
			if ( atom1_virtual || atom2_virtual ) {
				// NOOP! etable_energy.virtual_atom_pair_energy(emap);
			} else {
				if ( count_pair( i, j, weight, path_dist ) ) {
					//       std::cout << "Atom Pair Energy: " << i << " with " << j << "   ";
					etable_energy.atom_pair_energy( atom1, atom2, weight, emap, dsq );
					//     std::cout << "atr: " << emap[ coarse_fa_atr ] << " ";
					//     std::cout << "rep: " << emap[ coarse_fa_rep ] << " ";
					//     std::cout << "sol: " << emap[ coarse_fa_sol ] << " ";
					//     std::cout << std::endl;

				} else {
					dsq = atom1.xyz().distance_squared( atom2.xyz() );
				}
				if ( dsq < Hydrogen_interaction_cutoff2 ) {
					residue_fast_pair_energy_attached_H(
						res1, i, res2, j,
						r1hbegin[ i ], r1hend[ i ],
						r2hbegin[ j ], r2hend[ j ],
						count_pair, etable_energy , emap);
					//   std::cout << "atr: " << emap[ fa_atr ] << " ";
					//   std::cout << "rep: " << emap[ fa_rep ] << " ";
					//   std::cout << std::endl;

				}
			}
		}
	}
}

/// @brief intraresidue atom pair energy evaluations
template < class T, class T_Etable >
inline
void
inline_intraresidue_atom_pair_energy(
	conformation::Residue const & res,
	T_Etable const & etable_energy,
	T const & count_pair,
	EnergyMap & emap
)
{
	using conformation::Atom;

	//std::cout << "inline_residue_atom_pair_energy( res" << res1.seqpos() << ", res";
	//std::cout << res2.seqpos() << ");" << std::endl;

	DistanceSquared dsq;
	Weight weight;
	Size path_dist;

	// get hydrogen interaction cutoff
	Real const Hydrogen_interaction_cutoff2
		( etable_energy.hydrogen_interaction_cutoff2() );

	typedef utility::vector1< Size > const & vect;

	vect rhbegin( res.attached_H_begin() );
	vect rhend(   res.attached_H_end()   );

	int const resnheavyatoms = res.nheavyatoms();

	// Atom pairs
	for ( int i=1; i <= resnheavyatoms; ++i ) {
		Atom const & atom1( res.atom(i) );
		bool atom1_virtual(res.atom_type(i).is_virtual());
		for ( int j=i+1; j <= resnheavyatoms; ++j ) {
			bool atom2_virtual(res.atom_type(j).is_virtual());
			Atom const & atom2( res.atom(j) );
			weight = 1.0;
			path_dist = 0;
			if ( atom1_virtual || atom2_virtual ) {
				//etable_energy.virtual_atom_pair_energy(emap);
			} else {

				if ( count_pair( i, j, weight, path_dist ) ) {
					// std::cout << "Intra Res Atom Pair Energy: atom " << i << " with atom " << j << "   ";
					etable_energy.atom_pair_energy( atom1, atom2, weight, emap, dsq );
					//std::cout << "atr: " << emap[ coarse_fa_atr ] << " ";
					//std::cout << "rep: " << emap[ coarse_fa_rep ] << " ";
					//std::cout << "sol: " << emap[ coarse_fa_sol ] << " ";
					//std::cout << std::endl;
				} else {
					dsq = atom1.xyz().distance_squared( atom2.xyz() );
				}

				if ( dsq < Hydrogen_interaction_cutoff2 ) {
					residue_fast_pair_energy_attached_H(
						res, i, res, j,
						rhbegin[ i ], rhend[ i ],
						rhbegin[ j ], rhend[ j ],
						count_pair, etable_energy, emap);
					//std::cout << "atr: " << emap[ fa_atr ] << " ";
					//std::cout << "rep: " << emap[ fa_rep ] << " ";
					//std::cout << std::endl;

				}
			}
		}
	}
}


template < class T, class T_Etable >
inline
void
inline_residue_atom_pair_energy(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	T_Etable const & etable_energy,
	T const & count_pair,
	EnergyMap & emap
)
{
	inline_residue_atom_pair_energy(
		res1, res2, etable_energy, count_pair, emap,
		1, res1.nheavyatoms(), 1, res2.nheavyatoms() );
}

template < class T, class T_Etable>
inline
void
inline_residue_atom_pair_energy_sidechain_backbone(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	T_Etable const & etable_energy,
	T const & count_pair,
	EnergyMap & emap
)
{
	inline_residue_atom_pair_energy(
		res1, res2, etable_energy, count_pair, emap,
		res1.first_sidechain_atom(), res1.nheavyatoms(), 1, res2.last_backbone_atom() );
}

template < class T, class T_Etable >
inline
void
inline_residue_atom_pair_energy_sidechain_whole(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	T_Etable const & etable_energy,
	T const & count_pair,
	EnergyMap & emap
)
{
	inline_residue_atom_pair_energy(
		res1, res2, etable_energy, count_pair, emap,
		res1.first_sidechain_atom(), res1.nheavyatoms(), 1, res2.nheavyatoms() );
}

template < class T, class T_Etable>
inline
void
inline_residue_atom_pair_energy_backbone_backbone(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	T_Etable const & etable_energy,
	T const & count_pair,
	EnergyMap & emap
)
{
	inline_residue_atom_pair_energy(
		res1, res2, etable_energy, count_pair, emap,
		1, res1.last_backbone_atom(), 1, res2.last_backbone_atom() );
}

template < class T, class T_Etable>
inline
void
inline_residue_atom_pair_energy_sidechain_sidechain(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	T_Etable const & etable_energy,
	T const & count_pair,
	EnergyMap & emap
)
{
	inline_residue_atom_pair_energy(
		res1, res2, etable_energy, count_pair, emap,
		res1.first_sidechain_atom(), res1.nheavyatoms(), res2.first_sidechain_atom(), res2.nheavyatoms() );
}


} // namespace scoring
} // namespace core

#endif
