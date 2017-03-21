// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/scoring/dna/DNA_EnvPairPotential.hh
/// @brief  dna scoring
/// @author Phil Bradley

#ifndef INCLUDED_core_scoring_dna_DNA_EnvPairPotential_HH
#define INCLUDED_core_scoring_dna_DNA_EnvPairPotential_HH

#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/chemical/AA.hh>

#include <utility/vector1.hh>
#include <utility/vector0.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.fwd.hh>

#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/FArray4D.hh>

namespace core {
namespace scoring {
namespace dna {


class DNA_EnvPairPotential : public utility::pointer::ReferenceCount {
public:
	typedef chemical::AA AA;

public:
	/// ctor -- read the datafiles
	DNA_EnvPairPotential();

	///
	Real
	nbr_dis2_threshold() const
	{
		return nbr_dis2_threshold_;
	}

	///
	Vector
	centroid_xyz( conformation::Residue const & rsd ) const;

	///
	Real
	residue_pair_score( AA const & na, AA const & aa, Real const dis2 ) const;

	///
	Real
	residue_env_score( AA const & aa, Size const nbr_count ) const;

private:
	void
	load_score_tables();

	inline
	Size
	get_pair_disbin( Real const dis2 ) const;

	inline
	Size
	get_env_nbr_bin( Size const nbr_count ) const;

private:
	static Real const nbr_dis2_threshold_;

	static Real const min_env_enrichment_;
	static Real const max_env_enrichment_;
	static Real const min_pair_enrichment_;
	static Real const max_pair_enrichment_;

private:

	Size n_pair_disbins_;
	Size max_nbrs_;

	utility::vector1< utility::vector1< utility::vector1< Real > > > pair_stats_;
	utility::vector1< utility::vector1< Real > > env_stats_;

	utility::vector1< Real > pair_bin_lower_bounds_;
	utility::vector0< Size > nbr_bin_; // NOTE THIS IS A VECTOR0 !!!! SINCE WE CAN HAVE 0 DNA NEIGHBORS
};

inline
Size
DNA_EnvPairPotential::get_pair_disbin( Real const dis2 ) const
{
	debug_assert( dis2 <= nbr_dis2_threshold_ );
	for ( Size bin= n_pair_disbins_; bin >= 2; --bin ) {
		if ( pair_bin_lower_bounds_[ bin ] < dis2 ) return bin;
	}
	return 1;
}

inline
Size
DNA_EnvPairPotential::get_env_nbr_bin( Size const nbr_count ) const
{
	return nbr_bin_[ std::min( max_nbrs_, nbr_count ) ];
}



} // namespace dna
} // scoring
} // core

#endif
