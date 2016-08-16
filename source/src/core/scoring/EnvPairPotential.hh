// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/EnvPairPotential.hh
/// @brief  Statistically derived rotamer pair potential class implementation
/// @author Phil Bradley
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_core_scoring_EnvPairPotential_hh
#define INCLUDED_core_scoring_EnvPairPotential_hh

#include <core/types.hh>

// Unit headers
#include <core/scoring/EnvPairPotential.fwd.hh>

// Package headers
#include <core/conformation/Residue.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

#include <basic/datacache/CacheableData.hh>

// Utility headers

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>

#include <utility/vector1.hh>


// C++


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Keep track of the cenlist information

class CenListInfo : public basic::datacache::CacheableData {

public:

	CenListInfo(): calculated_(false) {};

	CenListInfo( CenListInfo const & src );

	basic::datacache::CacheableDataOP
	clone() const
	{
		return basic::datacache::CacheableDataOP( new CenListInfo( *this ) );
	}

	Size
	size() const {
		return fcen6_.size();
	}

	Real
	fcen6( Size const seqpos ) const
	{
		return fcen6_[ seqpos ];
	}

	Real &
	fcen6( Size const seqpos )
	{
		return fcen6_[ seqpos ];
	}

	Real
	fcen10( Size const seqpos ) const
	{
		return fcen10_[ seqpos ];
	}

	Real &
	fcen10( Size const seqpos )
	{
		return fcen10_[ seqpos ];
	}

	Real
	fcen12( Size const seqpos ) const
	{
		return fcen12_[ seqpos ];
	}

	Real &
	fcen12( Size const seqpos )
	{
		return fcen12_[ seqpos ];
	}

	bool
	calculated() const
	{
		return calculated_;
	}

	bool &
	calculated()
	{
		return calculated_;
	}

	void
	initialize( pose::Pose const & pose );

	// Setter functions
	void
	set_fcen6( Size const seqpos, Real value )
	{
		fcen6_[ seqpos ] = value;
	}

	void
	set_fcen10( Size const seqpos, Real value )
	{
		fcen10_[ seqpos ] = value;
	}

	void
	set_fcen12( Size const seqpos, Real value )
	{
		fcen12_[ seqpos ] = value;
	}

private:
	utility::vector1< Real > fcen6_;
	utility::vector1< Real > fcen10_;
	utility::vector1< Real > fcen12_;
	bool calculated_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

class EnvPairPotential : public utility::pointer::ReferenceCount {

public:
	EnvPairPotential();


	void
	compute_centroid_environment(
		pose::Pose & pose
	) const;

	void
	finalize( pose::Pose & pose ) const;


	void
	evaluate_env_and_cbeta_scores(
		pose::Pose const & pose,
		conformation::Residue const & rsd,
		Real & env_score,
		Real & cb_score6,
		Real & cb_score12
	) const;


	void
	evaluate_pair_and_cenpack_score(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		Real const cendist,
		Real & pair_contribution,
		Real & cenpack_contribution
	) const;

private:
	// this shouldn't be protected
	Real const cen_dist_cutoff2;

public:
	static CenListInfo const & cenlist_from_pose( pose::Pose const & );
	static CenListInfo & nonconst_cenlist_from_pose( pose::Pose & );

private:

	void
	fill_cenlist(
		CenListInfo & cenlist,
		Size const res1,
		Size const res2,
		Real const cendist
	) const;

	void
	truncate_cenlist_values( CenListInfo & cenlist ) const;

private: // data

	ObjexxFCL::FArray2D< Real > env_log_;
	ObjexxFCL::FArray1D< Real > cbeta_den6_;
	ObjexxFCL::FArray1D< Real > cbeta_den12_;
	ObjexxFCL::FArray3D< Real > pair_log_;
	ObjexxFCL::FArray1D< Real > cenpack_log_;

	// unused Real const cen_dist6sqr_;
	// unused Real const cen_dist10sqr_;
	// unused Real const cen_dist12sqr_;

	//cems transition regions between environment bins
	//cems transition is from +/- sqrt(36+pad6) +/- sqrt(100+pad10) etc
	Real const cen_dist5_pad;
	Real const cen_dist6_pad;
	Real const cen_dist7_pad;
	Real const cen_dist10_pad;
	Real const cen_dist12_pad;

	Real const cen_dist5_pad_plus ;
	Real const cen_dist6_pad_plus ;
	Real const cen_dist7_pad_plus ;
	Real const cen_dist10_pad_plus;
	Real const cen_dist12_pad_plus;

	Real const cen_dist5_pad_minus ;
	Real const cen_dist7_pad_minus ;
	Real const cen_dist10_pad_minus;
	Real const cen_dist12_pad_minus;

	Real const cen_dist5_pad_hinv ;
	Real const cen_dist6_pad_hinv ;
	Real const cen_dist7_pad_hinv ;
	Real const cen_dist10_pad_hinv;
	Real const cen_dist12_pad_hinv;

	Real const cen_dist_cutoff_12_pad;

};
} // ns scoring
} // ns core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_EnvPairPotential )
#endif // SERIALIZATION


#endif
