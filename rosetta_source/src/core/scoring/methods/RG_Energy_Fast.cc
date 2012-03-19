// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/RG_Energy_Fast.cc
/// @brief  Radius of gyration energy function definition.
/// @author James Thompson


// Unit headers
#include <core/scoring/methods/RG_Energy_Fast.hh>
#include <core/scoring/methods/RG_Energy_FastCreator.hh>

// Package headers
#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>
#include <core/chemical/VariantType.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/id/AtomID.hh>

// Utility headers
#include <basic/prof.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <core/pose/datacache/CacheableDataType.hh>

#include <core/scoring/EnergyMap.hh>
#include <utility/vector1.hh>



// C++


namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the RG_Energy_Fast class,
/// never an instance already in use
methods::EnergyMethodOP
RG_Energy_FastCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return new RG_Energy_Fast;
}

ScoreTypes
RG_Energy_FastCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( rg );
	return sts;
}


/// c-tor
RG_Energy_Fast::RG_Energy_Fast() :
	parent( new RG_Energy_FastCreator )
{}


/// clone
EnergyMethodOP
RG_Energy_Fast::clone() const
{
	return new RG_Energy_Fast;
}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

void
RG_Energy_Fast::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap & totals
) const {
	using namespace conformation;

	PROF_START( basic::RG );

	totals[ rg ] = calculate_rg_score( pose );

	PROF_STOP( basic::RG );
} // finalize_total_energy

core::Real
RG_Energy_Fast::calculate_rg_score( core::pose::Pose const & pose ) const
{

	Size const nres( pose.total_residue() );
	Size nres_counted=0;

	///////////////////////////////////////
	//
	// RG SCORE

	// calculate center of mass
	Vector center_of_mass( 0, 0, 0 );
	for ( Size i = 1; i <= nres; ++i ) {
		// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
		if ( pose.residue(i).has_variant_type( core::chemical::REPLONLY ) ) continue;
		if ( pose.residue(i).aa() == core::chemical::aa_vrt ) continue;

		Vector const v( pose.residue(i).nbr_atom_xyz() );
		center_of_mass += v;
		nres_counted++;
	}
	center_of_mass /= nres_counted;

	// calculate RG based on distance from center of mass
	Real rg_score = 0;
	for ( Size i = 1; i <= nres; ++i ) {
		// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
		if ( pose.residue(i).has_variant_type( core::chemical::REPLONLY ) ) continue;
		if ( pose.residue(i).aa() == core::chemical::aa_vrt ) continue;

		Vector const v( pose.residue(i).nbr_atom_xyz() );
		rg_score += v.distance_squared( center_of_mass );
	}

	// This definition of rg differs with the conventional definition which
	// divides by nres and not by nres-1.  For the sake of matching r++, it's
	// being left at nres-1 for now, but is a candidate for change in the near
	// future.
	rg_score /= (nres_counted - 1);

	return sqrt( rg_score );

}

core::Real
RG_Energy_Fast::calculate_rg_score(
	core::pose::Pose const & pose,
	utility::vector1< bool > const & relevant_residues) const
{

	Size const nres( pose.total_residue() );
	Size nres_counted=0;

	///////////////////////////////////////
	//
	// RG SCORE

	// calculate center of mass
	Vector center_of_mass( 0, 0, 0 );
	for ( Size i = 1; i <= nres; ++i ) {
		if (!relevant_residues[i]) continue;

		// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
		if ( pose.residue(i).has_variant_type( core::chemical::REPLONLY ) ) continue;
		if ( pose.residue(i).aa() == core::chemical::aa_vrt ) continue;

		Vector const v( pose.residue(i).nbr_atom_xyz() );
		center_of_mass += v;
		nres_counted++;
	}
	center_of_mass /= nres_counted;

	// calculate RG based on distance from center of mass
	Real rg_score = 0;
	for ( Size i = 1; i <= nres; ++i ) {
		if (!relevant_residues[i]) continue;

		// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
		if ( pose.residue(i).has_variant_type( core::chemical::REPLONLY ) ) continue;
		if ( pose.residue(i).aa() == core::chemical::aa_vrt ) continue;

		Vector const v( pose.residue(i).nbr_atom_xyz() );
		rg_score += v.distance_squared( center_of_mass );
	}

	// This definition of rg differs with the conventional definition which
	// divides by nres and not by nres-1.  For the sake of matching r++, it's
	// being left at nres-1 for now, but is a candidate for change in the near
	// future.
	rg_score /= (nres_counted - 1);

	return sqrt( rg_score );

}


void
RG_Energy_Fast::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & sf) const {
	RG_MinData &mindata = nonconst_mindata_from_pose( pose );

	// calculate center of mass
	Size const nres( pose.total_residue() );
	mindata.nres_scored = 0;
	mindata.com = Vector( 0, 0, 0 );
	for ( Size i = 1; i <= nres; ++i ) {
		// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
		if ( pose.residue(i).has_variant_type( core::chemical::REPLONLY ) ) continue;
		if ( pose.residue(i).aa() == core::chemical::aa_vrt ) continue;

		Vector const v( pose.residue(i).nbr_atom_xyz() );
		mindata.com += v;
		mindata.nres_scored++;
	}
	mindata.com /= mindata.nres_scored;

	// calulate score
	mindata.rg = 0;
	for ( Size i = 1; i <= nres; ++i ) {
		// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
		if ( pose.residue(i).has_variant_type( core::chemical::REPLONLY ) ) continue;
		if ( pose.residue(i).aa() == core::chemical::aa_vrt ) continue;

		Vector const v( pose.residue(i).nbr_atom_xyz() );
		mindata.rg += v.distance_squared( mindata.com );
	}
	mindata.rg = sqrt( mindata.rg / (mindata.nres_scored - 1) );
}


void
RG_Energy_Fast::eval_atom_derivative(
	id::AtomID const & id,
	pose::Pose const & pose,
	kinematics::DomainMap const &, // domain_map,
	ScoreFunction const & ,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const {
	RG_MinData const &mindata = mindata_from_pose( pose );

	int resid = id.rsd();
	int atmid = id.atomno();
	core::conformation::Residue const &rsd_i = pose.residue(resid);
	numeric::xyzVector<core::Real> X = pose.xyz(id);

	if (atmid != rsd_i.nbr_atom())
		return;

	numeric::xyzVector<core::Real> drg_dx = (X-mindata.com) / (mindata.rg*(mindata.nres_scored - 1));

	numeric::xyzVector<core::Real> atom_x = X;
	numeric::xyzVector<core::Real> const f2( drg_dx );
	numeric::xyzVector<core::Real> atom_y = -f2 + atom_x;
	Vector const f1( atom_x.cross( atom_y ) );

	F1 += weights[ rg ] * f1;
	F2 += weights[ rg ] * f2;
}


RG_MinData const & 
RG_Energy_Fast::mindata_from_pose( pose::Pose const & pose) const {
	using namespace core::pose::datacache;
	return *( static_cast< RG_MinData const * >( pose.data().get_const_ptr( CacheableDataType::RG_MINDATA )() ));

}

RG_MinData &
RG_Energy_Fast::nonconst_mindata_from_pose( pose::Pose & pose) const {
	if ( pose.data().has( core::pose::datacache::CacheableDataType::RG_MINDATA ) ) {
		return *( static_cast< RG_MinData* >( pose.data().get_ptr( core::pose::datacache::CacheableDataType::RG_MINDATA )() ));
	}
	// else
	RG_MinDataOP rgmindata = new RG_MinData;
	pose.data().set( core::pose::datacache::CacheableDataType::RG_MINDATA, rgmindata );
	return *rgmindata;
}



core::Size
RG_Energy_Fast::version() const
{
	return 1; // Initial versioning
}

} // methods
} // scoring
} // core
