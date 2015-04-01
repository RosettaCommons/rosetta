// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/RG_LocalEnergy.cc
/// @brief  Radius of gyration for local region
/// @author TJ Brunette


// Unit headers
#include <core/scoring/methods/RG_LocalEnergy.hh>
#include <core/scoring/methods/RG_LocalEnergyCreator.hh>

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

#include <basic/options/keys/score.OptionKeys.gen.hh>
// C++


namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the RG_Energy_Fast class,
/// never an instance already in use
methods::EnergyMethodOP
RG_LocalEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new RG_LocalEnergy );
}

ScoreTypes
RG_LocalEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( rg_local );
	return sts;
}


/// was going to have this read in the blueprint. But core can't relly on something in protocols.
//RG_LocalEnergy::RG_LocalEnergy():RG_Energy_Fast()
RG_LocalEnergy::RG_LocalEnergy():
		RG_Energy_Fast( EnergyMethodCreatorOP( new RG_LocalEnergyCreator ) ) 
{
using namespace basic::options;
using namespace basic::options::OptionKeys;
firstRes_= 0;
lastRes_=0;
utility::vector1<core::Size> res;
if(!option[OptionKeys::score::rg_local_span].user())
		utility_exit_with_message("must set rg_local_res <first_res> <last_res> when using rg_local");
res = option[OptionKeys::score::rg_local_span]();
if(res.size() != 2)
		utility_exit_with_message("wrong format for rg_local_res. Please use <first_res> <last_res> when using rg_local");	
firstRes_ = res[1];
lastRes_ = res[2];
if(firstRes_ == lastRes_)
		utility_exit_with_message("<first_res> should not equal <last_res> when using rg_local");	
assert(firstRes_ > 0);
assert(lastRes_ > 0);
}


/// clone
EnergyMethodOP
RG_LocalEnergy::clone() const
{
	return EnergyMethodOP( new RG_LocalEnergy() );
}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

void RG_LocalEnergy::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap & totals
) const {
	using namespace conformation;

	PROF_START( basic::RG_LOCAL );

	totals[ rg_local ] = calculate_rg_score( pose );

	PROF_STOP( basic::RG_LOCAL );
} // finalize_total_energy

core::Real
RG_LocalEnergy::calculate_rg_score( core::pose::Pose const & pose ) const
{	
	utility::vector1< bool > relevant_residues;
	for (uint ii = 1; ii <= pose.total_residue(); ++ii){
		if ((ii >= firstRes_) && (ii <= lastRes_))
			relevant_residues.push_back(false);
		else
			relevant_residues.push_back(true);
	}
	return(RG_Energy_Fast::calculate_rg_score(pose,relevant_residues));
}

core::Real
RG_LocalEnergy::calculate_rg_score(
	core::pose::Pose const & pose,
	utility::vector1< bool > const & relevant_residues) const
{
	utility::vector1< bool > updated_relevant_residues;
	for (uint ii = 1; ii <= pose.total_residue(); ++ii){
		if (((ii >= firstRes_) && (ii <= lastRes_)) || (relevant_residues [ii] == false))
			updated_relevant_residues.push_back(false);
		else
			updated_relevant_residues.push_back(true);
	}
	return(RG_Energy_Fast::calculate_rg_score(pose,updated_relevant_residues));
}

void
RG_LocalEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const {
	RG_Local_MinData &mindata = nonconst_mindata_from_pose( pose );

	// calculate center of mass
	mindata.nres_scored = 0;
	mindata.com = Vector( 0, 0, 0 );
	for ( Size i = firstRes_; i <= lastRes_; ++i ) {
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
	for ( Size i = firstRes_; i <= lastRes_; ++i ) {
		// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
		if ( pose.residue(i).has_variant_type( core::chemical::REPLONLY ) ) continue;
		if ( pose.residue(i).aa() == core::chemical::aa_vrt ) continue;

		Vector const v( pose.residue(i).nbr_atom_xyz() );
		mindata.rg += v.distance_squared( mindata.com );
	}
	mindata.rg = sqrt( mindata.rg / (mindata.nres_scored - 1) );
}

//not inherited bc data is stored in RG_Local_MinData 
void
RG_LocalEnergy::eval_atom_derivative(
	id::AtomID const & id,
	pose::Pose const & pose,
	kinematics::DomainMap const &, // domain_map,
	ScoreFunction const & ,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const {
	RG_Local_MinData const &mindata = mindata_from_pose( pose );

	Size resid = id.rsd();
	Size atmid = id.atomno();
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


RG_Local_MinData const &
RG_LocalEnergy::mindata_from_pose( pose::Pose const & pose) const {
	using namespace core::pose::datacache;
	return *( utility::pointer::static_pointer_cast< core::scoring::methods::RG_Local_MinData const > ( pose.data().get_const_ptr( CacheableDataType::RG_LOCAL_MINDATA ) ));

}

RG_Local_MinData & RG_LocalEnergy::nonconst_mindata_from_pose( pose::Pose & pose) const {
	if ( pose.data().has( core::pose::datacache::CacheableDataType::RG_LOCAL_MINDATA ) ) {
		return *( utility::pointer::static_pointer_cast< core::scoring::methods::RG_Local_MinData > ( pose.data().get_ptr( core::pose::datacache::CacheableDataType::RG_LOCAL_MINDATA ) ));
	}
	// else
	RG_Local_MinDataOP rgmindata( new RG_Local_MinData );
	pose.data().set( core::pose::datacache::CacheableDataType::RG_LOCAL_MINDATA, rgmindata );
	return *rgmindata;
}


} // methods
} // scoring
} // core
