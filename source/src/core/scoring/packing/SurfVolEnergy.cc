// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/packing/SurfVolEnergy.cc
/// @brief  Packing Score
/// @author Will Sheffler


//Unit headers
#include <core/scoring/packing/SurfVolEnergy.hh>
#include <core/scoring/packing/SurfVolEnergyCreator.hh>

// AUTO-REMOVED #include <core/scoring/packing/PoseBalls.hh>

#include <core/scoring/packing/surf_vol.hh>

//Package headers

//#include <core/scoring/ScoringManager.hh>
// AUTO-REMOVED #include <core/scoring/EnergyGraph.hh>
// AUTO-REMOVED #include <core/scoring/rms_util.hh>
#include <core/pose/Pose.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <core/pose/datacache/CacheableDataType.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
#include <core/id/CacheableAtomID_MapVector.hh>
// AUTO-REMOVED #include <basic/options/option.hh>

// AUTO-REMOVED #include <core/scoring/packing/compute_holes_score_res.hh>

// AUTO-REMOVED #include <basic/prof.hh>

//numeric headers
#include <numeric/numeric.functions.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
// AUTO-REMOVED #include <numeric/xyz.functions.hh>

//utility headers
// AUTO-REMOVED #include <utility/vector1.hh>
#include <utility/exit.hh>

//C++ headers
#include <iostream>

#include <core/chemical/AtomType.hh>
#include <core/pose/util.hh>
#include <core/scoring/EnergyMap.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>

namespace core {
namespace scoring {
namespace packing {


/// @details This must return a fresh instance of the SurfVolEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
SurfVolEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return new SurfVolEnergy;
}

ScoreTypes
SurfVolEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	// sts.push_back( dab_sasa );
	sts.push_back( dab_sev );
	return sts;
}


//////////////////////////////////////////////////////
//@brief
//////////////////////////////////////////////////////
SurfVolEnergy::SurfVolEnergy() :
	parent( methods::EnergyMethodCreatorOP( new SurfVolEnergyCreator ) )
{}


//////////////////////////////////////////////////////
//@brief
//////////////////////////////////////////////////////
void
SurfVolEnergy::finalize_total_energy(
  pose::Pose & pose,
  ScoreFunction const &,
  EnergyMap & totals
) const
{
	SurfVol sv = get_surf_vol(pose,1.4);
	totals[ dab_sasa ] = sv.tot_surf;
	totals[ dab_sev  ] = sv.tot_vol;
}

void
SurfVolEnergy::setup_for_derivatives(
	pose::Pose & pose,
	ScoreFunction const &
) const {
	std::cerr << "SurfVolEnergy::setup_for_derivatives" << std::endl;

	//using namespace core::pose::datacache::CacheableDataType;
	using namespace basic::datacache;
	using namespace id;
	using namespace numeric;
	using basic::datacache::DataCache_CacheableData;

	if( !pose.data().has( core::pose::datacache::CacheableDataType::DAB_SASA_POSE_INFO ) ) {
		pose.data().set( core::pose::datacache::CacheableDataType::DAB_SASA_POSE_INFO, DataCache_CacheableData::DataOP( new CacheableAtomID_MapVector ) );
		pose.data().set(  core::pose::datacache::CacheableDataType::DAB_SEV_POSE_INFO, DataCache_CacheableData::DataOP( new CacheableAtomID_MapVector ) );
	}

	CacheableDataOP dat1( pose.data().get_ptr( core::pose::datacache::CacheableDataType::DAB_SASA_POSE_INFO ) );
	CacheableAtomID_MapVectorOP cachemap1 = static_cast< CacheableAtomID_MapVector * >(dat1());
	AtomID_Map<xyzVector<Real> > & sasa_derivs(cachemap1->map());

	CacheableDataOP dat2( pose.data().get_ptr( core::pose::datacache::CacheableDataType::DAB_SASA_POSE_INFO ) );
	CacheableAtomID_MapVectorOP cachemap2 = static_cast< CacheableAtomID_MapVector * >(dat2());
	AtomID_Map<xyzVector<Real> > & sev_derivs(cachemap2->map());

	SurfVolDeriv svd = get_surf_vol_deriv( pose, 1.4 );

	core::pose::initialize_atomid_map_heavy_only(sasa_derivs,pose);
	core::pose::initialize_atomid_map_heavy_only( sev_derivs,pose);

	for( Size ir = 1; ir <= sasa_derivs.size(); ir++ ) {
		for( Size ia = 1; ia <= sasa_derivs.n_atom(ir); ia++ ) {
			AtomID const i(ia,ir);
			sasa_derivs[i] = svd.dsurf[i];
			 sev_derivs[i] = svd. dvol[i];
		}
	}

}


void
SurfVolEnergy::eval_atom_derivative(
	id::AtomID const & aid,
	pose::Pose const & pose,
	kinematics::DomainMap const &,
	ScoreFunction const &,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{
	// std::cerr << "SurfVolEnergy::eval_atom_derivative " << aid << std::endl;
	//using namespace core::pose::datacache::CacheableDataType;
	using namespace basic::datacache;
	using namespace id;
	using namespace numeric;

	// CacheableDataCOP dat( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::HOLES_POSE_INFO ) );
	// CacheableAtomID_MapVector const *cachemap = (CacheableAtomID_MapVector const *)dat();
	// AtomID_Map<xyzVector<Real> > const & derivs(cachemap->map());

	CacheableDataCOP dat1( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::DAB_SASA_POSE_INFO ) );
	CacheableAtomID_MapVectorCOP cachemap1 = static_cast< CacheableAtomID_MapVector const * > (dat1());
	AtomID_Map<xyzVector<Real> > const & sasa_derivs(cachemap1->map());

	CacheableDataCOP dat2( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::DAB_SASA_POSE_INFO ) );
	CacheableAtomID_MapVectorCOP cachemap2 = static_cast< CacheableAtomID_MapVector const * >(dat2());
	AtomID_Map<xyzVector<Real> > const & sev_derivs(cachemap2->map());

	if( aid.rsd() > sasa_derivs.n_residue() || aid.atomno() > sasa_derivs.n_atom(aid.rsd()) ) {
		return;
	}
	// std::cerr << "eval_atom_derivative " << aid << " " << derivs[aid].x() << " " << derivs[aid].y() << " " << derivs[aid].z() << std::endl;
	// F2 += weights * derivs[aid];
	// F1 += weights * derivs[aid].cross(xyzVector<Real>(1,0,0));

	{
  		numeric::xyzVector<core::Real> atom_x = pose.xyz(aid);
		numeric::xyzVector<core::Real> const f2( sasa_derivs[aid] );
		numeric::xyzVector<core::Real> const atom_y = atom_x - f2;   // a "fake" atom in the direcion of the gradient
		numeric::xyzVector<core::Real> const f1( atom_x.cross( atom_y ) );
		F1 += weights[ dab_sasa ] * f1;
		F2 += weights[ dab_sasa ] * f2;
	}
	{
  		numeric::xyzVector<core::Real> atom_x = pose.xyz(aid);
		numeric::xyzVector<core::Real> const f2( sev_derivs[aid] );
		numeric::xyzVector<core::Real> const atom_y = atom_x - f2;   // a "fake" atom in the direcion of the gradient
		numeric::xyzVector<core::Real> const f1( atom_x.cross( atom_y ) );
		F1 += weights[ dab_sev ] * f1;
		F2 += weights[ dab_sev ] * f2;
	}

}
core::Size
SurfVolEnergy::version() const
{
	return 1; // Initial versioning
}




} // packing
} // scoring
} // core
