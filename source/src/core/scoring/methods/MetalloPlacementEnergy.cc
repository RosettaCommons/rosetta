// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/MetalloPlacementEnergy.cc
/// @brief  Low-res placement score for metal sites
/// @author Will Sheffler


// Unit headers
#include <core/scoring/methods/MetalloPlacementEnergy.hh>

// Package headers
#include <core/scoring/ScoreFunction.hh>
//#include <core/scoring/ScoringManager.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <ObjexxFCL/FArray2D.hh>

#include <numeric/model_quality/rms.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/util.hh>




// Utility headers

#include <ObjexxFCL/format.hh>


namespace core {
namespace scoring {
namespace methods {

/// c-tor THIS CLASS REQUIRES AN EnergyMethodCreator to function properly
MetalloPlacementEnergy::MetalloPlacementEnergy()
{
	//add_score_type( metal_placement );
	collision_thresh2_ = 11.0*11.0;
	cb_cb_dis_thresh2_ = 3.0*3.0;
}



/// clone
EnergyMethodOP
MetalloPlacementEnergy::clone() const
{
	return new MetalloPlacementEnergy( *this );
}



/// all scoring happens here
void
MetalloPlacementEnergy::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const & /*scorefxn*/,
	EnergyMap & totals
) const
{
	using namespace core::conformation::symmetry;
	core::Size mmres,nsub;
	{
		assert( is_symmetric( pose ) );
		SymmetricConformation & symm_conf ( dynamic_cast<SymmetricConformation & > ( pose.conformation()) );
		mmres = symm_conf.Symmetry_Info().num_independent_residues();
		nsub  = symm_conf.Symmetry_Info().subunits();
	}
// std::cerr << "MetalloPlacementEnergy.cc:74 (" << mmres << " " << nsub << ")" << std::endl;


	totals[ metal_placement ] = 9e9;
	using core::id::AtomID;

	FArray2D<numeric::Real> sf4crd1(3,4);
	FArray2D<numeric::Real> sf4crd2(3,4);
	core::Real const R(3.94858713768);
	sf4crd1(1,1)= 0.000*R; sf4crd1(2,1)= 0.000*R; sf4crd1(3,1)= 1.000*R;
   sf4crd1(1,2)= 0.943*R; sf4crd1(2,2)= 0.000*R; sf4crd1(3,2)=-0.333*R;
   sf4crd1(1,3)=-0.471*R; sf4crd1(2,3)= 0.816*R; sf4crd1(3,3)=-0.333*R;
   sf4crd1(1,4)=-0.471*R; sf4crd1(2,4)=-0.816*R; sf4crd1(3,4)=-0.333*R;
   sf4crd2(1,1)= 0.943*R; sf4crd2(2,1)= 0.000*R; sf4crd2(3,1)=-0.333*R;
	sf4crd2(1,2)= 0.000*R; sf4crd2(2,2)= 0.000*R; sf4crd2(3,2)= 1.000*R;
   sf4crd2(1,3)=-0.471*R; sf4crd2(2,3)= 0.816*R; sf4crd2(3,3)=-0.333*R;
   sf4crd2(1,4)=-0.471*R; sf4crd2(2,4)=-0.816*R; sf4crd2(3,4)=-0.333*R;

	for(Size i =   1; i <= mmres     ; ++i ) {
		if( pose.residue_type(i).nheavyatoms() < 5 ) continue;
	for(Size j = i+1; j <= mmres*nsub; ++j ) {
		if( pose.residue_type(j).nheavyatoms() < 5 ) continue;
		if( i%mmres == j%mmres ) continue;
		if( cb_cb_dis_thresh2_ < pose.xyz(AtomID(5,i)).distance_squared(pose.xyz(AtomID(5,j))) ) continue;
	for(Size k = j+1; k <= mmres*nsub; ++k ) {
		if( pose.residue_type(k).nheavyatoms() < 5 ) continue;
		if( i%mmres == k%mmres ) continue;
		if( j%mmres == k%mmres ) continue;
		if( cb_cb_dis_thresh2_ < pose.xyz(AtomID(5,i)).distance_squared(pose.xyz(AtomID(5,k))) ) continue;
		if( cb_cb_dis_thresh2_ < pose.xyz(AtomID(5,j)).distance_squared(pose.xyz(AtomID(5,k))) ) continue;
	for(Size l = k+1; l <= mmres*nsub; ++l ) {
		if( pose.residue_type(k).nheavyatoms() < 5 ) continue;
		if( i%mmres == l%mmres ) continue;
		if( j%mmres == l%mmres ) continue;
		if( k%mmres == l%mmres ) continue;
		if( cb_cb_dis_thresh2_ < pose.xyz(AtomID(5,i)).distance_squared(pose.xyz(AtomID(5,l))) ) continue;
		if( cb_cb_dis_thresh2_ < pose.xyz(AtomID(5,j)).distance_squared(pose.xyz(AtomID(5,l))) ) continue;
		if( cb_cb_dis_thresh2_ < pose.xyz(AtomID(5,k)).distance_squared(pose.xyz(AtomID(5,l))) ) continue;
std::cerr << "MetalloPlacementEnergy.cc:112 (" << i << " " << j << " " << k << " " << l << ")" << std::endl;

		// check for space for SF4
		numeric::xyzVector<core::Real> const a = pose.xyz(AtomID(5,i));
		numeric::xyzVector<core::Real> const b = pose.xyz(AtomID(5,j));
		numeric::xyzVector<core::Real> const c = pose.xyz(AtomID(5,k));
		numeric::xyzVector<core::Real> const d = pose.xyz(AtomID(5,l));
		numeric::xyzVector<core::Real> const aa = pose.xyz(AtomID(2,i));
		numeric::xyzVector<core::Real> const ba = pose.xyz(AtomID(2,j));
		numeric::xyzVector<core::Real> const ca = pose.xyz(AtomID(2,k));
		numeric::xyzVector<core::Real> const da = pose.xyz(AtomID(2,l));
		numeric::xyzVector<core::Real> com = (a+b+c+d)/4.0;

		// check CA-CB pointing roughly right way
		if( (com-a).dot(a-aa) < 0.0 || (com-b).dot(b-ba) < 0.0 || (com-c).dot(c-ca) < 0.0 || (com-d).dot(d-da) < 0.0 ) {
			std::cerr << "MetalloPlacementEnergy.cc:127 ( fail dot prod )" << std::endl;
			continue;
		}

		// std::cout << "checking sf4 spot " << i << " " << j << " " << k << " " << l;

		bool collision = false;
		for(Size ir = 1; ir <= mmres*nsub; ++ir) {
			for(Size ia = 1; ia <= pose.residue_type(ir).nheavyatoms(); ++ia) {
				if( com.distance_squared(pose.xyz(AtomID(ia,ir))) < collision_thresh2_ ) {
					collision = true;
					break;
				}
			}
			if(collision) break;
		}
		if( collision ) {
			std::cerr << "MetalloPlacementEnergy.cc:144 (" << "fail collision" << ")" << std::endl;
			continue;
		}

		FArray2D<core::Real> crd(3,4);
		crd(1,1) = a.x(); crd(2,1) = a.y(); crd(3,1) = a.z();
		crd(1,2) = b.x(); crd(2,2) = b.y(); crd(3,2) = b.z();
		crd(1,3) = c.x(); crd(2,3) = c.y(); crd(3,3) = c.z();
		crd(1,4) = d.x(); crd(2,4) = d.y(); crd(3,4) = d.z();

		core::Real rms1 = numeric::model_quality::rms_wrapper(4,sf4crd1,crd);
		core::Real rms2 = numeric::model_quality::rms_wrapper(4,sf4crd2,crd);
		core::Real score = numeric::min(rms1,rms2);

		com.z(0);
		using numeric::min;
		using numeric::max;
		Size nt =         min((i-1)%mmres,min((j-1)%mmres,min((k-1)%mmres,(l-1)%mmres)))+1;
		Size ct = mmres - max((i-1)%mmres,max((j-1)%mmres,max((k-1)%mmres,(l-1)%mmres)))  ;
		score += com.length() / 2.0;
		Real tailpen = -0.5*( max(0.0,4.0-nt) + max(0.0,4.0-ct) );
		score += tailpen;

		if( score < totals[ metal_placement ] ) {
			std::cerr << "metal_placement " << score << " " << tailpen << " " << com.length() / 2.0 << std::endl;
			totals[ metal_placement ] = score;
		}

	} // l
	} // k
	} // j
	} // i

// std::cerr << "MetalloPlacementEnergy.cc:177 (DONE)" << std::endl;

}




}
}
}
