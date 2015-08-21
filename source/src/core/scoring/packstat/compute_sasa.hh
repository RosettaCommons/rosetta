// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/packstat/compute_sasa.hh
///
/// @brief
/// @author will sheffler


#ifndef INCLUDED_core_scoring_packstat_compute_sasa_hh
#define INCLUDED_core_scoring_packstat_compute_sasa_hh

#include <core/scoring/packstat/types.hh>
#include <core/scoring/packstat/CavityBall.hh>
#include <core/scoring/packstat/PackingScore.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <numeric/constants.hh>
#include <numeric/trig.functions.hh>

#include <ObjexxFCL/ubyte.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/Fmath.hh>


#include <core/id/AtomID_Map.fwd.hh>
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace packstat {

namespace old {

// this lookup table is used in sasa computation (also in void.cc)
extern short const bit_count[];

extern int const nbytes;
extern int const nphi;
extern int const ntheta;
extern int const nolp;
extern int const nori;
extern int const maskbits;

extern ObjexxFCL::FArray2D_int angles;
extern ObjexxFCL::FArray2D_ubyte masks;

///cj    Reads in sampling/SASA-angles.dat  sampling/SASA-masks.dat
void input_sasa_dats();

///cj    getting overlap from a to b (or i to j, see below)
///cj    this returns the degree of overlap between two atoms
///cj    adapted from erics code in area.c GetD2
///cj    returns value from 1 to 100
///cj    This calculation is based on the law of cosines,
///cj    see LeGrand and Merz, Journal of Computational
///cj    Chemistry 14(3):349-52 (1993).
///cj    Note that equation (4) is wrong, the denominator
///cj    should be 2r r   instead of 2r r
///cj                i iq              i q
inline
void
get_overlap(
	//Vector const & a,
	//FArray1_float const & a,
	PackstatReal const ra,
	//Vector const & b,
	PackstatReal const rb,
	PackstatReal const dist,
	int & olp
)
{
	PackstatReal epsilon,costh;

	//cj    min distance cutoff
	epsilon = 0.01;

	if ( dist < epsilon ) {
		//cj    atoms too close, causes round off error
		//cj    use this cutoff
		if ( ra < rb ) {
			olp = 100;
		} else {
			olp = 1;
		}
	} else if ( rb+dist <= ra ) {
		//cj    If atom a completely engulfs atom b, consider a to have
		//cj    no overlap due to atom b.
		olp = 1;
	} else if ( rb+dist <= ra ) {
		//cj    If atom a is completely engulfed by atom b, then turn it
		//cj    completely off (i.e. d2 = 99).
		olp = 100;
	} else {
		//cj    Otherwise, compute the amount of overlap using the law of cosines.
		//cj    "costh" is the angle of the cone of intersection that atom b
		//cj    imposes on atom a.  "ra" is the radius of atom a, and "rb" is
		//cj    the radius of atom b.  "sqrtd" is the actual distance between
		//cj    the a and b centers, while "dist" is the square of this distance.
		costh = (ra*ra+dist*dist-rb*rb)/(2*ra*dist);
		olp = static_cast< int >((1.0f-costh)*50)+1;
		if ( olp > 100 ) {
			olp = 100;
		} else if ( olp < 0 ) {
			//cj       We already hopefully accounted for this possibility by requiring that
			//cj       dist < epsilon, but in case not we don't want a potential bug to go
			//cj       unnoticed.
			// TRcs << "problem in calculating overlap between:" << std::endl;
			//    TRcs << "a  " <<
			//     F( 7, 3, a(1) ) << ' ' << F( 7, 3, a(2) ) << ' ' << F( 7, 3, a(3) ) <<
			//     std::endl;
			//    TRcs << "b  " <<
			//     F( 7, 3, b(1) ) << ' ' << F( 7, 3, b(2) ) << ' ' << F( 7, 3, b(3) ) <<
			//     std::endl;
			// TRcs << "ra=" << SS( ra ) << std::endl;
			// TRcs << "rb=" << SS( rb ) << std::endl;
			// TRcs << "dist=" << SS( dist ) << std::endl;
			// TRcs << "costh=" << SS( costh ) << std::endl;
			// TRcs << "Teminiating calculation" << std::endl;
			utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
		}
	}
}

///cj    gets the orientation of a to b (i to j, see below)
///cj    does this by calculating two angles, aphi and theta
inline
void
get_orientation(
	XYZ const & a,
	XYZ const & b,
	//FArray1_float const & a,
	//FArray1_float const & b,
	int & aphi,
	int & theta,
	PackstatReal dist
)
{
	using namespace numeric::constants::d;
	using numeric::sin_cos_range;

	// pb -- static can cause problems in multi-threading
	//apl allocate this once only
	//static FArray1D_float diff( 3 );
	XYZ diff( ( a - b ) / dist );

	//cj    figure the difference between a and b
	//apl - You've already computed the distance! reuse it.
	//diff(1) = (a(1)-b(1) ) / dist;
	//diff(2) = (a(2)-b(2) ) / dist;
	//diff(3) = (a(3)-b(3) ) / dist;

	//cj    now figure out polar values of aphi and theta
	//cj    Normalize the difference
	//cj    first get the length of the vector
	//vector_normalize(diff);

	//cj    figuring aphi

	PackstatReal p = std::acos( sin_cos_range( diff(3) ) );

	p *= nphi / pi_2;
	aphi = static_cast< int >( p );
	++aphi; // for fortran goes from 1 to n
	if ( aphi > nphi ) aphi = 1;

	//cj    figuring theta
	PackstatReal t = std::atan2(diff(2),diff(1));
	t *= ntheta / pi_2;
	theta = static_cast< int >( t );
	++theta; // for fortran goes from 1 to n
	if ( theta < 1 ) {
		theta += ntheta;
	} else if ( theta > ntheta ) {
		theta = 1;
	}

}

} // end namespace old


using core::pose::Pose;
using core::Real;
using utility::vector1;

struct SasaResult : public utility::pointer::ReferenceCount {
	SasaResult(std::size_t Nprobes, std::size_t Nspheres);
	ObjexxFCL::FArray2D<PackstatReal> sphere_sasa;
	CavBalls cavballs;
	// utility::vector1<numeric::xyzVector<PackstatReal> > sasa_centers;
};

struct SasaOptions : public utility::pointer::ReferenceCount {
	SasaOptions();
	Floats probe_radii;
	int   prune_max_iters;
	int   prune_max_delta;
	Floats prune_cavity_burial_probe_radii;
	PackstatReal min_hole_radius;
	int   num_cav_ball_layers;
	PackstatReal frac_cav_ball_required_exposed;
	PackstatReal area_cav_ball_required_exposed;
	PackstatReal min_cav_ball_radius;
	size_t num_surrounding_sasa_bins;
	size_t surrounding_sasa_smoothing_window;
	bool dont_compute_cav_balls;
	Real min_cluster_overlap;
	Real cluster_min_volume;
};

typedef utility::pointer::shared_ptr< SasaResult  > SasaResultOP;
typedef utility::pointer::shared_ptr< SasaOptions > SasaOptionsOP;

struct CavityBallCluster {
	numeric::xyzVector<core::Real> center;
	core::Real volume,surface_area,surface_accessibility;
	utility::vector1<CavityBall> cavballs;
};
struct OrderCBC {
	bool operator() ( CavityBallCluster const & a, CavityBallCluster const & b ) {
		return a.volume > b.volume;
	}
};


PosePackData
pose_to_pack_data(
	core::pose::Pose const & pose,
	int include_water = -1
);

SasaResultOP compute_sasa(
	Spheres & spheres,
	SasaOptions const & opts
);

CavBalls
prune_hidden_cavity_balls(
	CavBalls & cavballs,
	SasaOptions const & opts
);

CavBalls
prune_cavity_balls(
	Spheres & spheres,
	CavBalls & cavballs,
	SasaOptions const & opts
);

void
compute_cav_ball_volumes(
	CavBalls & cavballs,
	SasaOptions const & opts
);

void
compute_cav_ball_neighbor_count(
	Spheres & spheres,
	CavBalls & cavballs,
	PackstatReal dis
);

utility::vector1< CavityBallCluster >
compute_cav_ball_clusters(
	CavBalls & cavballs,
	SasaOptions const & opts
);


PackingScoreResDataOP
compute_surrounding_sasa(
	XYZ const & xyz,
	Spheres & spheres, // assumes spheres is sorted on x!
	SasaResultOP result,
	SasaOptions const & opts
);

CavBalls
select_cav_balls(
	CavBalls cavballs,
	PackstatReal spacing
);


core::Real
compute_packing_score(
	PosePackData & pd,
	core::Size oversample = 0
);

core::Real
compute_packing_score(
	Pose const & pose,
	core::Size oversample = 0
);

utility::vector1<core::Real>
compute_residue_packing_scores(
	PosePackData & pd,
	core::Size oversample = 0
);

utility::vector1<core::Real>
compute_residue_packing_scores(
	Pose const & pose,
	core::Size oversample = 0
);

core::Real
compute_residue_packing_score(
	PosePackData & pd,
	int const seqpos,
	core::Size oversample = 0
);

core::Real
compute_residue_packing_score(
	Pose const & pose,
	int const seqpos,
	core::Size oversample = 0
);

// std::pair<core::Real,core::Real>
// compute_packing_scores(
//  PosePackData & pd,
//  core::Size oversample = 0
// );
//
// std::pair<core::Real,core::Real>
// compute_packing_scores(
//  Pose const & pose,
//  core::Size oversample = 0
// );

utility::vector1<core::Real>
compute_atom_packing_scores(
	PosePackData & pd,
	core::Size oversample = 0
);

core::id::AtomID_Map<core::Real>
compute_atom_packing_scores(
	Pose const & pose,
	core::Size oversample = 0
);

//std::map<id::AtomID,Real>
//cavity_distance_constraint( Pose & pose, Size rsd );

template< class T >
PackstatReal compute_sasa_generic( utility::vector1< T > & S, PackstatReal probe, bool csa = true ) {
	using namespace core;
	using namespace numeric;
	using namespace utility;
	using namespace old;

	input_sasa_dats();

	int olp, aphi, theta, point, masknum;

	Size const Nspheres( S.size() );
	ObjexxFCL::FArray2D_ubyte atom_sasa_masks( old::nbytes, Nspheres, NULL );

	for ( std::size_t i = 1; i <= Nspheres; ++i ) {
		for ( std::size_t j = 1; j < i; ++j ) {
			PackstatReal const dist_sq = S[i].xyz.distance_squared(S[j].xyz);
			PackstatReal const dth = S[i].radius+S[j].radius+2*probe;
			if ( dist_sq > dth*dth ) continue;
			if ( dist_sq <= 0.0 ) continue;
			PackstatReal const dist = std::sqrt(dist_sq);
			PackstatReal const irad = S[i].radius + probe;
			PackstatReal const jrad = S[j].radius + probe;
			// account for j overlapping i:
			get_overlap(irad,jrad,dist,olp);
			get_orientation(S[i].xyz,S[j].xyz,aphi,theta,dist);
			point = angles(aphi,theta);
			masknum = point*100+olp;
			for ( int bb = 1, l = atom_sasa_masks.index(bb,i); bb <= nbytes; ++bb, ++l ) {
				atom_sasa_masks[ l ] = ObjexxFCL::bit::bit_or( atom_sasa_masks[ l ], masks(bb,masknum) );
			}
			// account for i overlapping j:
			get_overlap(jrad,irad,dist,olp);
			get_orientation(S[j].xyz,S[i].xyz,aphi,theta,dist);
			point = angles(aphi,theta);
			masknum = point*100+olp;
			for ( int bb = 1, l = atom_sasa_masks.index(bb,j); bb <= nbytes; ++bb, ++l ) {
				atom_sasa_masks[ l ] = ObjexxFCL::bit::bit_or( atom_sasa_masks[ l ], masks(bb,masknum) );
			}
		} // sphere j
	} // sphere i

	// compute sasas
	PackstatReal total = 0.0;
	PackstatReal fraction,total_sa,expose;
	for ( std::size_t is = 1; is <= Nspheres; ++is ) {
		PackstatReal const irad = S[is].radius;
		int ctr = 0;
		for ( std::size_t bb = 1, l = atom_sasa_masks.index(bb,is); (int)bb <= nbytes; ++bb, ++l ) {
			ctr += bit_count[atom_sasa_masks[ l ]];
		}
		fraction = static_cast< PackstatReal >( ctr ) / maskbits;
		if ( csa ) total_sa = 4.0 * numeric::constants::d::pi * ( irad * irad ); // 4*pi*r**2 -- this is M.S.A, not SASA
		else      total_sa = 4.0 * numeric::constants::d::pi * ( (irad+probe) * (irad+probe) ); // 4*pi*r**2 -- this is M.S.A, not SASA
		expose = ( 1.0f - fraction ) * total_sa;
		S[is].sasa = expose;
		total += expose;
	} // is

	return total;

}


void output_packstat_pdb( core::pose::Pose & pose, std::ostream & out );


} // namespace packstat
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_packstat_compute_sasa_HH
