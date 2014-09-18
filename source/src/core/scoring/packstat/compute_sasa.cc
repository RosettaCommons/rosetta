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



#include <cmath>
#include <basic/database/open.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <core/io/pdb/pose_io.hh>
#include <basic/options/keys/packstat.OptionKeys.gen.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/packstat/compute_sasa.hh>
#include <core/scoring/packstat/packing_score_params.hh>
#include <core/scoring/packstat/sasa_dot_locations.hh>
#include <core/scoring/packstat/SimplePDB.hh>
#include <core/scoring/packstat/util.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <numeric/numeric.functions.hh>
#include <numeric/NumericTraits.hh>
// AUTO-REMOVED #include <numeric/random/random.hh>
// AUTO-REMOVED #include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/ubyte.hh>
#include <set>
#include <sstream>
// AUTO-REMOVED #include <time.h>
// AUTO-REMOVED #include <utility/basic_sys_util.hh>
#include <utility/exit.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

#include <core/scoring/packstat/AtomRadiusMap.hh>
#include <utility/vector1.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.io.hh>
#include <numeric/random/random.fwd.hh>
#include <ObjexxFCL/format.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>
//Auto using namespaces
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS
//Auto using namespaces end

//Auto using namespaces
namespace ObjexxFCL { namespace format { } } using namespace ObjexxFCL::format; // AUTO USING NS
//Auto using namespaces end



namespace core {
namespace scoring {
namespace packstat {

static thread_local basic::Tracer TRcs( "protocols.packstat" );

typedef numeric::xyzMatrix<PackstatReal> Rot;
typedef std::pair< numeric::xyzMatrix<PackstatReal>, numeric::xyzMatrix<PackstatReal> > RotPair;

using utility::vector1;
using core::Real;
using core::Size;
using numeric::xyzVector;

namespace old {

	// this lookup table is used in sasa computation (also in void.cc)
	short const bit_count[] = { // lookup table for number of 1 bits in a ubyte
		0,1,1,2,1,2,2,3, 1,2,2,3,2,3,3,4,   1,2,2,3,2,3,3,4, 2,3,3,4,3,4,4,5,  // 0x 1x
		1,2,2,3,2,3,3,4, 2,3,3,4,3,4,4,5,   2,3,3,4,3,4,4,5, 3,4,4,5,4,5,5,6,  // 2x 3x
		1,2,2,3,2,3,3,4, 2,3,3,4,3,4,4,5,   2,3,3,4,3,4,4,5, 3,4,4,5,4,5,5,6,  // 4x 5x
		2,3,3,4,3,4,4,5, 3,4,4,5,4,5,5,6,   3,4,4,5,4,5,5,6, 4,5,5,6,5,6,6,7,  // 6x 7x
		1,2,2,3,2,3,3,4, 2,3,3,4,3,4,4,5,   2,3,3,4,3,4,4,5, 3,4,4,5,4,5,5,6,  // 8x 9x
		2,3,3,4,3,4,4,5, 3,4,4,5,4,5,5,6,   3,4,4,5,4,5,5,6, 4,5,5,6,5,6,6,7,  // Ax Bx
		2,3,3,4,3,4,4,5, 3,4,4,5,4,5,5,6,   3,4,4,5,4,5,5,6, 4,5,5,6,5,6,6,7,  // Cx Dx
		3,4,4,5,4,5,5,6, 4,5,5,6,5,6,6,7,   4,5,5,6,5,6,6,7, 5,6,6,7,6,7,7,8,  // Ex Fx
	};

	int const nbytes = { 21 };
	int const nphi = { 64 };
	int const ntheta = { 64 };
	int const nolp = { 100 };
	int const nori = { 162 };
	int const maskbits = { 162 };

	FArray2D_int angles( nphi, ntheta );
	FArray2D_ubyte masks( nbytes, nolp*nori );

	///cj    Reads in sampling/SASA-angles.dat  sampling/SASA-masks.dat
	void
	input_sasa_dats()
	{
		static bool done_init = false;
		if( done_init ) return;

		FArray1D_short tmp( nbytes );
		static bool init = { false };

		if ( init ) return;
		init = true;

	//cj    inputting the masks they are 21 ubytes long, 162x100 (see header)
	//cj    expects file to be complete
		utility::io::izstream masks_stream( basic::database::full_name("sampling/SASA-masks.dat" ) );

		for ( int i = 1; i <= nolp*nori; ++i ) {
			for ( int j = 1; j <= nbytes; ++j ) {
				masks_stream >> tmp(j);
			} masks_stream >> skip;

			for ( int j = 1; j <= nbytes; ++j ) {
				masks(j,i) = static_cast< ubyte>(tmp(j));
			}
		}
		masks_stream.close();

	//cj    inputting the angle lookup for the mask, need to add a 1 to each number
		utility::io::izstream angles_stream( basic::database::full_name( "sampling/SASA-angles.dat" ) );

	//cj    2D array is aphi by theta
		for ( int i = 1; i <= nphi; ++i ) {
			for ( int j = 1; j <= ntheta; ++j ) {
				angles_stream >> angles(i,j);
			} angles_stream >> skip;
	//cj       for ( j = 1; j <= ntheta; ++j ) {
	//cj          ++angles(i,j);
	//cj       }
		}
		angles_stream.close();

		done_init = true;
	}


} // end namespace old


SasaResult::SasaResult(size_t Nprobes, size_t Nspheres) :
	sphere_sasa( Nspheres, Nprobes, 0.0 )//,
	// sasa_centers( Nspheres, numeric::xyzVector<PackstatReal>(-12345.0f,-12345.0f,-12345.0f) )
{}

SasaOptions::SasaOptions() :
probe_radii(),
prune_max_iters(30),
prune_max_delta(5),
prune_cavity_burial_probe_radii(),
min_hole_radius(0.5),
num_cav_ball_layers(0),
frac_cav_ball_required_exposed(0.0f),
area_cav_ball_required_exposed(0.0f),
min_cav_ball_radius(0.5f),
num_surrounding_sasa_bins(7),
surrounding_sasa_smoothing_window(1),
dont_compute_cav_balls( false ),
min_cluster_overlap( 0.1 ),
cluster_min_volume(10)
{
}

RotPair rand_rot() {
	using namespace numeric;
	using namespace numeric::random;
	xyzVector<PackstatReal> axis(uniform(),uniform(),uniform());
	while( axis.length() > 1 ) axis = xyzVector<PackstatReal>(uniform(),uniform(),uniform());
	PackstatReal mag = uniform() * 2 * numeric::NumericTraits<Real>::pi();
	return RotPair( rotation_matrix<PackstatReal>( axis, mag ), rotation_matrix<PackstatReal>( axis, -mag ) );
}

struct OrderSphereOnX {
  bool operator()( Sphere const & a, Sphere const & b ) {
		return a.xyz.x() < b.xyz.x();
  }
};
struct OrderSphereOnID {
  bool operator()( Sphere const & a, Sphere const & b ) {
		return a.id < b.id;
  }
};
struct OrderCavBallOnX {
  bool operator()( CavityBall const & a, CavityBall const & b ) {
		return a.xyz().x() < b.xyz().x();
  }
};
struct OrderCavBallOnRmE {
  bool operator()( CavityBall const & a, CavityBall const & b ) {
		return (a.radius() - a.exposed_radius) > (b.radius() - b.exposed_radius);
  }
};
struct OrderCavBallOnR {
  bool operator()( CavityBall const & a, CavityBall const & b ) {
		return a.radius() > b.radius();
  }
};
struct OrderCavBallOnRmAnb {
  bool operator()( CavityBall const & a, CavityBall const & b ) {
		return ( (a.radius() - ((PackstatReal)a.anb)/50.0) > (b.radius() - ((PackstatReal)b.anb)/50.0 ) );
  }
};


PosePackData
pose_to_pack_data( Pose const & pose, int include_water ) {
	using namespace std;
	using namespace core::io::pdb;
	ostringstream oss;
	bool tmp = basic::options::option[basic::options::OptionKeys::out::file::output_virtual]();
	basic::options::option[basic::options::OptionKeys::out::file::output_virtual].value(true);
	dump_pdb(pose,oss);
	basic::options::option[basic::options::OptionKeys::out::file::output_virtual].value(tmp);
	istringstream iss( oss.str() );
	AtomRadiusMap arm( include_water );
  	SimplePDB pdb;
  	iss >> pdb;
	PosePackData p;
	p.spheres = pdb.get_spheres(arm);
	p.centers = pdb.get_res_centers();
	p.labels = pdb.res_labels();
	return p;
}



SasaResultOP
compute_sasa(
	Spheres & spheres,
	SasaOptions const & opts
) {
	using namespace old;
	using namespace utility;
	using namespace numeric;

	input_sasa_dats();

	//SphereIter i = spheres.begin();
	// PackstatReal prev_x = i->xyz.x();
	// bool must_sort = false;
	// for( i = i+1; i != spheres.end(); ++i ) {
	// 	if( i->xyz.x() < prev_x ) must_sort = true;
	// 	prev_x = i->xyz.x();
	// }
	// if( must_sort ) TRcs << "compute_sasa.cc: sorting your spheres ( input arg has changed! )" << std::endl;
	for( Size i = 1; i <= spheres.size(); ++i ) spheres[i].id = i;
	/*if( must_sort ) */std::sort( spheres.begin(), spheres.end(), OrderSphereOnX() );

	SasaResultOP result_op( new SasaResult( opts.probe_radii.size(), spheres.size() ) );
	SasaResult & result( *result_op );
	size_t Nspheres = spheres.size();
	size_t Nprobes  = opts.probe_radii.size();

	int /*in,it,      jn,jt,*/      olp,      aphi,theta,      point,masknum;  //,aa,lig_int_count;
	PackstatReal fraction,total_sa,expose;//,dist_sq,total_sasa5,total_sasa14,
	PackstatReal const largest_probe( opts.probe_radii[1] ); // first probe is largest

	ObjexxFCL::FArray3D_ubyte atom_sasa_masks( old::nbytes, Nspheres, Nprobes, NULL );

	vector1< RotPair > rotations;

	for ( size_t pr_bin = 1; pr_bin <= Nprobes; ++pr_bin ) {
  	PackstatReal const probe_radius( opts.probe_radii[pr_bin] );
		//if( pr_bin%10 == 0 ) TRcs << std::sqrt( ((PackstatReal)pr_bin)/((PackstatReal)Nprobes) ) << " ";

		RotPair rots = rand_rot();
		rotations.push_back( rots );
		Rot rot = rots.first;
		Spheres rspheres;
		for( SphereCIter i = spheres.begin(); i != spheres.end(); ++i ) {
			rspheres.push_back( Sphere( rot * i->xyz, i->radius ) );
		}

		for( size_t i = 1; i <= Nspheres; ++i ) {
			// TRcs << ((PackstatReal)i)/((PackstatReal)Nspheres) << std::endl;
			XYZ   const xyz1( rspheres[i].xyz );
			PackstatReal const rad1( rspheres[i].radius );

			for( size_t j = 1; j < i; ++j ) {
				XYZ   const xyz2( rspheres[j].xyz );
				PackstatReal const rad2( rspheres[j].radius );

				PackstatReal const dist_sq = xyz1.distance_squared(xyz2);
				PackstatReal const dth = rad1+rad2+2*largest_probe;

				if ( dist_sq > dth*dth ) continue;
	         if ( dist_sq <= 0.0 ) continue;
	         PackstatReal const dist = std::sqrt(dist_sq);

				PackstatReal const dth2 = rad1+rad2+2*probe_radius;
				if ( dist > dth2 ) continue;

            PackstatReal const irad = rad1 + probe_radius;
            PackstatReal const jrad = rad2 + probe_radius;

				// TRcs << i << " " << j << " " << dist_sq << " " << rad1 << " " << rad2 << " " << std::endl;

	      // account for j overlapping i:
	         get_overlap(irad,jrad,dist,olp);
	         get_orientation(xyz1,xyz2,aphi,theta, dist);
	         point = angles(aphi,theta);
	         masknum = point*100+olp;
            for ( int bb = 1, l = atom_sasa_masks.index(bb,i,pr_bin); bb <= nbytes; ++bb, ++l ) {
               atom_sasa_masks[ l ] = bit::bit_or( atom_sasa_masks[ l ], masks(bb,masknum) );
	         }

            // account for i overlapping j:
            get_overlap(jrad,irad,dist,olp);
            get_orientation(xyz2,xyz1,aphi,theta, dist);
            point = angles(aphi,theta);
            masknum = point*100+olp;
            for ( int bb = 1, l = atom_sasa_masks.index(bb,j,pr_bin); bb <= nbytes; ++bb, ++l ) {
              atom_sasa_masks[ l ] = bit::bit_or( atom_sasa_masks[ l ], masks(bb,masknum) );
            }

			} // pr_bin
		} // sphere j
	} // sphere i
	// TRcs <<  std::endl;
	// TRcs << "SASA work/all  " << ((PackstatReal)work) / ((PackstatReal)(clock() - loop)) << std::endl;

	// compute sasas
  for( size_t pr_bin = 1; pr_bin <= Nprobes; ++pr_bin ) {
		for( size_t is = 1; is <= Nspheres; ++is ) {
			PackstatReal const irad = spheres[is].radius;
      int ctr = 0;
      for ( size_t bb = 1, l = atom_sasa_masks.index(bb,is,pr_bin); (int)bb <= nbytes; ++bb, ++l ) {
        ctr += bit_count[atom_sasa_masks[ l ]];
      }
      fraction = static_cast< PackstatReal >( ctr ) / maskbits;
      total_sa = 12.56637 * ( irad * irad ); // 4*pi*r**2 -- this is M.S.A, not SASA
      expose = ( 1.0f - fraction ) * total_sa;
			result.sphere_sasa(spheres[is].id,pr_bin) = expose;
    } // is
  } // pr_bin

	if( opts.dont_compute_cav_balls ) {
		std::sort( spheres.begin(), spheres.end(), OrderSphereOnID() );
		return result_op;
	}

	///////////////////////////////////////////////////////////////////////
	/// compute cavity balls
	///////////////////////////////////////////////////////////////////////
	XYZs const sasa_dots( get_sasa_dot_locations() );

	size_t ball_count = 0;

  for ( size_t is = 1; is <= Nspheres; ++is ) {
    // get largest accessible probe radius for atom ir,ia
		numeric::xyzVector<PackstatReal> center(0.0f,0.0f,0.0f), median;
		size_t ii;
    for ( ii = opts.probe_radii.size(); ii >=  1; ii--) {
      if ( result.sphere_sasa(is,ii) < 0.0001f ) {
        break;
      }
    }
    if( ii < 1 || ii == opts.probe_radii.size() ) continue; // not big enough hole

    size_t const largest_pr_bin = ii+1;
		PackstatReal const probe_radius = opts.probe_radii[largest_pr_bin];
		if( probe_radius < opts.min_hole_radius ) {
			center = -12345.0f;
			continue; // next sphere
		}

		for ( size_t pr_bin = largest_pr_bin;
					pr_bin <= std::min(largest_pr_bin+opts.num_cav_ball_layers,opts.probe_radii.size());
					++pr_bin )
		{
			Rot rev_rot = rotations[pr_bin].second;
			PackstatReal tmp_probe_radius = opts.probe_radii[pr_bin];
			if( tmp_probe_radius < opts.min_cav_ball_radius ) continue;

	    for ( int bb = 1, l = atom_sasa_masks.index(bb,is,pr_bin); bb <= nbytes; ++bb, ++l ) {
        ubyte const mask = atom_sasa_masks[ l ]; // atom_sasa_masks_(bb,ia,ir,pr_bin)
        int NSHORT = 8; if ( bb == 21 ) NSHORT = 2; // 162 dots...
        for ( short k = 0; k < NSHORT; ++k ) {
          if ( !( bit::bit_test( static_cast<short>(mask), k ) ) ) {
            const size_t dot = 8*(bb-1)+k+1;

            median = (rev_rot * sasa_dots[dot]) * (spheres[is].radius+tmp_probe_radius) + spheres[is].xyz;
            CavityBall h( ++ball_count, is, median, tmp_probe_radius );
            result.cavballs.push_back( h );

          }
        }
	    } // end bit tests
		}

  } // is

	// std::sort( spheres.begin(), spheres.end(), OrderSphereOnID() );
	std::sort( result.cavballs.begin(), result.cavballs.end(), OrderCavBallOnX() );

	return result_op;
}

void
compute_cav_ball_volumes(
	CavBalls & cavballs,
	SasaOptions const &/*opts*/
) {
	using namespace ObjexxFCL;
	using namespace old;

	FArray2D_ubyte prune_sasa_masks( old::nbytes, cavballs.size(), NULL );

	std::sort( cavballs.begin(), cavballs.end(), OrderCavBallOnX() ); // already sorted

	PackstatReal maxrad = 0;
	for( CavBallIter i = cavballs.begin(); i != cavballs.end(); ++i ) {
		if( maxrad < i->radius() ) {
			maxrad = i->radius();
		}
		i->area = 0.0f;
		i->vol  = 0.0f;
	}

	for( PackstatReal dpr = 0.0f; dpr <= maxrad; dpr += 0.1 ) {
		PackstatReal dth = 2*(maxrad-dpr);
		// TRcs << dpr << " ";

		for( size_t ib = 1; ib <= cavballs.size(); ib++ )	{
			CavityBall const & cb1( cavballs[ib] );
			if( cb1.radius() - dpr < 0.1 ) continue;

			for( size_t jb = ib+1; jb <= cavballs.size(); jb++ ) {
				CavityBall const & cb2( cavballs[jb] );
				if( cb2.radius() - dpr < 0.1 ) continue;
				if( cb2.xyz().x() - cb1.xyz().x() > dth ) break;

				int olp,aphi,theta,point,masknum;//,aa,lig_int_count;
				PackstatReal const dis_thresh = (cb1.radius_+cb2.radius_-2*dpr );
				PackstatReal const dis2 = cb1.xyz().distance_squared( cb2.xyz() );
				if( dis2 > dis_thresh*dis_thresh || dis2 == 0 ) continue;
				PackstatReal const dist = sqrt( dis2 );

				get_overlap( cb1.radius()-dpr, cb2.radius()-dpr , dist, olp );
				get_orientation( cb1.xyz(), cb2.xyz(), aphi, theta, dist );
				point = angles(aphi,theta);
				masknum = point*100+olp;
				for ( int bb = 1, l = prune_sasa_masks.index(bb,ib); bb <= nbytes; ++bb, ++l ) {
					prune_sasa_masks[ l ] = bit::bit_or( prune_sasa_masks[ l ], masks(bb,masknum) );
				}
				get_overlap( cb2.radius()-dpr, cb1.radius()-dpr, dist, olp );
				get_orientation( cb2.xyz(), cb1.xyz() , aphi, theta, dist );
				point = angles(aphi,theta);
				masknum = point*100+olp;
				for ( int bb = 1, l = prune_sasa_masks.index(bb,jb); bb <= nbytes; ++bb, ++l ) {
					prune_sasa_masks[ l ] = bit::bit_or( prune_sasa_masks[ l ], masks(bb,masknum) );
				}

			}
		} // cav ball

		for( size_t ib = 1; ib <= cavballs.size(); ib++ )	{
			if( cavballs[ib].radius() - dpr < 0.1 ) continue;
			int ctr = 0;
			for ( int bb = 1, l = prune_sasa_masks.index(bb,ib); bb <= nbytes; ++bb, ++l ) {
				ctr += bit_count[prune_sasa_masks[ l ]];
			}
			PackstatReal expose = ((PackstatReal)(old::maskbits-ctr)) / ((PackstatReal)maskbits) ;
			expose *= 12.56637 * ( cavballs[ib].radius() - dpr ) * ( cavballs[ib].radius() - dpr );
			cavballs[ib].vol += 0.1 * expose;
			if( dpr == 0.0 ) {
			 	cavballs[ib].area = expose;
			}
		}

	}
	// TRcs << std::endl;
}

CavBalls
prune_hidden_cavity_balls(
	CavBalls & cavballs,
	SasaOptions const & opts
) {
	using namespace ObjexxFCL;
	using namespace old;

	FArray2D_ubyte prune_sasa_masks( old::nbytes, cavballs.size(), NULL );

	std::sort( cavballs.begin(), cavballs.end(), OrderCavBallOnX() );

	PackstatReal maxrad = max_rad( cavballs );

	for( size_t ib = 1; ib <= cavballs.size(); ib++ )	{
		CavityBall const & cb1( cavballs[ib] );

		for( size_t jb = ib+1; jb <= cavballs.size(); jb++ ) {
			CavityBall const & cb2( cavballs[jb] );

			if( cb2.xyz().x() - cb1.xyz().x() > maxrad ) break;

			int olp,aphi,theta,point,masknum;//,aa,lig_int_count;
			PackstatReal const dis_thresh = (cb1.radius_+cb2.radius_ );
			PackstatReal const dis2 = cb1.xyz().distance_squared( cb2.xyz() );
			if( dis2 > dis_thresh*dis_thresh || dis2 == 0 ) continue;
			PackstatReal const dist = sqrt( dis2 );

			get_overlap( cb1.radius(), cb2.radius() , dist, olp );
			get_orientation( cb1.xyz(), cb2.xyz() , aphi, theta, dist );
			point = angles(aphi,theta);
			masknum = point*100+olp;
			for ( int bb = 1, l = prune_sasa_masks.index(bb,ib); bb <= nbytes; ++bb, ++l ) {
				prune_sasa_masks[ l ] = bit::bit_or( prune_sasa_masks[ l ], masks(bb,masknum) );
			}

			get_overlap( cb2.radius(), cb1.radius() , dist, olp );
			get_orientation( cb2.xyz(), cb1.xyz() , aphi, theta, dist );
			point = angles(aphi,theta);
			masknum = point*100+olp;
			for ( int bb = 1, l = prune_sasa_masks.index(bb,jb); bb <= nbytes; ++bb, ++l ) {
				prune_sasa_masks[ l ] = bit::bit_or( prune_sasa_masks[ l ], masks(bb,masknum) );
			}

		}
	} // cav ball

	// make new list with only buried cav balls
	vector1< CavityBall > not_hidden_cav_balls;

	for( size_t ii = 1; ii <= cavballs.size(); ii++ ) {
		int ctr = 0;
		for ( int bb = 1, l = prune_sasa_masks.index(bb,ii); bb <= nbytes; ++bb, ++l ) {
			ctr += bit_count[prune_sasa_masks[ l ]];
		}
		PackstatReal expose = (PackstatReal)(old::maskbits-ctr) / (PackstatReal)maskbits;
		expose *= 12.56637 * cavballs[ii].radius() * cavballs[ii].radius();
		PackstatReal dot_th = (PackstatReal)(old::maskbits) * (1.0-opts.frac_cav_ball_required_exposed);
		if( (PackstatReal)ctr < dot_th && expose > opts.area_cav_ball_required_exposed ) {
			not_hidden_cav_balls.push_back( cavballs[ii] );
		}
	}
	//TRcs << "pruned " << cavballs.size() - not_hidden_cav_balls.size() << " hidden balls" << std::endl;
	return not_hidden_cav_balls;
}

size_t
prune_1pass(
	Spheres & spheres,
	CavBalls & cavballs,
	PackstatReal pr
) {
	using namespace ObjexxFCL;
	using namespace old;

	FArray2D_ubyte prune_sasa_masks( old::nbytes, cavballs.size(), 0 );

	PackstatReal const max_cbrad = max_rad( cavballs );
	PackstatReal const max_dis_sphere  = max_rad( spheres ) + max_cbrad + 2*pr;
	PackstatReal const max_dis_cavball = 2*max_cbrad+2*pr;

	for( size_t ib = 1; ib <= cavballs.size(); ib++ )	{
		CavityBall const & cb( cavballs[ib] );
		if( cb.exposed_radius >= pr ) continue;

		size_t begin = search_x( spheres, cb.xyz().x()-max_dis_sphere );
		for( size_t js = begin; js <= spheres.size(); ++js ) {
			Sphere const & s( spheres[js] );
			if( s.xyz.x() - cb.xyz().x() > max_dis_sphere ) break;

			int olp,aphi,theta,point,masknum;//,aa,lig_int_count;
			PackstatReal const dis_thresh = (cb.radius_+s.radius+2*pr );
			PackstatReal const dis2 = cb.xyz().distance_squared( s.xyz );
			if( dis2 > dis_thresh*dis_thresh ) continue;
			PackstatReal const dist = sqrt( dis2 );

			get_overlap( cb.radius()+pr, s.radius+pr , dist, olp );
			get_orientation( cb.xyz(), s.xyz , aphi, theta, dist );
			point = angles(aphi,theta);
			masknum = point*100+olp;
			for ( int bb = 1, l = prune_sasa_masks.index(bb,ib); bb <= nbytes; ++bb, ++l ) {
				prune_sasa_masks[ l ] = bit::bit_or( prune_sasa_masks[ l ], masks(bb,masknum) );
			}
	  } // sphere

		for( size_t jb = ib+1; jb < cavballs.size(); jb++ ) {
			CavityBall const & cb2( cavballs[jb] );
			if( cb2.xyz().x() - cb.xyz().x() > max_dis_cavball ) break;
			if( cb2.exposed_radius >= pr ) continue;

			int olp,aphi,theta,point,masknum;//,aa,lig_int_count;
			PackstatReal const dis_thresh = (cb.radius_+cb2.radius_+2*pr );
			PackstatReal const dis2 = cb.xyz().distance_squared( cb2.xyz() );
			if( dis2 > dis_thresh*dis_thresh || dis2 == 0 ) continue;
			PackstatReal const dist = sqrt( dis2 );

			get_overlap( cb.radius()+pr, cb2.radius_+pr , dist, olp );
			get_orientation( cb.xyz(), cb2.xyz() , aphi, theta, dist );
			point = angles(aphi,theta);
			masknum = point*100+olp;
			for ( int bb = 1, l = prune_sasa_masks.index(bb,ib); bb <= nbytes; ++bb, ++l ) {
				prune_sasa_masks[ l ] = bit::bit_or( prune_sasa_masks[ l ], masks(bb,masknum) );
			}

		}

	} // cav ball

	// make new list with only buried cav balls
	size_t num_exposed = 0;
	vector1< CavityBall > buried_cav_balls;
	for( size_t ii = 1; ii <= cavballs.size(); ii++ ) {
		if( cavballs[ii].exposed_radius >= pr ) {
			++num_exposed;
			continue;
		}
		int ctr = 0;
		for ( int bb = 1, l = prune_sasa_masks.index(bb,ii); bb <= nbytes; ++bb, ++l ) {
			ctr += bit_count[prune_sasa_masks[ l ]];
		}
		if( ctr != old::maskbits ) { // not buried
			cavballs[ii].exposed_radius = pr;
			++num_exposed;
		}
	}
	return num_exposed;
}

CavBalls
prune_cavity_balls(
	Spheres & spheres,
	CavBalls & cavballs,
	SasaOptions const & opts
) {
	using namespace std;
	//TRcs << "about to prune cavity balls" << std::endl;
	assert( opts.prune_cavity_burial_probe_radii.size() >= 1 );

	PackstatReal largest_probe_radius = opts.prune_cavity_burial_probe_radii[1];
	size_t delta_th = opts.prune_max_delta;
	int    iter_th  = opts.prune_max_iters;

	{
		TRcs << "prune raduis " << largest_probe_radius << " ";
		int num_prune_steps(0);
		size_t last_num = 0;
		size_t num_exposed;
		while( true ) {
		 	num_exposed = prune_1pass( spheres, cavballs, largest_probe_radius );
			++num_prune_steps;
			TRcs << cavballs.size() - num_exposed << " ";
			if ( num_prune_steps >= iter_th || (num_exposed-last_num) <= delta_th ) break;
			last_num = num_exposed;
		}
		TRcs << std::endl;
	}

	// make new list with exposed removed
	CavBalls buried_cav_balls;
	int id = 0;
	for( CavBallIter i = cavballs.begin(); i != cavballs.end(); ++i ) {
		if( i->exposed_radius < largest_probe_radius ) {
			buried_cav_balls.push_back( *i );
			buried_cav_balls[ buried_cav_balls.size() ].id_ = ++id;
		}
	}

	// mark balls for rest of radii
	for( size_t i = 2; i <= opts.prune_cavity_burial_probe_radii.size(); ++i ) {
		PackstatReal const pr = opts.prune_cavity_burial_probe_radii[i];
		TRcs << "prune raduis " << pr << " ";
		int num_prune_steps(0);
		size_t last_num = 0;
		size_t num_exposed;
		while( true ) {
		 	num_exposed = prune_1pass( spheres, buried_cav_balls, pr );
			++num_prune_steps;
			TRcs << buried_cav_balls.size() - num_exposed << " ";
			if ( num_prune_steps >= iter_th || (num_exposed-last_num) <= delta_th ) break;
			last_num = num_exposed;
		}
		TRcs << std::endl;
	}

	return buried_cav_balls;
}

void
compute_cav_ball_neighbor_count(
	Spheres & spheres,
	CavBalls & cavballs,
	PackstatReal dis
) {
	PackstatReal dis2 = dis*dis;
	for( CavBallIter c = cavballs.begin(); c != cavballs.end(); ++c ) {
		c->anb = 0;
		for( SphereIter s = spheres.begin(); s != spheres.end(); ++s ) {
			if( s->xyz.distance_squared( c->xyz() ) <= dis2 ) {
				c->anb++;
			}
		}
	}
}

CavBalls select_cav_balls(
	CavBalls cavballs,
	PackstatReal spacing
) {
	std::sort( cavballs.begin(), cavballs.end(), OrderCavBallOnR() );
	CavBalls selected_cavballs;
	PackstatReal const dist_th = spacing*spacing;
	for( CavBallIter i = cavballs.begin(); i != cavballs.end(); ++i ) {
		bool ok_to_add = true;
		for( CavBallIter j = selected_cavballs.begin(); j != selected_cavballs.end(); ++j ) {
			if( i->xyz().distance_squared(j->xyz()) < dist_th ) ok_to_add = false;
		}
		if( ok_to_add ) selected_cavballs.push_back( *i );
	}
	return selected_cavballs;
}


Real overlap(CavityBall const & cb1, CavityBall const & cb2) {
	core::Real const d2 = cb1.xyz().distance_squared(cb2.xyz());
	if( d2 > 36.0 ) return 0.0;
	core::Real const d = std::sqrt(d2);
	// return cb1.radius() + cb2.radius() - d;
	core::Real r1 = cb1.radius();
	core::Real r2 = cb2.radius();
	if( d > r1+r2 ) return 0.0;
	core::Real alpha = r2*sin( acos( (d*d + r2*r2 - r1*r1) / ( 2*d*r2 ) ) );
	// std::cerr << r1 << " " << r2 << " " << d << " " << alpha << std::endl;
	return alpha*alpha*numeric::NumericTraits<Real>::pi();
}

vector1< CavityBallCluster >
compute_cav_ball_clusters( CavBalls & cavballs, SasaOptions const & opts ) {
	using namespace numeric;

	for( Size i = 1; i <= cavballs.size(); i++ ) {
		cavballs[i].cluster_ = i;
	}

	for( Size i = 1; i <= cavballs.size(); i++ ) {
		CavityBall & cb1( cavballs[i] );
	// assume that cb1 has the right cluster before each 2nd loop
		for( Size j = i+1; j <= cavballs.size(); j++ ) {
			CavityBall & cb2( cavballs[j] );
			if( abs( cb1.xyz().x() - cb2.xyz().x() ) > 6.0 ) continue;
			// core::Real d = cb1.xyz().distance(cb2.xyz());
			if( overlap(cb1,cb2) >= opts.min_cluster_overlap ) {
				int c = numeric::min(cb1.cluster_,cb2.cluster_);
				cb1.cluster_ = c;
				cb2.cluster_ = c;
			}
		}
	}

	while(true) {
		bool fail = false;
		for( Size i = 1; i <= cavballs.size(); i++ ) {
			CavityBall & cb1( cavballs[i] );
			for( Size j = i+1; j <= cavballs.size(); j++ ) {
				CavityBall & cb2( cavballs[j] );
				if( abs( cb1.xyz().x() - cb2.xyz().x() ) > 6.0 ) continue;
				// core::Real d = cb1.xyz().distance(cb2.xyz());
				if( cb2.cluster_ != cb1.cluster_ && overlap(cb1,cb2) >= opts.min_cluster_overlap ) {
					fail = true;
					int from = numeric::max(cb1.cluster_,cb2.cluster_);
					int to   = numeric::min(cb1.cluster_,cb2.cluster_);
					for( Size k = 1; k <= cavballs.size(); k++ ) {
						if(cavballs[k].cluster_==from) cavballs[k].cluster_ = to;
					}
				}
				if(fail) break;
			}
			if(fail) break;
		}
		if(!fail) break;
	}

	// now renumber the clusters starting from 1
	std::map<Size,Size> perm;
	Size count = 1;
	for( Size i = 1; i <= cavballs.size(); i++ ) {
		if( perm.find(cavballs[i].cluster_) == perm.end() ) {
			perm[cavballs[i].cluster_] = count;
			count++;
		}
	}
	for( Size i = 1; i <= cavballs.size(); i++ ) {
		cavballs[i].cluster_ = perm[cavballs[i].cluster_];
	}

	vector1< xyzVector<core::Real> > cluster_centers;
	vector1< core::Real > cluster_volume;
	vector1< core::Real > cluster_surf;
	vector1< core::Real > cluster_counts;
	vector1< core::Real > cluster_surface_accessibility;
	for( Size i = 1; i <= cavballs.size(); i++ ) {
		CavityBall const & cb( cavballs[i] );
		if( (Size)cb.cluster_ > cluster_centers.size() ) {
			cluster_centers.resize(cb.cluster_);
			cluster_volume .resize(cb.cluster_);
			cluster_surf   .resize(cb.cluster_);
			cluster_counts .resize(cb.cluster_);
			cluster_surface_accessibility.resize(cb.cluster_);
		}
	}

	for( Size i = 1; i <= cluster_centers.size(); i++ ) {
		cluster_centers[i] = numeric::xyzVector<core::Real>(0.0,0.0,0.0);
		cluster_volume[i] = 0.0;
		cluster_counts[i] = 0.0;
		cluster_surf  [i] = 0.0;
		cluster_surface_accessibility[i] = 0.0;
	}

	for( Size i = 1; i <= cavballs.size(); i++ ) {
		CavityBall const & cb( cavballs[i] );
		if( (Size)cb.cluster_ > cluster_centers.size() ) {
			cluster_centers.resize(cb.cluster_);
			cluster_volume .resize(cb.cluster_);
			cluster_surf   .resize(cb.cluster_);
			cluster_counts .resize(cb.cluster_);
		}
		cluster_centers[cb.cluster_] += cb.xyz();
		cluster_surf   [cb.cluster_] += cb.area;
		cluster_volume [cb.cluster_] += cb.vol;
		cluster_counts [cb.cluster_] += 1.0;
		cluster_surface_accessibility[cb.cluster_] = numeric::max(cb.exposed_radius,cluster_surface_accessibility[cb.cluster_]);
		// std::cerr << cb.xyz().z() << " " << cluster_centers[cb.cluster_].z() << std::endl;
	}
	for( Size i = 1; i <= cluster_counts.size(); i++ ) {
		cluster_centers[i] /= cluster_counts[i];
		// std::cerr << cluster_counts[i] << " " << cluster_centers[i] << std::endl;
	}

	vector1< CavityBallCluster > result;
	for( Size i = 1; i <= cluster_counts.size(); i++ ) {
		// TRcs << "Cavity: " << i
		//   			  << " volume " << cluster_volume[i]
		//   			  << " surf "   << cluster_surf[i]
		//   			  << " Nballs " << cluster_counts[i]
		//   			  << " center " << cluster_centers[i].x()
		//   			  << ", "       << cluster_centers[i].y()
		//   			  << ", "       << cluster_centers[i].z()
		//   			  << std::endl;
  		if( cluster_volume[i] > opts.cluster_min_volume ) {
			CavityBallCluster cbc;
			cbc.volume = cluster_volume[i];
			cbc.surface_area = cluster_surf[i];
			cbc.center = cluster_centers[i];
			cbc.surface_accessibility = cluster_surface_accessibility[i];
			for( Size j = 1; j <= cavballs.size(); j++ ) {
				CavityBall const & cb( cavballs[j] );
				if( cb.cluster_ == (int)i ) {
					// std::cout << "cluster " << i << " add cav ball" << std::endl;
					cbc.cavballs.push_back(cb);
				}
			}
			result.push_back( cbc );
		}
	}

	std::sort( result.begin(), result.end(), OrderCBC() );

	for( Size i = 1; i <= result.size(); i++ ) {
		for( Size j = 1; j <= result[i].cavballs.size(); j++ ) {
			result[i].cavballs[j].cluster_ = i;
		}
	}

	return result;
}

PackingScoreResDataOP
compute_surrounding_sasa(
	XYZ const & xyz,
	Spheres & spheres, // assumes spheres is sorted on x!
	SasaResultOP result,
	SasaOptions const & opts
) {
	using namespace utility;
	using namespace numeric;


	// PackstatReal dist_th = (PackstatReal)opts.num_surrounding_sasa_bins;
	size_t Nprobes = opts.probe_radii.size() / opts.surrounding_sasa_smoothing_window;
	FArray2D<PackstatReal> tot_sasa( Nprobes, opts.num_surrounding_sasa_bins, 0.0f );
	size_t begin = 1;//search_x( spheres, xyz.x() - dist_th );
	// TRcs << "begin " << begin << std::endl;
	for( size_t is = begin; is <= spheres.size(); ++is ) {
		Sphere const & sphere( spheres[is] );
		//if( sphere.xyz.x() > xyz.x() + dist_th ) break;
		PackstatReal dist = xyz.distance( sphere.xyz );
		for( size_t id = 1; id <= opts.num_surrounding_sasa_bins; ++id ) {
			if( dist <= (PackstatReal)id ) {
				for( size_t pr = 1; pr <= Nprobes; ++pr ) {
					for( size_t io = 1; io <= opts.surrounding_sasa_smoothing_window; ++io ) {
					// std::cerr << "SASA " << is << " " << id << " " << pr << " " << result->sphere_sasa(is,pr) << std::endl;
						size_t sspr = (pr-1)*opts.surrounding_sasa_smoothing_window+io;
						// std::cerr << "pr " << pr << " sspr " << sspr << std::endl;
						tot_sasa(pr,id) += result->sphere_sasa(is,sspr);
					}
				}
				break;
			}
		}
	}
	if( opts.surrounding_sasa_smoothing_window > 1 ) {
		for( size_t id = 1; id <= opts.num_surrounding_sasa_bins; ++id ) {
			for( size_t pr = 1; pr <= Nprobes; ++pr ) {
				tot_sasa(pr,id) /= opts.surrounding_sasa_smoothing_window;
			}
		}
	}
	PackingScoreResDataOP psrdOP( new PackingScoreResData(opts.num_surrounding_sasa_bins,Nprobes-1) );
	for( size_t id = 1; id <= opts.num_surrounding_sasa_bins; ++id ) {
		for( size_t pr = 2; pr <= Nprobes; ++pr ) {
			psrdOP->msa(id,pr-1) = tot_sasa(pr,id) - tot_sasa(1,id); // stupid reverse indicies...
			// std::cout << psrdOP->msa(id,pr-1) << " ";
		}
	}
	// std::cout << std::endl;
	// std::exit(-1);
	return psrdOP;
}

core::Real
compute_packing_score(
	PosePackData & pd,
	core::Size oversample
) {

	assert( pd.spheres.size() > 0 );
	assert( pd.centers.size() > 0 );

	SasaOptions opts;
	opts.prune_max_iters = 0;
	opts.prune_max_delta = 0;
	opts.num_cav_ball_layers = 0;
	opts.frac_cav_ball_required_exposed = 0.00;
	opts.area_cav_ball_required_exposed = 0.00;
	opts.surrounding_sasa_smoothing_window = 1+2*oversample;
	opts.num_surrounding_sasa_bins = 7;
	for( core::Size ipr = 1; ipr <= 31; ++ipr ) {
		PackstatReal pr = 3.0 - ((double)(ipr-1))/10.0;
		PackstatReal ostep = 0.1 / (oversample*2.0+1.0);
		for( core::Size i = 1; i <= oversample; ++i )	opts.probe_radii.push_back( pr + i*ostep );
		opts.probe_radii.push_back( pr );
		for( core::Size i = 1; i <= oversample; ++i )	opts.probe_radii.push_back( pr - i*ostep );
	}
	// for( core::Size i = 1; i <= opts.probe_radii.size(); ++i )
	// 	std::cerr << "PR " << i << " " << opts.probe_radii[i] << std::endl;
	// opts.prune_cavity_burial_probe_radii.push_back( 1.6 );

	SasaResultOP result = compute_sasa( pd.spheres, opts );
	//std::cerr << "compute_sasa: num cav balls: " << result->cavballs.size() << std::endl;
	vector1< PackingScoreResDataCOP > psrds;
	for( core::Size i = 1; i <= pd.centers.size(); ++i ) {
		// std::cout << i << " ";
		psrds.push_back( compute_surrounding_sasa( pd.centers[i], pd.spheres, result, opts ) );
	}
	// std::exit(-1);
	PackingScore /*ps_respred(7,30,false),*/ ps_discrim(7,30,true);
	init_packing_score_discrim( ps_discrim );
	// init_packing_score_respred( ps_respred );

	// std::pair<core::Real,core::Real> score;
	// std::cerr << "DISCRIM" << std::endl;
	return ps_discrim.score( psrds );
	// std::cerr << "RESPRED" << std::endl;
	// score.second = ps_respred.score( psrds );

	// return score;

}

vector1<core::Real>
compute_residue_packing_scores(
  PosePackData & pd, core::Size oversample
) {
	return compute_atom_packing_scores( pd, oversample );
}

vector1<core::Real>
compute_atom_packing_scores(
	PosePackData & pd,
	core::Size oversample
) {
	assert( pd.spheres.size() > 0 );
	assert( pd.centers.size() > 0 );

	SasaOptions opts;
	opts.prune_max_iters = 0;
	opts.prune_max_delta = 0;
	opts.num_cav_ball_layers = 0;
	opts.frac_cav_ball_required_exposed = 0.00;
	opts.area_cav_ball_required_exposed = 0.00;
	opts.surrounding_sasa_smoothing_window = 1+2*oversample;
	opts.num_surrounding_sasa_bins = 7;
	for( core::Size ipr = 1; ipr <= 31; ++ipr ) {
		PackstatReal pr = 3.0 - ((double)(ipr-1))/10.0;
		PackstatReal ostep = 0.1 / (oversample*2.0+1.0);
		for( core::Size i = 1; i <= oversample; ++i )	opts.probe_radii.push_back( pr + i*ostep );
		opts.probe_radii.push_back( pr );
		for( core::Size i = 1; i <= oversample; ++i )	opts.probe_radii.push_back( pr - i*ostep );
	}

	PackingScore ps_discrim(7,30,true);
	init_packing_score_discrim( ps_discrim );
	vector1<Real> res_scores;
	SasaResultOP result = compute_sasa( pd.spheres, opts );
	for( core::Size i = 1; i <= pd.centers.size(); ++i ) {
		PackingScoreResDataCOP ss = compute_surrounding_sasa( pd.centers[i], pd.spheres, result, opts );
		res_scores.push_back( ps_discrim.score( ss ) );
	}

	return res_scores;
}

core::Real
compute_packing_score(
	Pose const & pose,
	core::Size oversample
) {
	PosePackData pd( pose_to_pack_data(pose) );
	return compute_packing_score( pd, oversample );
}


vector1<core::Real>
compute_residue_packing_scores(
	Pose const & pose,
	core::Size oversample
) {
	PosePackData pd( pose_to_pack_data(pose) );
	return compute_residue_packing_scores( pd, oversample );
}

core::Real
compute_residue_packing_score(
	Pose const & pose,
	int const seqpos,
	core::Size oversample
) {
	PosePackData pd( pose_to_pack_data(pose) );
	return compute_residue_packing_score( pd, seqpos, oversample );
}


core::Real
compute_residue_packing_score(
	PosePackData & pd,
	int const seqpos,
	core::Size oversample
) {
	assert( pd.spheres.size() > 0 );
	assert( pd.centers.size() > 0 );

	SasaOptions opts;
	opts.prune_max_iters = 0;
	opts.prune_max_delta = 0;
	opts.num_cav_ball_layers = 0;
	opts.frac_cav_ball_required_exposed = 0.00;
	opts.area_cav_ball_required_exposed = 0.00;
	opts.surrounding_sasa_smoothing_window = 1+2*oversample;
	opts.num_surrounding_sasa_bins = 7;
	for( core::Size ipr = 1; ipr <= 31; ++ipr ) {
		PackstatReal pr = 3.0 - ((double)(ipr-1))/10.0;
		PackstatReal ostep = 0.1 / (oversample*2.0+1.0);
		for( core::Size i = 1; i <= oversample; ++i )	opts.probe_radii.push_back( pr + i*ostep );
		opts.probe_radii.push_back( pr );
		for( core::Size i = 1; i <= oversample; ++i )	opts.probe_radii.push_back( pr - i*ostep );
	}

	PackingScore ps_discrim(7,30,true);
	init_packing_score_discrim( ps_discrim );

	SasaResultOP result = compute_sasa( pd.spheres, opts );

	PackingScoreResDataCOP ss = compute_surrounding_sasa( pd.centers[seqpos], pd.spheres, result, opts );

	for( int i = 1; i <= (int)ss->nrad(); ++i ) {
		for( int j = 1; j <= (int)ss->npr() ; ++j ) {
			std::cout << ss->msa(i,j) << " ";
		}
	}

	return ps_discrim.score( ss );
}

core::id::AtomID_Map<core::Real>
compute_atom_packing_scores(
	Pose const & pose,
	core::Size oversample
) {
	core::id::AtomID_Map<core::Real> atom_scores(0.0);
	core::pose::initialize_atomid_map( atom_scores, pose, 0.0 );
	PosePackData pd( pose_to_pack_data(pose) );
	vector1<Real> res_scores = compute_atom_packing_scores( pd, oversample );
	assert( res_scores.size() == pose.total_residue() );
	for( core::Size ir = 1; ir <= pose.total_residue(); ++ir ) {
		// numeric::xyzVector<Real> center(0,0,0);
		for( core::Size ia = 1; ia <= pose.residue(ir).nheavyatoms(); ++ia ) {
			// center += pose.residue(ir).xyz(ia);
			atom_scores[ core::id::AtomID(ia,ir) ] = res_scores[ir];
			// std::cerr << "RES SCORES " << ir << " " << ia << " " << res_scores[ir] << std::endl;
		}
		// center /= pose.residue(ir).nheavyatoms();
	}
	return atom_scores;
}

Real weight_func( Real d0, Real d ) {
	Real w = exp( -pow(d-d0,2.0)/10.0 ) + exp(-d/10.0) - 1.5/(d+1.0);
	if( w < 0 ) w = 0.0;
	return w;
}


vector1<std::map<id::AtomID,Real> >
cavity_distance_constraint(
   core::pose::Pose & pose,
   utility::vector1<core::Size> rois,
   core::pose::PoseOP native
) {

	using namespace numeric;

	PosePackData pd = pose_to_pack_data(pose);

	Spheres & spheres(pd.spheres);
	// vector1< xyzVector<PackstatReal> > & centers( pd.centers );

	SasaOptions opts;
	opts.prune_max_iters = 999;
	opts.prune_max_delta = 0;
	opts.num_cav_ball_layers = 10;
	opts.frac_cav_ball_required_exposed = 0.00;
	opts.area_cav_ball_required_exposed = 0.00;
	opts.surrounding_sasa_smoothing_window = 3;
	opts.min_cav_ball_radius = 0.7;
	opts.min_cluster_overlap = 0.1;   //how much overlap before considered to be in the same cluster
	opts.cluster_min_volume = 10.0;
	for( PackstatReal pr = 3.0; pr >= 0.4; pr -= 0.1 ) opts.probe_radii.push_back(pr);
	opts.prune_cavity_burial_probe_radii.push_back( 1.6 );
	// if( surface_accessibility ) {
	// 	for( PackstatReal pr = burial_radius-0.1; pr >= 0.1; pr -= 0.1 ) {
	// 		opts.prune_cavity_burial_probe_radii.push_back(pr);
	// 	}
	// }

	//TRcs << "compute MSAs" << std::endl;
	SasaResultOP sr = compute_sasa( spheres, opts );

	////////////////////////////////////////////////////////////////////////////////////////////////
	CavBalls cavballs = sr->cavballs;
	//TRps << "pruning hidden cav balls " << cavballs.size() << std::endl;
	cavballs = prune_hidden_cavity_balls( cavballs, opts );

	//TRps << "pruning exposed cav balls " << cavballs.size() << std::endl;
	cavballs = prune_cavity_balls( spheres, cavballs, opts );

	//TRps << "compute cav ball volumes	" << cavballs.size() << std::endl;
	compute_cav_ball_volumes( cavballs, opts );
	compute_cav_ball_neighbor_count( spheres,	cavballs, 10.0 );

	vector1< CavityBallCluster > clusters =
	compute_cav_ball_clusters( cavballs, opts );
	assert(clusters.size() > 0);

	vector1< std::map<id::AtomID,Real> > constraints_list;

	for( Size roi_i = 1; roi_i <= rois.size(); ++roi_i ) {
		Size roi = rois[roi_i];

		if( pose.residue(roi).nheavyatoms() <= 4 ) {
			utility_exit_with_message("residue of interest must have a CB!");
		}
		numeric::xyzVector<core::Real> roi_ca( pose.residue(roi).xyz(2) );
		core::Real dist_max = 0;
		numeric::xyzVector<core::Real> roi_max;
		numeric::xyzVector<core::Real> native_roi_max;
		for( Size ia = 5; ia <= pose.residue(roi).nheavyatoms(); ++ia ) {
			if( roi_ca.distance(pose.residue(roi).xyz(ia)) > dist_max ) {
				roi_max = pose.residue(roi).xyz(ia);
				if( native ) native_roi_max = native->residue(roi).xyz(ia);
				dist_max = roi_max.distance(roi_ca);
			}
		}
		Real const d0 = roi_max.distance(roi_ca);
		TRcs << "PACKSTAT_ROI " << roi << " d0 is " << d0 << std::endl;
		Real best_clust_score = 0.0;
		//Size best_clust = 123456789;
		numeric::xyzVector<core::Real> best_clust_wcen;
		for( Size i = 1; i <= clusters.size(); ++i ) {
			Real clust_score = 0.0, clust_wtot = 0.0;
			numeric::xyzVector<core::Real> clust_wcen(0,0,0);
			for( Size j = 1; j <= clusters[i].cavballs.size(); ++j ) {
				numeric::xyzVector<core::Real> xyz( clusters[i].cavballs[j].xyz() );
				Real const d = roi_ca.distance(xyz);
				Real const w = weight_func(d0,d);
				clust_score+= w*clusters[i].cavballs[j].vol*(clusters[i].cavballs[j].anb-50);
				clust_wtot += w*clusters[i].cavballs[j].vol*(clusters[i].cavballs[j].anb-50);
				clust_wcen += w*clusters[i].cavballs[j].vol*(clusters[i].cavballs[j].anb-50)*clusters[i].cavballs[j].xyz();
			}
			if( clust_score <= 0 ) continue;
			clust_wcen /= clust_wtot;

			TRcs << "PACKSTAT_ROI_CLUSTER " << roi << " clust " << i << " vol " << clusters[i].volume << " score " << clust_score;
			if( native ) {
				Real dnat = clust_wcen.distance(native_roi_max);
				TRcs << " dnat " << dnat;
			}
			TRcs << std::endl;

			if( clust_score > best_clust_score ) {
				best_clust_score = clust_score;
				best_clust_wcen = clust_wcen;
				//best_clust = i;  // set but never used ~Labonte
			}
		}

		std::set<id::AtomID> nbrs;
		for( Size ir = 1; ir <= pose.n_residue(); ++ir ) {
			for( Size ia = 1; ia <= pose.residue(ir).nheavyatoms(); ++ia ) {
				id::AtomID const aid(ia,ir);
				if( best_clust_wcen.distance_squared(pose.xyz(aid)) > 49.0 ) continue;
				// possibly some effort to pick only atoms lining cluster?
				nbrs.insert(aid);
			}
		}

		TRcs << "PACKSTAT_ROI " << roi << " cluster weighted center " << best_clust_wcen << std::endl;

		std::map<id::AtomID,Real> constraints;

		for( std::set<id::AtomID>::iterator i = nbrs.begin(); i != nbrs.end(); ++i ) {
			id::AtomID aid( *i );
			constraints[aid] = best_clust_wcen.distance(pose.xyz(aid));
			TRcs  << "PACKSTAT_ROI " << roi << " constraint: " << aid << " " << " " << std::endl;
		}

		constraints_list.push_back(constraints);

		/////////////////////////////// test output ///////////////////////////////////////////////////
		utility::io::ozstream out( ("PACKSTAT_ROI_"+string_of(roi)+".pdb").c_str() );

		pose.dump_pdb(out);

		for( Size i = 1; i <= clusters.size(); i++ ) {
			if( clusters[i].volume < 10 ) continue;
			for( Size j = 1; j <= clusters[i].cavballs.size(); j++ ) {
				out << clusters[i].cavballs[j].hetero_atom_line(i,i) << std::endl;
			}
		}

		out << "HETATM" + I( 5, 0 ) + "  C   CEN X"
				+ I( 4, 0 ) + "    "
				+ F( 8, 3, best_clust_wcen.x() ) + F( 8, 3, best_clust_wcen.y() ) + F( 8, 3, best_clust_wcen.z() )
				+ F( 6, 2, 0.0 ) + ' ' + F( 5, 2,2.0);


		out.close();

	}

	return constraints_list;
}



void output_packstat_pdb( core::pose::Pose & pose, std::ostream & out ) {
	using namespace core::scoring::packstat;
 	using namespace std;
	using namespace core;
	using namespace basic::options;
	using namespace ObjexxFCL::format;
	using namespace numeric;
	using namespace utility;

	bool surface_accessibility = option[ OptionKeys::packstat::surface_accessibility ]();
	core::Real burial_radius   = option[ OptionKeys::packstat::cavity_burial_probe_radius ]();

	PosePackData pd = pose_to_pack_data(pose);

	TRcs << "spheres len: " << pd.spheres.size() << std::endl;
	TRcs << "centers len: " << pd.centers.size() << std::endl;

	SasaOptions opts;
	opts.prune_max_iters = 999;
	opts.prune_max_delta = 0;
	opts.num_cav_ball_layers = 10;
	opts.frac_cav_ball_required_exposed = 0.00;
	opts.area_cav_ball_required_exposed = 0.00;
	opts.surrounding_sasa_smoothing_window = 3;
	opts.min_cav_ball_radius = option[ OptionKeys::packstat::min_cav_ball_radius ]();
	opts.min_cluster_overlap = option[ OptionKeys::packstat::min_cluster_overlap ]();
	opts.cluster_min_volume = option[ OptionKeys::packstat::cluster_min_volume ]();
	for( PackstatReal pr = 3.0; pr >  2.0; pr -= 0.2 ) opts.probe_radii.push_back(pr);
	for( PackstatReal pr = 2.0; pr >= 0.7; pr -= 0.1 ) opts.probe_radii.push_back(pr);
	opts.prune_cavity_burial_probe_radii.push_back( burial_radius );
	if( surface_accessibility ) {
		for( PackstatReal pr = burial_radius-0.1; pr >= 0.1; pr -= 0.1 ) {
			opts.prune_cavity_burial_probe_radii.push_back(pr);
		}
	}

	//TRps << "compute MSAs" << std::endl;
	SasaResultOP sr = compute_sasa( pd.spheres, opts );

	////////////////////////////////////////////////////////////////////////////////////////////////
	CavBalls cavballs = sr->cavballs;
	//TRps << "pruning hidden cav balls " << cavballs.size() << std::endl;
	cavballs = prune_hidden_cavity_balls( cavballs, opts );

	//TRps << "pruning exposed cav balls " << cavballs.size() << std::endl;
	cavballs = prune_cavity_balls( pd.spheres, cavballs, opts );

	//TRps << "compute cav ball volumes	" << cavballs.size() << std::endl;
	compute_cav_ball_volumes( cavballs, opts );

	vector1< CavityBallCluster > clusters = compute_cav_ball_clusters( cavballs, opts );

	for( Size i = 1; i <= clusters.size(); i++ ) {
		for( Size j = 1; j <= clusters[i].cavballs.size(); j++ ) {
			if( clusters[i].cavballs[j].radius() > 0.6 )
				out << clusters[i].cavballs[j].hetero_atom_line( pose.total_residue()+i, i, 0.0 ) << std::endl;
		}
	}

   // std::ofstream dbg("DEBUG.pdb");
   // pose.dump_pdb(dbg);
   // for( Size i = 1; i <= clusters.size(); i++ ) {
   //    for( Size j = 1; j <= clusters[i].cavballs.size(); j++ ) {
   //       dbg << clusters[i].cavballs[j].hetero_atom_line( pose.total_residue()+i, i ) << std::endl;
   //    }
   // }
   // dbg.close();

}


} // namespace packstat
} // namespace scoring
} // namespace core






