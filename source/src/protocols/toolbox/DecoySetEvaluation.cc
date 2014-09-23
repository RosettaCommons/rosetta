// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author smart money says oliver


#include <protocols/toolbox/DecoySetEvaluation.hh>
#include <protocols/toolbox/DecoySetEvaluation.impl.hh>
#include <protocols/toolbox/superimpose.hh>

#include <core/types.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/conformation/Residue.hh>

#include <core/scoring/constraints/ConstraintSet.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/BoundConstraint.hh>
#include <core/scoring/constraints/NamedAtomPairConstraint.hh>
#include <core/scoring/constraints/LocalCoordinateConstraint.hh>
#include <core/io/silent/SilentFileData.hh>

#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/LoopsFileIO.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>
#include <numeric/model_quality/rms.hh>
#include <basic/options/option_macros.hh>

//// C++ headers
#include <string>
#include <iostream>
#include <cstdio>

#include <core/id/NamedStubID.hh>
#include <utility/vector1.hh>

OPT_1GRP_KEY( Real, dist_cst, max_dist )
OPT_1GRP_KEY( Real, dist_cst, sd )
OPT_1GRP_KEY( Real, dist_cst, ub_fact )
OPT_1GRP_KEY( Real, dist_cst, lb_fact )
OPT_1GRP_KEY( Real, dist_cst, min_ub )
OPT_1GRP_KEY( Real, dist_cst, grow_fact )
OPT_1GRP_KEY( Real, dist_cst, grow_fact_lb )
OPT_1GRP_KEY( Real, dist_cst, max_spread )
OPT_1GRP_KEY( Real, dist_cst, min_spread )
OPT_1GRP_KEY( File, dist_cst, excl_rigid )
OPT_1GRP_KEY( Boolean, dist_cst, median )
OPT_1GRP_KEY( Integer, dist_cst, min_seq_sep )
OPT_1GRP_KEY( File, dist_cst, dump )


bool protocols::toolbox::DecoySetEvaluation::options_registered_( false );

void protocols::toolbox::DecoySetEvaluation::register_options() {
	if ( !options_registered_ ) {
		options_registered_ = true;
		NEW_OPT( dist_cst::max_dist, "only use distances whose average distance is less than", 10 );
		NEW_OPT( dist_cst::sd, "curvature for BoundFunc after ub is reached",1 );
		NEW_OPT( dist_cst::ub_fact, "multiply this to sqrt of distance variance to get upper bound", 0.05 );
		NEW_OPT( dist_cst::lb_fact, "multiply this to sqrt of distance variance to get lower bound", 0.02 );
		NEW_OPT( dist_cst::grow_fact, "grow IN (positive) or OUTWARD( negative ) by grow_fact*(ub-lb)", -0.05 );
		NEW_OPT( dist_cst::max_spread, "omit constraints with more than X between ub and lb (suggest: 5-7A) ", 5.0 );
		NEW_OPT( dist_cst::min_spread, "enforce minimum flat bottom of X ", 2.0 );
		NEW_OPT( dist_cst::min_ub, "make upper bound at least " , 0.0);
		NEW_OPT( dist_cst::grow_fact_lb, "on top of grow_fact also grow_fact_lb only for the lb side", 0.00 );
		NEW_OPT( dist_cst::min_seq_sep, "only constraints if residues are x apart", 4 );
		NEW_OPT( dist_cst::median, "base bounds on sorted distances", false );
		NEW_OPT( dist_cst::dump, "put all dists into this file to play around in matlab or R", "coords.dat" );
		NEW_OPT( dist_cst::excl_rigid, "exclude constraints between residues that are in rigid core", "core.rigid" );
	}
}

namespace protocols {
namespace toolbox {

using namespace ObjexxFCL;

static thread_local basic::Tracer tr( "protocols.toolbox.DecoySetEvaluation", basic::t_info );

using namespace core;
using namespace numeric::model_quality; //for rms functions

typedef core::Real matrix[3][3];
typedef core::Real rvec[3];


DecoySetEvaluation::DecoySetEvaluation() : COM( 3 ),n_decoys_( 0 ), n_atoms_( 0 ), n_decoys_max_( 0 )
{}

DecoySetEvaluation::~DecoySetEvaluation() {}

void DecoySetEvaluation::reserve( core::Size n_resize ) {
	if ( n_resize == n_decoys_max_ ) return;
	n_decoys_max_ = n_resize;
	if ( n_atoms_ ) coords_.redimension( 3, n_atoms_, n_decoys_max_ );
}

void DecoySetEvaluation::prepare_push_back( core::Size nres ) {
	if ( n_decoys_ == 0 ) {
		n_atoms_ = nres;
		coords_.dimension( 3, n_atoms_, n_decoys_max_ );
		weights_.dimension( n_atoms_, 1.0 );
		ref_structure_.dimension( 3, n_atoms_ );
	} else {
		if ( n_atoms_ != nres ) {
			utility_exit_with_message( "can't insert poses with varying sizes into DecoySetEvaluation" );
		}
	}

	if ( n_decoys_ >= n_decoys_max_ ) {
		// we could also just resize... but I want the user to know.. because he should keep track which decoy is which...
		throw utility::excn::EXCN_RangeError( "you can't add any more decoys to DecoySetEvaluation ");
	}
	++n_decoys_;
}

void DecoySetEvaluation::push_back( core::pose::Pose& pose ) {
	//count residues with CA
	Size nres=0;
	for ( core::Size i = 1; i <= pose.total_residue(); i++ ) {
		if ( !pose.residue_type( i ).is_protein() ) break;
		++nres;
	}

	prepare_push_back( nres );

	// fill coords
	for ( core::Size i = 1; i <= nres; i++ ) {
		id::NamedAtomID idCA( "CA", i );
		PointPosition const& xyz = pose.xyz( idCA );
		for ( core::Size d = 1; d<=3; ++d ) {
			coords_( d, i, n_decoys_ ) = xyz[ d-1 ];
		}
	}

	if ( n_decoys_ == 1 ) {
		ref_pose_ = pose;
	}

}


void DecoySetEvaluation::push_back_CA_xyz( ObjexxFCL::FArray2_double const& xyz, core::Size nres ) {
	prepare_push_back( nres ); //increments n_decoys_

	// fill coords
	for ( core::Size i = 1; i <= nres; i++ ) {
		for ( core::Size d = 1; d<=3; ++d ) {
			coords_( d, i, n_decoys_ ) = xyz( d, i );
		}
	}
}

/*
void DecoySetEvaluation::push_back_CA_xyz( ObjexxFCL::FArray2D< core::Real > const& xyz, core::Size nres ) {
	prepare_push_back( nres );

	// fill coords
	for ( core::Size i = 1; i <= n_atoms_; i++ ) {
		for ( core::Size d = 1; d<=3; ++d ) {
			coords_( d, i, n_decoys_ ) = xyz( d, i );
		}
	}
}
*/

void DecoySetEvaluation::pop_back_CA_xyz( ){
	//erases last decoy of decoy_set
	std::cerr << "DSE erasing the last decoy in the set: " << n_decoys_ << std::endl;
	n_decoys_--;
	std::cerr << "DSE now has " << n_decoys_ << std::endl;
	std::cerr << "DSE redimensioning: old dimension: " << coords_.u1() << " " << coords_.u2() << " " << coords_.u3() << std::endl;
	coords_.redimension( 3, (n_atoms_*n_decoys_) , n_decoys_ );
	std::cerr << "DSE new dimensions : " << coords_.u1() << " " << coords_.u2() << " " << coords_.u3() << std::endl;
}

void DecoySetEvaluation::push_back_CA_xyz_from_silent_file( io::silent::SilentFileData const& sfd, bool store_energies ) {
	Size const n_new_decoys( sfd.size() );
	push_back_CA_xyz_from_silent_file( n_new_decoys, sfd.begin(), sfd.end(), store_energies );
}

void DecoySetEvaluation::set_n_atom( core::Size natoms ) {
	if ( n_atoms() ) tr.Warning << "Overriding n_atom in DecoySetEvaluation " << std::endl;
	n_atoms_ = natoms;
}

void DecoySetEvaluation::superimpose( Size icenter ) {
	FArray1D_double const weights( n_atoms_, 1.0 );
	superimpose( weights, icenter );
}

void DecoySetEvaluation::set_weights( ObjexxFCL::FArray1_double const& weights ) {
	for ( Size i=1; i<=n_atoms_; ++i ) weights_( i ) = weights( i );
}

core::Real  DecoySetEvaluation::rmsd( FArray1_double const& weights, FArray2_double& xx_ref, FArray2_double& xx ) const {
	FArray1D_double transvec( 3 );
	reset_x( n_atoms_, xx_ref, weights, transvec );
	reset_x( n_atoms_, xx, weights, transvec );
	Matrix R;
	//	FArray2P_double xx_ref( coords_( 1, 1, n ), 3, n_atoms_ ); //proxy array provides view to coords_
	fit_centered_coords( n_atoms_, weights, xx_ref, xx, R );
	Real rmsd( 0.0 );
	Real tot_weight( 0.0 );
	for ( Size n = 1; n <= n_atoms(); n++) {
		tot_weight += weights_( n );
		for ( Size d = 1; d<=3; ++d ) {
			rmsd += ( xx( d, n ) -  xx_ref( d, n ) ) * ( xx( d, n ) - xx_ref( d, n ) ) * weights_( n );
		}
	}
	return rmsd = sqrt( rmsd/tot_weight );
}

void DecoySetEvaluation::superimpose( FArray1_double const& weights, Size icenter ) {
	if ( n_decoys() == 0 ) return;
	FArray2P_double xx_ref( coords_( 1, 1, icenter ), 3, n_atoms_ ); //proxy array provides view to coords_
	tr.Debug << "superimpose with " << n_decoys() << " with " << icenter << " as reference structure " << std::endl;
	Size go_around( icenter + n_decoys_ - 1 );
	Size offset( 0 );
	for ( Size ni = icenter; ni <= go_around; ++ni ) {
		Size n( ni - offset );
		if ( n > n_decoys_ ) {
			offset = n_decoys_;
			n = 1;
		}
 		FArray2P_double xx( coords_( 1, 1, n ), 3, n_atoms_ ); //proxy array provides view to coords_
// 		FArray1D_double transvec( 3 );
// 		reset_x( n_atoms_, xx, weights, transvec );
		center_structure( n, weights );

		//fit
		if ( n != icenter ) {
			Matrix R;
			//	FArray2P_double xx2( coords_( 1, 1, n ), 3, n_atoms_ ); //proxy array provides view to coords_
			fit_centered_coords( n_atoms_, weights, xx_ref, xx, R );
		}
	} // for ni
}

void DecoySetEvaluation::center_structure( core::Size i ) {
	FArray1D_double const weights( n_atoms_, 1.0 );
	center_structure( i, weights );
}

void DecoySetEvaluation::center_structure( core::Size i, FArray1_double const& weights ) {
	FArray2P_double xx( coords_( 1, 1, i ), 3, n_atoms_ ); //proxy array provides view to coords_
	FArray1D_double transvec( 3 );
	reset_x( n_atoms_, xx, weights, transvec );
}

void DecoySetEvaluation::center_all( FArray1_double const& weights ) {
	for ( Size n = 1; n <= n_decoys_; ++n ) {
		center_structure( n, weights );
	}
}

//only makes sense after "superimpose"
void DecoySetEvaluation::compute_average_structure( FArray2_double& average_structure ) const {
	for ( Size n = 1; n<= n_decoys_; ++n ) {
		for ( Size d = 1; d<=3; ++d ) {
			for ( Size iatom = 1; iatom<=n_atoms_; ++iatom ) {
				average_structure( d, iatom ) += coords_( d, iatom, n )/n_decoys_;
			}
		}
	}
}

Size DecoySetEvaluation::find_closest_to_average( FArray2_double& average_structure ) const {
	Real best_dist( 100000000 );
	Size closest_structure( 1 );

	for ( Size n = 1; n<=n_decoys_; ++n ) {
		Real dist2( 0 );
		for ( Size iatom = 1; iatom <=n_atoms_; ++iatom ) {
			for ( Size d=1; d<=3; ++d ) {
				dist2 += (average_structure( d, iatom ) - coords_( d, iatom, n ))*(average_structure( d, iatom ) - coords_( d, iatom, n ));
			}
		}
		tr.Trace << "dist to average " << dist2 << std::endl;
		if ( dist2 < (best_dist - 0.01) ) { //make sure it is a significant improvement in distance
			best_dist = dist2;
			closest_structure = n;
		}
	}
	return closest_structure;
}

Real DecoySetEvaluation::rmsf( core::Size pos ) {
	//Real rms( 0 );
	Vector rmsd_x( 0.0 );
	Vector rmsd_av( 0.0 );
	Real invn( 1.0 /n_decoys() );
	for ( Size idecoy = 1; idecoy <= n_decoys(); ++idecoy ) {
		for ( Size d = 1; d<=3; ++d ) {
			Real dx = coords_( d, pos, idecoy );
			rmsd_x( d ) += dx * dx * invn;
			rmsd_av( d ) += dx * invn;
		}
	}
	Real rms( 0 );
	for ( Size d = 1; d<=3; ++d ) {
		rms += rmsd_x( d ) - rmsd_av( d )*rmsd_av( d );
	}
	return sqrt( rms );
}

void DecoySetEvaluation::rmsf( utility::vector1< Real >& result ) {
	result.clear();
	result.reserve( n_atoms_ );

	for ( Size pos = 1; pos <= n_atoms_; ++pos ) {
		result.push_back( rmsf( pos ) );
	}
}

void DecoySetEvaluation::rmsf( FArray1_double& result ) {
	for ( Size pos = 1; pos <= n_atoms_; ++pos ) {
		result( pos ) = rmsf( pos );
	}
}

//return icenter
Size DecoySetEvaluation::wRMSD( Real sigma2, Real tolerance, ObjexxFCL::FArray1_double& weights ) {
	if ( n_decoys() == 0 || n_atoms_ == 0 ) return 0;
	//FArray1D_double weights( n_atoms_, 1.0 );
	Real wsum_old = 100;
	Real wsum ( 1.0 );
	Real invn ( 1.0/n_atoms_ );
	Size ct ( 0 );
	Size i_center( 1 );
	Size max_iter( 100 );
	tr.Debug << "run wRMSD iterations with " << n_decoys() << " decoys of " << n_atoms() << " atoms " << std::endl;
	tr.Info << "wRMSD: iter  wsum   wRMSD icenter " << std::endl;
	while ( ( std::abs( wsum_old - wsum ) > tolerance ) && ( --max_iter > 0 ) ) {
		superimpose( weights, i_center );
		FArray2D_double average_structure( 3, n_atoms_, 0.0 );
		compute_average_structure( average_structure );
		i_center = find_closest_to_average( average_structure );
		rmsf( weights );
		wsum_old = wsum;
		wsum = 0.0;
		Real wMSD = 0.0;
		for ( Size i = 1; i <= n_atoms_; ++i ) {
			Real di2 = weights( i )*weights( i );
			weights( i ) = exp( - di2 / sigma2 );
			wsum += weights( i )*invn;
			wMSD += weights( i ) * di2;
		}
		tr.Info << "wRMSD: " << ++ct << " " << wsum << " " << sqrt( wMSD ) << " " << i_center << std::endl;
	}
	return i_center;
}

void DecoySetEvaluation::compute_distance_matrix( ObjexxFCL::FArray2D_double& dist) const {
	dist.dimension( n_decoys(), n_decoys(), 0.0 );
	int count = 0;
	Real invn( 1.0 / n_atoms() );

	Real sum_w( 0.0 );
	for ( Size n = 1; n <= n_atoms(); n++) {
		sum_w+=weights_( n );
	}
	invn = 1.0 / sum_w;

	for ( Size i = 1; i <= n_decoys(); i++ ) {

		FArray2P_double xx( coords_( 1, 1, i ), 3, n_atoms_ ); //proxy array provides view to coords_
		FArray1D_double transvec( 3 );
		reset_x( n_atoms_, xx, weights_, transvec );

		for ( Size j = 1; j<i; j++ ) {
			//already centered since j<i
			Matrix R;
			FArray2P_double xx2( coords_( 1, 1, j ), 3, n_atoms_ ); //proxy array provides view to coords_
			fit_centered_coords( n_atoms_, weights_, xx, xx2, R );
			Real rmsd( 0.0 );
			for ( Size n = 1; n <= n_atoms(); n++) {
				for ( Size d = 1; d<=3; ++d ) {
					rmsd += ( xx( d, n ) -  xx2( d, n ) ) * ( xx( d, n ) - xx2( d, n ) ) * invn * weights_( n );
				}
			}
			rmsd = sqrt( rmsd );
			dist( i, j ) = rmsd;
			dist( j, i ) = rmsd;
			//			tr.Trace << i << " " << j << " " << rmsd << " " << n_atoms_ << std::endl;
			// print some stats of progress
			count ++;
			if ( count % 50000 == 0 ) {
				Real const percent_done ( 200.0 * static_cast< Real > ( count ) / ( (n_decoys() - 1) * n_decoys() ) );
				tr.Info << count
								<< "/" << ( n_decoys() - 1 )*( n_decoys() )/2
								<< " ( " << ObjexxFCL::format::F(8,1,percent_done) << "% )"
								<< std::endl;
			}
		}
		dist ( i,i )=0.0;
	}
}


void DecoySetEvaluation::create_dist_constraints(
	 scoring::constraints::ConstraintSet& cst_set
) const {

using namespace basic::options;
using namespace basic::options::OptionKeys;
 runtime_assert( options_registered_ );

 if ( option[ dist_cst::median ]() ) {
	 create_dist_constraints_median( cst_set );
	 return;
 }
 ObjexxFCL::FArray2D_double ivm( n_atoms(), n_atoms(), 0.0 );
 Real const invn( 1.0/n_decoys() );

	for ( Size i = 1; i <= n_atoms(); i++ ) {
		for ( Size j = i+1; j<=n_atoms(); j++ ) {
			Real var( 0.0 );
			Real av( 0.0 );

			for ( Size n = 1; n<=n_decoys(); n++ ) {
				core::Vector xi( coords_( 1, i, n ), coords_( 2, i, n ), coords_( 3, i, n ));
				core::Vector xj( coords_( 1, j, n ), coords_( 2, j, n ), coords_( 3, j, n ));
				Real dist = xi.distance(xj);
				var += dist * dist * invn;
				av += dist * invn;
			}
			ivm( j, i ) = var - av*av;
			ivm( i, j ) = ivm( j, i);
			tr.Debug << i << " " << j << " " << ivm( i, j ) << "\n";
			if ( av < option[ dist_cst::max_dist ] && ( (int) j - (int) i) >= option[ dist_cst::min_seq_sep ]) {
				using namespace scoring::constraints;
				Real const lb( sqrt( ivm( i,j ) )*option[ dist_cst::lb_fact ] );
				Real const ub( sqrt( ivm( i,j ) )*option[ dist_cst::ub_fact ] );
				Real const sd( option[ dist_cst::sd ] );
				core::scoring::func::FuncOP fx( new BoundFunc( av - lb, av + ub, sd, "CM_DECOYS" ) );
				NamedAtomPairConstraint cst(
						id::NamedAtomID( "CA", i ),
						id::NamedAtomID( "CA", j ),
						fx
				);
				cst_set.add_constraint( cst.clone() );
				pose::Pose dummy_pose;
				cst.show_def( std::cout, dummy_pose );
			}
		}
	}
	tr.Debug << std::endl;
}


void DecoySetEvaluation::create_dist_constraints_median(
	scoring::constraints::ConstraintSet& cst_set
) const {
using namespace basic::options;
using namespace basic::options::OptionKeys;


	Real const invn( 1.0/n_decoys() );
	Size const Nlb( static_cast< Size > ( std::ceil( n_decoys()*option[ dist_cst::lb_fact ] )));
	Size const Nub( static_cast< Size > ( std::ceil( n_decoys()*option[ dist_cst::ub_fact ] )));
	Real const grow_fact( option[ dist_cst::grow_fact ] );
	Real const grow_fact_lb( option[ dist_cst::grow_fact_lb ] );
	loops::Loops rigid;
	if ( option[ dist_cst::excl_rigid ].user() ) {
		std::ifstream is( option[ dist_cst::excl_rigid ]().name().c_str() );

		if (!is.good()) {
			utility_exit_with_message( "[ERROR] Error opening RBSeg file '" + option[ dist_cst::excl_rigid ]().name() + "'" );
		}

		loops::PoseNumberedLoopFileReader reader;
		reader.hijack_loop_reading_code_set_loop_line_begin_token( "RIGID" );
		loops::SerializedLoopList loops = reader.read_pose_numbered_loops_file(
			is, option[ dist_cst::excl_rigid ](), false );
		rigid = loops::Loops( loops );
	}
	//	utility::vector1< core::Real > dist_sorted;
	for ( Size i = 1; i <= n_atoms(); i++ ) {
		for ( Size j = i+1; j<=n_atoms(); j++ ) {
			if ( rigid.size() && rigid.is_loop_residue( i ) && rigid.is_loop_residue( j ) ) continue;
			Real av( 0.0 );
			//	dist_sorted.clear();
			//			dist_sorted.reserve(n_decoys());
			Real* dist_sorted = new Real[ n_decoys() ];
			for ( Size n = 1; n<=n_decoys(); n++ ) {
				core::Vector xi( coords_( 1, i, n ), coords_( 2, i, n ), coords_( 3, i, n ));
				core::Vector xj( coords_( 1, j, n ), coords_( 2, j, n ), coords_( 3, j, n ));
				Real dist = xi.distance(xj);
				av += dist * invn;
				//				dist_sorted.push_back( dist );
				dist_sorted[ n - 1 ]=dist;
			}
			if ( av < option[ dist_cst::max_dist ] && ( (int) j - (int) i) >= option[ dist_cst::min_seq_sep ]) {
				std::sort( dist_sorted, dist_sorted + n_decoys() );
				tr.Debug << "Nlb: " << Nlb << " "<< dist_sorted[ Nlb ] << " Nub: " << Nub << " " << dist_sorted[ n_decoys()-Nub ] << std::endl;
				Real  lb( dist_sorted[ Nlb ] );
				Real  ub( dist_sorted[ n_decoys() - Nub ] );
				lb = lb + (ub-lb)*(grow_fact + grow_fact_lb );
				ub = ub - (ub-lb)*grow_fact;
				if ( ub - lb < option[ dist_cst::max_spread ]() ) {
					if ( ub < option[ dist_cst::min_ub ] ) ub = option[ dist_cst::min_ub ];
					if ( (ub - lb) < option[ dist_cst::min_spread ]() ) {
						ub = (ub+lb)*0.5 + 0.5*option[ dist_cst::min_spread ]();
						lb = (ub+lb)*0.5 - 0.5*option[ dist_cst::min_spread ]();
					}
					Real const sd( option[ dist_cst::sd ] );
					using namespace scoring::constraints;
					core::scoring::func::FuncOP fx( new BoundFunc( lb, ub, sd, "CM_DECOYS" ) );
					NamedAtomPairConstraint cst(
																			id::NamedAtomID( "CA", i ),
																			id::NamedAtomID( "CA", j ),
																			fx
					);
					cst_set.add_constraint( cst.clone() );
					pose::Pose dummy_pose;
					cst.show_def( std::cout, dummy_pose );
				}
			}
			delete[] dist_sorted;
		}
	}
	tr.Debug << std::endl;
}


void DecoySetEvaluation::create_xyz_constraints_median(
	 scoring::constraints::ConstraintSet& cst_set,
	 core::pose::Pose const& ref_pose,
	 Size root
) const {
using namespace basic::options;
using namespace basic::options::OptionKeys;

	Real const invn( 1.0/n_decoys() );
	//Size const Nlb( std::ceil( n_decoys()*option[ dist_cst::lb_fact ] ));
	Size const Nub( static_cast< Size > ( std::ceil( n_decoys()*option[ dist_cst::ub_fact ] )));
	if ( Nub >= n_decoys() ) {
		utility_exit_with_message( "set ub_fact to 0.05 such that e.g., 5% of decoys are allowed to violate upper bound." );
	}
	Real const grow_fact( option[ dist_cst::grow_fact ] );
	Real const min_ub( option[ dist_cst::min_ub ] );
	//Size const grow_fact_lb( option[ dist_cst::grow_fact_lb ] );

	//	utility::vector1< core::Real > dist_sorted;
	for ( Size pos = 1; pos <= n_atoms(); pos++ ) {
		core::Vector xyz_av( 0.0 );
		//compute average position
		for ( Size n = 1; n<=n_decoys(); n++ ) {
			for ( Size d = 1; d<=3; ++d ) {
				core::Real dx = coords_( d, pos, n );
				xyz_av( d ) += dx * invn;
			}
		}

		//use this to throw outliers away.
		typedef std::list< std::pair< core::Real, Size > > PointList;
		PointList sorted_points;
		for ( Size n = 1; n<=n_decoys(); n++ ) {
			core::Vector xi( coords_( 1, pos, n ), coords_( 2, pos, n ), coords_( 3, pos, n ));
			Real dist = xyz_av.distance(xi);
			sorted_points.push_back( std::make_pair( dist, n ) );
		}
		sorted_points.sort();

		//now compute new average only on the 70% points closest to the average xyz
		Size ct = 1;
		Size const N_find_center( static_cast< Size >( std::ceil( n_decoys()*0.7 )));
		tr.Debug << "use " << N_find_center << " points to compute average " << std::endl;
		xyz_av = core::Vector( 0.0 );
		for ( PointList::const_iterator it = sorted_points.begin(); ct<=N_find_center; ++ct, ++it ) {
			for ( Size d = 1; d<=3; ++d ) {
				Real dx = coords_( d, pos, it->second );
				xyz_av( d ) += dx / N_find_center;
			}
		}
		char buf[300];
		sprintf( buf, "%8.3f%8.3f%8.3f",xyz_av.x(),xyz_av.y(),xyz_av.z() );
		std::cout << "ATOM  " << ObjexxFCL::format::RJ( 5, pos) << "  CA " << " ALA A " << ObjexxFCL::format::RJ( 3, pos) << "    " << std::setw( 3) <<
			buf << std::endl;


		Real* dist_sorted = new Real[ n_decoys() ];
		for ( Size n = 1; n<=n_decoys(); n++ ) {
			core::Vector xi( coords_( 1, pos, n ), coords_( 2, pos, n ), coords_( 3, pos, n ));
			Real dist = xyz_av.distance(xi);
			tr.Debug << "push_dist " << dist << std::endl;
			dist_sorted[ n - 1 ]=dist;
		}

		std::sort( dist_sorted, dist_sorted + n_decoys() );

		Real  lb( 0.0 );
		Real  ub( dist_sorted[ n_decoys() - Nub ] );
		delete[] dist_sorted;
		tr.Debug << "ub-raw " << ub << std::endl;
		ub = ub - ub*grow_fact;
		if ( ub <= min_ub ) ub = min_ub;
		using namespace scoring::constraints;
		core::scoring::func::FuncOP fx( new BoundFunc( lb, ub, 1.0, "CM_DECOYS" ) );
		LocalCoordinateConstraint cst(
				id::AtomID( ref_pose.residue(pos).atom_index("CA"), pos),
					core::pose::named_stub_id_to_stub_id( id::NamedStubID( "N", "CA", "C", root ), ref_pose ),
				xyz_av,
				fx
		);
		cst_set.add_constraint( cst.clone() );
	}
}


// void DecoySetEvaluation::dump_coords( utility::io::ozstream& out) const {

// 	// dump file with i, j dist1 dist2  ..... dist N
// 	for ( Size i = 1; i <= n_atoms(); i++ ) {
// 		for ( Size j = i+3; j<=n_atoms(); j++ ) {
// 			out << i << " " << j << " ";
// 			for ( Size n = 1; n<=n_decoys(); n++ ) {
// 				core::Vector xi( coords_( 1, i, n ), coords_( 2, i, n ), coords_( 3, i, n ));
// 				core::Vector xj( coords_( 1, j, n ), coords_( 2, j, n ), coords_( 3, j, n ));
// 				Real dist = distance(xi, xj);
// 				out << dist << " ";
// 			}
// 			out << std::endl;
// 		}
// 	}
// }



} //evaluation
} //protocols
