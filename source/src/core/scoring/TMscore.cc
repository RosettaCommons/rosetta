// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/TMscore.hh
/// @brief  reimplementation of TMscore sequence-based structure superposition
/// @author Hahnbeom Park
/// @detail original comment from TMscore.f code:
///
///*************************************************************************
///     This program is to compare two protein structures and identify the
///     best superposition that has the highest TM-score. Input structures
///     must be in the PDB format. By default, TM-score is normalized by
///     the second protein. Users can obtain a brief instruction by simply
///     running the program without arguments. For comments/suggestions,
///     please contact email: zhng@umich.edu.
///
///     Reference:
///     Yang Zhang, Jeffrey Skolnick, Proteins, 2004 57:702-10.
///
///     Permission to use, copy, modify, and distribute this program for
///     any purpose, with or without fee, is hereby granted, provided that
///     the notices on the head, the reference information, and this
///     copyright notice appear in all copies or substantial portions of
///     the Software. It is provided "as is" without express or implied
///     warranty.
///******************* Updating history ************************************

#include <core/types.hh>
#include <core/scoring/TMscore.hh>

#include <ObjexxFCL/FArray2D.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/model_quality/rms.hh>
#include <basic/Tracer.hh>

#include <cmath>
#include <cstdio>

namespace core {
namespace scoring {

static THREAD_LOCAL basic::Tracer TR( "core.scoring.TMscore" );

using namespace ObjexxFCL;

// Use FArray2D to be compatible with rms_util.cc
TMscore::TMscore( FArray2D< core::Real > const &p1 )
{
	nseq_ = p1.size()/3;
	convert_FArray2D_to_vector0( p1, xyza_ );
	set_default();
}

TMscore::~TMscore(){}

void
TMscore::set_default()
{
	d0_ = std::pow( 1.24*(nseq_-15), (1.0/3.0)-1.8 );
	if ( d0_<0.5 ) d0_ = 0.5;

	d0_search_ = d0_;
	if ( d0_search_ > 8 ) d0_search_ = 8;
	if ( d0_search_ < 4.5 ) d0_search_ = 4.5;
	TR.Debug << "d0_search set as: " << d0_search_ << std::endl;
}

void
TMscore::convert_FArray2D_to_vector0( FArray2D< core::Real > const & p1,
	utility::vector0< core::Vector > &xyz
)
{
	core::Size n( p1.size()/3 );
	xyz.resize( n , core::Vector(0,0,0));
	for ( core::Size i = 0; i < n; ++i ) {
		core::Vector v;
		v[0] = p1( 1, i+1 ); v[1] = p1( 2, i+1 ); v[2] = p1( 3, i+1 ); // Be careful of indexing
		xyz[i] = v;
	}
}

void
TMscore::get_ali_params()
{
	L_ini_.resize( 0 );

	// iA_ and iB_ are originally designed to map crd index to alignment index,
	// but becomes trivial once we get already aligned p1 & p2 as input...
	for ( core::Size i = 0; i < nseq_; ++i ) {
		iA_.push_back( i );
		iB_.push_back( i );
	}

	// L_ini_: N aligning residue window on each iter
	// for example, if nseq_ == 200:
	// L_ini_ = 200, 100, 50, 25, 12, 6
	core::Size const L_ini_min = 4;
	core::Size const n_init_max = 6;
	core::Size n_init( 0 );

	core::Size i( 0 );
	for ( i=1; i<=n_init_max-1; ++i ) {
		core::Size denom( 1 );
		for ( core::Size j=0; j<=i-1; ++j ) denom *= 2;

		L_ini_.push_back( core::Size(nseq_ / denom) );
		if ( L_ini_[i] <= L_ini_min ) {
			L_ini_[i] = L_ini_min;
			break;
		}
	}

	if ( i == n_init_max-1 ) {
		n_init++;
		L_ini_[n_init] = L_ini_min;
	}
}

// Use FArray2D to be compatible with rms_util.cc
void
TMscore::apply( FArray2D< core::Real > const &p2 )
{
	utility::vector0< core::Vector > u;
	core::Vector t;
	apply( p2, u, t, false );
}

void
TMscore::apply( FArray2D< core::Real > const &p2,
	utility::vector0< core::Vector > &u,
	core::Vector &t,
	bool const get_ut )

{

	// Minimize usage of FArray2D... could be easily converted into other array type
	convert_FArray2D_to_vector0( p2, xyzb_ );
	runtime_assert( xyza_.size() == xyzb_.size() );

	// Initialize scorestore
	score_.set_nseq( nseq_ );

	// Initialize alignment params
	get_ali_params();

	// Initialize u,t - we know they are not inputs
	core::Vector I( 0.0 );
	u.resize( 3 ); u[0] = I; u[1] = I; u[2] = I;
	t = I;

	utility::vector0< core::Size > k_ali0; //best continuous alignment
	core::Real d;

	// Iter by different window size
	// L_init: Window size, iL: starting residue num in window
	// I know that this variable names are weird, but let's keep as they are in the original code
	// to make easier to compare each other
	for ( core::Size i_init=1; i_init <= L_ini_.size(); ++i_init ) {
		core::Size L_init = L_ini_[i_init];
		core::Size iL_max = nseq_ - L_init + 1;

		// Iter by scanning window, iL: window starting residue
		for ( core::Size iL=1; iL <= iL_max; ++iL ) {
			TR.Debug << "i_init/L_init/iL/GDT: " << i_init << " " << L_init << " " << iL;
			TR.Debug << " " << get_GDTTS() << std::endl;

			// Get window residues within iL <= i < iL+L_init
			utility::vector0< core::Size > k_ali( L_init );
			for ( core::Size i=0; i<L_init; ++i ) {
				k_ali[i] = iL+i-1;
			}

			// Get whole transformed ref_crd aligned over k_ali
			utility::vector0< core::Vector > vt = get_transrot_ref( k_ali, u, t );

			d = d0_search_ - 1.0; //Search by smaller radii first
			// i_ali: matched residues, regardless of continuity
			utility::vector0< core::Size > i_ali = score_fun( d, vt, score_ );

			// update best alignment info
			if ( score_.TM_max < score_.TM ) {
				score_.TM_max = score_.TM;
				k_ali0 = k_ali;
			}
			score_.update();

			// Extend to larger distance, get best score given the alignment
			d = d0_search_ + 1.0; //extend cutoff
			extend( d, score_, u, t, i_ali, k_ali0 );

		} // iter over window
	} // iter over window size

	// All the calculations are done, but update u & t below
	// only if superposition matrix is required
	if ( get_ut ) {
		utility::vector0< core::Vector > vt = get_transrot_ref( k_ali0 , u, t );
	}

	return;

} /// end apply

/// 1, collect those residues with dis<d
/// 2. calculate score_GDT, score_maxsub, score_TM
utility::vector0< core::Size >
TMscore::score_fun( core::Real const d,
	utility::vector0< core::Vector > const vt,
	TMscoreStore &score ) const
{
	utility::vector0< core::Size > i_ali;
	core::Real d_tmp( d );

	while ( true ) {
		score.clear();
		i_ali.resize( 0 );

		for ( core::Size k=0; k< nseq_; ++k ) {
			core::Size i = iA(k);
			core::Size j = iB(k); // although i,j are always the same...
			core::Real const dis = vt[i].distance( xyzb_[j] );

			/// mark the residue-pairs in dis<d
			if ( dis < d_tmp ) i_ali.push_back( k );

			// update residue info
			score.add_residue_dis( d0_, dis );
		}

		if ( i_ali.size() < 3 && nseq_ > 3 ) {
			// If failed
			d_tmp += 0.5;
		} else {
			break;
		}
	}

	// Update finally
	score.update();
	score.apply();
	return i_ali;
}

void
TMscore::extend( core::Real const &d,
	TMscoreStore &score,
	utility::vector0< core::Vector > &u,
	core::Vector &t,
	utility::vector0< core::Size > &i_ali,
	utility::vector0< core::Size > &k_ali0
) const
{

	// iterative parameters;
	core::Size const n_it( 20 );
	// dummy

	// Iterate until
	// 1) all the residues are matched
	// 2) or iteration exceed
	core::Size it;
	for ( it=1; it<=n_it; ++it ) {
		// Pick residue matches that are within distance cut
		utility::vector0< core::Size > k_ali( i_ali );
		utility::vector0< Vector > vt = get_transrot_ref( k_ali, u, t );

		// get TMscore for given d & ref_crd(=vt)
		// update i_ali: The size will be changed!
		i_ali = score_fun( d, vt, score );

		//printf("iter/i_ali/k_ali: %4d %4d %4d\n", it, i_ali.size(), k_ali.size());

		// update current best: alignment first, followed by score
		if ( score.TM_max < score.TM ) k_ali0 = k_ali;
		score.update();

		// 1. If iteration exceed
		if ( it == n_it ) return;

		// 2. When the condition i_ali == k_ali be satisfied?
		// i_ali: segment alignment, k_ali: window thread
		// so if d<dis_cut becomes threaded alignment, it will be terminated
		if ( i_ali == k_ali ) return;
	}

	TR.Debug << "extension done after " << it << " with condition" << (it == n_it) << std::endl;
	return;
}

utility::vector0< core::Vector >
TMscore::get_transrot_ref( utility::vector0< core::Size > const &k_ali,
	utility::vector0< core::Vector > &u,
	core::Vector &t
) const
{
	// 1. Let's get u,t first
	// Kabsch using Rosetta built-in function
	{
		core::Size const n( k_ali.size() );

		// fill in r_1, r_2
		utility::vector1< core::Vector > r_1, r_2;
		r_1.resize( n );
		r_2.resize( n );
		for ( core::Size i = 0; i < n; ++i ) {
			r_1[i+1] = xyza_[k_ali[i]];
			r_2[i+1] = xyzb_[k_ali[i]];
		}

		// COM
		core::Vector t1( 0.0 ), t2( 0.0 );
		for ( core::Size j = 0; j < 3; ++j ) {
			for ( core::Size i = 1; i <= n; ++i ) {
				t1[j] += r_1[i][j];
				t2[j] += r_2[i][j];
			}
			t1[j] /= (core::Real)(n);
			t2[j] /= (core::Real)(n);
		}

		core::Real rmsd;
		utility::vector1< core::Real > ww( n, 1.0 );
		numeric::xyzMatrix< core::Real > uu;
		numeric::model_quality::findUU( r_1, r_2, ww, (int)(n), uu, rmsd ); // r_1,2 should be vector1< Vector >

		for ( core::Size i = 0; i < 3; ++i ) {
			for ( core::Size j = 0; j < 3; ++j ) {
				u[i][j] = uu(j+1, i+1);
			}
		}

		for ( core::Size j = 0; j < 3; ++j ) {
			t[j] = t2[j] - u[j][0]*t1[0] - u[j][1]*t1[1] - u[j][2]*t1[2];
		}
	}

	// apply transrot
	utility::vector0< core::Vector > vt( nseq_ );
	for ( core::Size j=0; j < nseq_; ++j ) {
		core::Vector const &xyzj( xyza_[j]);
		vt[j][0] = t[0] + u[0][0]*xyzj[0] + u[0][1]*xyzj[1] + u[0][2]*xyzj[2];
		vt[j][1] = t[1] + u[1][0]*xyzj[0] + u[1][1]*xyzj[1] + u[1][2]*xyzj[2];
		vt[j][2] = t[2] + u[2][0]*xyzj[0] + u[2][1]*xyzj[1] + u[2][2]*xyzj[2];
	}

	return vt;
}

} // namespace scoring
} // namespace core
