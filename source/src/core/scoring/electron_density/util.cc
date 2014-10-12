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
/// @author

#include <core/scoring/electron_density/util.hh>
#include <core/scoring/electron_density/SplineInterp.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/options/option.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>  // b factors

#include <iostream>

// option key includes
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/patterson.OptionKeys.gen.hh>

#include <utility/vector1.hh>


// C++ Headers

namespace core {
namespace scoring {
namespace electron_density {

/// @brief read density weights from the cmd line into the scorefunction
void add_dens_scores_from_cmdline_to_scorefxn( core::scoring::ScoreFunction &scorefxn ) {
	using namespace basic::options;

	if ( option[ OptionKeys::patterson::weight ].user() ) {
		scorefxn.set_weight( core::scoring::patterson_cc,
		                     option[ OptionKeys::patterson::weight ]() );
	}
	if ( option[ OptionKeys::edensity::fastdens_wt ].user() ) {
		scorefxn.set_weight( core::scoring::elec_dens_fast,
		                     option[ OptionKeys::edensity::fastdens_wt ]() );
	}
	if ( option[ OptionKeys::edensity::sliding_window_wt ].user() ) {
		scorefxn.set_weight( core::scoring::elec_dens_window,
		                     option[ OptionKeys::edensity::sliding_window_wt ]() );
	}
	if ( option[ OptionKeys::edensity::whole_structure_ca_wt ].user() ) {
		scorefxn.set_weight( core::scoring::elec_dens_whole_structure_ca,
		                     option[ OptionKeys::edensity::whole_structure_ca_wt ]() );
	}
	if ( option[ OptionKeys::edensity::whole_structure_allatom_wt ].user() ) {
		scorefxn.set_weight( core::scoring::elec_dens_whole_structure_allatom,
		                     option[ OptionKeys::edensity::whole_structure_allatom_wt ]() );
	}
}

bool pose_has_nonzero_Bs( core::pose::Pose const & pose ) {
	if (!pose.pdb_info()) return false;

	// to save time check only the first atom of each residue
	for (uint i = 1; i <= pose.total_residue(); ++i) {
		core::conformation::Residue const & rsd_i ( pose.residue(i) );
		if ( rsd_i.aa() == core::chemical::aa_vrt ) continue;
		Real B = pose.pdb_info()->temperature( i, 1 );
		if (B > 0) {
			return true;
		}
	}
	return false;
}

bool pose_has_nonzero_Bs( poseCoords const & pose ) {
	for (int i=1 ; i<=(int)pose.size(); ++i) {
		if( pose[i].B_ > 0) {
			return true;
		}
	}
	return false;
}


/// @brief spline interpolation with periodic boundaries
core::Real interp_spline(
                         ObjexxFCL::FArray3D< double > & coeffs ,
                         numeric::xyzVector< core::Real > const & idxX ) {
	int dims[3] = { coeffs.u3(), coeffs.u2(), coeffs.u1() };
	core::Real pt[3] = { idxX[2]-1.0 , idxX[1]-1.0, idxX[0]-1.0 };
	core::Real retval = SplineInterp::interp3(&coeffs[0], dims, pt);
	return retval;
}

/// @brief spline interpolation with periodic boundaries
numeric::xyzVector<core::Real> interp_dspline(
                         ObjexxFCL::FArray3D< double > & coeffs ,
                         numeric::xyzVector< core::Real > const & idxX ) {
	int dims[3] = { coeffs.u3(), coeffs.u2(), coeffs.u1() };
	core::Real pt[3] = { idxX[2]-1.0 , idxX[1]-1.0, idxX[0]-1.0 };
	core::Real grad[3] = { 0,0,0 };
	SplineInterp::grad3(&grad[0], &coeffs[0], dims, pt);
	return numeric::xyzVector<core::Real>(grad[2],grad[1],grad[0]);
}

void spline_coeffs(
           ObjexxFCL::FArray3D< double > const & data ,
           ObjexxFCL::FArray3D< double > & coeffs) {
	int dims[3] = { data.u3(), data.u2(), data.u1() };
	coeffs = data;
	SplineInterp::compute_coefficients3( const_cast<double*>(&coeffs[0]) , dims );  // external code wants nonconst even though array is unchanged
}

void spline_coeffs(
                         ObjexxFCL::FArray3D< float > const & data ,
                         ObjexxFCL::FArray3D< double > & coeffs) {
	int N = data.u3()*data.u2()*data.u1();
	ObjexxFCL::FArray3D< double > data_d(data.u1(),data.u2(),data.u3()) ;
	for (int i=0; i<N; ++i)
		data_d[i] = (double)data[i];
	spline_coeffs( data_d, coeffs );
}



void conj_map_times(ObjexxFCL::FArray3D< std::complex<double> > & map_product, ObjexxFCL::FArray3D< std::complex<double> > const & mapA, ObjexxFCL::FArray3D< std::complex<double> > const & mapB) {
	assert(mapA.u1() == mapB.u1());
	assert(mapA.u2() == mapB.u2());
	assert(mapA.u3() == mapB.u3());

	map_product.dimension(mapA.u1(), mapA.u2(), mapA.u3());
	for (Size i=0; i < mapA.size(); i++) {
		map_product[i] = std::conj(mapA[i]) * mapB[i];
	}
}

//∑[A(x) * B(y-x)] over y
ObjexxFCL::FArray3D< double > convolute_maps( ObjexxFCL::FArray3D< double > const & mapA, ObjexxFCL::FArray3D< double > const & mapB) {

	ObjexxFCL::FArray3D< std::complex<double> > FmapA;
	numeric::fourier::fft3(mapA, FmapA);

	ObjexxFCL::FArray3D< std::complex<double> > FmapB;
	numeric::fourier::fft3(mapB, FmapB);

	ObjexxFCL::FArray3D< std::complex<double> > Fconv_map;
	conj_map_times(Fconv_map, FmapB, FmapA );

	ObjexxFCL::FArray3D< double > conv_map;
	numeric::fourier::ifft3(Fconv_map , conv_map);

	return conv_map;
}


///
/// 4D interpolants

/// @brief spline interpolation with periodic boundaries
core::Real interp_spline(
                         ObjexxFCL::FArray4D< double > & coeffs ,
                         core::Real slab,
                         numeric::xyzVector< core::Real > const & idxX )
{
	int dims[4] = { coeffs.u4(), coeffs.u3(), coeffs.u2(), coeffs.u1() };
	core::Real pt[4] = { slab-1.0, idxX[2]-1.0 , idxX[1]-1.0, idxX[0]-1.0 };
	core::Real retval = SplineInterp::interp4(&coeffs[0], dims, pt);

	return retval;
}

/// @brief spline interpolation with periodic boundaries
void interp_dspline(
		ObjexxFCL::FArray4D< double > & coeffs ,
		numeric::xyzVector< core::Real > const & idxX ,
		core::Real slab,
		numeric::xyzVector< core::Real > & gradX,
		core::Real & gradSlab )
{
	int dims[4] = { coeffs.u4(), coeffs.u3(), coeffs.u2(), coeffs.u1() };
	core::Real pt[4] = { slab-1.0, idxX[2]-1.0 , idxX[1]-1.0, idxX[0]-1.0 };
	core::Real grad[4] = { 0,0,0,0 };
	SplineInterp::grad4(&grad[0], &coeffs[0], dims, pt);
	gradX = numeric::xyzVector<core::Real>(grad[3],grad[2],grad[1]);
	gradSlab = grad[0];
}

void spline_coeffs(
           ObjexxFCL::FArray4D< double > const & data ,
           ObjexxFCL::FArray4D< double > & coeffs)
{
	int dims[4] = { coeffs.u4(), coeffs.u3(), coeffs.u2(), coeffs.u1() };
	coeffs = data;
	SplineInterp::compute_coefficients4( const_cast<double*>(&coeffs[0]) , dims );
}

void spline_coeffs(
           ObjexxFCL::FArray4D< float > const & data ,
           ObjexxFCL::FArray4D< double > & coeffs)
{
 	int N = data.u4()*data.u3()*data.u2()*data.u1();
	ObjexxFCL::FArray4D< double > data_d(data.u1(),data.u2(),data.u3(),data.u4()) ;
	for (int i=0; i<N; ++i)
		data_d[i] = (double)data[i];
	spline_coeffs( data_d, coeffs );
}


} // namespace constraints
} // namespace scoring
} // namespace core
