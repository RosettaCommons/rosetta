// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file

#include <protocols/electron_density/DensitySymmInfo.hh>


#include <ObjexxFCL/format.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/string_util.hh>

#include <protocols/electron_density/util.hh>

#include <utility/tag/Tag.hh>

#include <numeric/random/random.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyz.functions.hh>

#include <basic/Tracer.hh>

#include <ObjexxFCL/format.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.electron_density.DensitySymmInfo" );

namespace protocols {
namespace electron_density {

// auto-detect symmetry axes from an electron density map
void
DensitySymmInfo::detect_axes( core::scoring::electron_density::ElectronDensity const &e ) {
	if ( count_primary == 1 ) return; // make sure a spacegroup is assigned

	TR << "Detecting symmetry axes." << std::endl;

	utility::vector1< numeric::xyzVector<core::Real> > axes_to_test;
	axes_to_test.push_back(numeric::xyzVector<core::Real>(1,0,0));
	axes_to_test.push_back(numeric::xyzVector<core::Real>(0,1,0));
	axes_to_test.push_back(numeric::xyzVector<core::Real>(0,0,1));

	utility::vector1< core::Real > correls( axes_to_test.size() );
	utility::vector1< numeric::xyzVector< core::Real > > centers( axes_to_test.size() );

	for ( int i=1; i<=(int)axes_to_test.size(); ++i ) {
		autocorrelate( e, 2*numeric::constants::d::pi / count_primary,  axes_to_test[i], centers[i], correls[i] );
	}

	// find best autocorrel
	core::Real best_correl=-1e30;
	core::Size best_correl_idx=0;
	for ( int i=1; i<=(int)axes_to_test.size(); ++i ) {
		if ( correls[i] > best_correl ) {
			best_correl = correls[i];
			best_correl_idx = i;
		}
	}

	if ( best_correl < 0.5 ) {
		TR << "ERROR!  C symmetry not detected in the density." << std::endl;
		utility_exit();
	}

	symm_center = centers[best_correl_idx];
	axis_primary = axes_to_test[best_correl_idx];
	core::Size old_best_correl_idx = best_correl_idx;

	TR << "Found " << count_primary << "-fold symm axis: " << axis_primary[0] << ","<< axis_primary[1] << ","<< axis_primary[2]
		<< "  center: " << symm_center[0] << ","<< symm_center[1] << ","<< symm_center[2]
		<< "  autocorrelation=" << best_correl << std::endl;

	if ( type == 'D' ) {
		// for D2, also find 2nd best autocorrel
		best_correl=-1e30;
		best_correl_idx=0;
		if ( count_primary != 2 ) {
			// for Dn, n>2, rerun correlation using twofold
			for ( int i=1; i<=(int)axes_to_test.size(); ++i ) {
				if ( i==(int)old_best_correl_idx ) continue;
				autocorrelate( e, numeric::constants::d::pi, axes_to_test[i], centers[i], correls[i] );
			}
		}

		for ( int i=1; i<=(int)axes_to_test.size(); ++i ) {
			if ( i==(int)old_best_correl_idx ) continue;
			if ( correls[i] > best_correl ) {
				best_correl = correls[i];
				best_correl_idx = i;
			}
		}

		if ( best_correl < 0.5 ) {
			TR << "ERROR!  D symmetry not detected in the density." << std::endl;
			utility_exit();
		}

		// resolve symm center if we are D symmetric
		// intersection of two symm axes
		axis_secondary = axes_to_test[best_correl_idx];
		core::Real error=0.0;

		symm_center = resolve_symm_axes( symm_center, axis_primary, centers[best_correl_idx], axis_secondary, error );

		TR << "Found 2-fold symm axis: " << axis_secondary[0] << ","<< axis_secondary[1] << ","<< axis_secondary[2]
			<< "  resolved center: " << symm_center[0] << ","<< symm_center[1] << ","<< symm_center[2]
			<< "  with error=" << error
			<< "  autocorrelation=" << best_correl << std::endl;

		if ( error > 1.0 ) {
			TR << "ERROR!  D symmetry center could not be resolved." << std::endl;
			utility_exit();
		}
	}
}


numeric::xyzVector< core::Real >
DensitySymmInfo::resolve_symm_axes(
	numeric::xyzVector< core::Real > const &p1, numeric::xyzVector< core::Real > const &p21,
	numeric::xyzVector< core::Real > const &p3, numeric::xyzVector< core::Real > const &p43,
	core::Real &error )
{
	numeric::xyzVector< core::Real > p2 = p1+p21, p4 = p3+p43;
	numeric::xyzVector< core::Real > p13 = p1-p3;

	core::Real d1343 = p13[0] * p43[0] + p13[1] * p43[1] + p13[2] * p43[2];
	core::Real d4321 = p43[0] * p21[0] + p43[1] * p21[1] + p43[2] * p21[2];
	core::Real d1321 = p13[0] * p21[0] + p13[1] * p21[1] + p13[2] * p21[2];
	core::Real d4343 = p43[0] * p43[0] + p43[1] * p43[1] + p43[2] * p43[2];
	core::Real d2121 = p21[0] * p21[0] + p21[1] * p21[1] + p21[2] * p21[2];

	core::Real denom = d2121 * d4343 - d4321 * d4321;
	core::Real numer = d1343 * d4321 - d1321 * d4343;

	core::Real mua, mub;

	if ( std::fabs( denom ) < 1e-6 ) { // lines are parallel OR p1 & p3 are coincident
		mua = mub = 0;
	} else {
		mua = numer / denom;
		mub = (d1343 + d4321 * mua) / d4343;
	}

	error = ( (p1+mua*p21) - (p3+mub*p43) ).length();

	return (( (p1+mua*p21) + (p3+mub*p43) )/2.0 );
}

void
DensitySymmInfo::autocorrelate(
	core::scoring::electron_density::ElectronDensity const &e,
	core::Real angle,
	numeric::xyzVector< core::Real > const &axis,
	numeric::xyzVector< core::Real > &symm_center,
	core::Real &autocorrelation
) {
	core::Size undersample = 1; // don't need to use fullmap

	// resample
	ObjexxFCL::FArray3D< float > const & densdata = e.get_data();
	numeric::xyzVector<core::Size> grid(
		densdata.u1()/undersample, densdata.u2()/undersample, densdata.u3()/undersample );
	ObjexxFCL::FArray3D< double > dens, rot;
	rot.dimension( grid[0], grid[1], grid[2] );
	dens.dimension( grid[0], grid[1], grid[2] );
	dens=0;

	for ( int k=1; k<= dens.u3(); ++k ) {
		for ( int j=1; j<= dens.u2(); ++j ) {
			for ( int i=1; i<= dens.u1(); ++i ) {
				for ( int dk=1; dk<=(int)undersample; ++dk ) {
					for ( int dj=1; dj<=(int)undersample; ++dj ) {
						for ( int di=1; di<=(int)undersample; ++di ) {
							dens(i,j,k) +=
								densdata( (i-1)*undersample+di, (j-1)*undersample+dj, (k-1)*undersample+dk );
						}
					}
				}
			}
		}
	}

	// rotate density map about center of mass by angle (no wrapping)
	ObjexxFCL::FArray3D< std::complex<double> > Fdens, Frot;
	rot.dimension( grid[0], grid[1], grid[2] );
	numeric::xyzMatrix<core::Real> R = numeric::rotation_matrix( axis, angle );

	ObjexxFCL::FArray3D< double > coeffs_density;
	core::scoring::electron_density::spline_coeffs( dens , coeffs_density );

	numeric::xyzVector<core::Real> idx_center ( grid[0]/2.0+1 , grid[1]/2.0+1 , grid[2]/2.0+1 );

	for ( int k=1; k<= dens.u3(); ++k ) {
		for ( int j=1; j<= dens.u2(); ++j ) {
			for ( int i=1; i<= dens.u1(); ++i ) {
				numeric::xyzVector<core::Real> cart_ijk =
					e.get_f2c()* (numeric::xyzVector<core::Real>(
					(i-idx_center[0])/grid[0], (j-idx_center[1])/grid[1], (k-idx_center[2])/grid[2] ));
				numeric::xyzVector<core::Real> Rfrac_ijk = e.get_c2f()*R*cart_ijk;
				numeric::xyzVector<core::Real> Ridx_ijk =
					numeric::xyzVector<core::Real>( Rfrac_ijk[0]*grid[0], Rfrac_ijk[1]*grid[1], Rfrac_ijk[2]*grid[2] ) + idx_center;
				rot(i,j,k) = core::scoring::electron_density::interp_spline( coeffs_density , Ridx_ijk );
			}
		}
	}

	// debug
	{
		core::scoring::electron_density::ElectronDensity ecopy = e;
		ecopy.set_data( rot );
		ecopy.writeMRC( "rho_rot.mrc" );
	}

	// FFT correlation
	numeric::fourier::fft3(dens, Fdens);
	numeric::fourier::fft3(rot, Frot);
	core::Size npoints = grid[0]*grid[1]*grid[2];
	for ( int i=0; i<(int)npoints ; ++i ) {
		Frot[i] = Fdens[i] * std::conj( Frot[i] );  // reuse Frot storage
	}
	numeric::fourier::ifft3(Frot, rot);

	// parse results
	core::Real maxcorrel=-1e30, dens2=0;
	numeric::xyzVector<core::Real> best_idx(0,0,0);
	for ( int k=1; k<= dens.u3(); ++k ) {
		for ( int j=1; j<= dens.u2(); ++j ) {
			for ( int i=1; i<= dens.u1(); ++i ) {

				if ( rot(i,j,k) > maxcorrel ) {
					maxcorrel = rot(i,j,k);
					best_idx = numeric::xyzVector<core::Real>(i,j,k);
				}
				dens2 += dens(i,j,k)*dens(i,j,k);
			}
		}
	}

	// subvoxel sampling...
	core::scoring::electron_density::spline_coeffs( rot , coeffs_density );
	numeric::xyzVector<core::Real> currX=best_idx, prevX;
	core::Real func = core::scoring::electron_density::interp_spline( coeffs_density , currX ), func_prev;
	core::Real step = 0.0001;
	for ( int ITER=1; ITER<=50; ++ITER ) {
		prevX = currX;
		func_prev = func;
		numeric::xyzVector<core::Real> dX = core::scoring::electron_density::interp_dspline( coeffs_density , currX );
		if ( dX.length() < 0.0001*dens2 ) break;  //?
		currX += step*dX;
		func = core::scoring::electron_density::interp_spline( coeffs_density , currX );
		if ( func <= func_prev+1e-4 ) {
			step /= 2;
			currX = prevX;
			func = func_prev;
		}
	}
	best_idx = currX;
	maxcorrel = func;

	// currently we have 'B': the translation between the original and rotated center
	// we need to get 'T': the translation from the rotation center to the true center
	best_idx = best_idx-1.0;
	if ( best_idx[0]>=dens.u1()/2 ) best_idx[0]-=dens.u1();
	if ( best_idx[1]>=dens.u2()/2 ) best_idx[1]-=dens.u2();
	if ( best_idx[2]>=dens.u3()/2 ) best_idx[2]-=dens.u3();

	// first get 'T_hat', a vector from rot axis to midpoint of orig and rotated center
	core::Real lenT = best_idx.length() / (2 * tan(angle/2));

	// direction of T
	numeric::xyzVector<core::Real> dirT;
	if ( std::fabs( angle - numeric::constants::d::pi ) < 1e-4 ) {
		dirT = -best_idx;
	} else {
		dirT = axis.cross( best_idx );
	}
	if ( dirT.length() > 1e-4 ) {
		dirT.normalize();
	}

	best_idx = best_idx/2.0-lenT*dirT;

	e.idx2cart( ((core::Real)undersample)*(best_idx+idx_center), symm_center );

	autocorrelation = maxcorrel /(dens2);

	TR << "Found autocorrelation of " << autocorrelation << " for axis [" << axis[0] << ","<< axis[1] << ","<< axis[2] << "]" << std::endl;
}

void
DensitySymmInfo::mask_asu( ObjexxFCL::FArray3D< double > &vol, core::scoring::electron_density::ElectronDensity const &e, double value ) {
	numeric::xyzVector< core::Real > grid(vol.u1(), vol.u2(), vol.u3()), idx_center;

	if ( type=='C' || type == 'D' ) {
		// get grid coords of symm center
		e.cart2idx( symm_center, idx_center );

		// build our coord system
		numeric::xyzVector< core::Real > X,Y,Z=axis_primary;
		if ( type == 'D' ) {
			X = axis_secondary;
		} else {
			if ( axis_primary[1]!=0 || axis_primary[2]!=0 ) {
				X = axis_primary.cross( numeric::xyzVector< core::Real >(1,0,0) );
			} else {
				X = axis_primary.cross( numeric::xyzVector< core::Real >(0,1,0) );
			}
			X.normalize();
		}
		Y = Z.cross(X);

		numeric::xyzMatrix< core::Real > R_canonic = numeric::xyzMatrix< core::Real >::rows(X,Y,Z);
		for ( int k=1; k<= (int)vol.u3(); ++k ) {
			for ( int j=1; j<= (int)vol.u2(); ++j ) {
				for ( int i=1; i<= (int)vol.u1(); ++i ) {
					numeric::xyzVector<core::Real> cart_ijk_canonic =
						R_canonic * e.get_f2c()* (numeric::xyzVector<core::Real>(
						(i-idx_center[0])/grid[0], (j-idx_center[1])/grid[1], (k-idx_center[2])/grid[2] ));

					// find the angle w.r.t. canonic

					core::Real theta = atan2( cart_ijk_canonic[1], cart_ijk_canonic[0] );
					bool theta_good = (theta >= 0) && (theta <= (2.0 * numeric::constants::d::pi / count_primary));
					if ( type == 'C' ) {
						if ( theta_good ) continue;
					} else if ( type == 'D' ) {
						core::Real phi = asin( cart_ijk_canonic[2] / cart_ijk_canonic.length() );
						bool phi_good = (phi >= 0);
						if ( theta_good && phi_good ) {
							continue;
						}
					}

					vol(i,j,k) = value;
				}
			}
		}
	} // type = C or D
}

// min symm dist between X and _any symmetric copy_ of Y
core::Real
DensitySymmInfo::min_symm_dist2(
	numeric::xyzVector< core::Real > const &X,
	numeric::xyzVector< core::Real > const &Y
) const {
	if ( !enabled() ) {
		return (X - Y).length_squared();
	}

	core::Real mindist = (X - Y).length_squared();
	if ( type == 'C' ) {
		for ( int i=2; i<=(int)count_primary; ++i ) {
			core::Real angle = (i-1) * 2*numeric::constants::d::pi / count_primary;
			numeric::xyzMatrix<core::Real> R_i = numeric::rotation_matrix( axis_primary, angle );
			numeric::xyzVector<core::Real> Y_i = R_i*(Y-symm_center) + symm_center;
			mindist = std::min( (X - Y_i).length_squared() , mindist );
		}
	} else if ( type == 'D' ) {
		for ( int i=1; i<=(int)count_primary; ++i ) {
			core::Real angle_i = (i-1) * 2*numeric::constants::d::pi / count_primary;
			numeric::xyzMatrix<core::Real> R_i = numeric::rotation_matrix( axis_primary, angle_i );
			for ( int j=1; j<=2; ++j ) {
				if ( i==1 && j==1 ) continue; // identity
				core::Real angle_j = (j==1) ? 0 : numeric::constants::d::pi;
				numeric::xyzMatrix<core::Real> R_j = numeric::rotation_matrix( axis_secondary, angle_j );
				numeric::xyzVector<core::Real> Y_i = R_j*R_i*(Y-symm_center) + symm_center;
				mindist = std::min( (X - Y_i).length_squared() , mindist );
			}
		}
	}

	return mindist;
}

}
}
