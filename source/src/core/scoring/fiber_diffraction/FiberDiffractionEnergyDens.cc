// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/fiber_diffractiobn/FiberDiffractionEnergyDens.cc
/// @brief  FiberDiffractionDens scoring (low resolution, electron-density based)
/// @author Wojciech Potrzebowski and Ingemar Andre

// Unit headers
#include <core/scoring/fiber_diffraction/FiberDiffractionEnergyDens.hh>
#include <core/scoring/fiber_diffraction/FiberDiffraction.hh>
#include <core/scoring/fiber_diffraction/util.hh>
#include <core/scoring/fiber_diffraction/util_ndft.hh>
#include <core/scoring/fiber_diffraction/hankel_kiss_fft.hh>
#include <core/scoring/fiber_diffraction/xray_scattering.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/conformation/Residue.hh>

#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/io/pdb/file_data.hh>

#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/symmetry/util.hh>

// Options
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>

#include <basic/datacache/BasicDataCache.hh>
#include <core/pose/datacache/CacheableDataType.hh>


// ObjexxFCL headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/FArray2D.hh>

#include <core/id/AtomID.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>

// Utility headers
#include <basic/prof.hh>
#include <basic/database/open.hh>
#include <basic/Tracer.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/basic_sys_util.hh>
#include <boost/math_fwd.hpp>

// C++
#include <iomanip>
#include <string>
//#include <sys/time.h>

#ifdef WIN32
#define _USE_MATH_DEFINES
#include <math.h>
#endif

namespace core {
namespace scoring {
namespace fiber_diffraction {

// tracer
static basic::Tracer TR("core.scoring.fiber_diffraction.FiberDiffractionEnergyDens");

ScoreTypes FiberDiffractionEnergyDensCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( fiberdiffractiondens );
	return sts;
}

methods::EnergyMethodOP FiberDiffractionEnergyDensCreator::create_energy_method( methods::EnergyMethodOptions const &) const {
	return methods::EnergyMethodOP( new FiberDiffractionEnergyDens() );
}

FiberDiffractionEnergyDens::FiberDiffractionEnergyDens() :
	parent( methods::EnergyMethodCreatorOP( new FiberDiffractionEnergyDensCreator ) ) {
	chi2_=0;
	dchi2_d.clear();
	dchi2_d_cross_R.clear();
}

/// clone
methods::EnergyMethodOP FiberDiffractionEnergyDens::clone() const {
	return methods::EnergyMethodOP( new FiberDiffractionEnergyDens( *this ) );
}

void FiberDiffractionEnergyDens::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const {

	if ( !core::pose::symmetry::is_symmetric(pose) ) {
		utility_exit_with_message("Structure needs to be symmetric! Aborting...");
	}

	// load fiber diffraction data
	utility::vector0< utility::vector1< core::Real > >::iterator layer_lines_I;
	utility::vector0< utility::vector1< core::Real > >::iterator layer_lines_R;
	utility::vector0 < utility::vector0 < int > >::iterator nvals;

	core::Size lmax, Rmax;

	if ( basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::a ].user() ) {
		a_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::a ]();
	}

	if ( basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::b ].user() ) {
		b_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::b ]();
	}

	if ( !(a_ > 0 ) ) {
		utility_exit_with_message("The number of subunits per repeat, score::fiber_diffraction::a, must be set!");
	}

	if ( !(b_ > 0 ) ) {
		utility_exit_with_message("The number of turns, score::fiber_diffraction::b, must be set!");
	}

	if ( basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::p ].user() ) {
		p_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::p ]();
	} else {
		find_pitch( pose, p_ );
	}
	c_ = p_*a_;

	core::Real  res_cutoff_low_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::resolution_cutoff_low ]();

	core::Real res_cutoff_high_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::resolution_cutoff_high ]();

	getFiberDiffractionData(c_, res_cutoff_high_, res_cutoff_low_).getAllFiberData(layer_lines_I, layer_lines_R, nvals, lmax, Rmax);
	TR << "Rmax : lmax " << "( " << Rmax << " : " << lmax << " )" << std::endl;

	core::Size grid_r_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::grid_r ]();
	core::Size grid_phi_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::grid_phi ]();
	core::Size grid_z_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::grid_z ]();

	core::Real qfht_K1_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::qfht_K1 ]();
	core::Real qfht_K2_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::qfht_K2 ]();

	Size max_bessel_orders(0);
	for ( Size l=0; l <= lmax; ++l ) {
		Size max_b_order( nvals[l].size() );
		if ( max_b_order > max_bessel_orders ) max_bessel_orders = max_b_order;
	}

	ObjexxFCL::FArray1D< float > rc, phic, zc, Rinv_ht;
	phic.dimension( grid_phi_ );
	zc.dimension( grid_z_ );

	rc.dimension( grid_r_ );
	Rinv_ht.dimension( grid_r_ );

	double max_r_value;
	find_max_r( pose, max_r_value );
	double PAD_R (5 );
	max_r_value += PAD_R;

	//Hankel transform functions:
	set_r_array( grid_r_, qfht_K1_, qfht_K2_, max_r_value, rc );
	set_r_inv_array( grid_r_, qfht_K1_, qfht_K2_, max_r_value, Rinv_ht );

	// calculate rho
	//Setting up cylindrical coordinates
	for ( Size i=0; i<grid_phi_; i++ ) {
		phic(i+1)=float(i)/float(grid_phi_)*2*M_PI;
	}

	for ( Size i=0; i<grid_z_; i++ ) {
		zc(i+1)=i*c_/grid_z_;
	}

	calculate_rho_fast2( pose, rc, phic, zc, c_ );

	ObjexxFCL::FArray3D< std::complex<float> > Gnl;
	Gnl.dimension( grid_r_, lmax+1, max_bessel_orders );

	gnl_R_qfht( rho_cylindrical_, nvals, lmax, max_r_value, qfht_K1_, qfht_K2_, Gnl );

	I.resize(lmax+1);
	for ( Size l=0; l<= lmax; ++l ) {
		I[l].resize(Rmax);
		for ( Size i=1; i<= Rmax; ++i ) I[l][i] = 0.0;
	}


	utility::vector0< utility::pointer::shared_ptr< numeric::interpolation::spline::Interpolator > > bessel_order_splines;
	Size lsum=0;
	for ( Size l=0; l <= lmax; ++l ) {
		std::complex< float > GNLR(0,0);
		Size max_b_order( nvals[l].size() );
		for ( Size b_order=1; b_order <= max_b_order; ++b_order ) {
			ObjexxFCL::FArray1D < float > Ia, Ra;
			Ia.dimension(grid_r_);
			Ra.dimension(grid_r_);
			for ( Size Rinv=1; Rinv <= grid_r_; ++Rinv ) {
				std::complex< float > GNLR = Gnl( Rinv,l+1,b_order);
				double absval = abs(GNLR);
				Ra(Rinv) = Rinv_ht(Rinv);
				Ia(Rinv) = absval*absval;
			}
			utility::pointer::shared_ptr< numeric::interpolation::spline::Interpolator > inter = fit_layer_lines_with_splines( Ra, Ia );
			bessel_order_splines.push_back(inter);
		}
		lsum +=max_b_order;
	}
	// outdebug.close();

	core::Size total_b_order(0);
	for ( Size l=0; l <= lmax; ++l ) {
		Size max_b_order( nvals[l].size() );
		for ( Size b_order=1; b_order <= max_b_order; ++b_order ) {
			for ( Size R=1; R <= layer_lines_R[l].size(); ++R ) {
				utility::pointer::shared_ptr< numeric::interpolation::spline::Interpolator > inter = bessel_order_splines[total_b_order];
				Real Ival;
				Real dy = 0;
				inter->interpolate( layer_lines_R[l][R], Ival, dy );
				I[l][R]  +=  fabs(Ival);
			}
			++total_b_order;
		}
	}

	Real prod(0);
	square_obs_ = 0;
	sum_obs_ = 0;
	for ( Size l=0; l <= lmax; ++l ) {
		for ( Size R=1; R<=layer_lines_R[l].size(); ++R ) {
			Real F_obs_square ( layer_lines_I[l][R]*layer_lines_I[l][R] );
			prod += I[l][R]*F_obs_square;
			square_obs_ +=  F_obs_square*F_obs_square;
			sum_obs_ +=  F_obs_square;
		}
	}

	scale_factor_ = square_obs_/prod;
	TR << scale_factor_ << std::endl;

	bool output_fiber_spectra_(false);
	if ( basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::output_fiber_spectra ].user() ) {
		output_fiber_spectra_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::output_fiber_spectra ]();
	}
	std::ofstream out;
	if ( output_fiber_spectra_ ) {
		std::string outfile = "IntensityDens.txt";
		out.open(outfile.c_str(), std::ios::out);
	}

	chi2_=0;
	for ( Size l=0; l <= lmax; ++l ) {
		for ( Size R=1; R<=layer_lines_R[l].size(); ++R ) {
			chi2_ +=  (scale_factor_*I[l][R]-layer_lines_I[l][R]*layer_lines_I[l][R])*(scale_factor_*I[l][R]-layer_lines_I[l][R]*layer_lines_I[l][R]);
			if ( output_fiber_spectra_ ) {
				out << I[l][R] << " " << layer_lines_R[l][R] <<" "<< l << std::endl;
			}
		}
	}
	chi2_ /=square_obs_;

	TR << "chi2 " << chi2_ << " sum_obs_ " << sum_obs_ << std::endl;
	if ( output_fiber_spectra_ ) {
		out.close();
	}
}

void FiberDiffractionEnergyDens::finalize_total_energy(pose::Pose & /*pose*/, ScoreFunction const &, EnergyMap & emap) const {
	// all the work is done in setup for scoring
	// just return results here
	emap[ fiberdiffractiondens ] += chi2_;
}


void
FiberDiffractionEnergyDens::calculate_rho_fast2(
	pose::Pose & pose,
	ObjexxFCL::FArray1D< float > & rc,
	ObjexxFCL::FArray1D< float > & phic,
	ObjexxFCL::FArray1D< float > & zc,
	Real const c_
) const
{

	// Are we symmetric?
	const core::conformation::symmetry::SymmetryInfo *symminfo=NULL;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		symminfo = dynamic_cast<const core::conformation::symmetry::SymmetricConformation & >(
			pose.conformation()).Symmetry_Info().get();
	}

	if ( !symminfo ) {
		utility_exit_with_message("Structure needs to be symmetric! Aborting...");
	}

	utility::vector1< OneGaussianScattering > sig_centroid_(setup_centroid_scatter( pose ));
	utility::vector1< OneGaussianScattering >::iterator sig_centroid( sig_centroid_.begin() );

	core::Real minX, minY, minZ, maxX, maxY, maxZ;

	find_min_xyz(pose, minX, minY, minZ, maxX, maxY, maxZ);

	if ( fabs(maxZ-minZ) > c_ ) maxZ = minZ+c_;
	core::Real minZZ(minZ), maxZZ(maxZ);

	core::Real pad(5);

	minX -= pad;
	minY -= pad;
	minZ -= pad;
	maxX += pad;
	maxY += pad;
	maxZ += pad;


	// initialize grid and density
	core::Real grid_reso = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::grid_reso ]();
	//core::Real grid_reso(0.5);

	core::Real maxDelta( fabs(maxX-minX) );
	if ( fabs(maxY-minY)  > maxDelta ) maxDelta = fabs(maxY-minY);
	if ( fabs(maxZ-minZ)  > maxDelta ) maxDelta = fabs(maxZ-minZ);
	int grid_points(maxDelta/grid_reso);

	ObjexxFCL::FArray3D< float > rho_cartesian;
	rho_cartesian.dimension(grid_points , grid_points, grid_points );
	rho_cylindrical_.dimension( rc.size(), phic.size() , zc.size() );
	Size nscatterers(0);
	find_num_scattering_atoms( pose, nscatterers );

	for ( int i=0; i< rho_cartesian.u1()*rho_cartesian.u2()*rho_cartesian.u3(); ++i ) rho_cartesian[i]=0.0;
	for ( int i=0; i< rho_cylindrical_.u1()*rho_cylindrical_.u2()*rho_cylindrical_.u3(); ++i ) rho_cylindrical_[i]=0.0;

	// B factor...
	core::Real b_factor_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::b_factor ]();

	// Sidechain weight
	core::Real SC_scaling = 0.92;
	if ( basic::options::option[ basic::options::OptionKeys::edensity::sc_scaling ].user() ) {
		SC_scaling = basic::options::option[ basic::options::OptionKeys::edensity::sc_scaling ]();
	}
	if ( pose.is_fullatom() ) SC_scaling = 1.0;

	TR.Debug << "SC_scaling: " << SC_scaling << std::endl;

	// ATOM_CUTOFF
	core::Real ATOM_CUTOFF = 2*3.2;
	int MAX_GRID_DEVIATON = int( 3.2/grid_reso ) + 1;
	utility::vector0 < int > grid_pos_x, grid_pos_y,  grid_pos_z;
	for ( int x=-MAX_GRID_DEVIATON; x<= MAX_GRID_DEVIATON; ++x ) {
		for ( int y=-MAX_GRID_DEVIATON; y<= MAX_GRID_DEVIATON; ++y ) {
			for ( int z=-MAX_GRID_DEVIATON; z<= MAX_GRID_DEVIATON; ++z ) {
				core::Real d2 = grid_reso*grid_reso*( x*x + y*y + z*z );
				if ( d2 < ATOM_CUTOFF*ATOM_CUTOFF ) {
					grid_pos_x.push_back(x);
					grid_pos_y.push_back(y);
					grid_pos_z.push_back(z);
				}
			}
		}
	}

	int nres = pose.total_residue();
	numeric::xyzVector< core::Real > atmi_xyz_(0.0,0.0,0.0);
	int zr(0);

	for ( int i=1 ; i<=nres; ++i ) {
		if ( ! symminfo->bb_is_independent( i ) ) continue;
		conformation::Residue const &rsd_i (pose.residue(i));
		// skip vrts & masked reses
		if ( rsd_i.aa() == core::chemical::aa_vrt ) continue;
		//    if ( scoring_mask_.find(i) != scoring_mask_.end() ) continue;
		int nheavyatoms = rsd_i.nheavyatoms();
		for ( int j=1 ; j<=nheavyatoms; ++j ) {
			//   id::AtomID id( j, i ); //deriv
			conformation::Atom const &atm_i( rsd_i.atom(j) ); //(pose.residue(i).atom("CA"));
			chemical::AtomTypeSet const & atom_type_set( rsd_i.atom_type_set() );
			std::string elt_i = atom_type_set[ rsd_i.atom_type_index( j ) ].element();
			OneGaussianScattering sig_j;
			if ( elt_i == "X" ) {
				int aa_num ( rsd_i.aa() -1 );
				sig_j = sig_centroid[ aa_num  ];
			} else {
				sig_j = core::scoring::fiber_diffraction::get_A( elt_i );
			}
			core::Real k = sig_j.k( b_factor_, grid_reso);
			core::Real C = sig_j.C( k );
			// sidechain weight
			if ( (Size) j > rsd_i.last_backbone_atom() ) {
				C *= SC_scaling;
			}

			if ( atm_i.xyz()[2]>maxZZ ) {
				zr = int(fabs(atm_i.xyz()[2]-minZZ)/c_);

				atmi_xyz_[0] =  atm_i.xyz()[0];
				atmi_xyz_[1] =  atm_i.xyz()[1];
				atmi_xyz_[2] =  atm_i.xyz()[2]-zr*c_;
			} else {
				atmi_xyz_[0] =  atm_i.xyz()[0];
				atmi_xyz_[1] =  atm_i.xyz()[1];
				atmi_xyz_[2] =  atm_i.xyz()[2];
			}
			int gridX =  int( ( atmi_xyz_[0] - minX )/grid_reso + 0.5 );
			int gridY =  int( ( atmi_xyz_[1] - minY )/grid_reso + 0.5 );
			int gridZ =  int( ( atmi_xyz_[2] - minZ )/grid_reso + 0.5 );


			for ( core::Size pos = 0; pos < grid_pos_x.size(); ++pos ) {
				int x = grid_pos_x[pos];
				int y = grid_pos_y[pos];
				int z = grid_pos_z[pos];
				if ( gridX + x  > grid_points || gridX + x < 0 ) continue;
				if ( gridY + y  > grid_points || gridY + y < 0 ) continue;
				if ( gridZ + z  > grid_points || gridZ + z < 0 ) continue;
				numeric::xyzVector< core::Real > gridpos( (gridX+x)*grid_reso + minX, (gridY+y)*grid_reso + minY, (gridZ+z)*grid_reso + minZ );
				//numeric::xyzVector< core::Real > dist (gridpos - atm_i.xyz());
				numeric::xyzVector< core::Real > dist (gridpos - atmi_xyz_);
				double d2 = dist.length_squared();
				core::Real atm = C*exp(-k*d2);
				rho_cartesian( gridX + x, gridY + y, gridZ + z) += atm;
			}
		}
	}

	if ( basic::options::option[ basic::options::OptionKeys::edensity::debug ]() ) {
		core::scoring::electron_density::ElectronDensity(rho_cartesian,grid_reso, numeric::xyzVector< core::Real >(0,0,0), false).writeMRC( "rho_calc.mrc" );
	}

	for ( int i=1; i<= rc.u1(); ++i ) {
		for ( int j=1; j<= phic.u1(); ++j ) {
			for ( int k=1; k<= zc.u1(); ++k ) {
				core::Real r = rc(i);
				core::Real phi = phic(j);

				//TODO: Resolve this minZZ issue
				core::Real x = cos(phi)*r;
				core::Real y = sin(phi)*r;
				core::Real z = zc(k) + minZZ;

				int xindex = int(( x - minX )/grid_reso + 0.5 );
				int yindex = int(( y - minY )/grid_reso + 0.5 );
				int zindex = int(( z - minZ )/grid_reso + 0.5 );

				if ( xindex > grid_points ||
						yindex > grid_points ||
						zindex > grid_points ||
						xindex < 1 ||
						yindex < 1 ||
						zindex < 1
						) { rho_cylindrical_(i,j,k) = 0.0;}
				else {
					rho_cylindrical_(i,j,k) = rho_cartesian(xindex,yindex,zindex);
				}
			}
		}
	}
}


core::Size
FiberDiffractionEnergyDens::version() const
{
	return 1; // Initial versioning
}


void gnl_R_qfht(
	ObjexxFCL::FArray3D< float > & fourier_in,
	utility::vector0 < utility::vector0 < int > >::iterator & nvals,
	core::Size & lmax,
	core::Real & max_r_value,
	double & qfht_K1,
	double & qfht_K2,
	ObjexxFCL::FArray3D< std::complex<float> > & Gnl
)
{
	Size max_n_val(0);
	for ( Size i=0; i<= lmax; ++i ) {
		if ( nvals[i].size() > max_n_val ) max_n_val = nvals[i].size();
	}
	Size total_rvals ( fourier_in.size1() );
	ObjexxFCL::FArray3D< std::complex<float> > fourier_out(total_rvals, lmax+1, max_n_val );

	//This is  NFFT part /////////////////////////////////
	ft_nfft(fourier_in, nvals, lmax, total_rvals, fourier_out);
	///////////////////////////////////////////////////////////

	double fder[] = {0,0};
	double *f;   /* f traverses the Hankel data array. */
	for ( Size lindex=1; lindex<=lmax; ++lindex ) {
		for ( Size nindex=1; nindex<=nvals[lindex-1].size() ; ++nindex ) {
			int n ( nvals[lindex-1][nindex-1] );
			Hankel *p_hankel;
			p_hankel = hankel_make_input( total_rvals , qfht_K1 , qfht_K2 , max_r_value , 0 , fder , fourier_out , lindex, nindex, n );
			hankel_trans_no_lec( p_hankel );
			f = p_hankel->f;
			for ( core::Size pt = 0 ; pt < p_hankel->n ; pt++ ) {
				std::complex<float> val( *f, *(f + 1)  );
				Gnl(pt+1,lindex,nindex) = val;
				f += 2;
			}
			hankel_free( p_hankel );
		}
	}
}


utility::pointer::shared_ptr< numeric::interpolation::spline::Interpolator >
fit_layer_lines_with_splines(
	ObjexxFCL::FArray1D < float > xvals,
	ObjexxFCL::FArray1D < float > yvals
)
{
	Real minX,minY;
	Real maxX,maxY;

	minX = maxX = xvals(1);
	maxY = minY = yvals(1);
	for ( Size i=2; i<=xvals.size(); i++ ) {
		if ( xvals(i) < minX ) minX = xvals(i);
		if ( xvals(i) > maxX ) maxX = xvals(i);
		if ( yvals(i) < minY ) minY = yvals(i);
		if ( yvals(i) > maxY ) maxY = yvals(i);
	}

	Real deltaR = 1;
	Real deltaD = minY/10;
	numeric::interpolation::spline::SplineGenerator gen( minX-deltaR, minY+deltaD, 0, maxX+deltaR, maxY-deltaD, 0 );
	for ( Size i = 1; i <= xvals.size(); ++i ) {
		gen.add_known_value( xvals(i),yvals(i) );
	}
	utility::pointer::shared_ptr< numeric::interpolation::spline::Interpolator > spline_interpolator = gen.get_interpolator();

	return spline_interpolator;
}

} // fiber_diffraction
} // scoring
} // core

