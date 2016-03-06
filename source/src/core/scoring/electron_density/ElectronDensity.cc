// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/electron_density/ElectronDensity.cc
/// @brief  Scoring a structure against an electron density map
/// @author Frank DiMaio

// Unit Headers
#include <core/scoring/electron_density/ElectronDensity.hh>
#include <core/scoring/electron_density/util.hh>

#ifdef WIN32
#define _USE_MATH_DEFINES
#include <math.h>
#ifndef WIN_PYROSETTA
#include <windows.h>
#endif
#endif

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>  // b factors
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/scoring/Energies.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/Atom.hh>
#include <basic/options/option.hh>
#include <core/scoring/electron_density/xray_scattering.hh>
#include <basic/Tracer.hh>

#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.io.hh>
#include <numeric/statistics/functions.hh>
#include <numeric/fourier/FFT.hh>
#include <numeric/fourier/SHT.hh>


#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/resource_manager/ResourceManager.hh>
#include <basic/resource_manager/util.hh>

//Auto Headers
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/id/AtomID.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>
#include <core/scoring/EnergyGraph.hh>

// Utility headers
#include <utility/string_util.hh>
#include <ObjexxFCL/format.hh>

#include <boost/math/special_functions/bessel.hpp>

// C++ headers
#include <fstream>
#include <limits>
#include <complex>

#ifndef WIN32
#include <pthread.h>
#endif

namespace core {
namespace scoring {
namespace electron_density {

/// @details Auto-generated virtual destructor
ElectronDensity::~ElectronDensity() {}

using basic::T;
using basic::Tracer;
static THREAD_LOCAL basic::Tracer TR( "core.scoring.electron_density.ElectronDensity" );

#ifdef GL_GRAPHICS
// protect access to density from viewer thread
pthread_mutex_t density_map_db_mut_ = PTHREAD_MUTEX_INITIALIZER;
#endif

using namespace core;
using namespace basic::options;

///
///  SHORT HELPER FUNCTIONS
inline float d2r(float d) { return (d*M_PI/180.0); }
inline double d2r(double d) { return (d*M_PI/180.0); }
inline float  square(float  x) { return (x*x); }
inline double square(double x) { return (x*x); }


// x mod y, returns z in [-y/2,y/2]
inline int min_mod(int x,int y) {
	int r=x%y; if ( r<-y/2 ) { r+=y; } if ( r>=y/2 ) { r-=y; }
	return r;
}
inline float min_mod(float x,float y) {
	float r=std::fmod(x,y); if ( r<-0.5*y ) { r+=y; } if ( r>=0.5*y ) { r-=y; }
	return r;
}
inline double min_mod(double x,double y) {
	double r=std::fmod(x,y); if ( r<-0.5*y ) { r+=y; } if ( r>=0.5*y ) { r-=y; }
	return r;
}

// Endianness swap
// Only works with aligned 4-byte quantities
static void swap4_aligned(void *v, long ndata) {
	int *data = (int *) v;
	long i;
	for ( i=0; i<ndata; i++ ) {
		int *N;
		N = data + i;
		*N=(((*N>>24)&0xff) | ((*N&0xff)<<24) | ((*N>>8)&0xff00) | ((*N&0xff00)<<8));
	}
}

///
///  ELECTRON DENSITY CLASS DEFINITIONS
///
ElectronDensity& getDensityMap(std::string filename, bool force_reload) {
	if ( basic::resource_manager::ResourceManager::get_instance()->
			has_resource_with_description("electron_density") ) {

		ElectronDensityOP electron_density(
			basic::resource_manager::get_resource< ElectronDensity >(
			"electron_density"));

		return *electron_density;
	} else {
		return getDensityMap_legacy(filename, force_reload);
	}
}


ElectronDensity& getDensityMap_legacy(std::string filename, bool force_reload) {
	static ElectronDensity theDensityMap;

#ifdef GL_GRAPHICS
	pthread_mutex_lock(&density_map_db_mut_);
#endif

	if ( !theDensityMap.isMapLoaded() || force_reload ) {
		// load map from disk
		TR << "Loading Density Map" << std::endl;
		if ( !basic::options::option[ basic::options::OptionKeys::edensity::mapfile ].user() && filename.length()==0 ) {
			TR.Warning << "[ Warning ] No density map specified" << std::endl;
		} else {
			bool map_loaded=false;
			std::string mapfile;

			if ( filename.length()==0 ) {  // use CMD line args
				mapfile = basic::options::option[ basic::options::OptionKeys::edensity::mapfile ]();
				core::Real mapreso = basic::options::option[ basic::options::OptionKeys::edensity::mapreso ]();
				core::Real mapsampling = basic::options::option[ basic::options::OptionKeys::edensity::grid_spacing ]();

				// Initialize ElectronDensity object
				TR << "Loading density map" << mapfile << std::endl;
				map_loaded = theDensityMap.readMRCandResize( mapfile , mapreso , mapsampling );
			} else {
				mapfile = filename;
				core::Real mapreso = basic::options::option[ basic::options::OptionKeys::edensity::mapreso ]();
				core::Real mapsampling = basic::options::option[ basic::options::OptionKeys::edensity::grid_spacing ]();

				TR << "Loading density map" << mapfile << std::endl;
				map_loaded = theDensityMap.readMRCandResize( mapfile , mapreso , mapsampling );
			}
			if ( !map_loaded ) {
				TR << "[ ERROR ] Error loading density map named '" << mapfile << "'" << std::endl;
				utility_exit();
			}
		}
	}

#ifdef GL_GRAPHICS
	pthread_mutex_unlock(&density_map_db_mut_);
#endif

	return theDensityMap;
}


core::Real ElectronDensity::NUM_DERIV_H = 0.01;

/// null constructor
ElectronDensity::ElectronDensity() {
	init();
}


/// @brief "rho_calc" constructor: make a map from a vector of poses
ElectronDensity::ElectronDensity( utility::vector1< core::pose::PoseOP > poses, core::Real reso, core::Real apix ) {
	init();

	core::Size nposes = poses.size();
	this->reso = reso;

	//1  get bounds
	numeric::xyzVector< core::Real > d_min(0,0,0), d_max(0,0,0);
	bool is_set = false;
	const core::Real FLUFF = 10.0; // add a bounding box
	for ( core::Size n=1; n<=nposes; ++n ) {
		core::pose::Pose &pose = *(poses[n]);
		int nres = pose.total_residue();

		for ( int i=1 ; i<=nres; ++i ) {
			conformation::Residue const &rsd_i (pose.residue(i));
			if ( (rsd_i.aa() == core::chemical::aa_vrt) || (scoring_mask_.find(i) != scoring_mask_.end()) ) continue;
			int nheavyatoms = rsd_i.nheavyatoms();
			for ( int j=1 ; j<=nheavyatoms; ++j ) {
				numeric::xyzVector< core::Real > const &xyz_ij = rsd_i.atom(j).xyz();
				if ( !is_set ) {
					d_min = d_max = xyz_ij;
					is_set = true;
				}
				d_min[0] = std::min(d_min[0],xyz_ij[0]); d_min[1] = std::min(d_min[1],xyz_ij[1]); d_min[2] = std::min(d_min[2],xyz_ij[2]);
				d_max[0] = std::max(d_max[0],xyz_ij[0]); d_max[1] = std::max(d_max[1],xyz_ij[1]); d_max[2] = std::max(d_max[2],xyz_ij[2]);
			}
		}
	}

	// figure out our grid
	numeric::xyzVector< core::Real > extent = ( (d_max - d_min) + 2*FLUFF)/apix;
	grid[0] = findSampling5(extent[0], 2);
	grid[1] = findSampling5(extent[1], 2);
	grid[2] = findSampling5(extent[2], 2);
	numeric::xyzVector< core::Real > real_apix;
	real_apix[0] = ( (d_max[0] - d_min[0]) + 2*FLUFF) / ((core::Real)grid[0]);
	real_apix[1] = ( (d_max[1] - d_min[1]) + 2*FLUFF) / ((core::Real)grid[1]);
	real_apix[2] = ( (d_max[2] - d_min[2]) + 2*FLUFF) / ((core::Real)grid[2]);
	TR << "    grid: " << grid[0] << " x " << grid[1] << " x " << grid[2] << std::endl;
	TR << "    real_apix: " << real_apix[0] << " x " << real_apix[1] << " x " << real_apix[2] << std::endl;

	// make fake crystal data
	cellAngles[0] = cellAngles[1] = cellAngles[2] = 90;
	cellDimensions[0] = grid[0]*real_apix[0];
	cellDimensions[1] = grid[1]*real_apix[1];
	cellDimensions[2] = grid[2]*real_apix[2];
	computeCrystParams();
	TR << "    celldim: " << cellDimensions[0] << " x " << cellDimensions[1] << " x " << cellDimensions[2] << std::endl;
	TR << " cellangles: " << cellAngles[0] << " x " << cellAngles[1] << " x " << cellAngles[2] << std::endl;

	// encode input resolution as minimum B factor
	core::Real max_del_grid = std::max( real_apix[0] , real_apix[1] );
	max_del_grid = std::max( max_del_grid , real_apix[2] );
	minimumB = 16*max_del_grid*max_del_grid;
	max_del_grid *= 1.5;
	effectiveB = 16*max_del_grid*max_del_grid;

	minimumB = std::max( minimumB, 4*reso*reso );
	effectiveB = std::max( minimumB, 4*reso*reso );

	TR << "    minimum B   (min B if input pose contains B factors) = " << minimumB << std::endl;
	TR << "    effective B (min B if input pose does not contain B factors) = " << effectiveB << std::endl;

	// find the origin
	numeric::xyzVector< core::Real > frac_dmin = c2f*(d_min - FLUFF);
	origin[0] = std::floor(frac_dmin[0]*grid[0]);
	origin[1] = std::floor(frac_dmin[1]*grid[1]);
	origin[2] = std::floor(frac_dmin[2]*grid[2]);
	TR << "    origin: " << origin[0] << " x " << origin[1] << " x " << origin[2] << std::endl;

	// atom_mask
	core::Real mask_min = 3.0 * sqrt( effectiveB / (2*M_PI*M_PI) );
	ATOM_MASK = mask_min;
	ATOM_MASK_PADDING = 2;

	// 2 rho_calc
	density.dimension(grid[0],grid[1],grid[2]);
	for ( int i=0; i<density.u1()*density.u2()*density.u3(); ++i ) density[i]=0.0;
	for ( Size n=1; n<=nposes; ++n ) {
		bool use_Bs = pose_has_nonzero_Bs( *(poses[n]) );
		if ( !use_Bs ) {
			TR << "Input pose has no nonzero B factors ... setting to " << effectiveB << std::endl;
		}

		// pose->poseCoords
		poseCoords litePose;
		for ( core::Size i = 1; i <= poses[n]->total_residue(); ++i ) {
			core::conformation::Residue const & rsd_i ( poses[n]->residue(i) );
			if ( rsd_i.aa() == core::chemical::aa_vrt ) continue;

			core::Size natoms = rsd_i.nheavyatoms();
			for ( core::Size j = 1; j <= natoms; ++j ) {
				core::chemical::AtomTypeSet const & atom_type_set( rsd_i.atom_type_set() );

				poseCoord coord_j;
				coord_j.x_ = rsd_i.xyz( j );

				if ( use_Bs && poses[n]->pdb_info() ) {
					coord_j.B_ = std::max( poses[n]->pdb_info()->temperature( i, j ), minimumB);
				}

				coord_j.elt_ = atom_type_set[ rsd_i.atom_type_index( j ) ].element();

				litePose.push_back( coord_j );
			}
		}

		ObjexxFCL::FArray3D< double > rhoC, rhoMask;
		calcRhoC( litePose, 0, rhoC, rhoMask, -1, 1e6 );
		for ( int i=0; i<density.u1()*density.u2()*density.u3(); ++i ) density[i] += rhoC[i];
	}
}

void
ElectronDensity::init() {
	isLoaded = false;

	grid = numeric::xyzVector< int >(0,0,0);
	origin = numeric::xyzVector< int >(0,0,0);
	cellDimensions = numeric::xyzVector< float >(1,1,1);
	cellAngles = numeric::xyzVector< float >(90,90,90);
	use_altorigin =  false;

	// command line overrides defaults
	reso = basic::options::option[ basic::options::OptionKeys::edensity::mapreso ]();
	ATOM_MASK = basic::options::option[ basic::options::OptionKeys::edensity::atom_mask ]();
	ATOM_MASK_PADDING = 2;
	CA_MASK = basic::options::option[ basic::options::OptionKeys::edensity::ca_mask ]();
	WINDOW_ = basic::options::option[ basic::options::OptionKeys::edensity::sliding_window ]();
	score_window_context_ = basic::options::option[ basic::options::OptionKeys::edensity::score_sliding_window_context ]();
	remap_symm_ = basic::options::option[ basic::options::OptionKeys::edensity::score_symm_complex ]();
	force_apix_on_map_load_ = basic::options::option[ basic::options::OptionKeys::edensity::force_apix ]();
	nkbins_ = basic::options::option[ basic::options::OptionKeys::edensity::n_kbins ]();

	// use-specified B factor (may be overridden)
	effectiveB = 0;
}

// gradient of density
numeric::xyzVector<core::Real> ElectronDensity::dens_grad (
	numeric::xyzVector<core::Real> const & idxX ) const {
	numeric::xyzVector< core::Real > dx;
	dx[0] = interp_spline( coeff_grad_x, idxX );
	dx[1] = interp_spline( coeff_grad_y, idxX );
	dx[2] = interp_spline( coeff_grad_z, idxX );
	return dx;
}

//////////////////////
// mapSphericalSamples
// - resample the map in spherical slices around a pose
void ElectronDensity::mapSphericalSamples (
	ObjexxFCL::FArray3D< double > &mapShellR,
	core::Size nRsteps, core::Real delR, core::Size B,
	numeric::xyzVector< core::Real > center
) {

	// make sure map is loaded
	if ( !isLoaded ) {
		TR << "![ ERROR ]  ElectronDensity::mapSHT called but no map is loaded!" << std::endl;
		utility_exit();
	}

	numeric::xyzVector< core::Real > cartOffset, idxX;
	core::Real theta, phi, r;
	if ( coeffs_density_.u1()*coeffs_density_.u2()*coeffs_density_.u3() == 0 ) {
		spline_coeffs( density , coeffs_density_ );
	}

	mapShellR.dimension(2*B,2*B,nRsteps);

	// idx_com1: index coord
	numeric::xyzVector<core::Real> idx_com1 ( center[0] , center[1] , center[2] );
	numeric::xyzVector<core::Real> idx_com2 ( center[0] + origin[0] - 1 ,
		center[1] + origin[1] - 1 ,
		center[2] + origin[2] - 1 );

	// resample
	for ( Size r_idx=1; r_idx<=nRsteps; ++r_idx ) {
		r = (core::Real) delR*(r_idx);

		for ( Size th_idx=1; th_idx<=2*B; ++th_idx ) {
			theta = (2.0*th_idx - 1.0) * M_PI / (4.0*B);

			for ( Size phi_idx=1; phi_idx<=2*B; ++phi_idx ) {
				phi = (2.0*phi_idx - 2.0) * M_PI / (2.0*B);
				// reverse X/Y -- needed for spharm xform
				cartOffset[1] = r*cos(phi)*sin(theta);
				cartOffset[0] = r*sin(phi)*sin(theta);
				cartOffset[2] = r*cos(theta);

				numeric::xyzVector<core::Real> fracX = c2f*cartOffset;
				idxX = numeric::xyzVector<core::Real>( fracX[0]*grid[0],
					fracX[1]*grid[1],
					fracX[2]*grid[2] );

				idxX = idxX + idx_com1;

				mapShellR(phi_idx, th_idx, r_idx) = interp_spline( coeffs_density_ , idxX );
			}
		}
	}
}

/////////////////////////////////////
// Match a centroid pose to the density map, returning correlation coefficient between
//    map and pose.
core::Real ElectronDensity::matchCentroidPose(
	core::pose::Pose const &pose,
	core::conformation::symmetry::SymmetryInfoCOP symmInfo /*=NULL*/,
	bool cacheCCs /* = false */)
{
	using namespace numeric::statistics;

	// make sure map is loaded
	if ( !isLoaded ) {
		TR << "[ ERROR ]  ElectronDensity::matchCentroidPose called but no map is loaded!\n";
		return 0.0;
	}

	ObjexxFCL::FArray3D< double >  rho_calc, inv_rho_mask;
	rho_calc.dimension(density.u1() , density.u2() , density.u3());
	inv_rho_mask.dimension(density.u1() , density.u2() , density.u3());
	for ( int i=0; i<density.u1()*density.u2()*density.u3(); ++i ) {
		rho_calc[i]=0.0;
		inv_rho_mask[i]=1.0;
	}

	int nres = pose.total_residue(); //reses.size();
	numeric::xyzVector< core::Real > cartX, fracX;
	numeric::xyzVector< core::Real > atm_i, atm_j, del_ij;

	// compute RHO_C --> a gaussian at each CA
	core::Real effReso = std::max( 2.4+0.8*reso , reso );
	core::Real k=square(M_PI/effReso);
	core::Real a=33.0;  // treat everything as ALA
	core::Real C=a*pow(k/M_PI,1.5);

	// per-atom derivs
	utility::vector1< numeric::xyzVector<core::Real> >                     atm_idx(nres);
	utility::vector1< utility::vector1< int > >                            rho_dx_pt(nres);
	utility::vector1< utility::vector1< numeric::xyzVector<core::Real> > > rho_dx_mask(nres), rho_dx_atm(nres);

	// symmetry
	bool isSymm = (symmInfo.get() != NULL);
	bool remapSymm = remap_symm_;

	///////////////////////////
	/// 1 COMPUTE RHO_C, MASK
	for ( int i=1 ; i<=nres; ++i ) {
		conformation::Residue const &rsd_i (pose.residue(i));

		// skip non-protein residues & masked reses
		if ( !pose.residue_type(i).is_protein() ) continue;
		if ( scoring_mask_.find(i) != scoring_mask_.end() ) continue;

		// symm
		if ( isSymm && !symmInfo->bb_is_independent(i) && !remapSymm ) {  // should this be fa_...??
			continue; // only score the independent monomer
		}

		conformation::Atom const &atm_i( rsd_i.atom("CA") );

		cartX = atm_i.xyz(); // - getTransform();
		fracX = c2f*cartX;
		atm_idx[i][0] = pos_mod (fracX[0]*grid[0] - origin[0] + 1 , (double)grid[0]);
		atm_idx[i][1] = pos_mod (fracX[1]*grid[1] - origin[1] + 1 , (double)grid[1]);
		atm_idx[i][2] = pos_mod (fracX[2]*grid[2] - origin[2] + 1 , (double)grid[2]);


		for ( int z=1; z<=density.u3(); ++z ) {
			atm_j[2] = z;
			del_ij[2] = (atm_idx[i][2] - atm_j[2]) / grid[2];
			// wrap-around??
			if ( del_ij[2] > 0.5 ) del_ij[2]-=1.0;
			if ( del_ij[2] < -0.5 ) del_ij[2]+=1.0;

			del_ij[0] = del_ij[1] = 0.0;
			if ( (f2c*del_ij).length_squared() > (CA_MASK+ATOM_MASK_PADDING)*(CA_MASK+ATOM_MASK_PADDING) ) continue;

			for ( int y=1; y<=density.u2(); ++y ) {
				atm_j[1] = y;

				// early exit?
				del_ij[1] = (atm_idx[i][1] - atm_j[1]) / grid[1] ;
				// wrap-around??
				if ( del_ij[1] > 0.5 ) del_ij[1]-=1.0;
				if ( del_ij[1] < -0.5 ) del_ij[1]+=1.0;
				del_ij[0] = 0.0;
				if ( (f2c*del_ij).length_squared() > (CA_MASK+ATOM_MASK_PADDING)*(CA_MASK+ATOM_MASK_PADDING) ) continue;

				for ( int x=1; x<=density.u1(); ++x ) {
					atm_j[0] = x;

					// early exit?
					del_ij[0] = (atm_idx[i][0] - atm_j[0]) / grid[0];
					// wrap-around??
					if ( del_ij[0] > 0.5 ) del_ij[0]-=1.0;
					if ( del_ij[0] < -0.5 ) del_ij[0]+=1.0;

					numeric::xyzVector< core::Real > cart_del_ij = (f2c*del_ij);  // cartesian offset from (x,y,z) to atom_i
					core::Real d2 = (cart_del_ij).length_squared();

					if ( d2 > (CA_MASK+ATOM_MASK_PADDING)*(CA_MASK+ATOM_MASK_PADDING) )  continue;

					core::Real atm = C*exp(-k*d2);
					core::Real sigmoid_msk = exp( d2 - (ATOM_MASK)*(ATOM_MASK)  );
					core::Real inv_msk = 1/(1+sigmoid_msk);

					rho_calc(x,y,z) += atm;
					inv_rho_mask(x,y,z) *= (1 - inv_msk);

					if ( !cacheCCs )  continue;

					int idx = (z-1)*density.u2()*density.u1() + (y-1)*density.u1() + x-1;
					rho_dx_pt[i].push_back  ( idx );
					rho_dx_atm[i].push_back ( (-2*k*atm)*cart_del_ij );

					core::Real eps_i = (1-inv_msk), inv_eps_i;
					if ( eps_i == 0 ) { // divide-by-zero
						inv_eps_i = sigmoid_msk;
					} else {
						inv_eps_i = 1/eps_i;
					}

					rho_dx_mask[i].push_back( (-2*sigmoid_msk*inv_msk*inv_msk*inv_eps_i)*cart_del_ij );
				}
			}
		}
	}

	//////////////////////////
	/// 2 COMPUTE SUMMARY STATISTICS
	core::Real sumC_i=0, sumO_i=0, sumCO_i=0, vol_i=0, CC_i=0;
	core::Real sumO2_i=0.0, sumC2_i=0.0, varC_i=0, varO_i=0;
	core::Real clc_x, obs_x, eps_x;

	for ( int x=0; x<density.u1()*density.u2()*density.u3(); ++x ) {
		// fetch this point
		clc_x = rho_calc[x];
		obs_x = density[x];
		eps_x = 1-inv_rho_mask[x]; //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal

		// SMOOTHED
		sumCO_i += eps_x*clc_x*obs_x;
		sumO_i  += eps_x*obs_x;
		sumO2_i += eps_x*obs_x*obs_x;
		sumC_i  += eps_x*clc_x;
		sumC2_i += eps_x*clc_x*clc_x;
		vol_i   += eps_x;
	}
	varC_i = (sumC2_i - sumC_i*sumC_i / vol_i );
	varO_i = (sumO2_i - sumO_i*sumO_i / vol_i ) ;
	if ( varC_i == 0 || varO_i == 0 ) {
		CC_i = 0;
	} else {
		CC_i = (sumCO_i - sumC_i*sumO_i/ vol_i) / sqrt( varC_i * varO_i );
	}

	if ( cacheCCs ) {
		CC_cen = CC_i;
	}


	///////////////////////////
	/// 4  CALCULATE PER-CA DERIVATIVES
	if ( ! cacheCCs ) return CC_i;

	//std::map< core::Size , numeric::xyzMatrix< core::Real > > symmRots;
	for ( int i=1 ; i<=nres; ++i ) {
		if ( isSymm && !symmInfo->bb_is_independent(i) && !remapSymm ) {  // should this be fa_...??
			continue; // only score the monomer
		}

		if ( !pose.residue_type(i).is_protein() ) continue;
		if ( scoring_mask_.find(i) != scoring_mask_.end() ) continue;

		numeric::xyzVector< core::Real > dVdx_ij(0,0,0), dOdx_ij(0,0,0), dO2dx_ij(0,0,0), dCOdx_ij(0,0,0), dC2dx_ij(0,0,0);

		utility::vector1< int > const &rho_dx_pt_ij   = rho_dx_pt[i];
		utility::vector1< numeric::xyzVector<core::Real> > const &rho_dx_mask_ij = rho_dx_mask[i];
		utility::vector1< numeric::xyzVector<core::Real> > const &rho_dx_atm_ij  = rho_dx_atm[i];

		int npoints = rho_dx_pt_ij.size();
		for ( int n=1; n<=npoints; ++n ) {
			const int x(rho_dx_pt_ij[n]);
			clc_x = rho_calc[x];
			obs_x = density[x];
			core::Real inv_eps_x = inv_rho_mask[x];

			numeric::xyzVector<double> del_mask = inv_eps_x*rho_dx_mask_ij[n];
			numeric::xyzVector<double> del_rhoc = rho_dx_atm_ij[n];

			dVdx_ij  += del_mask;
			dOdx_ij  += del_mask*obs_x;
			dO2dx_ij += del_mask*obs_x*obs_x;
			dCOdx_ij += del_rhoc*obs_x;
			dC2dx_ij += 2.0*del_rhoc*clc_x;
		}

		// finally compute dCC/dx_ij
		core::Real f = ( sumCO_i - sumC_i*sumO_i / vol_i );
		core::Real g = sqrt ( varO_i * varC_i );

		numeric::xyzVector<core::Real> fprime = dCOdx_ij - 1/(vol_i*vol_i) * ( dOdx_ij*sumC_i*vol_i - sumO_i*sumC_i*dVdx_ij);
		numeric::xyzVector<core::Real> gprime = 0.5 * (
			sqrt(varO_i)/sqrt(varC_i) * ( dC2dx_ij + ( sumC_i*sumC_i*dVdx_ij/(vol_i*vol_i) ) )  +
			sqrt(varC_i)/sqrt(varO_i) * ( dO2dx_ij - ( 1/(vol_i*vol_i) * ( 2*vol_i*sumO_i*dOdx_ij - sumO_i*sumO_i*dVdx_ij ) ) ) );

		dCCdxs_cen[i] = (g*fprime - f*gprime) / (g*g);
	}

	return CC_i;
}


/////////////////////////////////////
/// Match a residue to the density map, returning correlation coefficient between
///    map and pose
core::Real ElectronDensity::matchPose(
	core::pose::Pose const &pose,
	core::conformation::symmetry::SymmetryInfoCOP symmInfo /*=NULL*/,
	bool cacheCCs/*=false*/ )
{
	using namespace numeric::statistics;

	// make sure map is loaded
	if ( !isLoaded ) {
		TR << "[ ERROR ]  ElectronDensity::matchPose called but no map is loaded!\n";
		return 0.0;
	}

	ObjexxFCL::FArray3D< double >  rho_calc, inv_rho_mask;
	rho_calc.dimension(density.u1() , density.u2() , density.u3());
	inv_rho_mask.dimension(density.u1() , density.u2() , density.u3());
	for ( int i=0; i<density.u1()*density.u2()*density.u3(); ++i ) {
		rho_calc[i]=0.0;
		inv_rho_mask[i]=1.0;
	}

	int nres = pose.total_residue(); //reses.size();
	numeric::xyzVector< core::Real > cartX, fracX;
	numeric::xyzVector< core::Real > atm_i, atm_j, del_ij;

	core::Real SC_scaling = basic::options::option[ basic::options::OptionKeys::edensity::sc_scaling ]();

	// per-atom derivs
	utility::vector1< utility::vector1< numeric::xyzVector<core::Real> > >                     atm_idx(nres);
	utility::vector1< utility::vector1< utility::vector1< int > > >                            rho_dx_pt(nres);
	utility::vector1< utility::vector1< utility::vector1< numeric::xyzVector<core::Real> > > > rho_dx_mask(nres), rho_dx_atm(nres);

	// symmetry
	bool isSymm = (symmInfo.get() != NULL);
	bool remapSymm = remap_symm_;

	///////////////////////////
	/// 1 COMPUTE RHO_C, MASK
	for ( int i=1 ; i<=nres; ++i ) {
		conformation::Residue const &rsd_i (pose.residue(i)); //( *reses[i] );

		// skip vrts & masked reses
		if ( rsd_i.aa() == core::chemical::aa_vrt ) continue;
		if ( scoring_mask_.find(i) != scoring_mask_.end() ) continue;

		// symm
		if ( isSymm && !symmInfo->bb_is_independent(i) && !remapSymm ) {
			continue; // only score the independent monomer
		}

		int nheavyatoms = rsd_i.nheavyatoms();
		atm_idx[i].resize(nheavyatoms);
		rho_dx_pt[i].resize(nheavyatoms);
		rho_dx_mask[i].resize(nheavyatoms);
		rho_dx_atm[i].resize(nheavyatoms);

		for ( int j=1 ; j<=nheavyatoms; ++j ) {
			conformation::Atom const &atm_i( rsd_i.atom(j) );

			chemical::AtomTypeSet const & atom_type_set( rsd_i.atom_type_set() );
			std::string elt_i = atom_type_set[ rsd_i.atom_type_index( j ) ].element();
			OneGaussianScattering sig_j = get_A( elt_i );
			core::Real k = sig_j.k( effectiveB );
			core::Real C = sig_j.C( k );

			// sidechain weight
			if ( (Size) j > rsd_i.last_backbone_atom() ) {
				C *= SC_scaling;
			}

			// if this atom's weight is 0 continue
			if ( C < 1e-6 ) continue;

			cartX = atm_i.xyz(); // - getTransform();
			fracX = c2f*cartX;
			atm_idx[i][j][0] = pos_mod (fracX[0]*grid[0] - origin[0] + 1 , (double)grid[0]);
			atm_idx[i][j][1] = pos_mod (fracX[1]*grid[1] - origin[1] + 1 , (double)grid[1]);
			atm_idx[i][j][2] = pos_mod (fracX[2]*grid[2] - origin[2] + 1 , (double)grid[2]);


			for ( int z=1; z<=density.u3(); ++z ) {
				atm_j[2] = z;
				del_ij[2] = (atm_idx[i][j][2] - atm_j[2]) / grid[2];
				// wrap-around??
				if ( del_ij[2] > 0.5 ) del_ij[2]-=1.0;
				if ( del_ij[2] < -0.5 ) del_ij[2]+=1.0;

				del_ij[0] = del_ij[1] = 0.0;
				if ( (f2c*del_ij).length_squared() > (ATOM_MASK+ATOM_MASK_PADDING)*(ATOM_MASK+ATOM_MASK_PADDING) ) continue;

				for ( int y=1; y<=density.u2(); ++y ) {
					atm_j[1] = y;

					// early exit?
					del_ij[1] = (atm_idx[i][j][1] - atm_j[1]) / grid[1] ;
					// wrap-around??
					if ( del_ij[1] > 0.5 ) del_ij[1]-=1.0;
					if ( del_ij[1] < -0.5 ) del_ij[1]+=1.0;
					del_ij[0] = 0.0;
					if ( (f2c*del_ij).length_squared() > (ATOM_MASK+ATOM_MASK_PADDING)*(ATOM_MASK+ATOM_MASK_PADDING) ) continue;

					for ( int x=1; x<=density.u1(); ++x ) {
						atm_j[0] = x;

						// early exit?
						del_ij[0] = (atm_idx[i][j][0] - atm_j[0]) / grid[0];
						// wrap-around??
						if ( del_ij[0] > 0.5 ) del_ij[0]-=1.0;
						if ( del_ij[0] < -0.5 ) del_ij[0]+=1.0;

						numeric::xyzVector< core::Real > cart_del_ij = (f2c*del_ij);  // cartesian offset from (x,y,z) to atom_i
						core::Real d2 = (cart_del_ij).length_squared();

						if ( d2 > (ATOM_MASK+ATOM_MASK_PADDING)*(ATOM_MASK+ATOM_MASK_PADDING) )  continue;

						core::Real atm = C*exp(-k*d2);
						core::Real sigmoid_msk = exp( d2 - (ATOM_MASK)*(ATOM_MASK)  );
						core::Real inv_msk = 1/(1+sigmoid_msk);

						rho_calc(x,y,z) += atm;
						inv_rho_mask(x,y,z) *= (1 - inv_msk);

						if ( ! cacheCCs )  continue;

						int idx = (z-1)*density.u2()*density.u1() + (y-1)*density.u1() + x-1;

						core::Real eps_i = (1-inv_msk), inv_eps_i;
						if ( eps_i == 0 ) { // divide-by-zero
							inv_eps_i = sigmoid_msk;
						} else {
							inv_eps_i = 1/eps_i;
						}

						rho_dx_pt[i][j].push_back  ( idx );
						rho_dx_atm[i][j].push_back ( (-2*k*atm)*cart_del_ij );
						rho_dx_mask[i][j].push_back( (-2*sigmoid_msk*inv_msk*inv_msk*inv_eps_i)*cart_del_ij );
					}
				}
			}
		}
	}

	//////////////////////////
	/// 2 COMPUTE SUMMARY STATISTICS
	core::Real sumC_i=0, sumO_i=0, sumCO_i=0, vol_i=0, CC_i=0;
	core::Real sumO2_i=0.0, sumC2_i=0.0, varC_i=0, varO_i=0;
	core::Real clc_x, obs_x, eps_x;

	for ( int x=0; x<density.u1()*density.u2()*density.u3(); ++x ) {
		// fetch this point
		clc_x = rho_calc[x];
		obs_x = density[x];
		eps_x = 1-inv_rho_mask[x]; //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal

		// SMOOTHED
		sumCO_i += eps_x*clc_x*obs_x;
		sumO_i  += eps_x*obs_x;
		sumO2_i += eps_x*obs_x*obs_x;
		sumC_i  += eps_x*clc_x;
		sumC2_i += eps_x*clc_x*clc_x;
		vol_i   += eps_x;
	}
	varC_i = (sumC2_i - sumC_i*sumC_i / vol_i );
	varO_i = (sumO2_i - sumO_i*sumO_i / vol_i ) ;
	if ( varC_i == 0 || varO_i == 0 ) {
		CC_i = 0;
	} else {
		CC_i = (sumCO_i - sumC_i*sumO_i/ vol_i) / sqrt( varC_i * varO_i );
	}

	if ( cacheCCs ) {
		CC_aacen = CC_i;
	}


	///////////////////////////
	/// 4  CALCULATE PER-ATOM DERIVATIVES
	if ( ! cacheCCs )  return CC_i;

	//std::map< core::Size , numeric::xyzMatrix< core::Real > > symmRots;
	for ( int i=1 ; i<=nres; ++i ) {
		if ( isSymm && !symmInfo->bb_is_independent(i) && !remapSymm ) {  // should this be fa_...??
			continue; // only score the monomer
		}

		conformation::Residue const &rsd_i (pose.residue(i)); //( *reses[i] );

		if ( rsd_i.aa() == core::chemical::aa_vrt ) continue;
		if ( scoring_mask_.find(i) != scoring_mask_.end() ) continue;

		int nheavyatoms = atm_idx[i].size();
		dCCdxs_aacen[i].resize( nheavyatoms, numeric::xyzVector< core::Real >(0,0,0) );

		for ( int j=1 ; j<=nheavyatoms; ++j ) {
			numeric::xyzVector< core::Real > dVdx_ij(0,0,0), dOdx_ij(0,0,0), dO2dx_ij(0,0,0), dCOdx_ij(0,0,0), dC2dx_ij(0,0,0);

			utility::vector1< int > const &rho_dx_pt_ij   = rho_dx_pt[i][j];
			utility::vector1< numeric::xyzVector<core::Real> > const &rho_dx_mask_ij = rho_dx_mask[i][j];
			utility::vector1< numeric::xyzVector<core::Real> > const &rho_dx_atm_ij  = rho_dx_atm[i][j];

			int npoints = rho_dx_pt_ij.size();
			for ( int n=1; n<=npoints; ++n ) {
				const int x(rho_dx_pt_ij[n]);
				clc_x = rho_calc[x];
				obs_x = density[x];
				core::Real inv_eps_x = inv_rho_mask[x];

				numeric::xyzVector<double> del_mask = inv_eps_x*rho_dx_mask_ij[n];
				numeric::xyzVector<double> del_rhoc = rho_dx_atm_ij[n];

				dVdx_ij  += del_mask;
				dOdx_ij  += del_mask*obs_x;
				dO2dx_ij += del_mask*obs_x*obs_x;
				dCOdx_ij += del_rhoc*obs_x;
				dC2dx_ij += 2.0*del_rhoc*clc_x;
			}

			// finally compute dCC/dx_ij
			core::Real f = ( sumCO_i - sumC_i*sumO_i / vol_i );
			core::Real g = sqrt ( varO_i * varC_i );

			numeric::xyzVector<core::Real> fprime = dCOdx_ij - 1/(vol_i*vol_i) * ( dOdx_ij*sumC_i*vol_i - sumO_i*sumC_i*dVdx_ij);
			numeric::xyzVector<core::Real> gprime = 0.5 * (
				sqrt(varO_i)/sqrt(varC_i) * ( dC2dx_ij + ( sumC_i*sumC_i*dVdx_ij/(vol_i*vol_i) ) )  +
				sqrt(varC_i)/sqrt(varO_i) * ( dO2dx_ij - ( 1/(vol_i*vol_i) * ( 2*vol_i*sumO_i*dOdx_ij - sumO_i*sumO_i*dVdx_ij ) ) ) );

			dCCdxs_aacen[i][j] = (g*fprime - f*gprime) / (g*g);
		}
	}

	return CC_i;
}

void
ElectronDensity::getResolutionBins(
	core::Size nbuckets, core::Real maxreso, core::Real minreso,
	utility::vector1< core::Real > &resbins,
	utility::vector1< core::Size > &counts,
	bool S2_bin/*=false*/ ) {
	Real min_allowed = sqrt(S2( density.u1()/2, density.u2()/2, density.u3()/2 ));
	Real max_allowed = std::min( sqrt(S2( 1,0,0 )), sqrt(S2( 0,1,0 )) );
	max_allowed = std::min( max_allowed, sqrt(S2( 0,0,1 )) );
	if ( minreso>min_allowed || minreso==0 ) {
		TR.Debug << "Forcing res min to " << 1/min_allowed << std::endl;
		minreso = min_allowed;
	}
	if ( maxreso<max_allowed || maxreso==0 ) {
		TR.Debug << "Forcing res max to " << 1/max_allowed << std::endl;
		maxreso = max_allowed;
	}

	if ( S2_bin ) {
		minreso = minreso*minreso;
		maxreso = maxreso*maxreso;
	}

	// get bins
	Real step = (minreso-maxreso)/nbuckets;
	resbins.resize(nbuckets);
	counts.resize(nbuckets);
	for ( Size i=1; i<=nbuckets; ++i ) {
		counts[i] = 0;
		if ( S2_bin ) {
			resbins[i] = std::sqrt( maxreso + (i-0.5)*step );
		} else {
			resbins[i] =  maxreso + (i-0.5)*step;
		}
	}

	// get bin counts
	//int H,K,L;
	for ( int z=1; z<=(int)density.u3(); ++z ) {
		int H = (z < (int)density.u3()/2) ? z-1 : z-density.u3() - 1;
		for ( int y=1; y<=(int)density.u2(); ++y ) {
			int K = (y < (int)density.u2()/2) ? y-1 : y-density.u2()-1;
			for ( int x=1; x<=(int)density.u1(); ++x ) {
				int L = (x < (int)density.u1()/2) ? x-1 : x-density.u1()-1;
				Real s_i = (S2(H,K,L));
				if ( !S2_bin ) s_i=sqrt(s_i);
				int bucket_i = 1+(int)std::floor( (s_i-maxreso) / step );
				if ( bucket_i > 0 && bucket_i <= (int)nbuckets ) {
					counts[bucket_i]++;
				}
			}
		}
	}
}


void
ElectronDensity::getIntensities(
	ObjexxFCL::FArray3D< std::complex<double> > const &Fdensity,
	core::Size nbuckets,
	core::Real maxreso,
	core::Real minreso,
	utility::vector1< core::Real > &Imap,
	bool S2_bin/*=false*/ ) {

	Real min_allowed = sqrt(S2( density.u1()/2, density.u2()/2, density.u3()/2 ));
	Real max_allowed = std::min( sqrt(S2( 1,0,0 )), sqrt(S2( 0,1,0 )) );
	max_allowed = std::min( max_allowed, sqrt(S2( 0,0,1 )) );
	if ( minreso>min_allowed || minreso==0 ) {
		TR.Debug << "Forcing res min to " << 1/min_allowed << std::endl;
		minreso = min_allowed;
	}
	if ( maxreso<max_allowed || maxreso==0 ) {
		TR.Debug << "Forcing res max to " << 1/max_allowed << std::endl;
		maxreso = max_allowed;
	}

	if ( S2_bin ) {
		minreso = minreso*minreso;
		maxreso = maxreso*maxreso;
	}
	Real step = (minreso-maxreso)/nbuckets;

	utility::vector1< core::Real > sum_I2(nbuckets, 0.0);
	utility::vector1< core::Size > counts(nbuckets, 0);

	//int H,K,L;
	for ( int z=1; z<=(int)density.u3(); ++z ) {
		int H = (z < (int)density.u3()/2) ? z-1 : z-density.u3() - 1;
		for ( int y=1; y<=(int)density.u2(); ++y ) {
			int K = (y < (int)density.u2()/2) ? y-1 : y-density.u2()-1;
			for ( int x=1; x<=(int)density.u1(); ++x ) {
				int L = (x < (int)density.u1()/2) ? x-1 : x-density.u1()-1;
				Real s_i = (S2(H,K,L));
				if ( !S2_bin ) s_i=sqrt(s_i);
				int bucket_i = 1+(int)std::floor( (s_i-maxreso) / step );
				if ( bucket_i > 0 && bucket_i <= (int)nbuckets ) {
					sum_I2[bucket_i] += std::real( Fdensity(x,y,z)*std::conj(Fdensity(x,y,z)) );
					counts[bucket_i]++;
				}
			}
		}
	}
	for ( Size i=1; i<=nbuckets; ++i ) {
		sum_I2[i] /= counts[i];
	}

	smooth_intensities(sum_I2);
	Imap = sum_I2;
}

void
ElectronDensity::smooth_intensities(utility::vector1< core::Real > &Is) const {
	utility::vector1< core::Real > Is_in = Is;

	Is[1] = 0.65*Is_in[1]+0.23*Is_in[2]+0.12*Is_in[3];
	Is[2] = 0.3*Is_in[2]+0.35*Is_in[1]+0.23*Is_in[3]+0.12*Is_in[4];
	for ( uint i = 3; i <= Is.size() - 2; ++i ) {
		Is[i] = 0.3*Is_in[i]+0.23*Is_in[i-1]+0.23*Is_in[i+1]+0.12*Is_in[i-2]+0.12*Is_in[i+2];
	}
	Is[Is.size()-1] = 0.3*Is_in[Is.size()-1]+0.35*Is_in[Is.size()]+0.23*Is_in[Is.size()-2]+0.12*Is_in[Is.size()-3];
	Is[Is.size()] = 0.65*Is_in[Is.size()]+0.23*Is_in[Is.size()-1]+0.12*Is_in[Is.size()-2];
}


/// @brief Compute model-map FSC & errors
void
ElectronDensity::getFSC(
	ObjexxFCL::FArray3D< std::complex<double> > const &Fdensity,
	ObjexxFCL::FArray3D< std::complex<double> > const &Fdensity2,
	core::Size nbuckets, core::Real maxreso, core::Real minreso,
	utility::vector1< core::Real >& FSC,
	bool S2_bin/*=false*/) {
	Real min_allowed = sqrt(S2( density.u1()/2, density.u2()/2, density.u3()/2 ));
	Real max_allowed = std::min( sqrt(S2( 1,0,0 )), sqrt(S2( 0,1,0 )) );
	max_allowed = std::min( max_allowed, sqrt(S2( 0,0,1 )) );
	if ( minreso>min_allowed || minreso==0 ) {
		TR.Debug << "Forcing res min to " << 1/min_allowed << std::endl;
		minreso = min_allowed;
	}
	if ( maxreso<max_allowed || maxreso==0 ) {
		TR.Debug << "Forcing res max to " << 1/max_allowed << std::endl;
		maxreso = max_allowed;
	}

	if ( S2_bin ) {
		minreso = minreso*minreso;
		maxreso = maxreso*maxreso;
	}
	Real step = (minreso-maxreso)/nbuckets;

	runtime_assert( Fdensity.u1()==Fdensity2.u1() && Fdensity.u2()==Fdensity2.u2() && Fdensity.u3()==Fdensity2.u3() );

	// correl
	utility::vector1<core::Real> num(nbuckets, 0.0), denom1(nbuckets, 0.0), denom2(nbuckets, 0.0), counter(nbuckets, 0.0);
	FSC.resize(nbuckets);

	for ( Size i=1; i<=nbuckets; ++i ) {
		FSC[i] = 0;
	}

	// PASS 1: FSC, phase error
	//int H,K,L;
	for ( int z=1; z<=(int)density.u3(); ++z ) {
		int H = (z < (int)density.u3()/2) ? z-1 : z-density.u3() - 1;
		for ( int y=1; y<=(int)density.u2(); ++y ) {
			int K = (y < (int)density.u2()/2) ? y-1 : y-density.u2()-1;
			for ( int x=1; x<=(int)density.u1()/2; ++x ) {
				int L = x-1;
				Real s_i = (S2(H,K,L));
				if ( !S2_bin ) s_i=sqrt(s_i);
				int bucket_i = 1+(int)std::floor( (s_i-maxreso) / step );
				if ( bucket_i > 0 && bucket_i <= (int)nbuckets ) {
					Real num_i = std::real( Fdensity2(x,y,z) * std::conj( Fdensity(x,y,z) ) );

					num[bucket_i] += num_i;    // |E_1|*|E_2|*cos(a_1-a_2)
					denom1[bucket_i] += std::abs( Fdensity(x,y,z) ) * std::abs( Fdensity(x,y,z) );
					denom2[bucket_i] += std::abs( Fdensity2(x,y,z) ) * std::abs( Fdensity2(x,y,z) );
				}
			}
		}
	}

	for ( Size i=1; i<=nbuckets; ++i ) {
		denom1[i] = sqrt(denom1[i]);
		denom2[i] = sqrt(denom2[i]);
		FSC[i] = num[i] / (denom1[i]*denom2[i]);
	}
}

void
ElectronDensity::getPhaseError(
	ObjexxFCL::FArray3D< std::complex<double> > const &Fdensity,
	ObjexxFCL::FArray3D< std::complex<double> > const &Fdensity2,
	core::Size nbuckets, core::Real maxreso, core::Real minreso,
	utility::vector1< core::Real >& phaseError,
	bool S2_bin/*=false*/) {
	Real min_allowed = sqrt(S2( density.u1()/2, density.u2()/2, density.u3()/2 ));
	Real max_allowed = std::min( sqrt(S2( 1,0,0 )), sqrt(S2( 0,1,0 )) );
	max_allowed = std::min( max_allowed, sqrt(S2( 0,0,1 )) );
	if ( minreso>min_allowed || minreso==0 ) {
		TR.Debug << "Forcing res min to " << 1/min_allowed << std::endl;
		minreso = min_allowed;
	}
	if ( maxreso<max_allowed || maxreso==0 ) {
		TR.Debug << "Forcing res max to " << 1/max_allowed << std::endl;
		maxreso = max_allowed;
	}

	if ( S2_bin ) {
		minreso = minreso*minreso;
		maxreso = maxreso*maxreso;
	}
	Real step = (minreso-maxreso)/nbuckets;

	runtime_assert( Fdensity.u1()==Fdensity2.u1() && Fdensity.u2()==Fdensity2.u2() && Fdensity.u3()==Fdensity2.u3() );

	// correl
	utility::vector1<core::Real> num(nbuckets, 0.0), denom1(nbuckets, 0.0), denom2(nbuckets, 0.0), counter(nbuckets, 0.0);
	phaseError.resize(nbuckets);

	for ( Size i=1; i<=nbuckets; ++i ) {
		phaseError[i] = 0;
	}

	// PASS 1: FSC, phase error
	//int H,K,L;
	for ( int z=1; z<=(int)density.u3(); ++z ) {
		int H = (z < (int)density.u3()/2) ? z-1 : z-density.u3() - 1;
		for ( int y=1; y<=(int)density.u2(); ++y ) {
			int K = (y < (int)density.u2()/2) ? y-1 : y-density.u2()-1;
			for ( int x=1; x<=(int)density.u1()/2; ++x ) {
				int L = x-1;
				Real s_i = (S2(H,K,L));
				if ( !S2_bin ) s_i=sqrt(s_i);
				int bucket_i = 1+(int)std::floor( (s_i-maxreso) / step );
				if ( bucket_i > 0 && bucket_i <= (int)nbuckets ) {
					Real num_i = std::real( Fdensity2(x,y,z) * std::conj( Fdensity(x,y,z) ) );

					Real denom_i = std::abs(Fdensity(x,y,z)) * std::abs( Fdensity2(x,y,z) );
					if ( denom_i <= 1e-6 ) continue;

					Real ratio_i = std::max( std::min( num_i/denom_i,1.0) , -1.0 );
					Real fX = ratio_i;

					if ( H==0 || K==0 || L==0 ) fX = std::abs(fX);

					Real scale = (H==0 || K==0 || L==0) ? 0.5 : 1.0;
					counter[bucket_i] += scale;

					phaseError[bucket_i] += scale*fX;
				}
			}
		}
	}

	for ( Size i=1; i<=nbuckets; ++i ) {
		phaseError[i] = std::abs( phaseError[i]) / counter[i] ;
		phaseError[i] = phaseError[i] * (2-phaseError[i]*phaseError[i]) / (1-phaseError[i]*phaseError[i]);  // approx von mises Kappa
		phaseError[i] = sqrt(1/phaseError[i]);  // convert Kappa to sigma^2
	}
}


void
ElectronDensity::scaleIntensities(
	utility::vector1< core::Real > scale_i,
	core::Real maxreso, core::Real minreso,
	bool S2_bin/*=false*/ ) {
	if ( Fdensity.u1() == 0 ) numeric::fourier::fft3(density, Fdensity);
	Size nbuckets = scale_i.size();

	Real min_allowed = sqrt(S2( density.u1()/2, density.u2()/2, density.u3()/2 ));
	Real max_allowed = std::min( sqrt(S2( 1,0,0 )), sqrt(S2( 0,1,0 )) );
	max_allowed = std::min( max_allowed, sqrt(S2( 0,0,1 )) );
	if ( minreso>min_allowed || minreso==0 ) {
		TR.Debug << "Forcing res min to " << 1/min_allowed << std::endl;
		minreso = min_allowed;
	}
	if ( maxreso<max_allowed || maxreso==0 ) {
		TR.Debug << "Forcing res max to " << 1/max_allowed << std::endl;
		maxreso = max_allowed;
	}

	if ( S2_bin ) {
		minreso = minreso*minreso;
		maxreso = maxreso*maxreso;
	}
	Real step = (minreso-maxreso)/nbuckets;

	//int H,K,L;
	for ( int z=1; z<=(int)density.u3(); ++z ) {
		int H = (z < (int)density.u3()/2) ? z-1 : z-density.u3() - 1;
		for ( int y=1; y<=(int)density.u2(); ++y ) {
			int K = (y < (int)density.u2()/2) ? y-1 : y-density.u2()-1;
			for ( int x=1; x<=(int)density.u1(); ++x ) {
				int L = (x < (int)density.u1()/2) ? x-1 : x-density.u1()-1;

				Real s_i = (S2(H,K,L));
				if ( !S2_bin ) s_i=sqrt(s_i);

				//fpd smooth interpolate between buckets
				Real bucket = 0.5+((s_i-maxreso) / step);
				int bucket_i = (int)std::floor(bucket);
				Real bucket_offset0 = bucket-bucket_i;
				Real bucket_offset1 = 1.0-bucket_offset0;

				// fpd: no longer truncate here, if we want truncation, apply that separately
				if ( bucket_i >= (int)nbuckets ) {
					Fdensity(x,y,z) *= scale_i[nbuckets];
				} else if ( bucket_i <= 0 ) {
					Fdensity(x,y,z) *= scale_i[1];
				} else {
					Fdensity(x,y,z) *= bucket_offset1*scale_i[bucket_i] + bucket_offset0*scale_i[bucket_i+1];
				}

			}
		}
	}
	numeric::fourier::ifft3(Fdensity, density);

	// clear derived data
	density_change_trigger();
}


void
ElectronDensity::reciprocalSpaceFilter( core::Real maxreso, core::Real minreso, core::Real fadewidth ) {
	numeric::fourier::fft3(density, Fdensity);

	//int H,K,L;
	for ( int z=1; z<=(int)density.u3(); ++z ) {
		int H = (z < (int)density.u3()/2) ? z-1 : z-density.u3() - 1;
		for ( int y=1; y<=(int)density.u2(); ++y ) {
			int K = (y < (int)density.u2()/2) ? y-1 : y-density.u2()-1;
			for ( int x=1; x<=(int)density.u1(); ++x ) {
				int L = (x < (int)density.u1()/2) ? x-1 : x-density.u1()-1;

				Real r_i = 1.0 / sqrt(S2(H,K,L));

				Real fade = 1.0;
				if ( r_i < minreso-fadewidth/2.0 || r_i > maxreso+fadewidth/2.0 ) {
					fade = 0.0;
				} else if ( r_i < minreso+fadewidth/2.0 ) {
					Real del = (r_i-(minreso+fadewidth/2.0))/fadewidth;
					fade = (1-del*del);
					fade = fade*fade;
				} else if ( r_i > maxreso-fadewidth/2.0 ) {
					Real del = (r_i-(maxreso-fadewidth/2.0))/fadewidth;
					fade = (del*del-1);
					fade = fade*fade;
				}

				if ( H!=0 || K!=0 || L!=0 ) {
					Fdensity(x,y,z) *= fade;
				}
			}
		}
	}
	numeric::fourier::ifft3(Fdensity, density);

	// clear derived data
	density_change_trigger();
}


core::Real
ElectronDensity::getRSCC( ObjexxFCL::FArray3D< double > const &density2,  ObjexxFCL::FArray3D< double > const &mask) {
	runtime_assert( density.u1()==density2.u1() && density.u2()==density2.u2() && density.u3()==density2.u3() );

	core::Real sumC_i=0, sumO_i=0, sumCO_i=0, vol_i=0, CC_i=0;
	core::Real sumO2_i=0.0, sumC2_i=0.0, varC_i=0, varO_i=0;
	core::Real clc_x, obs_x, mask_x;
	for ( int x=0; x<density.u1()*density.u2()*density.u3(); ++x ) {
		clc_x = density2[x];
		obs_x = density[x];
		mask_x = (mask.u1() == 0) ? 1.0 : mask[x];

		sumCO_i += mask_x*clc_x*obs_x; sumO_i  += mask_x*obs_x;
		sumO2_i += mask_x*obs_x*obs_x; sumC_i  += mask_x*clc_x;
		sumC2_i += mask_x*clc_x*clc_x;
		vol_i   += mask_x;
	}
	varC_i = (sumC2_i - sumC_i*sumC_i / vol_i );
	varO_i = (sumO2_i - sumO_i*sumO_i / vol_i ) ;
	if ( varC_i == 0 || varO_i == 0 ) {
		CC_i = 0;
	} else {
		CC_i = (sumCO_i - sumC_i*sumO_i/ vol_i) / sqrt( varC_i * varO_i );
	}

	return CC_i;
}

core::Real
ElectronDensity::maxNominalRes() {
	Real S = (1/sqrt(3.)) * sqrt(S2( density.u1()/2, density.u2()/2, density.u3()/2 ));
	return 1.0/S;
}


void
ElectronDensity::calcRhoC(
	poseCoords const &pose,
	core::Real highres_limit,
	ObjexxFCL::FArray3D< double > &rhoC,
	ObjexxFCL::FArray3D< double > &mask,
	core::Real fixed_mask_B /* = -1 */,
	core::Real B_upper_limit /* = 600 */,
	core::Real force_mask /*=-1*/ ) {

	// get rho_c
	rhoC.dimension(density.u1() , density.u2() , density.u3());
	mask.dimension(density.u1() , density.u2() , density.u3());
	rhoC = 0.0;
	mask = 0.0;

	bool use_Bs = pose_has_nonzero_Bs( pose );
	if ( !use_Bs ) {
		TR << "Input pose has no nonzero B factors ... setting to baseline" << std::endl;
	}

	utility::vector1< core::Real > per_atom_ks(pose.size(), 0.0);
	utility::vector1< core::Real > per_atom_Cs(pose.size(), 0.0);
	utility::vector1< core::Real > ATOM_MASK_SQS(pose.size(), 0.0);
	poseCoords pose_grid = pose;
	for ( int i=1 ; i<=(int)pose.size(); ++i ) {
		std::string elt_i = pose[i].elt_;
		OneGaussianScattering sig_j = get_A( elt_i );
		per_atom_ks[i] = sig_j.k( pose[i].B_, B_upper_limit );
		if ( use_Bs ) {
			per_atom_ks[i] = std::min ( per_atom_ks[i], 4*M_PI*M_PI/minimumB );
		} else {
			per_atom_ks[i] = std::min ( per_atom_ks[i], 4*M_PI*M_PI/effectiveB );
		}
		per_atom_Cs[i] = sig_j.C( per_atom_ks[i] );

		core::Real C = per_atom_Cs[i], k = per_atom_ks[i];
		if ( C < 1e-6 ) continue;
		ATOM_MASK_SQS[i] = (1.0/k) * ( std::log( C ) - std::log(1e-4) );  // 1e-4 is the density value where the mask goes to 0

		// for b factor minimization, we don't want the mask to change.
		// so we fix the mask at a high B value
		if ( fixed_mask_B>0 ) {
			core::Real mK = sig_j.k( fixed_mask_B, B_upper_limit );
			core::Real mC = sig_j.C( mK );
			ATOM_MASK_SQS[i] = (1.0/mK) * (std::log( mC ) - std::log(1e-4));
		}

		// force mask
		if ( force_mask>0 ) {
			ATOM_MASK_SQS[i] = force_mask*force_mask;
		}

		// TO DO: OPTIONALLY control periodic/nonperiodic boundaries
		//atm_idx[0] = pos_mod (fracX[0]*grid[0] - origin[0] + 1 , (double)grid[0]);
		//atm_idx[1] = pos_mod (fracX[1]*grid[1] - origin[1] + 1 , (double)grid[1]);
		//atm_idx[2] = pos_mod (fracX[2]*grid[2] - origin[2] + 1 , (double)grid[2]);
		numeric::xyzVector< core::Real> cartX = pose[i].x_; // - getTransform();
		numeric::xyzVector< core::Real> fracX = c2f*cartX;
		numeric::xyzVector< core::Real> atm_idx (
			fracX[0]*grid[0] - origin[0] + 1 ,
			fracX[1]*grid[1] - origin[1] + 1 ,
			fracX[2]*grid[2] - origin[2] + 1 );
		pose_grid[i].x_ = atm_idx;
	}

	// PASS 1: rho_calc
	for ( int i=1 ; i<=(int)pose.size(); ++i ) {
		std::cout << "\rrho_calc " << i << " of " << pose.size() << std::flush;

		core::Real C = per_atom_Cs[i], k = per_atom_ks[i];
		if ( C < 1e-6 ) continue;

		// dist cutoff for mask
		core::Real ATOM_MASK_SQ = ATOM_MASK_SQS[i];

		// dist cutoff for density
		core::Real ATOM_DENS_SQ = (1.0/k) * (std::log( C ) - std::log(1e-4));

		numeric::xyzVector< core::Real> atm_j, del_ij, atm_idx;
		atm_idx =  pose_grid[i].x_;

		if ( atm_idx[0]<1.0 || atm_idx[0]>grid[0] ) continue;
		if ( atm_idx[1]<1.0 || atm_idx[1]>grid[1] ) continue;
		if ( atm_idx[2]<1.0 || atm_idx[2]>grid[2] ) continue;

		for ( int z=1; z<=density.u3(); ++z ) {
			atm_j[2] = z;
			del_ij[2] = (atm_idx[2] - atm_j[2]) / grid[2];
			del_ij[0] = del_ij[1] = 0.0;
			if ( (f2c*del_ij).length_squared() > (ATOM_MASK_SQ) ) continue;
			for ( int y=1; y<=density.u2(); ++y ) {
				atm_j[1] = y;
				del_ij[1] = (atm_idx[1] - atm_j[1]) / grid[1] ;
				del_ij[0] = 0.0;
				if ( (f2c*del_ij).length_squared() > (ATOM_MASK_SQ) ) continue;
				for ( int x=1; x<=density.u1(); ++x ) {
					atm_j[0] = x;
					del_ij[0] = (atm_idx[0] - atm_j[0]) / grid[0];
					numeric::xyzVector< core::Real > cart_del_ij = (f2c*del_ij);  // cartesian offset from (x,y,z) to atom_i
					core::Real d2 = (cart_del_ij).length_squared();

					if ( d2 <= (ATOM_MASK_SQ) ) {
						mask(x,y,z) = 1.0; // problem?
						if ( d2 <= (ATOM_DENS_SQ) ) {
							core::Real atm = C*exp(-k*d2);
							rhoC(x,y,z) += atm;
						}
					}
				}
			}
		}
	}
	std::cout << std::endl;

	if ( highres_limit == 0 ) return;

	// bandlimit mask at 'radius'
	ObjexxFCL::FArray3D< std::complex<double> > Fmask;
	numeric::fourier::fft3(mask, Fmask);
	//int H,K,L;
	for ( int z=1; z<=(int)grid[2]; ++z ) {
		int H = (z < (int)grid[2]/2) ? z-1 : z-grid[2] - 1;
		for ( int y=1; y<=(int)grid[1]; ++y ) {
			int K = (y < (int)grid[1]/2) ? y-1 : y-grid[1]-1;
			for ( int x=1; x<=(int)grid[0]; ++x ) {
				int L = (x < (int)grid[0]/2) ? x-1 : x-grid[0]-1;
				core::Real S2c =  S2(H,K,L);

				// exp fade
				core::Real scale = exp(-S2c*(highres_limit*highres_limit));
				Fmask(x,y,z) *= scale;
			}
		}
	}
	numeric::fourier::ifft3(Fmask, mask);
}


/////////////////////////////////////
/// setup oversampled maps for fast density scoring
void ElectronDensity::setup_fastscoring_first_time(core::pose::Pose const &pose) {
	// atom count
	core::Size natms=0,nres=0;
	for ( core::Size i=1; i<=pose.total_residue(); ++i ) {
		if ( pose.residue(i).aa() == core::chemical::aa_vrt ) { continue; } nres++;
		core::conformation::Residue const& rsd_i = pose.residue(i);
		natms += rsd_i.nheavyatoms();
	}
	core::Real scalefactor = ((core::Real)natms/(core::Real)nres);
	setup_fastscoring_first_time(scalefactor);
}


void ElectronDensity::setup_fastscoring_first_time(Real scalefactor) {
	fastgrid = grid;
	fastorigin = origin;

	utility::vector1< core::Real > fastdens_params;
	if ( basic::options::option[ basic::options::OptionKeys::edensity::fastdens_params ].user() ) {
		fastdens_params = basic::options::option[ basic::options::OptionKeys::edensity::fastdens_params ]();
		runtime_assert( fastdens_params.size() == 2 );
	} else {
		fastdens_params.push_back( 0.4 );
		fastdens_params.push_back( 0.6 );
	}

	ObjexxFCL::FArray3D< double > rhoc, fastdens_score_i;
	ObjexxFCL::FArray3D< std::complex<double> > Frhoo, Frhoc;
	rhoc.dimension(fastgrid[0], fastgrid[1], fastgrid[2]);

	fastdens_score.dimension( density.u1(), density.u2(), density.u3(), nkbins_);

	core::Real max_val = 0.0, min_val = 0.0;

	// compute k limits (corresponding to B=0 and B=1000)
	OneGaussianScattering S_C = get_A( "C" );
	OneGaussianScattering S_O = get_A( "O" );
	OneGaussianScattering S_N = get_A( "N" );
	OneGaussianScattering S_S = get_A( "S" );

	kmin_ = S_C.k( 1000 ); // "infinite"
	kmax_ = std::max( std::max( std::max (S_S.k( 0 ), S_O.k( 0 )) , S_N.k( 0 )), S_C.k( 0 ) );
	kmax_ = std::min( kmax_, 4*M_PI*M_PI/(minimumB) );
	TR << "Setting [kmin_,kmax_] to [" << kmin_ << "," << kmax_ << "]" << std::endl;
	if ( nkbins_ > 1 ) {
		kstep_ = (kmax_-kmin_) / (nkbins_-1);
	} else {
		kstep_ = 0.0;
	}

	numeric::fourier::fft3(density, Frhoo);

	for ( uint kbin = nkbins_; kbin >= 1; --kbin ) {
		rhoc = 0.0;
		numeric::xyzVector< core::Real > del_ij;
		core::Real k = (kbin-1)*kstep_ + kmin_;

		if ( nkbins_ == 1 ) {
			OneGaussianScattering S = get_A( "C" );
			numeric::xyzVector< core::Real > del_ij;
			k = std::min ( S.k(0), 4*M_PI*M_PI/effectiveB );
		}

		if ( k>0 ) {
			k = std::min ( k, 4*M_PI*M_PI/minimumB );

			core::Real C = pow(k, 1.5);
			core::Real ATOM_MASK_SQ = (1.0/k) * (std::log( C ) - std::log(1e-4));  // very generous mask (rho<1e-5)
			core::Real VV = voxel_volume();

			for ( int z=1; z<=(int)fastgrid[2]; ++z ) {
				if ( z < (int)fastgrid[2]/2 ) {
					del_ij[2] = ((core::Real)z - 1.0) / fastgrid[2];
				} else {
					del_ij[2] = (((core::Real)z - fastgrid[2] - 1.0)) / fastgrid[2];
				}
				del_ij[0] = del_ij[1] = 0.0;
				if ( (f2c*del_ij).length_squared() > ATOM_MASK_SQ ) continue;  // early exit
				for ( int y=1; y<=(int)fastgrid[1]; ++y ) {
					if ( y < (int)fastgrid[1]/2 ) {
						del_ij[1] = ((core::Real)y - 1.0) / fastgrid[1] ;
					} else {
						del_ij[1] = (((core::Real)y - fastgrid[1] - 1.0)) / fastgrid[1];
					}
					del_ij[0] = 0.0;
					if ( (f2c*del_ij).length_squared() > ATOM_MASK_SQ ) continue;  // early exit
					for ( int x=1; x<=(int)fastgrid[0]; ++x ) {
						if ( x < (int)fastgrid[0]/2 ) {
							del_ij[0] = ((core::Real)x - 1.0) / fastgrid[0];
						} else {
							del_ij[0] = (((core::Real)x - fastgrid[0] - 1.0)) / fastgrid[0];
						}
						numeric::xyzVector< core::Real > cart_del_ij = (f2c*del_ij);
						core::Real d2 = (cart_del_ij).length_squared();
						if ( d2 <= ATOM_MASK_SQ ) {
							rhoc(x,y,z) = C*VV*exp(-k*d2);
						}
					}
				}
			}

			// standardize density & convert to an approximate CC
			// this has been minimally tuned but might be improved
			for ( int i=0; i<fastgrid[0]*fastgrid[1]*fastgrid[2]; ++i ) {
				rhoc[i] = rhoc[i] * (std::pow( k, -fastdens_params[1] ) - fastdens_params[2]);
			}
		}

		// ffts
		numeric::fourier::fft3(rhoc, Frhoc);
		Frhoc(1,1,1) = Frhoo(1,1,1) = 0.0;

		TR << "Bin " << kbin << ":  B(C/N/O/S)=" << S_C.B(k) << " / " << S_N.B(k) << " / " << S_O.B(k) << " / " <<
			S_S.B(k) << "  sum=" << Frhoc(1,1,1) << std::endl;

		// convolute
		Frhoo(1,1,1) = 0.0;
		for ( int i=1; i<=density.u1(); i++ ) {
			for ( int j=1; j<=density.u2(); j++ ) {
				for ( int k=1; k<=density.u3(); k++ ) {
					Frhoc(i,j,k) *= ( Frhoo(i,j,k) );
				}
			}
		}
		numeric::fourier::ifft3( Frhoc , fastdens_score_i );

		// copy to big array
		for ( int i=1; i<=density.u1(); i++ ) {
			for ( int j=1; j<=density.u2(); j++ ) {
				for ( int k=1; k<=density.u3(); k++ ) {
					fastdens_score(i,j,k,kbin) = fastdens_score_i(i,j,k);
				}
			}
		}

		// normalization
		//    [0.5% -- 99.5%] in min temp bin
		if ( kbin == nkbins_ ) {
			std::vector<core::Real> scores_i( density.u1()*density.u2()*density.u3(),0.0 );
			for ( int i=0; i<density.u1()*density.u2()*density.u3(); ++i ) {
				scores_i[i] = fastdens_score_i[i];
			}

			core::Size cutat = 5;
			std::nth_element (scores_i.begin(), scores_i.begin()+cutat-1, scores_i.end());
			std::nth_element (scores_i.begin(), scores_i.end()-cutat, scores_i.end());
			min_val = scores_i[cutat-1];
			max_val = scores_i[scores_i.size()-cutat];
		}

		if ( basic::options::option[ basic::options::OptionKeys::edensity::debug ]() ) {
			core::Real mu=0.5*(max_val+min_val), sigma=0.5*scalefactor*(max_val-min_val);
			for ( int i=1; i<=density.u1(); i++ ) {
				for ( int j=1; j<=density.u2(); j++ ) {
					for ( int k=1; k<=density.u3(); k++ ) {
						fastdens_score_i(i,j,k) = (fastdens_score(i,j,k,kbin)-mu)/sigma;
					}
				}
			}
			std::ostringstream oss; oss << "fastdens" << kbin << ".mrc";
			ElectronDensity(fastdens_score_i,1.0, numeric::xyzVector< core::Real >(0,0,0), false ).writeMRC( oss.str() );
		}
	}

	core::Real mu=0.0, sigma=0.5*scalefactor*(max_val-min_val);
	for ( uint i = 0; i < nkbins_ * density.u1() * density.u2() * density.u3(); ++i ) {
		fastdens_score[i] = (fastdens_score[i]-mu)/sigma;
	}

	// spline coeffs
	ObjexxFCL::FArray4D< double > temp_coeffs;
	spline_coeffs( fastdens_score, temp_coeffs);
	fastdens_score = temp_coeffs;
}


// use a pose to rescale the fastscoring bins so they correspond to RSCC
void ElectronDensity::rescale_fastscoring_temp_bins(core::pose::Pose const &pose, bool initBs) {
	if ( fastdens_score.u1()*fastdens_score.u2()*fastdens_score.u3()*fastdens_score.u4() == 0 ) {
		setup_fastscoring_first_time(pose);
	}

	if ( nkbins_ == 1 ) return;

	// pose->poseCoords
	poseCoords litePose;
	utility::vector1< core::Real> allBs;
	core::conformation::symmetry::SymmetryInfoCOP symm_info;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		core::conformation::symmetry::SymmetricConformation const & SymmConf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation const &> ( pose.conformation()) );
		symm_info = SymmConf.Symmetry_Info();
	}
	for ( uint i = 1; i <= pose.total_residue(); ++i ) {
		if ( !remap_symm_ && symm_info && !symm_info->bb_is_independent( i ) ) continue;
		core::conformation::Residue const & rsd_i ( pose.residue(i) );
		if ( rsd_i.aa() == core::chemical::aa_vrt ) continue;

		core::Size natoms = rsd_i.nheavyatoms();
		for ( uint j = 1; j <= natoms; ++j ) {
			core::chemical::AtomTypeSet const & atom_type_set( rsd_i.atom_type_set() );

			poseCoord coord_j;
			coord_j.x_ = rsd_i.xyz( j );
			coord_j.B_ = 0;
			coord_j.elt_ = atom_type_set[ rsd_i.atom_type_index( j ) ].element();
			litePose.push_back( coord_j );

			core::Real B = 0.0;
			if ( !initBs ) B = pose.pdb_info()->temperature( i, j );
			allBs.push_back( B );
		}
	}

	// get avg K
	core::Real meanK=0.0;
	if ( !initBs ) {
		for ( uint i = 1; i <= litePose.size(); ++i ) {
			std::string &elt_i = litePose[i].elt_;
			OneGaussianScattering sig_j = get_A( elt_i );
			meanK += sig_j.k( allBs[i] );
		}
		meanK /= litePose.size();
	}

	for ( uint RPT=1; RPT<=5; ++RPT ) {
		core::Size ref=1;
		core::Real RSCC_ref=-999;
		utility::vector1< core::Real > RSCC(nkbins_), overlap(nkbins_);
		numeric::xyzVector< core::Real > fracX, idxX;
		for ( uint kbin = nkbins_; kbin >= 1; --kbin ) {
			core::Real k = (kbin-1)*kstep_ + kmin_;
			if ( k==0 ) continue;

			overlap[kbin] = 0.0;
			for ( uint i = 1; i <= litePose.size(); ++i ) {
				std::string &elt_i = litePose[i].elt_;
				OneGaussianScattering sig_j = get_A( elt_i );

				core::Real k_tgt = k;
				if ( !initBs ) k_tgt += sig_j.k( allBs[i] ) - meanK;
				litePose[i].B_ = std::max(sig_j.B(k_tgt),0.01);

				// compute overlap
				fracX = c2f*litePose[i].x_;
				idxX[0] = pos_mod (fracX[0]*fastgrid[0] - fastorigin[0] + 1 , (double)fastgrid[0]);
				idxX[1] = pos_mod (fracX[1]*fastgrid[1] - fastorigin[1] + 1 , (double)fastgrid[1]);
				idxX[2] = pos_mod (fracX[2]*fastgrid[2] - fastorigin[2] + 1 , (double)fastgrid[2]);
				idxX = numeric::xyzVector<core::Real>( fracX[0]*fastgrid[0] - fastorigin[0] + 1,
					fracX[1]*fastgrid[1] - fastorigin[1] + 1,
					fracX[2]*fastgrid[2] - fastorigin[2] + 1);
				core::Real score_i = interp_spline( fastdens_score , kbin, idxX );
				core::Real W = sig_j.a(  ) / 6.0;
				overlap[kbin] += W*score_i;
			}

			// compute RSCC
			ObjexxFCL::FArray3D< double > rhoC, rhoMask;
			calcRhoC( litePose, 0.0, rhoC, rhoMask );
			RSCC[kbin] = getRSCC(rhoC);
			if ( RSCC[kbin] > RSCC_ref ) {
				RSCC_ref = RSCC[kbin];
				ref = kbin;
			}
		}

		// rescale
		for ( int kbin=nkbins_; kbin>=1; --kbin ) {
			core::Real k = (kbin-1)*kstep_ + kmin_;
			if ( k==0 ) continue;

			core::Real scale_i = RSCC[kbin]*overlap[ref] / (RSCC[ref]*overlap[kbin]);
			if ( RPT == 5 ) TR << "BIN" << kbin << " : [" << k << "/" << meanK << "] " << RSCC[kbin] << "  " << overlap[kbin] << "  " << scale_i << std::endl;

			for ( int i=1; i<=density.u1(); i++ ) {
				for ( int j=1; j<=density.u2(); j++ ) {
					for ( int k=1; k<=density.u3(); k++ ) {
						fastdens_score(i,j,k,kbin) *= scale_i;
					}
				}
			}
		}
	}
}


///@brief  Compute symmetric transformations needed by -score_symm_complex
void ElectronDensity::compute_symm_rotations(
	core::pose::Pose const &pose,
	core::conformation::symmetry::SymmetryInfoCOP symmInfo /*=NULL*/
) {
	// symmetry
	bool isSymm = (symmInfo.get() != NULL);
	bool remapSymm = remap_symm_;

	if ( !isSymm || !remapSymm ) return;

	int nsubunits = symmInfo->subunits();
	int nres_per = symmInfo->num_independent_residues();

	// mapping at the (non-virtual) leaves
	// MUST BE DONE EVERY TIME THE POSE CHANGES
	for ( int subunit_i=1 ;  subunit_i<=nsubunits; ++subunit_i ) {
		// symm
		int startRes = (subunit_i-1)*nres_per+1;
		int source_subunit = subunit_i;
		numeric::xyzMatrix< core::Real > R_i = numeric::xyzMatrix<core::Real>::rows(1,0,0, 0,1,0, 0,0,1);

		if ( !symmInfo->bb_is_independent(startRes) ) {
			core::Size sourceRes = symmInfo->bb_follows( startRes );
			source_subunit = symmInfo->subunit_index( sourceRes );
			R_i = numeric::alignVectorSets(
				pose.residue(startRes).atom(1).xyz() - pose.residue(startRes).atom(2).xyz(),
				pose.residue(startRes).atom(3).xyz() - pose.residue(startRes).atom(2).xyz(),
				pose.residue(sourceRes).atom(1).xyz() - pose.residue(sourceRes).atom(2).xyz(),
				pose.residue(sourceRes).atom(3).xyz() - pose.residue(sourceRes).atom(2).xyz());
		}
		// vrt of 0 => from the non-vrt point of view

		utility::vector1<int> mapping_i(nsubunits,0);
		mapping_i[ subunit_i ] = source_subunit;
		symmap[ -subunit_i ] = make_pair( mapping_i, R_i );
	}

	// at the internal VRTs traverse up the tree to get rotations
	// mapping only needs to be calculated once for each fold tree, however,
	//    Rs must still be percolated up the tree.
	// because of that, we'll just recompute the mapping each time;
	//    however, correcting this could save some running time.
	int nres = pose.total_residue();
	int vrtStart = nres_per*nsubunits;
	int nvrts = nres - vrtStart;
	utility::vector1<bool> vrts_mapped( nvrts, false );
	bool fully_mapped = false;
	while ( !fully_mapped ) {
		fully_mapped = true;
		for ( int i=1 ; i<=nvrts; ++i ) {
			int resid = i+vrtStart;
			if ( vrts_mapped[i] ) {
				continue;
			}
			utility::vector1< core::kinematics::Edge > edges_i = pose.fold_tree().get_outgoing_edges(resid);
			int nchildren = edges_i.size();

			if ( nchildren == 0 ) {
				//TR.Debug << "[ WARNING ] VRT (" << resid << ") at leaf ... ignoring" << std::endl;
				symmap[ resid ] = make_pair( utility::vector1<int>(nsubunits,0), numeric::xyzMatrix<core::Real>::rows(1,0,0, 0,1,0, 0,0,1) );
				vrts_mapped[i] = true;
			} else if ( nchildren == 1 ) {
				// if 1 child, steal mapping from child
				int downstream = edges_i[1].stop();
				if ( downstream <= vrtStart || vrts_mapped[downstream-vrtStart] ) {
					if ( downstream <= vrtStart ) {   // do we jump to a subunit?
						symmap[ resid ] = symmap[ -symmInfo->subunit_index( downstream ) ];
					} else {
						symmap[ resid ] = symmap[ downstream ];
					}
					vrts_mapped[i] = true;
				} else {
					fully_mapped = false;
				}
			} else {
				// if >= 2 children, merge
				bool allChildrenMapped = true;
				for ( int j=1; j<=nchildren; ++j ) {
					int downstream = edges_i[j].stop();
					if ( downstream <= vrtStart ) {
						TR.Error << "[ Error ] VRT (" << resid << ") contains multiple jumps to subunits!  Exiting." << std::endl;
						utility_exit();
					}
					allChildrenMapped &= vrts_mapped[ downstream-vrtStart ];
				}

				if ( allChildrenMapped ) {
					int firstChild = edges_i[1].stop();

					// use rot from 1st child
					numeric::xyzMatrix< core::Real > R_i = symmap[ firstChild ].second;

					// the first child keeps its mapping
					utility::vector1<int> mapping_i = symmap[ firstChild ].first;

					// for each subunit of each child
					for ( int j=2; j<=nchildren; ++j ) {
						utility::vector1<int> const &mapping_j = symmap[ edges_i[j].stop() ].first;
						for ( int k=1; k<=nsubunits; ++k ) {
							if ( mapping_j[k] == 0 ) continue;   // subunit k is not under child j

							// find the subunit to which R_i maps k
							numeric::xyzMatrix< core::Real > R_ik = (symmap[ -k ].second) * R_i;  // R_i * R_k
							core::Real bestFit = 999999;
							int bestL=-1;
							for ( int l=1; l<=nsubunits; ++l ) {
								numeric::xyzMatrix< core::Real > R_diff = R_ik - (symmap[ -l ].second);
								core::Real thisErr = R_diff.col_x().length_squared() +  R_diff.col_y().length_squared() +  R_diff.col_z().length_squared();
								if ( thisErr < bestFit ) {
									bestFit = thisErr;
									bestL = l;
								}
							}
							//std::cerr << "at vrt (" << resid << ") ... mapping " << k << " to " << bestL << " / " << bestFit << std::endl;
							mapping_i[ k ] = bestL;
						}
					}

					symmap[ resid ] = make_pair( mapping_i, R_i );
					vrts_mapped[i] = true;
				} else {
					fully_mapped = false;
				}
			}
		}
	}
}

//FPD >>this should really live in pose datacache
void ElectronDensity::clear_dCCdx_res_cache( core::pose::Pose const &pose ) {
	for ( int i=1, iend=pose.total_residue(); i<=iend; ++i ) {
		if ( pose.residue(i).aa() == core::chemical::aa_vrt ) continue;
		core::Size nAtms = pose.residue(i).nheavyatoms();
		dCCdxs_res[i].resize( nAtms );
		std::fill (dCCdxs_res[i].begin(), dCCdxs_res[i].end(), numeric::xyzVector< core::Real >(0.0,0.0,0.0));
	}
}

/// Match a residue to the density map, returning correlation coefficient between
///    map and residue.
/// Caches information about scoring, to be used in derivative computation
core::Real ElectronDensity::matchRes(
	int resid,
	core::conformation::Residue const &rsd,
	core::pose::Pose const &pose,
	core::conformation::symmetry::SymmetryInfoCOP symmInfo /*=NULL*/,
	bool cacheCCs /* = false */
) {
	// make sure map is loaded
	if ( !isLoaded ) {
		TR << "[ ERROR ]  ElectronDensity::matchRes called but no map is loaded!\n";
		return 0.0;
	}

	if ( scoring_mask_.find(resid) != scoring_mask_.end() ) return 0.0;

	// symmetry
	bool isSymm = (symmInfo.get() != NULL);
	bool remapSymm = remap_symm_;

	// symm
	if ( isSymm && !symmInfo->bb_is_independent(resid) && !remapSymm ) return 0.0; // only score monomer

	//// grab atoms to be included in scoring
	////     context atoms - within the window; mask/derivatives computed
	////     neighbor atoms - atoms whose density _may_ fall within the mask; no mask/derivs computed
	int nres = (int)pose.total_residue();
	utility::vector1< conformation::Atom > neighborAtoms, contextAtoms;
	utility::vector1< std::pair<core::Size,core::Size> > neighborAtomIds, contextAtomIds;
	utility::vector1< OneGaussianScattering > neighborAtomAs, contextAtomAs;

	EnergyGraph const & energy_graph( pose.energies().energy_graph() );
	std::set< core::Size > neighborResids;

	Size HALFWINDOW = WINDOW_/2;
	Size win_start= resid-HALFWINDOW, win_stop = resid+HALFWINDOW;
	for ( int i=HALFWINDOW-1; i>=0; --i ) {
		if ( (resid-i)>1 && pose.fold_tree().is_cutpoint( resid-i-1 ) ) win_start = resid-i;
		if ( (resid+i)<=nres && pose.fold_tree().is_cutpoint( resid+i ) ) win_stop = resid+i;
	}
	win_start = std::max(1, (int)win_start);
	win_stop = std::min((int)win_stop, nres);

	// 1: grab context atom ids
	//    only extend left/right to the next chainbreak
	for ( Size i=win_start; i<=win_stop; ++i ) {
		core::conformation::Residue const &rsd_i( pose.residue(i) );

		// calc neighbor residues
		if ( score_window_context_ ) {
			for ( graph::Graph::EdgeListConstIter
					iru  = energy_graph.get_node(i)->const_edge_list_begin(),
					irue = energy_graph.get_node(i)->const_edge_list_end();
					iru != irue; ++iru ) {
				EnergyEdge const * edge( static_cast< EnergyEdge const *> (*iru) );
				Size const e1( edge->get_first_node_ind() );
				Size const e2( edge->get_second_node_ind() );
				Size const j = (e1==i) ? e2 : e1;
				if ( j<win_start || j>win_stop ) {
					neighborResids.insert( j );
				}
			}
		}

		if ( i==(Size)resid ) continue;  // already added these atoms
		if ( !rsd_i.is_polymer() ) continue;

		chemical::AtomTypeSet const & atom_type_set( rsd_i.atom_type_set() );

		// BB atoms
		for ( core::Size j=1; j<=rsd_i.last_backbone_atom(); ++j ) {
			std::string elt_i = atom_type_set[ rsd_i.atom_type_index( j ) ].element();
			OneGaussianScattering sig_j = get_A( elt_i );
			contextAtoms.push_back( rsd_i.atom( j ) );
			contextAtomIds.push_back( std::pair<core::Size,core::Size>(i,j) );
			contextAtomAs.push_back( sig_j );
		}
		// SC atoms
		for ( core::Size j=rsd_i.last_backbone_atom()+1; j<=rsd_i.nheavyatoms(); ++j ) {
			std::string elt_i = atom_type_set[ rsd_i.atom_type_index( j ) ].element();
			OneGaussianScattering sig_j = get_A( elt_i );
			neighborAtoms.push_back( rsd_i.atom( j ) );
			neighborAtomIds.push_back( std::pair<core::Size,core::Size>(i,j) );
			neighborAtomAs.push_back( sig_j );
		}
	}

	// 2: grab neighbor atom ids
	// all atoms from neighborgraph
	for ( std::set< core::Size >::iterator it=neighborResids.begin(), end=neighborResids.end(); it!=end; ++it ) {
		core::conformation::Residue const &rsd_i( pose.residue(*it) );
		if ( !rsd_i.is_polymer() ) continue;
		chemical::AtomTypeSet const & atom_type_set( rsd_i.atom_type_set() );
		for ( core::Size j=1; j<=rsd_i.nheavyatoms(); ++j ) {
			std::string elt_i = atom_type_set[ rsd_i.atom_type_index( j ) ].element();
			OneGaussianScattering sig_j = get_A( elt_i );
			neighborAtoms.push_back( rsd_i.atom( j ) );
			neighborAtomIds.push_back( std::pair<core::Size,core::Size>(*it,j) );
			neighborAtomAs.push_back( sig_j );
		}
	}

	int nResAtms = rsd.nheavyatoms(), nContextAtms=contextAtoms.size(), nNeighborAtoms=neighborAtoms.size();
	int nTotalAtms = nResAtms+nContextAtms;

	///////////////////////////
	/// 1 COMPUTE BOUNDING BOX
	numeric::xyzVector< core::Real > idxX_low(MAX_FLT,MAX_FLT,MAX_FLT), idxX_high(-MAX_FLT,-MAX_FLT,-MAX_FLT);
	utility::vector1< numeric::xyzVector< Real > > atmList;

	for ( int j=1; j<=nTotalAtms; ++j ) {
		conformation::Atom const &atom (j<=nResAtms? rsd.atom(j) : contextAtoms[j-nResAtms]);
		numeric::xyzVector< core::Real > cartX, idxX, fracX;

		cartX = atom.xyz(); // - getTransform();
		fracX = c2f*cartX;
		idxX = numeric::xyzVector< core::Real >( fracX[0]*grid[0] - origin[0] + 1 ,
			fracX[1]*grid[1] - origin[1] + 1  ,
			fracX[2]*grid[2] - origin[2] + 1  );
		atmList.push_back( idxX );

		fracX = c2f*( cartX-(ATOM_MASK+2*ATOM_MASK_PADDING) );
		idxX = numeric::xyzVector< core::Real >( fracX[0]*grid[0] - origin[0] + 1  ,
			fracX[1]*grid[1] - origin[1] + 1  ,
			fracX[2]*grid[2] - origin[2] + 1  );
		idxX_low[0] = std::min(idxX[0],idxX_low[0]);
		idxX_low[1] = std::min(idxX[1],idxX_low[1]);
		idxX_low[2] = std::min(idxX[2],idxX_low[2]);

		fracX = c2f*( cartX+(ATOM_MASK+2*ATOM_MASK_PADDING) );
		idxX = numeric::xyzVector< core::Real >(
			fracX[0]*grid[0] - origin[0] + 1  ,
			fracX[1]*grid[1] - origin[1] + 1  ,
			fracX[2]*grid[2] - origin[2] + 1  );
		idxX_high[0] = std::max(idxX[0],idxX_high[0]);
		idxX_high[1] = std::max(idxX[1],idxX_high[1]);
		idxX_high[2] = std::max(idxX[2],idxX_high[2]);
	}

	int firstMaskedAtom = 1;
	int lastMaskedAtom = atmList.size();
	if ( basic::options::option[ basic::options::OptionKeys::edensity::unmask_bb ]() ) {
		firstMaskedAtom = rsd.last_backbone_atom()+1;
	}

	// add in neighbor atoms
	// don't let them modify size of bounding box
	for ( int j=1; j<=nNeighborAtoms; ++j ) {
		conformation::Atom const &atom (neighborAtoms[j]);
		numeric::xyzVector< core::Real > cartX, idxX, fracX;
		cartX = atom.xyz(); // - getTransform();
		fracX = c2f*cartX;
		idxX = numeric::xyzVector< core::Real >(
			fracX[0]*grid[0] - origin[0] + 1 ,
			fracX[1]*grid[1] - origin[1] + 1  ,
			fracX[2]*grid[2] - origin[2] + 1  );
		atmList.push_back( idxX );
	}

	numeric::xyzVector< int > bbox_min, bbox_max, bbox_dims;
	bbox_min[0] = (int)floor(idxX_low[0])-1; bbox_max[0] = (int)ceil(idxX_high[0]); bbox_dims[0] = bbox_max[0]-bbox_min[0];
	bbox_min[1] = (int)floor(idxX_low[1])-1; bbox_max[1] = (int)ceil(idxX_high[1]); bbox_dims[1] = bbox_max[1]-bbox_min[1];
	bbox_min[2] = (int)floor(idxX_low[2])-1; bbox_max[2] = (int)ceil(idxX_high[2]); bbox_dims[2] = bbox_max[2]-bbox_min[2];

	ObjexxFCL::FArray3D< double >  rho_obs, inv_rho_mask, rho_calc_fg, rho_calc_bg;
	rho_obs.dimension(bbox_dims[0],bbox_dims[1],bbox_dims[2]);
	rho_calc_fg.dimension(bbox_dims[0],bbox_dims[1],bbox_dims[2]);
	rho_calc_bg.dimension(bbox_dims[0],bbox_dims[1],bbox_dims[2]);
	inv_rho_mask.dimension(bbox_dims[0],bbox_dims[1],bbox_dims[2]);

	for ( int x=0; x<bbox_dims[0]*bbox_dims[1]*bbox_dims[2]; ++x ) {
		rho_obs[x] = 0.0;
		rho_calc_fg[x] = 0.0;
		rho_calc_bg[x] = 0.0;
		inv_rho_mask[x] = 1.0;
	}

	utility::vector1< numeric::xyzVector<core::Real> >                     atm_idx(nTotalAtms);
	utility::vector1< utility::vector1< int > >                            rho_dx_pt(nTotalAtms);
	utility::vector1< utility::vector1< numeric::xyzVector<core::Real> > > rho_dx_mask(nTotalAtms), rho_dx_atm(nTotalAtms);

	///////////////////////////
	/// 2 COMPUTE RHO_C, MASK
	chemical::AtomTypeSet const & atom_type_set( rsd.atom_type_set() );
	core::Real clc_x, obs_x;
	int mapX,mapY,mapZ;

	for ( int i=1; i<=(int)atmList.size(); ++i ) {
		numeric::xyzVector< core::Real > & atm_i = atmList[i];
		atm_i[0] -= bbox_min[0];
		atm_i[1] -= bbox_min[1];
		atm_i[2] -= bbox_min[2];

		numeric::xyzVector< core::Real > atm_j, del_ij;

		core::Real k,C;
		if ( i <= nResAtms ) {
			std::string elt_i = atom_type_set[ rsd.atom_type_index( i ) ].element();
			OneGaussianScattering sig_j = get_A( elt_i );
			k = sig_j.k( effectiveB );
			C = sig_j.C( k );
		} else if ( i<= lastMaskedAtom ) {
			OneGaussianScattering const &sig_j = contextAtomAs[i-nResAtms];
			k = sig_j.k( effectiveB );
			C = sig_j.C( k );
		} else {
			OneGaussianScattering const &sig_j = neighborAtomAs[i-lastMaskedAtom];
			k = sig_j.k( effectiveB );
			C = sig_j.C( k );
		}

		for ( int z=1; z<=bbox_dims[2]; ++z ) {
			atm_j[2] = z;
			del_ij[2] = (atm_i[2] - atm_j[2]) / grid[2];
			// wrap-around??
			if ( del_ij[2] > 0.5 ) del_ij[2]-=1.0;
			if ( del_ij[2] < -0.5 ) del_ij[2]+=1.0;

			del_ij[0] = del_ij[1] = 0.0;
			if ( (f2c*del_ij).length_squared() > (ATOM_MASK+ATOM_MASK_PADDING)*(ATOM_MASK+ATOM_MASK_PADDING) ) continue;

			mapZ = (z+bbox_min[2]) % grid[2];
			if ( mapZ <= 0 ) mapZ += grid[2];
			if ( mapZ > density.u3() ) continue;

			for ( int y=1; y<=bbox_dims[1]; ++y ) {
				atm_j[1] = y;

				// early exit?
				del_ij[1] = (atm_i[1] - atm_j[1]) / grid[1] ;
				// wrap-around??
				if ( del_ij[1] > 0.5 ) del_ij[1]-=1.0;
				if ( del_ij[1] < -0.5 ) del_ij[1]+=1.0;
				del_ij[0] = 0.0;
				if ( (f2c*del_ij).length_squared() > (ATOM_MASK+ATOM_MASK_PADDING)*(ATOM_MASK+ATOM_MASK_PADDING) ) continue;

				mapY = (y+bbox_min[1]) % grid[1];
				if ( mapY <= 0 ) mapY += grid[1];
				if ( mapY > density.u2() ) continue;

				for ( int x=1; x<=bbox_dims[0]; ++x ) {
					atm_j[0] = x;

					// early exit?
					del_ij[0] = (atm_i[0] - atm_j[0]) / grid[0];
					// wrap-around??
					if ( del_ij[0] > 0.5 ) del_ij[0]-=1.0;
					if ( del_ij[0] < -0.5 ) del_ij[0]+=1.0;

					mapX = (x+bbox_min[0]) % grid[0];
					if ( mapX <= 0 ) mapX += grid[0];
					if ( mapX > density.u1() ) continue;

					numeric::xyzVector< core::Real > cart_del_ij = (f2c*del_ij);  // cartesian offset from atom_i to (x,y,z)
					core::Real d2 = (cart_del_ij).length_squared();
					if ( d2 > (ATOM_MASK+ATOM_MASK_PADDING)*(ATOM_MASK+ATOM_MASK_PADDING) )  continue;

					core::Real atm = C*exp(-k*d2);
					core::Real sigmoid_msk = exp( d2 - (ATOM_MASK)*(ATOM_MASK)  );
					core::Real inv_msk = 1/(1+sigmoid_msk);

					rho_obs(x,y,z) = density(mapX,mapY,mapZ);
					if ( i>=firstMaskedAtom && i<=lastMaskedAtom ) {
						rho_calc_fg(x,y,z) += atm;
						inv_rho_mask(x,y,z) *= (1 - inv_msk);
						if ( ! cacheCCs )  continue;

						int idx = (z-1)*rho_calc_fg.u2()*rho_calc_fg.u1() + (y-1)*rho_calc_fg.u1() + x-1;

						core::Real eps_i = (1-inv_msk), inv_eps_i;
						if ( eps_i == 0 ) { // divide-by-zero
							inv_eps_i = sigmoid_msk;
						} else {
							inv_eps_i = 1/eps_i;
						}

						rho_dx_pt[i].push_back  ( idx );
						rho_dx_atm[i].push_back ( (-2*k*atm)*cart_del_ij );
						rho_dx_mask[i].push_back( (-2*sigmoid_msk*inv_msk*inv_msk*inv_eps_i)*cart_del_ij );

					} else {
						rho_calc_bg(x,y,z) += atm;
					}
				}
			}
		}
	}


	if ( basic::options::option[ basic::options::OptionKeys::edensity::debug ]() && resid == 1 ) {
		ElectronDensity(rho_obs,1.0, numeric::xyzVector< core::Real >(0,0,0), false ).writeMRC( "rho_obs.mrc");
		ElectronDensity(inv_rho_mask,1.0, numeric::xyzVector< core::Real >(0,0,0), false ).writeMRC( "rho_mask.mrc");
		ElectronDensity(rho_calc_bg,1.0, numeric::xyzVector< core::Real >(0,0,0), false ).writeMRC( "rho_calc_bg.mrc");
		ElectronDensity(rho_calc_fg,1.0, numeric::xyzVector< core::Real >(0,0,0), false ).writeMRC( "rho_calc_fg.mrc");
	}

	// dumps map for every residue
	/*if (basic::options::option[ basic::options::OptionKeys::edensity::debug ]() && resid ) {
	ElectronDensity(rho_obs,1.0, numeric::xyzVector< core::Real >(0,0,0), false ).writeMRC( "rho_obs.mrc." + utility::to_string(resid) );
	ElectronDensity(inv_rho_mask,1.0, numeric::xyzVector< core::Real >(0,0,0), false ).writeMRC( "rho_mask.mrc." + utility::to_string(resid) );
	ElectronDensity(rho_calc_bg,1.0, numeric::xyzVector< core::Real >(0,0,0), false ).writeMRC( "rho_calc_bg.mrc." + utility::to_string(resid) );
	ElectronDensity(rho_calc_fg,1.0, numeric::xyzVector< core::Real >(0,0,0), false ).writeMRC( "rho_calc_fg.mrc." + utility::to_string(resid) );
	}*/

	//////////////////////////
	/// 2 COMPUTE SUMMARY STATISTICS
	core::Real sumC_i=0.0, sumO_i=0.0, sumCO_i=0.0, vol_i=0.0, CC_i=0.0;
	core::Real sumO2_i=0.0, sumC2_i=0.0, varC_i=0.0, varO_i=0.0;

	for ( int x=0; x<bbox_dims[0]*bbox_dims[1]*bbox_dims[2]; ++x ) {
		// fetch this point
		clc_x = rho_calc_bg[x] + rho_calc_fg[x];
		obs_x = rho_obs[x];

		core::Real wt = 1-inv_rho_mask[x];
		vol_i   += wt;
		sumC_i  += wt*clc_x;
		sumC2_i += wt*clc_x*clc_x;
		sumO_i  += wt*obs_x;
		sumO2_i += wt*obs_x*obs_x;
		sumCO_i += wt*clc_x*obs_x;
	}
	varC_i = (sumC2_i - sumC_i*sumC_i / vol_i );
	varO_i = (sumO2_i - sumO_i*sumO_i / vol_i ) ;
	if ( varC_i == 0 || varO_i == 0 ) {
		CC_i = 0;
	} else {
		CC_i = (sumCO_i - sumC_i*sumO_i/ vol_i) / sqrt( varC_i * varO_i );
	}

	if ( cacheCCs ) {
		CCs[resid] = CC_i;
	}

	///////////////////////////
	/// 4  CALCULATE PER-ATOM DERIVATIVES
	if ( cacheCCs ) {
		for ( int j=1 ; j<=(int)lastMaskedAtom; ++j ) {
			utility::vector1< int > const &rho_dx_pt_ij   = rho_dx_pt[j];
			utility::vector1< numeric::xyzVector<core::Real> > const &rho_dx_mask_ij = rho_dx_mask[j];
			utility::vector1< numeric::xyzVector<core::Real> > const &rho_dx_atm_ij  = rho_dx_atm[j];

			numeric::xyzVector< core::Real > dVdx_ij(0,0,0), dOdx_ij(0,0,0), dO2dx_ij(0,0,0), dCOdx_ij(0,0,0), dC2dx_ij(0,0,0), dCdx_ij(0,0,0);

			int npoints = rho_dx_pt_ij.size();
			for ( int n=1; n<=npoints; ++n ) {
				const int x(rho_dx_pt_ij[n]);
				clc_x = rho_calc_bg[x] + rho_calc_fg[x];
				obs_x = rho_obs[x];
				core::Real inv_eps_x = inv_rho_mask[x];
				core::Real eps_x = 1-inv_eps_x;

				numeric::xyzVector<double> del_mask = inv_eps_x*rho_dx_mask_ij[n];
				numeric::xyzVector<double> del_rhoc = rho_dx_atm_ij[n];

				dVdx_ij  += del_mask;
				// dC_dx = rhoc.*2.*(1-eps).*(x-xi).*sig_i.*eps_i.^2./(1-eps_i)  +  2 * rho_i .* (x-xi) .* (eps);
				dCdx_ij  += del_mask*clc_x + del_rhoc*eps_x;
				dOdx_ij  += del_mask*obs_x;
				dO2dx_ij += del_mask*obs_x*obs_x;
				dC2dx_ij += clc_x*(del_mask*clc_x + 2.0*del_rhoc*eps_x);
				dCOdx_ij += obs_x*(del_mask*clc_x + del_rhoc*eps_x);
			}

			// finally compute dCC/dx_ij
			core::Real f = ( sumCO_i - sumC_i*sumO_i / vol_i );
			core::Real g = sqrt ( varO_i * varC_i );

			numeric::xyzVector<core::Real> fprime = dCOdx_ij - 1/(vol_i*vol_i) * (
				(dOdx_ij*sumC_i + dCdx_ij*sumO_i)*vol_i - sumO_i*sumC_i*dVdx_ij);
			numeric::xyzVector<core::Real> gprime = 0.5 * (
				sqrt(varO_i)/sqrt(varC_i) * ( dC2dx_ij - ( 1/(vol_i*vol_i) * ( 2*vol_i*sumC_i*dCdx_ij - sumC_i*sumC_i*dVdx_ij ) ) ) +
				sqrt(varC_i)/sqrt(varO_i) * ( dO2dx_ij - ( 1/(vol_i*vol_i) * ( 2*vol_i*sumO_i*dOdx_ij - sumO_i*sumO_i*dVdx_ij ) ) ) );

			std::pair<core::Size,core::Size> const &atomid (j<=nResAtms?
				std::pair<core::Size,core::Size>(resid,j) : contextAtomIds[j-nResAtms]);
			dCCdxs_res[atomid.first][atomid.second] += (g*fprime - f*gprime) / (g*g);
		}
	}

	return( CC_i );
}

/////////////////////////////////////
/// Match a residue to the density map.  Use the fast version of the scoring function
core::Real
ElectronDensity::matchResFast(
	int resid,
	core::conformation::Residue const &rsd,
	core::pose::Pose const &pose,
	core::conformation::symmetry::SymmetryInfoCOP symmInfo /*=NULL*/,
	core::Real sc_scale /*=1.0*/
) {
	// make sure map is loaded
	if ( !isLoaded ) {
		TR << "[ ERROR ]  ElectronDensity::matchResFast called but no map is loaded!\n";
		return 0.0;
	}

	if ( fastdens_score.u1()*fastdens_score.u2()*fastdens_score.u3()*fastdens_score.u4() == 0 ) {
		setup_fastscoring_first_time(pose);
	}

	if ( scoring_mask_.find(resid) != scoring_mask_.end() ) return 0.0;

	// symmetry
	bool isSymm = (symmInfo.get() != NULL);
	bool remapSymm = remap_symm_;

	// symm
	if ( isSymm && !symmInfo->bb_is_independent(resid) && !remapSymm ) return 0.0; // only score monomer

	core::Real score = 0;
	numeric::xyzVector< core::Real > fracX, idxX;
	for ( Size i=1; i<=rsd.nheavyatoms(); ++i ) {
		chemical::AtomTypeSet const & atom_type_set( rsd.atom_type_set() );
		std::string elt_i = atom_type_set[ rsd.atom_type_index( i ) ].element();
		OneGaussianScattering sig_j = get_A( elt_i );
		core::Real B = pose.pdb_info() ? pose.pdb_info()->temperature( rsd.seqpos(), i ) : effectiveB;
		core::Real k = sig_j.k( B );

		core::Real kbin = 1;
		if ( nkbins_>1 ) {
			kbin = (k - kmin_)/kstep_ + 1;
		}
		if ( kbin<1 ) kbin=1;
		if ( kbin>nkbins_ ) kbin=nkbins_;

		fracX = c2f*rsd.atom(i).xyz();

		idxX[0] = pos_mod (fracX[0]*fastgrid[0] - fastorigin[0] + 1 , (double)fastgrid[0]);
		idxX[1] = pos_mod (fracX[1]*fastgrid[1] - fastorigin[1] + 1 , (double)fastgrid[1]);
		idxX[2] = pos_mod (fracX[2]*fastgrid[2] - fastorigin[2] + 1 , (double)fastgrid[2]);
		idxX = numeric::xyzVector<core::Real>( fracX[0]*fastgrid[0] - fastorigin[0] + 1,
			fracX[1]*fastgrid[1] - fastorigin[1] + 1,
			fracX[2]*fastgrid[2] - fastorigin[2] + 1);
		core::Real score_i = interp_spline( fastdens_score , kbin, idxX );
		core::Real W = sig_j.a(  ) / 6.0;
		if ( i>rsd.last_backbone_atom() ) {
			W *= sc_scale;
		}
		score += W*score_i;
	}

	return score;
}

/////////////////////////////////////
/// Match an atom to the density map.  Use the fast version of the scoring function and default atom type
core::Real
ElectronDensity::matchPointFast(
	numeric::xyzVector< core::Real > X
) {
	// make sure map is loaded
	if ( !isLoaded ) {
		TR << "[ ERROR ]  ElectronDensity::matchPointFast called but no map is loaded!\n";
		return 0.0;
	}

	if ( fastdens_score.u1()*fastdens_score.u2()*fastdens_score.u3()*fastdens_score.u4() == 0 ) {
		setup_fastscoring_first_time(0.1);
	}

	numeric::xyzVector<core::Real> fracX = c2f*X;
	numeric::xyzVector< core::Real > idxX;
	core::Real kbin = 1;
	idxX[0] = pos_mod (fracX[0]*fastgrid[0] - fastorigin[0] + 1 , (double)fastgrid[0]);
	idxX[1] = pos_mod (fracX[1]*fastgrid[1] - fastorigin[1] + 1 , (double)fastgrid[1]);
	idxX[2] = pos_mod (fracX[2]*fastgrid[2] - fastorigin[2] + 1 , (double)fastgrid[2]);
	core::Real score_i = interp_spline( fastdens_score , kbin, idxX );

	return score_i;
}


void
ElectronDensity::dCCdx_PointFast(
	numeric::xyzVector< core::Real > X,
	numeric::xyzVector< core::Real > & dCCdX
) {
	// make sure map is loaded
	if ( !isLoaded ) {
		TR << "[ ERROR ]  ElectronDensity::matchResFast called but no map is loaded!\n";
		return;
	}

	if ( fastdens_score.u1()*fastdens_score.u2()*fastdens_score.u3()*fastdens_score.u4() == 0 ) {
		setup_fastscoring_first_time(0.1);
	}

	numeric::xyzVector<core::Real> fracX = c2f*X;
	numeric::xyzVector< core::Real > idxX;
	core::Real kbin = 1;
	numeric::xyzVector<core::Real> dCCdX_grid;
	core::Real dkbin;
	idxX[0] = pos_mod (fracX[0]*fastgrid[0] - fastorigin[0] + 1 , (double)fastgrid[0]);
	idxX[1] = pos_mod (fracX[1]*fastgrid[1] - fastorigin[1] + 1 , (double)fastgrid[1]);
	idxX[2] = pos_mod (fracX[2]*fastgrid[2] - fastorigin[2] + 1 , (double)fastgrid[2]);
	interp_dspline( fastdens_score , idxX, kbin,  dCCdX_grid, dkbin );

	dCCdX[0] = dCCdX_grid[0]*c2f(1,1)*fastgrid[0] + dCCdX_grid[1]*c2f(2,1)*fastgrid[1] + dCCdX_grid[2]*c2f(3,1)*fastgrid[2];
	dCCdX[1] = dCCdX_grid[0]*c2f(1,2)*fastgrid[0] + dCCdX_grid[1]*c2f(2,2)*fastgrid[1] + dCCdX_grid[2]*c2f(3,2)*fastgrid[2];
	dCCdX[2] = dCCdX_grid[0]*c2f(1,3)*fastgrid[0] + dCCdX_grid[1]*c2f(2,3)*fastgrid[1] + dCCdX_grid[2]*c2f(3,3)*fastgrid[2];
}


//  Compute the gradient (fast density score)
void ElectronDensity::dCCdx_fastRes(
	int atmid, int resid,
	numeric::xyzVector<core::Real> const &X,
	core::conformation::Residue const &rsd,
	core::pose::Pose const &pose,
	numeric::xyzVector<core::Real> &dCCdX){

	// make sure map is loaded
	if ( !isLoaded ) {
		TR.Error << "[ ERROR ]  ElectronDensity::dCCdx_fastRes called but no map is loaded!\n";
		return;
	}

	if ( fastdens_score.u1()*fastdens_score.u2()*fastdens_score.u3()*fastdens_score.u4() == 0 ) {
		setup_fastscoring_first_time(pose);
	}

	if ( scoring_mask_.find(resid) != scoring_mask_.end() ) return;
	if ( pose.residue(resid).aa() == core::chemical::aa_vrt ) return;

	numeric::xyzVector<core::Real> fracX = c2f*X;
	numeric::xyzVector<core::Real> idxX;
	idxX[0] = pos_mod (fracX[0]*fastgrid[0] - fastorigin[0] + 1 , (double)fastgrid[0]);
	idxX[1] = pos_mod (fracX[1]*fastgrid[1] - fastorigin[1] + 1 , (double)fastgrid[1]);
	idxX[2] = pos_mod (fracX[2]*fastgrid[2] - fastorigin[2] + 1 , (double)fastgrid[2]);

	chemical::AtomTypeSet const & atom_type_set( rsd.atom_type_set() );
	std::string elt_i = atom_type_set[ rsd.atom_type_index( atmid ) ].element();
	OneGaussianScattering sig_j = get_A( elt_i );
	core::Real B = pose.pdb_info() ? pose.pdb_info()->temperature( rsd.seqpos(), atmid ) : effectiveB;

	core::Real k = sig_j.k( B );
	core::Real kbin = 1;
	if ( nkbins_>1 ) {
		kbin = (k - kmin_)/kstep_ + 1;
	}
	if ( kbin<1 ) kbin=1;
	if ( kbin>nkbins_ ) kbin=nkbins_;

	numeric::xyzVector<core::Real> dCCdX_grid;
	core::Real dkbin;

	core::Real W = sig_j.a(  ) / 6.0;
	interp_dspline( fastdens_score , idxX, kbin,  dCCdX_grid, dkbin );
	dCCdX[0] = W*dCCdX_grid[0]*c2f(1,1)*fastgrid[0] + W*dCCdX_grid[1]*c2f(2,1)*fastgrid[1] + W*dCCdX_grid[2]*c2f(3,1)*fastgrid[2];
	dCCdX[1] = W*dCCdX_grid[0]*c2f(1,2)*fastgrid[0] + W*dCCdX_grid[1]*c2f(2,2)*fastgrid[1] + W*dCCdX_grid[2]*c2f(3,2)*fastgrid[2];
	dCCdX[2] = W*dCCdX_grid[0]*c2f(1,3)*fastgrid[0] + W*dCCdX_grid[1]*c2f(2,3)*fastgrid[1] + W*dCCdX_grid[2]*c2f(3,3)*fastgrid[2];
}

// Compute the gradient (fast density score) w.r.t B factors
Real
ElectronDensity::dCCdB_fastRes(
	int atmid, int resid,
	core::conformation::Residue const &rsd,
	core::pose::Pose const &pose
)
{
	// make sure map is loaded
	if ( !isLoaded ) {
		TR << "[ ERROR ]  ElectronDensity::dCCdB_fast called but no map is loaded!\n";
		return 0;
	}

	if ( fastdens_score.u1()*fastdens_score.u2()*fastdens_score.u3()*fastdens_score.u4() == 0 ) {
		setup_fastscoring_first_time(pose);
	}

	if ( scoring_mask_.find(resid) != scoring_mask_.end() ) return 0;
	if ( pose.residue(resid).aa() == core::chemical::aa_vrt ) return 0;

	numeric::xyzVector<core::Real> const &X = rsd.xyz(atmid);
	numeric::xyzVector<core::Real> fracX = c2f*X;
	numeric::xyzVector<core::Real> idxX;
	idxX[0] = pos_mod (fracX[0]*fastgrid[0] - fastorigin[0] + 1 , (double)fastgrid[0]);
	idxX[1] = pos_mod (fracX[1]*fastgrid[1] - fastorigin[1] + 1 , (double)fastgrid[1]);
	idxX[2] = pos_mod (fracX[2]*fastgrid[2] - fastorigin[2] + 1 , (double)fastgrid[2]);

	chemical::AtomTypeSet const & atom_type_set( rsd.atom_type_set() );
	std::string elt_i = atom_type_set[ rsd.atom_type_index( atmid ) ].element();
	OneGaussianScattering sig_j = get_A( elt_i );
	core::Real B = pose.pdb_info()->temperature( rsd.seqpos(), atmid );

	core::Real k = sig_j.k( B );
	if ( nkbins_ == 1 ) {
		return 0;
	}
	core::Real kbin = (k - kmin_)/kstep_ + 1;
	if ( kbin<1 || kbin>nkbins_ ) {
		return 0;
	}

	numeric::xyzVector<core::Real> dCCdX_grid;
	core::Real dscore_dkbin=0.0, dscore_dk;
	interp_dspline( fastdens_score , idxX, kbin,  dCCdX_grid, dscore_dkbin );
	dscore_dk = dscore_dkbin / kstep_;

	core::Real W = sig_j.a(  ) / 6.0;
	Real dCCdb = W*dscore_dk*sig_j.dk( B );

	return dCCdb;
}

// Compute the gradient (fast density score) w.r.t B factors
void
ElectronDensity::dCCdBs(
	core::pose::Pose const &pose,
	utility::vector1< core::Real>  & dE_dvars,
	ObjexxFCL::FArray3D< double > &maskC
) {
	// 1 get rhoc
	core::scoring::electron_density::poseCoords litePose;

	core::conformation::symmetry::SymmetryInfoCOP symm_info;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		core::conformation::symmetry::SymmetricConformation const & SymmConf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation const &> ( pose.conformation()) );
		symm_info = SymmConf.Symmetry_Info();
	}

	for ( uint i = 1; i <= pose.total_residue(); ++i ) {
		if ( symm_info && !symm_info->bb_is_independent( i ) ) continue;
		core::conformation::Residue const & rsd_i ( pose.residue(i) );
		if ( rsd_i.aa() == core::chemical::aa_vrt ) continue;

		core::Size natoms = rsd_i.nheavyatoms();
		for ( uint j = 1; j <= natoms; ++j ) {
			if ( rsd_i.atom_type(j).is_virtual() ) continue;
			core::chemical::AtomTypeSet const & atom_type_set( rsd_i.atom_type_set() );

			core::scoring::electron_density::poseCoord coord_j;
			coord_j.x_ = rsd_i.xyz( j );
			coord_j.B_ = pose.pdb_info()->temperature( i, j );
			coord_j.elt_ = atom_type_set[ rsd_i.atom_type_index( j ) ].element();
			litePose.push_back( coord_j );
		}
	}
	Size natoms = litePose.size();
	dE_dvars.resize( natoms );

	ObjexxFCL::FArray3D< double > rhoC;
	rhoC.dimension( density.u1(), density.u2(), density.u3() );
	rhoC = 0;

	utility::vector1< utility::vector1< int > > rho_dx_pt(natoms);
	utility::vector1< utility::vector1< core::Real > > rho_dx_bs(natoms);

	utility::vector1< core::Real > per_atom_ks(natoms, 0.0);
	utility::vector1< core::Real > per_atom_dks(natoms, 0.0);
	utility::vector1< core::Real > per_atom_Cs(natoms, 0.0);

	poseCoords pose_grid = litePose;
	for ( int i=1 ; i<=(int)litePose.size(); ++i ) {
		std::string elt_i = litePose[i].elt_;
		OneGaussianScattering sig_j = get_A( elt_i );
		per_atom_ks[i] = sig_j.k( litePose[i].B_ );
		per_atom_ks[i] = std::min ( per_atom_ks[i], 4*M_PI*M_PI/minimumB );
		per_atom_dks[i] = sig_j.dk( litePose[i].B_ );
		per_atom_Cs[i] = sig_j.C( per_atom_ks[i] );

		numeric::xyzVector< core::Real> cartX = litePose[i].x_; // - getTransform();
		numeric::xyzVector< core::Real> fracX = c2f*cartX;
		numeric::xyzVector< core::Real> atm_idx (
			fracX[0]*grid[0] - origin[0] + 1 ,
			fracX[1]*grid[1] - origin[1] + 1 ,
			fracX[2]*grid[2] - origin[2] + 1 );
		pose_grid[i].x_ = atm_idx;
	}

	for ( uint i = 1 ; i <= natoms; ++i ) {
		std::cout << "\rrho_calc " << i << " of " << litePose.size() << std::flush;

		core::Real C = per_atom_Cs[i], k = per_atom_ks[i], dk = per_atom_dks[i];
		if ( C < 1e-6 ) continue;
		core::Real ATOM_MASK_SQ = (1.0/k) * (std::log( C ) - std::log(1e-4));

		numeric::xyzVector< core::Real> atm_j, del_ij, atm_idx;
		atm_idx =  pose_grid[i].x_;

		for ( int z=1; z<=density.u3(); ++z ) {
			atm_j[2] = z;
			del_ij[2] = (atm_idx[2] - atm_j[2]) / grid[2];
			del_ij[0] = del_ij[1] = 0.0;
			if ( (f2c*del_ij).length_squared() > (ATOM_MASK_SQ) ) continue;
			for ( int y=1; y<=density.u2(); ++y ) {
				atm_j[1] = y;
				del_ij[1] = (atm_idx[1] - atm_j[1]) / grid[1] ;
				del_ij[0] = 0.0;
				if ( (f2c*del_ij).length_squared() > (ATOM_MASK_SQ) ) continue;
				for ( int x=1; x<=density.u1(); ++x ) {
					atm_j[0] = x;
					del_ij[0] = (atm_idx[0] - atm_j[0]) / grid[0];
					numeric::xyzVector< core::Real > cart_del_ij = (f2c*del_ij);  // cartesian offset from (x,y,z) to atom_i
					core::Real d2 = (cart_del_ij).length_squared();

					if ( d2 <= ATOM_MASK_SQ ) {
						core::Real atm = C*exp(-k*d2);
						rhoC(x,y,z) += atm;
						int idx = (z-1)*density.u2()*density.u1() + (y-1)*density.u1() + x-1;
						rho_dx_pt[i].push_back ( idx );
						rho_dx_bs[i].push_back ( atm*dk*(1.5/k - d2) );
					}
				}
			}
		}
	}

	// CCs
	core::Real sumC_i=0, sumO_i=0, sumCO_i=0, vol_i=0, /*CC_i=0,*/ sumO2_i=0.0, sumC2_i=0.0, varC_i=0, varO_i=0;
	for ( int x=0; x<density.u1()*density.u2()*density.u3(); ++x ) {
		core::Real clc_x = rhoC[x];
		core::Real obs_x = density[x];
		core::Real eps_x = maskC[x];

		sumCO_i += eps_x*clc_x*obs_x;
		sumO_i  += eps_x*obs_x;
		sumO2_i += eps_x*obs_x*obs_x;
		sumC_i  += eps_x*clc_x;
		sumC2_i += eps_x*clc_x*clc_x;
		vol_i   += eps_x;
	}
	varC_i = (sumC2_i - sumC_i*sumC_i / vol_i );
	varO_i = (sumO2_i - sumO_i*sumO_i / vol_i ) ;

	// dCCdbs
	for ( core::uint i = 1; i <= natoms; ++i ) {
		utility::vector1< int > const &rho_dx_pt_i = rho_dx_pt[i];
		utility::vector1< core::Real > const &rho_dx_bs_i = rho_dx_bs[i];
		core::Size npoints = rho_dx_pt_i.size();

		core::Real dCOdx_i=0, dC2dx_i=0;
		for ( core::uint n = 1; n <= npoints; ++n ) {
			const int x(rho_dx_pt_i[n]);
			core::Real clc_x = rhoC[x];
			core::Real obs_x = density[x];
			core::Real del_rhoc = rho_dx_bs_i[n];
			dCOdx_i += del_rhoc*obs_x;
			dC2dx_i += 2.0*del_rhoc*clc_x;
		}

		// finally compute dCC/dx_ij
		core::Real f = ( sumCO_i - sumC_i*sumO_i / vol_i );
		core::Real g = sqrt ( varO_i * varC_i );

		core::Real fprime = dCOdx_i;
		core::Real gprime = 0.5 * sqrt(varO_i)/sqrt(varC_i) * ( dC2dx_i );

		dE_dvars[i] = (g*fprime - f*gprime) / (g*g);
	}
}

/////////////////////////////////////
//  Compute the gradient (sliding_window density score)
void ElectronDensity::dCCdx_res(
	int atmid,
	int resid,
	numeric::xyzVector< core::Real > const & /*X*/,
	core::conformation::Residue const & /*rsd*/,
	core::pose::Pose const & /*pose*/,
	numeric::xyzVector<core::Real> &dCCdX
)
{
	// make sure map is loaded
	if ( !isLoaded ) {
		TR << "[ ERROR ]  ElectronDensity::dCCdx_res called but no map is loaded!\n";
		dCCdX = numeric::xyzVector<core::Real>(0.0,0.0,0.0);
		exit(1);
	}

	// if we didn't score this residue
	//      (because it was missing density or masked or a symm copy)
	// then don't compute its derivative
	if ( (int)dCCdxs_res[resid].size() < atmid ) {
		dCCdX = numeric::xyzVector<core::Real>(0.0,0.0,0.0);
		return;
	}

	dCCdX = dCCdxs_res[resid][atmid];
}


/////////////////////////////////////
//  Compute the gradient (whole-structure CA density score)
void ElectronDensity::dCCdx_cen( int resid,
	numeric::xyzVector<core::Real> const & /*X*/,
	core::pose::Pose const & /*pose*/,
	numeric::xyzVector<core::Real> &dCCdX ) {

	// make sure map is loaded
	if ( !isLoaded ) {
		TR << "[ ERROR ]  ElectronDensity::dCCdx_cen called but no map is loaded!\n";
		dCCdX = numeric::xyzVector<core::Real>(0.0,0.0,0.0);
		utility_exit();
	}

	dCCdX = dCCdxs_cen[resid];
}


/////////////////////////////////////
//  Compute the gradient (whole-structure allatom density score)
void ElectronDensity::dCCdx_aacen( int atmid, int resid,
	numeric::xyzVector<core::Real> const & /*X*/,
	core::pose::Pose const & /*pose*/,
	numeric::xyzVector<core::Real> &dCCdX ) {
	// make sure map is loaded
	if ( !isLoaded ) {
		TR << "[ ERROR ]  ElectronDensity::dCCdx_aacen called but no map is loaded!\n";
		dCCdX = numeric::xyzVector<core::Real>(0.0,0.0,0.0);
		utility_exit();
	}


	// if we didn't score this residue
	//      (because it was missing density or masked or a symm copy)
	// then don't compute its derivative
	if ( (int)dCCdxs_aacen[resid].size() < atmid ) {
		dCCdX = numeric::xyzVector<core::Real>(0.0,0.0,0.0);
		return;
	}

	dCCdX = dCCdxs_aacen[resid][atmid];
}


// ElectronDensity::readMRC(std::string mapfile)
//      read an MRC/CCP4 density map
bool
ElectronDensity::readMRCandResize(
	std::string mapfile,
	core::Real reso,
	core::Real gridSpacing
) {
	std::ifstream mapin(mapfile.c_str() , std::ios::binary | std::ios::in );
	bool isLoaded(readMRCandResize(mapin, mapfile, reso, gridSpacing));
	mapin.close();

	return isLoaded;
}

// ElectronDensity::readMRC(std::istream mapfile)
//      read an MRC/CCP4 density map
bool
ElectronDensity::readMRCandResize(
	std::istream & mapin,
	std::string mapfile,
	core::Real reso,
	core::Real gridSpacing
) {
	char mapString[4], symData[81];

	int  crs2xyz[3], extent[3], mode, symBytes, grid[3], origin[3];
	int  xyz2crs[3], vol_xsize, vol_ysize, vol_zsize;
	int xIndex, yIndex, zIndex, vol_xySize, coord[3];
	unsigned long long dataOffset, filesize;
	float *rowdata;

	bool swap=false;

	if ( !mapin ) {
		TR << "[ ERROR ]  Error opening MRC map " << mapfile << ".  Not loading map." << std::endl;
		return false;
	}

	if (    !mapin.read(reinterpret_cast <char*> (extent), 3*sizeof(int))
			|| !mapin.read(reinterpret_cast <char*> (&mode), 1*sizeof(int))
			|| !mapin.read(reinterpret_cast <char*> (&origin[0]), 3*sizeof(int))
			|| !mapin.read(reinterpret_cast <char*> (&grid[0]), 3*sizeof(int))
			|| !mapin.read(reinterpret_cast <char*> (&cellDimensions[0]), 3*sizeof(float))
			|| !mapin.read(reinterpret_cast <char*> (&cellAngles[0]), 3*sizeof(float))
			|| !mapin.read(reinterpret_cast <char*> (crs2xyz), 3*sizeof(int)) )  {
		TR << "[ ERROR ]   Improperly formatted line in MRC map.  Not loading map." << std::endl;
		return false;
	}

	// Check the number of bytes used for storing symmetry operators
	mapin.seekg(92, std::ios::beg);
	if ( !mapin.read(reinterpret_cast <char*> (&symBytes), 1*sizeof(int)) ) {
		TR << "[ ERROR ]   Failed reading symmetry bytes record.  Not loading map." << "\n";
		return false;
	}

	// alt: MRC files have floating-point origin info at byte 196
	// read this and try to figure out if it is used
	float altorigin[3];
	mapin.seekg(196, std::ios::beg);
	if ( !mapin.read(reinterpret_cast <char*> (altorigin), 3*sizeof(float)) ) {
		TR << "[ ERROR ]   Improperly formatted line in MRC map.  Not loading map." << std::endl;
		return false;
	}

	// Check for the string "MAP" at byte 208, indicating a CCP4 file.
	mapin.seekg(208, std::ios::beg);
	mapString[3] = '\0';
	if ( !mapin.read(mapString, 3) || (std::string(mapString) != "MAP") ) {
		TR << "[ ERROR ]  'MAP' string missing, not a valid MRC map.  Not loading map." << std::endl;
		return false;
	}
	// Check the file endianness
	if ( mode != 2 ) {
		swap4_aligned(&mode, 1);
		if ( mode != 2 ) {
			TR << "[ ERROR ]   Non-real (32-bit float) data types are unsupported.  Not loading map." << std::endl;
			return false;
		} else {
			swap = true; // enable byte swapping
		}
	}

	// Swap all the information obtained from the header
	if ( swap ) {
		swap4_aligned(extent, 3);
		swap4_aligned(&origin[0], 3);
		swap4_aligned(&altorigin[0], 3);
		swap4_aligned(&grid[0], 3);
		swap4_aligned(&cellDimensions[0], 3);
		swap4_aligned(&cellAngles[0], 3);
		swap4_aligned(crs2xyz, 3);
		swap4_aligned(&symBytes, 1);
	}

	if ( reso != 0 ) TR << " Setting resolution to " << reso << "A" << std::endl;
	else TR << " Setting resolution to AUTO" << std::endl;
	TR << "          atom mask to " << ATOM_MASK << "A" << std::endl;
	TR << "            CA mask to " << CA_MASK << "A" << std::endl;
	TR << " Read density map'" << mapfile << "'" << std::endl;
	TR << "     extent: " << extent[0] << " x " << extent[1] << " x " << extent[2] << std::endl;
	TR << "     origin: " << origin[0] << " x " << origin[1] << " x " << origin[2] << std::endl;
	TR << "  altorigin: " << altorigin[0] << " x " << altorigin[1] << " x " << altorigin[2] << std::endl;
	TR << "       grid: " << grid[0] << " x " << grid[1] << " x " << grid[2] << std::endl;
	TR << "    celldim: " << cellDimensions[0] << " x " << cellDimensions[1] << " x " << cellDimensions[2] << std::endl;
	TR << " cellangles: " << cellAngles[0] << " x " << cellAngles[1] << " x " << cellAngles[2] << std::endl;

	// Check the dataOffset: this fixes the problem caused by files claiming
	// to have symmetry records when they do not.
	mapin.seekg(0, std::ios::end);
	filesize = mapin.tellg();
	dataOffset = filesize - 4L *
		(static_cast<long long>(extent[0]) * static_cast<long long>(extent[1]) * static_cast<long long>(extent[2]));
	if ( dataOffset != static_cast<unsigned long long>(1024 + symBytes) ) {
		if ( dataOffset == static_cast<unsigned long long>(1024) ) {
			// Bogus symmetry record information
			TR.Warning << "[ WARNING ] File contains bogus symmetry record.  Continuing." << std::endl;
			symBytes = 0;
		} else if ( dataOffset < static_cast<unsigned long long>(1024) ) {
			TR.Error << "[ ERROR ] File appears truncated and doesn't match header.  Not loading map." << std::endl;
			return false;
		} else if ( (dataOffset > static_cast<unsigned long long>(1024)) && (dataOffset < (1024*1024)) ) {
			// Fix for loading SPIDER files which are larger than usual
			// In this specific case, we must absolutely trust the symBytes record
			dataOffset = 1024 + symBytes;
			TR.Warning << "[ WARNING ]  File is larger than expected and doesn't match header.  Reading anyway." <<
				std::endl;
		} else {
			TR.Error << "[ ERROR ] File is MUCH larger than expected and doesn't match header.  Not loading map." <<
				std::endl;
			TR.Error << dataOffset  << std::endl;
			return false;
		}
	}

	// Read symmetry records -- organized as 80-byte lines of text.
	utility::vector1< std::string > symList;
	symData[80]='\0';
	if ( symBytes != 0 ) {
		TR << "Symmetry records found:" << std::endl;
		mapin.seekg(1024, std::ios::beg);
		for ( int i = 0; i < symBytes/80; i++ ) {
			mapin.read(symData, 80);
			symList.push_back(symData);
			TR << symData << std::endl;
		}
	} else {
		// no symm info; assume P 1
		symList.push_back("X,  Y,  Z");
	}
	initializeSymmOps( symList );

	// check extent and grid interval counts
	if ( grid[0] == 0 && extent[0] > 0 ) {
		grid[0] = extent[0] - 1;
		TR << "[ WARNING ] Fixed X interval count.  Continuing." << std::endl;
	}
	if ( grid[1] == 0 && extent[1] > 0 ) {
		grid[1] = extent[1] - 1;
		TR << "[ WARNING ]  Fixed Y interval count.  Continuing." << std::endl;
	}
	if ( grid[2] == 0 && extent[2] > 0 ) {
		grid[2] = extent[2] - 1;
		TR << "[ WARNING ]  Fixed Z interval count.  Continuing." << std::endl;
	}

	// Mapping between CCP4 column, row, section and Cartesian x, y, z.
	if ( crs2xyz[0] == 0 && crs2xyz[1] == 0 && crs2xyz[2] == 0 ) {
		TR << "[ WARNING ]  All crs2xyz records are zero." << std::endl;
		TR << "[ WARNING ]  Setting crs2xyz to 1, 2, 3 and continuing." << std::endl;
		crs2xyz[0] = 1; crs2xyz[1] = 2; crs2xyz[2] = 3;
	}

	xyz2crs[crs2xyz[0]-1] = 0; xyz2crs[crs2xyz[1]-1] = 1; xyz2crs[crs2xyz[2]-1] = 2;
	xIndex = xyz2crs[0]; yIndex = xyz2crs[1]; zIndex = xyz2crs[2];

	vol_xsize = extent[xIndex];
	vol_ysize = extent[yIndex];
	vol_zsize = extent[zIndex];
	vol_xySize = vol_xsize * vol_ysize;

	// coord = <col, row, sec>
	// extent = <colSize, rowSize, secSize>
	rowdata = new float[extent[0]];
	mapin.seekg(dataOffset, std::ios::beg);

	// 'alloc' changes ordering of "extent"
	density.dimension( vol_xsize,vol_ysize,vol_zsize );

	for ( coord[2] = 1; coord[2] <= extent[2]; coord[2]++ ) {
		for ( coord[1] = 1; coord[1] <= extent[1]; coord[1]++ ) {
			// Read an entire row of data from the file, then write it into the
			// datablock with the correct slice ordering.
			if ( mapin.eof() ) {
				TR << "[ ERROR ] Unexpected end-of-file. Not loading map." << std::endl;
				delete [] rowdata;
				return false;
			}
			if ( mapin.fail() ) {
				TR << "[ ERROR ] Problem reading the file. Not loading map." << std::endl;
				delete [] rowdata;
				return false;
			}
			if ( !mapin.read( reinterpret_cast< char* >(rowdata), sizeof(float)*extent[0]) ) {
				TR << "[ ERROR ] Error reading data row. Not loading map." << std::endl;
				delete [] rowdata;
				return false;
			}

			for ( coord[0] = 1; coord[0] <= extent[0]; coord[0]++ ) {
				density( coord[xyz2crs[0]], coord[xyz2crs[1]], coord[xyz2crs[2]]) = rowdata[coord[0]-1];
			}
		}
	}

	if ( swap == 1 ) {
		swap4_aligned( &density[0], vol_xySize * vol_zsize);
	}
	delete [] rowdata;

	this->origin[0] = origin[xyz2crs[0]];
	this->origin[1] = origin[xyz2crs[1]];
	this->origin[2] = origin[xyz2crs[2]];

	// grid doesnt seemed to get remapped in ccp4 maps
	this->grid[0] = grid[0];
	this->grid[1] = grid[1];
	this->grid[2] = grid[2];

	///////////////////////////////////
	/// POST PROCESSING
	// expand to unit cell
	this->computeCrystParams();
	TR << " voxel vol.: " << voxel_volume() << std::endl;

	// mrc format maps occasionally specify a real-valued origin in a different spot in the header
	if (  altorigin[0]!=0 &&  altorigin[1]!=0 &&  altorigin[2]!=0 &&
			( altorigin[0] > -10000 && altorigin[0] < 10000) &&
			( altorigin[1] > -10000 && altorigin[1] < 10000) &&
			( altorigin[2] > -10000 && altorigin[2] < 10000)
			) {
		this->origin[0] = altorigin[xyz2crs[0]];
		this->origin[1] = altorigin[xyz2crs[1]];
		this->origin[2] = altorigin[xyz2crs[2]];
		numeric::xyzVector<core::Real> fracX = c2f*(this->origin);
		this->origin = numeric::xyzVector<core::Real>( fracX[0]*grid[0] , fracX[1]*grid[1] , fracX[2]*grid[2] );

		TR << "Using ALTERNATE origin\n";
		TR << "     origin =" << this->origin[0] << " x " << this->origin[1] << " x " << this->origin[2] << std::endl;

		use_altorigin = true;
	} else {
		use_altorigin = false;
	}

	this->expandToUnitCell();

	// optionally, if force_apix_on_map_load_ is specified
	// override input apix with supplied value
	if ( force_apix_on_map_load_ > 0 ) {
		set_voxel_spacing( force_apix_on_map_load_ );
	}

	// resample the map
	if ( gridSpacing > 0 ) {
		this->resize( gridSpacing );
	}

	// Resolution logic is a bit convoluted
	//   We define a "minimum" resolution and an "effective" resolution
	//   The minimum resolution is the absolute minimum resolution component we can represent (grid_sampling/2)
	//   Effective resolution is the "actual" resolution of the data, and comes from -map_reso (if provided)
	//      or is guessed at grid_sampling
	core::Real max_del_grid = std::max( cellDimensions[0]/((double)this->grid[0]) , cellDimensions[1]/((double)this->grid[1]) );
	max_del_grid = std::max( max_del_grid , cellDimensions[2]/((double)this->grid[2]) );
	minimumB = 16*max_del_grid*max_del_grid;

	if ( reso == 0 ) {
		max_del_grid *= 1.5;
	} else {
		max_del_grid = std::max( max_del_grid, reso/2 );
	}

	TR << "Effective resolution = " << 2*max_del_grid << std::endl;
	effectiveB = 16*max_del_grid*max_del_grid;

	// mask width logic is also convoluted
	// make sure mask extends >= 3 carbon stdevs at effectiveB
	core::Real mask_min = 3.0 * sqrt( effectiveB / (2*M_PI*M_PI) );
	if ( basic::options::option[ basic::options::OptionKeys::edensity::atom_mask_min ].user() ) {
		mask_min = basic::options::option[ basic::options::OptionKeys::edensity::atom_mask_min ];
	}

	if ( ATOM_MASK < mask_min ) {
		ATOM_MASK = mask_min;
	}
	mask_min = 3.0 * sqrt(2.0) * std::max( 2.4+1.6*reso , reso ) / M_PI;
	if ( CA_MASK < mask_min ) {
		CA_MASK = mask_min;
	}

	// set map resolution
	this->reso = reso;

	// update/clear derived data
	density_change_trigger();

	// we're done!
	isLoaded = true;
	return isLoaded;
}


// expand density to cover complete unit cell
// maintain origin
void ElectronDensity::expandToUnitCell() {
	numeric::xyzVector< int > extent( density.u1(), density.u2(), density.u3() );

	// if it already covers unit cell do nothing
	if ( grid[0] == extent[0] && grid[1] == extent[1] && grid[2] == extent[2] ) {
		return;
	}

	ObjexxFCL::FArray3D< float > newDensity( grid[0],grid[1],grid[2], 0.0 );

	// copy the block
	int limX=std::min(extent[0],grid[0]),
		limY=std::min(extent[1],grid[1]),
		limZ=std::min(extent[2],grid[2]);
	for ( int x=1; x<=limX; ++x ) {
		for ( int y=1; y<=limY; ++y ) {
			for ( int z=1; z<=limZ; ++z ) {
				newDensity( x,y,z ) = density( x,y,z );
			}
		}
	}

	// apply symmetry
	// why backwards? it is a mystery
	for ( int x=grid[0]; x>=1; --x ) {
		for ( int y=grid[1]; y>=1; --y ) {
			for ( int z=grid[2]; z>=1; --z ) {
				if ( x <= limX && y <= limY && z <= limZ ) {
					continue;
				}

				numeric::xyzVector<core::Real> fracX(
					( (core::Real)x + origin[0] - 1 ) / grid[0],
					( (core::Real)y + origin[1] - 1 ) / grid[1],
					( (core::Real)z + origin[2] - 1 ) / grid[2] );

				for ( int symm_idx=1; symm_idx<=(int)symmOps.size(); symm_idx++ ) {
					numeric::xyzVector<core::Real> SfracX =
						symmOps[symm_idx].get_rotation() * fracX +  symmOps[symm_idx].get_translation();

					// indices of symm copy
					int Sx = pos_mod((int)floor(SfracX[0]*grid[0]+0.5 - origin[0]) , grid[0]) + 1;
					int Sy = pos_mod((int)floor(SfracX[1]*grid[1]+0.5 - origin[1]) , grid[1]) + 1 ;
					int Sz = pos_mod((int)floor(SfracX[2]*grid[2]+0.5 - origin[2]) , grid[2]) + 1 ;

					if ( Sx <= limX && Sy <= limY && Sz <= limZ ) {
						newDensity( x,y,z ) = density(Sx,Sy,Sz);
					}
				}
			}
		}
	}

	if ( basic::options::option[ basic::options::OptionKeys::edensity::debug ]() ) {
		ElectronDensity(density,1.0, numeric::xyzVector< core::Real >(0,0,0), false ).writeMRC( "rho_before_expand.mrc");
		ElectronDensity(newDensity,1.0, numeric::xyzVector< core::Real >(0,0,0), false ).writeMRC( "rho_after_expand.mrc");
	}

	// new map!
	density = newDensity;
}

// parse symmops from ccp4 map header
/// also sets MINMULT(XYZ)
/// eventually replace with cctbx
void ElectronDensity::initializeSymmOps( utility::vector1< std::string > const & symList ) {
	using core::kinematics::RT;

	symmOps.clear();

	if ( symList.size() == 0 ) { // no symminfo in header, assume P 1
		symmOps.push_back( RT( numeric::xyzMatrix< core::Real >::identity(),
			numeric::xyzVector< core::Real >(0.0,0.0,0.0) ) );
	}

	MINMULT[0] = MINMULT[1] = MINMULT[2] = 2;
	for ( int i=1; i<=(int)symList.size(); ++i ) {
		std::string line = symList[i];
		utility::vector1< std::string > rows = utility::string_split(line, ',');
		if ( rows.size() != 3 ) {
			TR.Error << "[ ERROR ] invalid symmop in map file" << std::endl;
			TR.Error << line << std::endl;
			TR.Error << "Setting symmetry to P1 and continuing!" << line << std::endl;

			// should we throw an exception here????  nah, just set symm to P1 and continue
			symmOps.clear();
			symmOps.push_back( RT( numeric::xyzMatrix< core::Real >::identity(),
				numeric::xyzVector< core::Real >(0.0,0.0,0.0) ) );

			return;
		}

		// _REALLY_ simple parser
		//int k;
		for ( int j=1; j<=3; ++j ) {
			int k = rows[j].find('/');
			if ( k != (int)std::string::npos ) {
				// numerator
				int startDenom = k+1;
				float denom = std::atof( &rows[j][startDenom]);

				// make sure this shift corresponds to a point in the map
				core::Size oldMinMult = MINMULT[j-1];
				while ( std::fmod(MINMULT[j-1] , denom) > 1e-6 ) MINMULT[j-1] += oldMinMult;
			}
		}
	}

	return;
}


bool ElectronDensity::writeMRC(std::string mapfilename) {
	std::fstream outx( (mapfilename).c_str() , std::ios::binary | std::ios::out );

	float buff_f;
	int buff_i;
	float buff_vf[3];
	int buff_vi[3];
	int symBytes = 0;

	if ( !outx ) {
		TR.Error << "[ ERROR ]  Error opening MRC map for writing." << std::endl;
		return false;
	}

	// extent
	buff_vi[0] = density.u1(); buff_vi[1] = density.u2(); buff_vi[2] = density.u3();
	outx.write(reinterpret_cast <char*>(buff_vi), sizeof(int)*3);

	// mode
	buff_i = 2;
	outx.write(reinterpret_cast <char*>(&buff_i), sizeof(int));

	// origin
	int ori_int[3];
	ori_int[0] = use_altorigin? 0 : (int)std::floor( origin[0] );
	ori_int[1] = use_altorigin? 0 : (int)std::floor( origin[1] );
	ori_int[2] = use_altorigin? 0 : (int)std::floor( origin[2] );
	outx.write(reinterpret_cast <char*>(ori_int), sizeof(int)*3);

	// grid
	outx.write(reinterpret_cast <char*>(&grid[0]), sizeof(int)*3);

	// cell params
	outx.write(reinterpret_cast <char*>(&cellDimensions), sizeof(float)*3);
	outx.write(reinterpret_cast <char*>(&cellAngles), sizeof(float)*3);

	// crs2xyz
	buff_vi[0] = 1; buff_vi[1] = 2; buff_vi[2] = 3;
	outx.write(reinterpret_cast <char*>(buff_vi), sizeof(int)*3);

	// min, max, mean dens
	buff_vf[0] = -100.0; buff_vf[1] = 100.0; buff_vf[2] = 0.0;
	outx.write(reinterpret_cast <char*>(buff_vf), sizeof(float)*3);

	// 4 bytes junk
	buff_i = 0;
	outx.write(reinterpret_cast <char*>(&buff_i), sizeof(int));

	// symmops (to do!)
	outx.write(reinterpret_cast <char*>(&symBytes), sizeof(int));

	// 104 bytes junk
	buff_i = 0;
	for ( int i=0; i<25; ++i ) {
		outx.write(reinterpret_cast <char*>(&buff_i), sizeof(int));
	}

	// alt origin (MRC)
	float ori_float[3]={0,0,0};
	numeric::xyzVector<core::Real> origin_realspace(0,0,0);
	if ( use_altorigin ) {
		idx2cart ( numeric::xyzVector<core::Real>(1,1,1), origin_realspace );
		ori_float[0] = (float)( origin_realspace[0] );
		ori_float[1] = (float)( origin_realspace[1] );
		ori_float[2] = (float)( origin_realspace[2] );
	}
	outx.write(reinterpret_cast <char*>(ori_float), sizeof(float)*3);

	// Write "MAP" at byte 208, indicating a CCP4 file.
	char buff_s[80]; strcpy(buff_s, "MAP DD  ");
	outx.write(reinterpret_cast <char*>(buff_s), 8);

	// fill remainder of head with junk
	int nJunkWords = (1024 - 216) /4;
	buff_i = 0;
	for ( int i=0; i<nJunkWords; ++i ) {
		outx.write(reinterpret_cast <char*>(&buff_i), sizeof(int));
	}

	// data
	int coord[3];
	for ( coord[2] = 1; coord[2] <= density.u3(); coord[2]++ ) {
		for ( coord[1] = 1; coord[1] <= density.u2(); coord[1]++ ) {
			for ( coord[0] = 1; coord[0] <= density.u1(); coord[0]++ ) {
				buff_f = (float) density(coord[0],coord[1],coord[2]);
				outx.write(reinterpret_cast <char*>(&buff_f), sizeof(float));
			}
		}
	}

	return true;
}

core::Real
ElectronDensity::get(numeric::xyzVector<core::Real> X) {
	if ( !isLoaded ) {
		utility_exit_with_message("ElectronDensity::get_interp_dens called but no map is loaded!");
	}

	if ( coeffs_density_.u1()*coeffs_density_.u2()*coeffs_density_.u3() == 0 ) {
		spline_coeffs( density , coeffs_density_ );
	}

	return (interp_spline( coeffs_density_ , X ));
}


void
ElectronDensity::set_voxel_spacing( numeric::xyzVector<core::Real> apix ) {
	// 1) scale origin
	//origin[0] *= (apix[0] * this->grid[0]) / cellDimensions[0];
	//origin[1] *= (apix[1] * this->grid[1]) / cellDimensions[1];
	//origin[2] *= (apix[2] * this->grid[2]) / cellDimensions[2];

	// 2) scale cell dimensions
	cellDimensions[0] = apix[0] * this->grid[0];
	cellDimensions[1] = apix[1] * this->grid[1];
	cellDimensions[2] = apix[2] * this->grid[2];

	// 3) update f2c,c2f
	this->computeCrystParams();

	TR << "Forcing apix to " << apix[0] << "," << apix[1] << "," << apix[2] << std::endl;
	//TR << "new origin:" << this->origin[0] << "," << this->origin[1] << "," << this->origin[2] << std::endl;

	//fpd no need for full density change trigger, just invalidate a few things
	// Fdensity and coeffs_density_ are still valid
	fastdens_score.clear();
	fastgrid = numeric::xyzVector< core::Real >(0,0,0);
}


void ElectronDensity::density_change_trigger() {
	fastdens_score.clear();
	fastgrid = numeric::xyzVector< core::Real >(0,0,0);

	Fdensity.clear();
	coeffs_density_.clear();

	computeGradients();  // visualization only
}


/////////////////////////////////////
// resize a map (using FFT-interpolation)
void ElectronDensity::resize( core::Real approxGridSpacing ) {
	// potentially expand map to cover entire unit cell
	if ( grid[0] != density.u1() || grid[1] != density.u2() || grid[2] != density.u3() ) {
		TR << "[ ERROR ] resize() not supported for maps not covering the entire unit cell."<< std::endl;
		TR << "   " << grid[0] << " != " << density.u1()
			<< " || " << grid[1] << " != " << density.u2()
			<< " || " << grid[2] << " != " << density.u3() << std::endl;
		utility_exit();
	}

	// compute new dimensions & origin
	numeric::xyzVector<int> newDims,  newGrid;
	numeric::xyzVector<double> newOri;

	//fpd since we're doing a bunch of FFTs now resize this to something with no large prime factors
	newDims[0] = findSampling( cellDimensions[0] / approxGridSpacing , MINMULT[0] );
	newDims[1] = findSampling( cellDimensions[1] / approxGridSpacing , MINMULT[1] );
	newDims[2] = findSampling( cellDimensions[2] / approxGridSpacing , MINMULT[2] );

	newOri[0] = newDims[0]*origin[0] / ((core::Real)grid[0]);
	newOri[1] = newDims[1]*origin[1] / ((core::Real)grid[1]);
	newOri[2] = newDims[2]*origin[2] / ((core::Real)grid[2]);
	newGrid = newDims;

	ObjexxFCL::FArray3D< double > newDensity;

	resample( density, newDensity, newDims );
	TR << "Resizing " << density.u1() << "x" << density.u2() << "x" << density.u3() << " to "
		<< newDensity.u1() << "x" << newDensity.u2() << "x" << newDensity.u3() << std::endl;

	// update density
	density.dimension( newDims[0], newDims[1], newDims[2] );
	for ( int i=0; i< newDims[0]*newDims[1]*newDims[2] ; ++i ) {
		density[i] = (float)newDensity[i];
	}
	this->grid = newGrid;
	this->origin = newOri;

	TR << " new extent: " << density.u1() << " x " << density.u2() << " x " << density.u3() << std::endl;
	TR << " new origin: " << origin[0] << " x " << origin[1] << " x " << origin[2] << std::endl;
	TR << "   new grid: " << grid[0] << " x " << grid[1] << " x " << grid[2] << std::endl;
}


/////////////////////////////////////
//  compute gradients of density --- used to calculate surface normals in visualizer
void ElectronDensity::computeGradients() {
#ifdef GL_GRAPHICS
	numeric::xyzVector<int> dims( density.u1(), density.u2(), density.u3() );

	// Allocate arrays
	ObjexxFCL::FArray3D< float > grad_x, grad_y, grad_z;
	grad_x.dimension( dims[0] , dims[1] , dims[2] ); grad_x = 0;
	grad_y.dimension( dims[0] , dims[1] , dims[2] ); grad_y = 0;
	grad_z.dimension( dims[0] , dims[1] , dims[2] ); grad_z = 0;

	for (int x=2; x<dims[0]; ++x) {
		for (int y=2; y<dims[1]; ++y) {
			for (int z=2; z<dims[2]; ++z) {
				grad_x(x,y,z) = 0.5*(density(x+1,y,z) - density(x-1,y,z));
				grad_x(x,y,z) = 0.125*(density(x+1,y+1,z) - density(x-1,y+1,z));
				grad_x(x,y,z) = 0.125*(density(x+1,y-1,z) - density(x-1,y-1,z));
				grad_x(x,y,z) = 0.125*(density(x+1,y,z+1) - density(x-1,y,z+1));
				grad_x(x,y,z) = 0.125*(density(x+1,y,z-1) - density(x-1,y,z-1));

				grad_y(x,y,z) = 0.5*(density(x,y+1,z) - density(x,y-1,z));
				grad_y(x,y,z) = 0.125*(density(x+1,y+1,z) - density(x+1,y-1,z));
				grad_y(x,y,z) = 0.125*(density(x-1,y+1,z) - density(x-1,y-1,z));
				grad_y(x,y,z) = 0.125*(density(x,y+1,z+1) - density(x,y-1,z+1));
				grad_y(x,y,z) = 0.125*(density(x,y+1,z-1) - density(x,y-1,z-1));

				grad_z(x,y,z) = 0.5*(density(x,y,z+1) - density(x,y,z-1));
				grad_z(x,y,z) = 0.125*(density(x+1,y,z+1) - density(x+1,y,z-1));
				grad_z(x,y,z) = 0.125*(density(x-1,y,z+1) - density(x-1,y,z-1));
				grad_z(x,y,z) = 0.125*(density(x,y+1,z+1) - density(x,y+1,z-1));
				grad_z(x,y,z) = 0.125*(density(x,y-1,z+1) - density(x,y-1,z-1));
			}
		}
	}
	// spline coeff of dCOdx's
	spline_coeffs( grad_x , coeff_grad_x );
	spline_coeffs( grad_y , coeff_grad_y );
	spline_coeffs( grad_z , coeff_grad_z );

	TR << "Finished computing gradient maps" << std::endl;
#endif
}

/////////////////////////////////////
// compute map statistics
void ElectronDensity::computeStats() {
	core::Real sum=0, sum2=0;
	dens_max = -std::numeric_limits< core::Real >::max();
	dens_min = std::numeric_limits< core::Real >::max();

	int N = density.u1()*density.u2()*density.u3();

	for ( int i=1; i<=density.u1(); i++ ) {
		for ( int j=1; j<=density.u2(); j++ ) {
			for ( int k=1; k<=density.u3(); k++ ) {
				sum += density(i,j,k);
				sum2 += density(i,j,k)*density(i,j,k);
				dens_max = std::max(dens_max, (core::Real)density(i,j,k));
				dens_min = std::min(dens_min, (core::Real)density(i,j,k));
			}
		}
	}
	dens_mean = sum/N;
	dens_stdev = sqrt( sum2/N-dens_mean*dens_mean );
}

/////////////////////////////////////
// void ElectronDensity::computeCrystParams()
void ElectronDensity::computeCrystParams() {
	// recompute reciprocal cell
	// f2c, c2f
	core::Real ca = cos(d2r(cellAngles[0])), cb = cos(d2r(cellAngles[1])), cg = cos(d2r(cellAngles[2]));
	core::Real sa = sin(d2r(cellAngles[0])), sb = sin(d2r(cellAngles[1])), sg = sin(d2r(cellAngles[2]));

	// conversion from fractional cell coords to cartesian coords
	this->f2c = numeric::xyzMatrix<core::Real>::rows(
		cellDimensions[0]  , cellDimensions[1] * cg, cellDimensions[2] * cb,
		0.0, cellDimensions[1] * sg, cellDimensions[2] * (ca - cb*cg) / sg,
		0.0, 0.0   , cellDimensions[2] * sb * sqrt(1.0 - square((cb*cg - ca)/(sb*sg))) );
	core::Real D = this->f2c.det();
	if ( D == 0 ) {
		TR << "[ WARNING ] Invalid crystal cell dimensions." << std::endl;
		return;
	}
	// c2f is inverse of f2c
	this->c2f = numeric::xyzMatrix<core::Real>::rows(
		(f2c(2,2)*f2c(3,3)-f2c(2,3)*f2c(3,2))/D,
		-(f2c(1,2)*f2c(3,3)-f2c(1,3)*f2c(3,2))/D,
		(f2c(1,2)*f2c(2,3)-f2c(1,3)*f2c(2,2))/D,
		-(f2c(2,1)*f2c(3,3)-f2c(3,1)*f2c(2,3))/D,
		(f2c(1,1)*f2c(3,3)-f2c(1,3)*f2c(3,1))/D,
		-(f2c(1,1)*f2c(2,3)-f2c(1,3)*f2c(2,1))/D,
		(f2c(2,1)*f2c(3,2)-f2c(3,1)*f2c(2,2))/D,
		-(f2c(1,1)*f2c(3,2)-f2c(1,2)*f2c(3,1))/D,
		(f2c(1,1)*f2c(2,2)-f2c(1,2)*f2c(2,1))/D );
	this->V = cellDimensions[0]*cellDimensions[1]*cellDimensions[2]* sqrt(1-square(ca)-square(cb)-square(cg)+2*ca*cb*cg);

	// reciprocal space cell dimensions
	this->RcellDimensions[0] = cellDimensions[1]*cellDimensions[2]*sa/V;
	this->RcellDimensions[1] = cellDimensions[0]*cellDimensions[2]*sb/V;
	this->RcellDimensions[2] = cellDimensions[0]*cellDimensions[1]*sg/V;
	this->cosRcellAngles[0] = cos(  asin( std::min( std::max( V/(cellDimensions[0]*cellDimensions[1]*cellDimensions[2]*sb*sg) , -1.0) , 1.0) )  );
	this->cosRcellAngles[1] = cos(  asin( std::min( std::max( V/(cellDimensions[0]*cellDimensions[1]*cellDimensions[2]*sa*sg) , -1.0) , 1.0) )  );
	this->cosRcellAngles[2] = cos(  asin( std::min( std::max( V/(cellDimensions[0]*cellDimensions[1]*cellDimensions[2]*sa*sb) , -1.0) , 1.0) )  );
	this->RV = 1.0/V;
}

}
}
}
