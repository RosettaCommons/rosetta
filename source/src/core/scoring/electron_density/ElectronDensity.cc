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
#include <numeric/statistics.functions.hh>
#include <numeric/fourier/FFT.hh>

//
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/patterson.OptionKeys.gen.hh>

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
basic::Tracer TR("core.scoring.electron_density.ElectronDensity");

#ifdef GL_GRAPHICS
// protect access to density from viewer thread
pthread_mutex_t density_map_db_mut_ = PTHREAD_MUTEX_INITIALIZER;
#endif

using namespace core;
using namespace basic::options;

const int CCP4HDSIZE = 1024;  // size of CCP4/MRC header


///////////////////////////////  ///////////////////////////////
///////////////////////////////  ///////////////////////////////
///////////////////////////////  ///////////////////////////////
//
// one-liners
inline float d2r(float d) { return (d*M_PI/180.0); }
inline double d2r(double d) { return (d*M_PI/180.0); }
inline float  square(float  x) { return (x*x); }
inline double square(double x) { return (x*x); }


// x mod y, returns z in [-y/2,y/2]
inline int min_mod(int x,int y) {
	int r=x%y; if (r<-y/2) r+=y;if (r>=y/2) r-=y;
	return r;
}
inline float min_mod(float x,float y) {
	float r=std::fmod(x,y); if (r<-0.5*y) r+=y;if (r>=0.5*y) r-=y;
	return r;
}
inline double min_mod(double x,double y) {
	double r=std::fmod(x,y); if (r<-0.5*y) r+=y;if (r>=0.5*y) r-=y;
	return r;
}

// missing density test
// pose_from_pdb randomizing missing density with:
// 		ai.x = ai.x + 900.000 + RG.uniform()*100.000;
inline bool is_missing_density( numeric::xyzVector< core::Real > const &X ) {
	if ( X.length() >= 1500 ) {
	  	return true;
	}
	return false;
}

// Endianness swap
// Only works with aligned 4-byte quantities
static void swap4_aligned(void *v, long ndata) {
	int *data = (int *) v;
	long i;
	int *N;
	for (i=0; i<ndata; i++) {
		N = data + i;
		*N=(((*N>>24)&0xff) | ((*N&0xff)<<24) | ((*N>>8)&0xff00) | ((*N&0xff00)<<8));
	}
}

///////////////////////////////  ///////////////////////////////
///////////////////////////////  ///////////////////////////////
///////////////////////////////  ///////////////////////////////

ElectronDensity& getDensityMap(std::string filename, bool force_reload) {
	/*
	if (basic::options::option[ basic::options::OptionKeys::edensity::force_legacy ]()) {
		return getDensityMap_legacy();
	}
	 */

	if(basic::resource_manager::ResourceManager::get_instance()->
		has_resource_with_description("electron_density")){

		ElectronDensityOP electron_density(
			basic::resource_manager::get_resource< ElectronDensity >(
				"electron_density"));

		return *electron_density;
	} else {
		return getDensityMap_legacy(filename, force_reload);
	}
}

//
ElectronDensity& getDensityMap_legacy(std::string filename, bool force_reload) {
	static ElectronDensity theDensityMap;

#ifdef GL_GRAPHICS
	pthread_mutex_lock(&density_map_db_mut_);
#endif

	if ( !theDensityMap.isMapLoaded() || force_reload ) {
		// load map from disk
		TR << "Loading Density Map" << std::endl;
		if (!basic::options::option[ basic::options::OptionKeys::edensity::mapfile ].user() && filename.length()==0 ) {
			TR.Warning << "[ Warning ] No density map specified" << std::endl;
		} else {
			bool map_loaded=false;
			std::string mapfile;

			if (filename.length()==0) {  // use CMD line args
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
			if (!map_loaded) {
				TR << "[ ERROR ] Error loading density map named '" << mapfile << "'" << std::endl;
				exit(1);
			}
		}
	}

#ifdef GL_GRAPHICS
	pthread_mutex_unlock(&density_map_db_mut_);
#endif

	return theDensityMap;
}



/// null constructor
ElectronDensity::ElectronDensity() {
	init();
}


// "rho_calc" constructor
ElectronDensity::ElectronDensity( utility::vector1< core::pose::PoseOP > poses, core::Real reso, core::Real apix ) {
	core::Size nposes = poses.size();
	this->reso = reso;

	//1  get bounds
	numeric::xyzVector< core::Real > d_min(0,0,0), d_max(0,0,0);
	bool is_set = false;
	const core::Real FLUFF = 10.0; // add a bounding box
	for (core::Size n=1; n<=nposes; ++n) {
		core::pose::Pose &pose = *(poses[n]);
		int nres = pose.total_residue();

		for (int i=1 ; i<=nres; ++i) {
			conformation::Residue const &rsd_i (pose.residue(i));
			if ( (rsd_i.aa() == core::chemical::aa_vrt) || (scoring_mask_.find(i) != scoring_mask_.end()) ) continue;
			int nheavyatoms = rsd_i.nheavyatoms();
			for (int j=1 ; j<=nheavyatoms; ++j) {
				numeric::xyzVector< core::Real > const &xyz_ij = rsd_i.atom(j).xyz();
				if (is_missing_density( xyz_ij )) continue;
				if (!is_set) {
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

	// make fake crystal data
	cellAngles[0] = cellAngles[1] = cellAngles[2] = 90;
	cellDimensions[0] = grid[0]*real_apix[0];
	cellDimensions[1] = grid[1]*real_apix[1];
	cellDimensions[2] = grid[2]*real_apix[2];
	computeCrystParams();
	TR << "    celldim: " << cellDimensions[0] << " x " << cellDimensions[1] << " x " << cellDimensions[2] << std::endl;
	TR << " cellangles: " << cellAngles[0] << " x " << cellAngles[1] << " x " << cellAngles[2] << std::endl;


	// figure out effective B
	core::Real max_del_grid = std::max( real_apix[0] , real_apix[1] );
	max_del_grid = std::max( max_del_grid , real_apix[2] );
	if (reso == 0) max_del_grid *= 1.5;
	else max_del_grid = std::max( max_del_grid, reso/2 );

	TR << "Effective resolution = " << 2*max_del_grid << std::endl;
	effectiveB = 16*max_del_grid*max_del_grid;
	TR << "Effective B factor = " << effectiveB << std::endl;

	// find the origin
	numeric::xyzVector< core::Real > frac_dmin = c2f*(d_min - FLUFF);
	origin[0] = std::floor(frac_dmin[0]*grid[0]);
	origin[1] = std::floor(frac_dmin[1]*grid[1]);
	origin[2] = std::floor(frac_dmin[2]*grid[2]);
	efforigin = origin;

	// atom_mask
	OneGaussianScattering cscat = get_A( "C" );
	core::Real mask_min = 3.0 * sqrt( 2.0 / cscat.k(effectiveB) );
	ATOM_MASK = mask_min;
	ATOM_MASK_PADDING = 2;

	// 2 rho_calc
	density.dimension(grid[0],grid[1],grid[2]);
	for (int i=0; i<density.u1()*density.u2()*density.u3(); ++i) density[i]=0.0;
	numeric::xyzVector< core::Real > cartX, fracX;
	numeric::xyzVector< core::Real > atm_i, atm_j, del_ij;
	for (Size n=1; n<=nposes; ++n) {
		core::pose::Pose &pose = *(poses[n]);
		int nres = pose.total_residue();
		for (int i=1 ; i<=nres; ++i) {
			conformation::Residue const &rsd_i (pose.residue(i));

			// skip vrts & masked reses
			if ( rsd_i.aa() == core::chemical::aa_vrt ) continue;
			if ( scoring_mask_.find(i) != scoring_mask_.end() ) continue;
			int nheavyatoms = rsd_i.nheavyatoms();
			for (int j=1 ; j<=nheavyatoms; ++j) {
				conformation::Atom const &atom_i( rsd_i.atom(j) );
				chemical::AtomTypeSet const & atom_type_set( rsd_i.atom_type_set() );
				std::string elt_i = atom_type_set[ rsd_i.atom_type_index( j ) ].element();
				OneGaussianScattering sig_j = get_A( elt_i );
				core::Real k = sig_j.k( effectiveB );
				core::Real C = sig_j.C( k );

				if ( is_missing_density( atom_i.xyz() ) ) continue;
				if ( C < 1e-6 ) continue;

				cartX = atom_i.xyz() - getTransform();
				fracX = c2f*cartX;
				atm_i[0] = pos_mod (fracX[0]*grid[0] - origin[0] + 1 , (double)grid[0]);
				atm_i[1] = pos_mod (fracX[1]*grid[1] - origin[1] + 1 , (double)grid[1]);
				atm_i[2] = pos_mod (fracX[2]*grid[2] - origin[2] + 1 , (double)grid[2]);

				for (int z=1; z<=density.u3(); ++z) {
					atm_j[2] = z;
					del_ij[2] = (atm_i[2] - atm_j[2]) / grid[2];
					// wrap-around??
					if (del_ij[2] > 0.5) del_ij[2]-=1.0;
					if (del_ij[2] < -0.5) del_ij[2]+=1.0;

					del_ij[0] = del_ij[1] = 0.0;
					if ((f2c*del_ij).length_squared() > (ATOM_MASK+ATOM_MASK_PADDING)*(ATOM_MASK+ATOM_MASK_PADDING)) continue;

					for (int y=1; y<=density.u2(); ++y) {
						atm_j[1] = y;

						// early exit?
						del_ij[1] = (atm_i[1] - atm_j[1]) / grid[1] ;
						// wrap-around??
						if (del_ij[1] > 0.5) del_ij[1]-=1.0;
						if (del_ij[1] < -0.5) del_ij[1]+=1.0;
						del_ij[0] = 0.0;
						if ((f2c*del_ij).length_squared() > (ATOM_MASK+ATOM_MASK_PADDING)*(ATOM_MASK+ATOM_MASK_PADDING)) continue;

						for (int x=1; x<=density.u1(); ++x) {
							atm_j[0] = x;

							// early exit?
							del_ij[0] = (atm_i[0] - atm_j[0]) / grid[0];
							// wrap-around??
							if (del_ij[0] > 0.5) del_ij[0]-=1.0;
							if (del_ij[0] < -0.5) del_ij[0]+=1.0;

							numeric::xyzVector< core::Real > cart_del_ij = (f2c*del_ij);  // cartesian offset from (x,y,z) to atom_i
							core::Real d2 = (cart_del_ij).length_squared();

							if (d2 <= (ATOM_MASK+ATOM_MASK_PADDING)*(ATOM_MASK+ATOM_MASK_PADDING)) {
								core::Real atm = C*exp(-k*d2);
								density(x,y,z) += atm;
							}
						}
					}
				}
			}
		}
	}
}

void
ElectronDensity::init() {
	isLoaded = false;

	grid = numeric::xyzVector< int >(0,0,0);
	efforigin = origin = numeric::xyzVector< int >(0,0,0);
	cellDimensions = numeric::xyzVector< float >(1,1,1);
	cellAngles = numeric::xyzVector< float >(90,90,90);
	use_altorigin =  false;

	legacy_ = basic::options::option[ basic::options::OptionKeys::edensity::legacy_fastdens_score ]();

	// command line overrides defaults
	reso = basic::options::option[ basic::options::OptionKeys::edensity::mapreso ]();
	ATOM_MASK = basic::options::option[ basic::options::OptionKeys::edensity::atom_mask ]();
	ATOM_MASK_PADDING = 2;
	CA_MASK = basic::options::option[ basic::options::OptionKeys::edensity::ca_mask ]();
	WINDOW_ = basic::options::option[ basic::options::OptionKeys::edensity::sliding_window ]();
	score_window_context_ = basic::options::option[ basic::options::OptionKeys::edensity::score_sliding_window_context ]();
	remap_symm_ = basic::options::option[ basic::options::OptionKeys::edensity::score_symm_complex ]();
	force_apix_ = basic::options::option[ basic::options::OptionKeys::edensity::force_apix ]();
	nkbins_ = basic::options::option[ basic::options::OptionKeys::edensity::n_kbins ]();

	// more defaults
	DensScoreInMinimizer = true;
	ExactDerivatives = basic::options::option[ basic::options::OptionKeys::edensity::debug_derivatives ]();

	// use-specified B factor (may be overridden)
	effectiveB = 0;
	if (basic::options::option[ basic::options::OptionKeys::patterson::model_B ].user())
		effectiveB = basic::options::option[ basic::options::OptionKeys::patterson::model_B ]();
	PattersonMinR = 3.0;
	PattersonMaxR = 20.0;
	if (basic::options::option[ basic::options::OptionKeys::patterson::radius_cutoffs ].user() ) {
		utility::vector1< core::Real > radius_cuts = basic::options::option[ basic::options::OptionKeys::patterson::radius_cutoffs ]();
		if (radius_cuts.size() == 1) {
			PattersonMinR = radius_cuts[1];
		} else if (radius_cuts.size() >= 2) {
			PattersonMinR = std::min( radius_cuts[1] , radius_cuts[2] );
			PattersonMaxR = std::max( radius_cuts[1] , radius_cuts[2] );
		}
	}

	p_extent = p_origin = numeric::xyzVector< core::Real >(0,0,0);
	p_grid = numeric::xyzVector< core::Size >(0,0,0);

	// only used when computing numeric derivatives
	NUM_DERIV_H = 0.1;
	NUM_DERIV_H_CEN = NUM_DERIV_H;
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


/////////////////////////////////////
/// Match a residue to the density map, returning correlation coefficient between
///    map and pose
numeric::xyzMatrix< core::Real > ElectronDensity::rotAlign2DPose(
		core::pose::Pose const &pose,
		std::string axis )
{
	// make sure map is loaded
	if (!isLoaded) {
		TR << "[ ERROR ]  ElectronDensity::matchPose called but no map is loaded!\n";
		return numeric::xyzMatrix< core::Real >(0);
	}

	core::Size axis_Z = 2, axis_X = 0, axis_Y = 1;
	if (axis == "X") {
		axis_Z = 0; axis_X = 1; axis_Y = 2;
	}
	if (axis == "Y") {
		axis_Z = 1; axis_X = 2; axis_Y = 0;
	}

	rho_calc.dimension(density.u1() , density.u2() , density.u3());
	for (int i=0; i<density.u1()*density.u2()*density.u3(); ++i) rho_calc[i]=0.0;

	int nres = pose.total_residue();
	numeric::xyzVector< core::Real > cartX, fracX;
	numeric::xyzVector< core::Real > atm_i, atm_j, del_ij, atm_idx_ij;
	core::Real SC_scaling = basic::options::option[ basic::options::OptionKeys::edensity::sc_scaling ]();

	// stats
	core::Real maxRadius = 0;
	core::Real minHeight = 0;
	core::Real maxHeight = 0;
	numeric::xyzVector< core::Real > poseCoM(0,0,0);
	core::Size sum_atoms = 0;
	for (int i=1 ; i<=nres; ++i) {
		conformation::Residue const &rsd_i (pose.residue(i));
		if ( (rsd_i.aa() == core::chemical::aa_vrt) || (scoring_mask_.find(i) != scoring_mask_.end()) ) continue;
		int nheavyatoms = rsd_i.nheavyatoms();
		for (int j=1 ; j<=nheavyatoms; ++j) {
			numeric::xyzVector< core::Real > const &xyz_ij = rsd_i.atom(j).xyz();
			if (is_missing_density( xyz_ij )) continue;
			poseCoM = poseCoM + xyz_ij;
			sum_atoms+=1;
		}
	}
	poseCoM /= sum_atoms;

	/// 1: rho_c
	for (int i=1 ; i<=nres; ++i) {
		conformation::Residue const &rsd_i (pose.residue(i));

		// skip vrts & masked reses
		if ( rsd_i.aa() == core::chemical::aa_vrt ) continue;
		if ( scoring_mask_.find(i) != scoring_mask_.end() ) continue;
		int nheavyatoms = rsd_i.nheavyatoms();

		for (int j=1 ; j<=nheavyatoms; ++j) {
			conformation::Atom const &atm_i( rsd_i.atom(j) );
			chemical::AtomTypeSet const & atom_type_set( rsd_i.atom_type_set() );
			std::string elt_i = atom_type_set[ rsd_i.atom_type_index( j ) ].element();
			OneGaussianScattering sig_j = get_A( elt_i );
			core::Real k = sig_j.k( effectiveB );
			core::Real C = sig_j.C( k );

			if ( (Size) j > rsd_i.last_backbone_atom())
				C *= SC_scaling;
			if ( is_missing_density( atm_i.xyz() ) ) continue;
			if ( C < 1e-6 ) continue;

			// stats
			core::Real thisH = (atm_i.xyz()[axis_Z] - poseCoM[axis_Z]);
			minHeight = std::min( minHeight , thisH );
			maxHeight = std::max( maxHeight , thisH );

			core::Real thisR2 = square(atm_i.xyz()[axis_X] - poseCoM[axis_X]) + square(atm_i.xyz()[axis_Y] - poseCoM[axis_Y]);
			maxRadius = std::max( maxRadius , thisR2 );

			cartX = atm_i.xyz() - getTransform();
			fracX = c2f*cartX;
			atm_idx_ij[0] = pos_mod (fracX[0]*grid[0] - origin[0] + 1 , (double)grid[0]);
			atm_idx_ij[1] = pos_mod (fracX[1]*grid[1] - origin[1] + 1 , (double)grid[1]);
			atm_idx_ij[2] = pos_mod (fracX[2]*grid[2] - origin[2] + 1 , (double)grid[2]);

			for (int z=1; z<=density.u3(); ++z) {
				atm_j[2] = z;
				del_ij[2] = (atm_idx_ij[2] - atm_j[2]) / grid[2];
				if (del_ij[2] > 0.5) del_ij[2]-=1.0;
				if (del_ij[2] < -0.5) del_ij[2]+=1.0;

				del_ij[0] = del_ij[1] = 0.0;
				if ((f2c*del_ij).length_squared() > (ATOM_MASK+1)*(ATOM_MASK+1)) continue;  // early exit

				for (int y=1; y<=density.u2(); ++y) {
					atm_j[1] = y;

					del_ij[1] = (atm_idx_ij[1] - atm_j[1]) / grid[1] ;
					if (del_ij[1] > 0.5) del_ij[1]-=1.0;
					if (del_ij[1] < -0.5) del_ij[1]+=1.0;
					del_ij[0] = 0.0;
					if ((f2c*del_ij).length_squared() > (ATOM_MASK+1)*(ATOM_MASK+1)) continue;  // early exit

					for (int x=1; x<=density.u1(); ++x) {
						atm_j[0] = x;
						del_ij[0] = (atm_idx_ij[0] - atm_j[0]) / grid[0];
						if (del_ij[0] > 0.5) del_ij[0]-=1.0;
						if (del_ij[0] < -0.5) del_ij[0]+=1.0;

						numeric::xyzVector< core::Real > cart_del_ij = (f2c*del_ij);  // cartesian offset from (x,y,z) to atom_i
						core::Real d2 = (cart_del_ij).length_squared();
						if (d2 <= (ATOM_MASK+1)*(ATOM_MASK+1)) {
							core::Real atm = C*exp(-k*d2);
							rho_calc(x,y,z) += atm;
						}
					}
				}
			}
		}
	}
	maxRadius = sqrt( maxRadius );

	if (basic::options::option[ basic::options::OptionKeys::edensity::debug ]()) {
		ElectronDensity(rho_calc,1.0, numeric::xyzVector< core::Real >(0,0,0), false).writeMRC( "rho_calc.mrc" );
		ElectronDensity(density,1.0, numeric::xyzVector< core::Real >(0,0,0), false).writeMRC( "rho_obs.mrc" );
	}

	/// 2: resample in cylindrical shells
	core::Real shell_spacing = std::max( cellDimensions[0]/((double)this->grid[0]) , cellDimensions[1]/((double)this->grid[1]) );
	shell_spacing = 2*std::max( shell_spacing , cellDimensions[2]/((double)this->grid[2]) );
	core::Size nshells = (core::Size) std::ceil( maxRadius / shell_spacing );
	core::Size height_samples = (core::Size) std::ceil( (maxHeight-minHeight) / shell_spacing );
	core::Real height_spacing = (maxHeight-minHeight) / (height_samples-1);
	core::Size radial_samples = (core::Size) std::ceil( 4*M_PI*maxRadius / shell_spacing );  // 2x oversample
	core::Real radial_spacing = 2*M_PI/radial_samples;
	TR << "Realignment with max radius = " << maxRadius << std::endl;
	TR << "     #shells = " << nshells << " ( spacing = " << shell_spacing << ")" << std::endl;
	TR << "     #hsteps = " << height_samples << " ( spacing = " << height_spacing << ")" << std::endl;
	TR << "     #sampls = " << radial_samples << " ( spacing = " << radial_spacing*180/M_PI << "deg )" << std::endl;

	utility::vector1< ObjexxFCL::FArray1D< double > > rho_o_cyl(nshells*height_samples);
	utility::vector1< ObjexxFCL::FArray1D< double > > rho_c_cyl(nshells*height_samples);
	for (int r=0; r<(int)nshells ; ++r ) {
		for (int h=0; h<(int)height_samples ; ++h ) {
			rho_o_cyl[r*height_samples+h+1].dimension( radial_samples );
			rho_c_cyl[r*height_samples+h+1].dimension( radial_samples );
		}
	}

	core::Real sumN=0.0;
	core::Real sumRhoC=0.0, sumRhoO=0.0;
	core::Real sumRho2C=0.0, sumRho2O=0.0;
	core::Real Rwt=0.0;
	for (int r=0; r<(int)nshells ; ++r ){
		Rwt = 2*M_PI*(r+1)*shell_spacing;

		for (int h=0; h<(int)height_samples ; ++h ){
			cartX[axis_Z] = poseCoM[axis_Z] + minHeight + h*height_spacing;

			for (int theta=0; theta<(int)radial_samples ; ++theta ){
				cartX[axis_X] = poseCoM[axis_X] + (r+1)*shell_spacing*cos(theta*radial_spacing);
				cartX[axis_Y] = poseCoM[axis_Y] + (r+1)*shell_spacing*sin(theta*radial_spacing);

				fracX = c2f*cartX;
				atm_idx_ij[0] = pos_mod (fracX[0]*grid[0] - origin[0] + 1 , (double)grid[0]);
				atm_idx_ij[1] = pos_mod (fracX[1]*grid[1] - origin[1] + 1 , (double)grid[1]);
				atm_idx_ij[2] = pos_mod (fracX[2]*grid[2] - origin[2] + 1 , (double)grid[2]);

				core::Real thisRhoO = interp_linear( density , atm_idx_ij );
				core::Real thisRhoC = interp_linear( rho_calc , atm_idx_ij );
				rho_o_cyl[r*height_samples+h+1][theta] = thisRhoO;
				rho_c_cyl[r*height_samples+h+1][theta] = thisRhoC;

				// stats
				sumN += Rwt;
				sumRhoC += Rwt*thisRhoC;
				sumRhoO += Rwt*thisRhoO;
				sumRho2C += Rwt*thisRhoC*thisRhoC;
				sumRho2O += Rwt*thisRhoO*thisRhoO;
			}
		}
	}
	sumRhoC /= sumN;
	sumRhoO /= sumN;
	sumRho2C = std::sqrt( sumRho2C/sumN - sumRhoC*sumRhoC );
	sumRho2O = std::sqrt( sumRho2O/sumN - sumRhoO*sumRhoO );

	//TR << "rho_c:  mean = " << sumRhoC << "  sigma = " << sumRho2C << std::endl;
	//TR << "rho_o:  mean = " << sumRhoO << "  sigma = " << sumRho2O << std::endl;

	/// 3: bunch of 1D ffts to get overlap
	///    weighted sum over shells
	ObjexxFCL::FArray1D< double > correl_sum, correl_flip_sum, correl_i;
	ObjexxFCL::FArray1D< std::complex<double> > FrhoO, FrhoC, Fcorrel_i;

	correl_sum.dimension( radial_samples, 0.0 );
	correl_flip_sum.dimension( radial_samples, 0.0 );

	// standardize
	for (int r=0; r<(int)nshells ; ++r ) {
		for (int h=0; h<(int)height_samples ; ++h ){
			for (int theta=0; theta<(int)radial_samples ; ++theta ){
				rho_o_cyl[r*height_samples+h+1][theta] = (rho_o_cyl[r*height_samples+h+1][theta]-sumRhoO)/sumRho2O;
				rho_c_cyl[r*height_samples+h+1][theta] = (rho_c_cyl[r*height_samples+h+1][theta]-sumRhoC)/sumRho2C;
			}
		}
	}

	Fcorrel_i.dimension( radial_samples );
	for (int r=0; r<(int)nshells ; ++r ) {
		Rwt = 2*M_PI*(r+1)*shell_spacing;

		for (int h=0; h<(int)height_samples ; ++h ){
			/// 1 check non-flipped
			// fft convolution
			numeric::fourier::fft(rho_o_cyl[r*height_samples+h+1], FrhoO);
			numeric::fourier::fft(rho_c_cyl[r*height_samples+h+1], FrhoC);

			// rotate Fc to Fo
			for (int theta=0; theta<(int)radial_samples ; ++theta )
				Fcorrel_i[theta] = FrhoO[theta] * std::conj( FrhoC[theta] );

			numeric::fourier::ifft(Fcorrel_i, correl_i);

			for (int theta=0; theta<(int)radial_samples ; ++theta )
				correl_sum[theta] += Rwt*correl_i[theta];

			/// 2 check flipped
			// fft convolution
			numeric::fourier::fft(rho_c_cyl[r*height_samples+(height_samples-h)], FrhoC);

			// rotate Fc to Fo
			for (int theta=0; theta<(int)radial_samples ; ++theta )
				Fcorrel_i[theta] = FrhoO[theta] * FrhoC[theta];

			numeric::fourier::ifft(Fcorrel_i, correl_i);

			for (int theta=0; theta<(int)radial_samples ; ++theta )
				correl_flip_sum[theta] += Rwt*correl_i[theta];
		}
	}

	int maxTheta = 0;
	core::Real flipped=1;
	core::Real maxCC = -sumN;
	for (int theta=0; theta<(int)radial_samples ; ++theta ) {
		if (correl_sum[theta] > maxCC) {
			maxCC = correl_sum[theta];
			maxTheta = theta;
			flipped = 1;
		}
		if (correl_flip_sum[theta] > maxCC) {
			maxCC = correl_flip_sum[theta];
			maxTheta = theta;
			flipped = -1;
		}
	}
	TR << "Best alignment at theta = " << maxTheta*radial_spacing*180/M_PI << std::endl;
	TR << "                   flip = " << flipped << std::endl;
	TR << "                     cc = " << maxCC/sumN << " ( " << maxCC << " / " << sumN << " ) " << std::endl;

	// finally apply rotation to pose
	numeric::xyzMatrix< core::Real > rotation;
	rotation(axis_Z+1,axis_X+1) = 0; rotation(axis_Z+1,axis_Y+1) = 0;
	rotation(axis_X+1,axis_Z+1) = 0; rotation(axis_Y+1,axis_Z+1) = 0;
	rotation(axis_Z+1,axis_Z+1) = flipped;

	rotation(axis_X+1,axis_X+1) = cos(maxTheta*radial_spacing);
	rotation(axis_X+1,axis_Y+1) = -flipped*sin(maxTheta*radial_spacing);
	rotation(axis_Y+1,axis_X+1) = sin(maxTheta*radial_spacing);
	rotation(axis_Y+1,axis_Y+1) = flipped*cos(maxTheta*radial_spacing);

	return rotation;
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
	if (!isLoaded) {
		TR << "[ ERROR ]  ElectronDensity::matchCentroidPose called but no map is loaded!\n";
		return 0.0;
	}

	if (!DensScoreInMinimizer) cacheCCs = false;

	//ObjexxFCL::FArray3D< double >  rho_calc, inv_rho_mask;
	//ObjexxFCL::FArray3D< double >  inv_rho_mask;
	rho_calc.dimension(density.u1() , density.u2() , density.u3());
	inv_rho_mask.dimension(density.u1() , density.u2() , density.u3());
	for (int i=0; i<density.u1()*density.u2()*density.u3(); ++i) {
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
	for (int i=1 ; i<=nres; ++i) {
		conformation::Residue const &rsd_i (pose.residue(i));

		// skip non-protein residues & masked reses
		if ( !pose.residue_type(i).is_protein() ) continue;
		if ( scoring_mask_.find(i) != scoring_mask_.end() ) continue;

		// symm
		if (isSymm && !symmInfo->bb_is_independent(i) && !remapSymm) {  // should this be fa_...??
			continue; // only score the independent monomer
		}

		conformation::Atom const &atm_i( rsd_i.atom("CA") );

		// skip randomized residues
		if ( is_missing_density( atm_i.xyz() ) ) continue;

		cartX = atm_i.xyz() - getTransform();
		fracX = c2f*cartX;
		atm_idx[i][0] = pos_mod (fracX[0]*grid[0] - origin[0] + 1 , (double)grid[0]);
		atm_idx[i][1] = pos_mod (fracX[1]*grid[1] - origin[1] + 1 , (double)grid[1]);
		atm_idx[i][2] = pos_mod (fracX[2]*grid[2] - origin[2] + 1 , (double)grid[2]);


		for (int z=1; z<=density.u3(); ++z) {
			atm_j[2] = z;
			del_ij[2] = (atm_idx[i][2] - atm_j[2]) / grid[2];
			// wrap-around??
			if (del_ij[2] > 0.5) del_ij[2]-=1.0;
			if (del_ij[2] < -0.5) del_ij[2]+=1.0;

			del_ij[0] = del_ij[1] = 0.0;
			if ((f2c*del_ij).length_squared() > (CA_MASK+ATOM_MASK_PADDING)*(CA_MASK+ATOM_MASK_PADDING)) continue;

			for (int y=1; y<=density.u2(); ++y) {
				atm_j[1] = y;

				// early exit?
				del_ij[1] = (atm_idx[i][1] - atm_j[1]) / grid[1] ;
				// wrap-around??
				if (del_ij[1] > 0.5) del_ij[1]-=1.0;
				if (del_ij[1] < -0.5) del_ij[1]+=1.0;
				del_ij[0] = 0.0;
				if ((f2c*del_ij).length_squared() > (CA_MASK+ATOM_MASK_PADDING)*(CA_MASK+ATOM_MASK_PADDING)) continue;

				for (int x=1; x<=density.u1(); ++x) {
					atm_j[0] = x;

					// early exit?
					del_ij[0] = (atm_idx[i][0] - atm_j[0]) / grid[0];
					// wrap-around??
					if (del_ij[0] > 0.5) del_ij[0]-=1.0;
					if (del_ij[0] < -0.5) del_ij[0]+=1.0;

					numeric::xyzVector< core::Real > cart_del_ij = (f2c*del_ij);  // cartesian offset from (x,y,z) to atom_i
					core::Real d2 = (cart_del_ij).length_squared();

					if (d2 <= (CA_MASK+ATOM_MASK_PADDING)*(CA_MASK+ATOM_MASK_PADDING)) {
						core::Real atm = C*exp(-k*d2);
						core::Real sigmoid_msk = exp( d2 - (ATOM_MASK)*(ATOM_MASK)  );
						core::Real inv_msk = 1/(1+sigmoid_msk);

						rho_calc(x,y,z) += atm;
						inv_rho_mask(x,y,z) *= (1 - inv_msk);

						if (cacheCCs) {
							int idx = (z-1)*density.u2()*density.u1() + (y-1)*density.u1() + x-1;
							rho_dx_pt[i].push_back  ( idx );
							rho_dx_atm[i].push_back ( (-2*k*atm)*cart_del_ij );

							core::Real eps_i = (1-inv_msk), inv_eps_i;
							if (eps_i == 0) // divide-by-zero
								inv_eps_i = sigmoid_msk;
							else
								inv_eps_i = 1/eps_i;

							rho_dx_mask[i].push_back( (-2*sigmoid_msk*inv_msk*inv_msk*inv_eps_i)*cart_del_ij );
						}
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

	for (int x=0; x<density.u1()*density.u2()*density.u3(); ++x) {
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
	if (varC_i == 0 || varO_i == 0)
		CC_i = 0;
	else
		CC_i = (sumCO_i - sumC_i*sumO_i/ vol_i) / sqrt( varC_i * varO_i );

	if (cacheCCs)
		CC_cen = CC_i;


	///////////////////////////
	/// 3  CALCULATE SYMMETRIC ROTATION MATRICES + SYMM MAPPING at each level
	if (isSymm && remapSymm && cacheCCs) {
		compute_symm_rotations( pose, symmInfo );
	}

	///////////////////////////
	/// 4  CALCULATE PER-CA DERIVATIVES
	if (cacheCCs) {
		std::map< core::Size , numeric::xyzMatrix< core::Real > > symmRots;
		for (int i=1 ; i<=nres; ++i) {
			if (isSymm && !symmInfo->bb_is_independent(i) && !remapSymm) {  // should this be fa_...??
				continue; // only score the monomer
			}

			conformation::Residue const &rsd_i (pose.residue(i)); //( *reses[i] );

			//if ( rsd_i.aa() == core::chemical::aa_vrt ) continue;
			if ( !pose.residue_type(i).is_protein() ) continue;
			if ( scoring_mask_.find(i) != scoring_mask_.end() ) continue;

			numeric::xyzVector< core::Real > dVdx_ij(0,0,0), dOdx_ij(0,0,0), dO2dx_ij(0,0,0), dCOdx_ij(0,0,0), dC2dx_ij(0,0,0);

			conformation::Atom const &atm_i( rsd_i.atom("CA") );
	 		if ( is_missing_density( atm_i.xyz() ) ) continue;

			utility::vector1< int > const &rho_dx_pt_ij   = rho_dx_pt[i];
			utility::vector1< numeric::xyzVector<core::Real> > const &rho_dx_mask_ij = rho_dx_mask[i];
			utility::vector1< numeric::xyzVector<core::Real> > const &rho_dx_atm_ij  = rho_dx_atm[i];

			int npoints = rho_dx_pt_ij.size();
			for (int n=1; n<=npoints; ++n) {
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
	}
	// >> debugging <<
	//ElectronDensity(rho_calc, 2.0).writeMRC( "rho_calc.mrc" );
	//ElectronDensity(inv_rho_mask, 2.0).writeMRC( "inv_rho_mask.mrc" );
	//exit(1);

	//std::cerr << "ElectronDensity::matchCentroidPose() returning CC = " << CC_i << "   vol = " << vol_i << std::endl;

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
	if (!isLoaded) {
		TR << "[ ERROR ]  ElectronDensity::matchPose called but no map is loaded!\n";
		return 0.0;
	}

	if (!DensScoreInMinimizer) cacheCCs = false;

	//ObjexxFCL::FArray3D< double >  inv_rho_mask;
	rho_calc.dimension(density.u1() , density.u2() , density.u3());
	inv_rho_mask.dimension(density.u1() , density.u2() , density.u3());
	for (int i=0; i<density.u1()*density.u2()*density.u3(); ++i) {
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
	for (int i=1 ; i<=nres; ++i) {
		conformation::Residue const &rsd_i (pose.residue(i)); //( *reses[i] );

		// skip vrts & masked reses
		if ( rsd_i.aa() == core::chemical::aa_vrt ) continue;
		if ( scoring_mask_.find(i) != scoring_mask_.end() ) continue;

		// symm
		if (isSymm && !symmInfo->bb_is_independent(i) && !remapSymm) {
			continue; // only score the independent monomer
		}

		int nheavyatoms = rsd_i.nheavyatoms();
		atm_idx[i].resize(nheavyatoms);
		rho_dx_pt[i].resize(nheavyatoms);
		rho_dx_mask[i].resize(nheavyatoms);
		rho_dx_atm[i].resize(nheavyatoms);

		for (int j=1 ; j<=nheavyatoms; ++j) {
			conformation::Atom const &atm_i( rsd_i.atom(j) );

			chemical::AtomTypeSet const & atom_type_set( rsd_i.atom_type_set() );
			std::string elt_i = atom_type_set[ rsd_i.atom_type_index( j ) ].element();
			OneGaussianScattering sig_j = get_A( elt_i );
			core::Real k = sig_j.k( effectiveB );
			core::Real C = sig_j.C( k );

			// sidechain weight
			if ( (Size) j > rsd_i.last_backbone_atom())
				C *= SC_scaling;

			// skip randomized residues
			if ( is_missing_density( atm_i.xyz() ) ) continue;

			// if this atom's weight is 0 continue
			if ( C < 1e-6 ) continue;

			cartX = atm_i.xyz() - getTransform();
			fracX = c2f*cartX;
			atm_idx[i][j][0] = pos_mod (fracX[0]*grid[0] - origin[0] + 1 , (double)grid[0]);
			atm_idx[i][j][1] = pos_mod (fracX[1]*grid[1] - origin[1] + 1 , (double)grid[1]);
			atm_idx[i][j][2] = pos_mod (fracX[2]*grid[2] - origin[2] + 1 , (double)grid[2]);


			for (int z=1; z<=density.u3(); ++z) {
				atm_j[2] = z;
				del_ij[2] = (atm_idx[i][j][2] - atm_j[2]) / grid[2];
				// wrap-around??
				if (del_ij[2] > 0.5) del_ij[2]-=1.0;
				if (del_ij[2] < -0.5) del_ij[2]+=1.0;

				del_ij[0] = del_ij[1] = 0.0;
				if ((f2c*del_ij).length_squared() > (ATOM_MASK+ATOM_MASK_PADDING)*(ATOM_MASK+ATOM_MASK_PADDING)) continue;

				for (int y=1; y<=density.u2(); ++y) {
					atm_j[1] = y;

					// early exit?
					del_ij[1] = (atm_idx[i][j][1] - atm_j[1]) / grid[1] ;
					// wrap-around??
					if (del_ij[1] > 0.5) del_ij[1]-=1.0;
					if (del_ij[1] < -0.5) del_ij[1]+=1.0;
					del_ij[0] = 0.0;
					if ((f2c*del_ij).length_squared() > (ATOM_MASK+ATOM_MASK_PADDING)*(ATOM_MASK+ATOM_MASK_PADDING)) continue;

					for (int x=1; x<=density.u1(); ++x) {
						atm_j[0] = x;

						// early exit?
						del_ij[0] = (atm_idx[i][j][0] - atm_j[0]) / grid[0];
						// wrap-around??
						if (del_ij[0] > 0.5) del_ij[0]-=1.0;
						if (del_ij[0] < -0.5) del_ij[0]+=1.0;

						numeric::xyzVector< core::Real > cart_del_ij = (f2c*del_ij);  // cartesian offset from (x,y,z) to atom_i
						core::Real d2 = (cart_del_ij).length_squared();

						if (d2 <= (ATOM_MASK+ATOM_MASK_PADDING)*(ATOM_MASK+ATOM_MASK_PADDING)) {
							core::Real atm = C*exp(-k*d2);
							core::Real sigmoid_msk = exp( d2 - (ATOM_MASK)*(ATOM_MASK)  );
							core::Real inv_msk = 1/(1+sigmoid_msk);

							rho_calc(x,y,z) += atm;
							inv_rho_mask(x,y,z) *= (1 - inv_msk);

							if (cacheCCs) {
								int idx = (z-1)*density.u2()*density.u1() + (y-1)*density.u1() + x-1;

								core::Real eps_i = (1-inv_msk), inv_eps_i;
								if (eps_i == 0) // divide-by-zero
									inv_eps_i = sigmoid_msk;
								else
									inv_eps_i = 1/eps_i;

								rho_dx_pt[i][j].push_back  ( idx );
								rho_dx_atm[i][j].push_back ( (-2*k*atm)*cart_del_ij );
								rho_dx_mask[i][j].push_back( (-2*sigmoid_msk*inv_msk*inv_msk*inv_eps_i)*cart_del_ij );
							}
						}
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

	for (int x=0; x<density.u1()*density.u2()*density.u3(); ++x) {
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
	if (varC_i == 0 || varO_i == 0)
		CC_i = 0;
	else
		CC_i = (sumCO_i - sumC_i*sumO_i/ vol_i) / sqrt( varC_i * varO_i );

	if (cacheCCs)
		CC_aacen = CC_i;


	///////////////////////////
	/// 3  CALCULATE SYMMETRIC ROTATION MATRICES + SYMM MAPPING at each level
	if (isSymm && remapSymm && cacheCCs) {
		compute_symm_rotations( pose, symmInfo );
	}


	///////////////////////////
	/// 4  CALCULATE PER-ATOM DERIVATIVES
	if (cacheCCs) {
		std::map< core::Size , numeric::xyzMatrix< core::Real > > symmRots;
		for (int i=1 ; i<=nres; ++i) {
			if (isSymm && !symmInfo->bb_is_independent(i) && !remapSymm) {  // should this be fa_...??
				continue; // only score the monomer
			}

			conformation::Residue const &rsd_i (pose.residue(i)); //( *reses[i] );

			if ( rsd_i.aa() == core::chemical::aa_vrt ) continue;
			if ( scoring_mask_.find(i) != scoring_mask_.end() ) continue;

			int nheavyatoms = atm_idx[i].size();
			dCCdxs_aacen[i].resize( nheavyatoms, numeric::xyzVector< core::Real >(0,0,0) );

			for (int j=1 ; j<=nheavyatoms; ++j) {
				numeric::xyzVector< core::Real > dVdx_ij(0,0,0), dOdx_ij(0,0,0), dO2dx_ij(0,0,0), dCOdx_ij(0,0,0), dC2dx_ij(0,0,0);

				conformation::Atom const &atm_i( rsd_i.atom(j) );
		 		if ( is_missing_density( atm_i.xyz() ) ) continue;

				chemical::AtomTypeSet const & atom_type_set( rsd_i.atom_type_set() );
				std::string elt_i = atom_type_set[ rsd_i.atom_type_index( j ) ].element();

				utility::vector1< int > const &rho_dx_pt_ij   = rho_dx_pt[i][j];
				utility::vector1< numeric::xyzVector<core::Real> > const &rho_dx_mask_ij = rho_dx_mask[i][j];
				utility::vector1< numeric::xyzVector<core::Real> > const &rho_dx_atm_ij  = rho_dx_atm[i][j];

				int npoints = rho_dx_pt_ij.size();
				for (int n=1; n<=npoints; ++n) {
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
	if (minreso>min_allowed) {
		TR << "Forcing res min to " << 1/min_allowed << std::endl;
		minreso = min_allowed;
	}
	if (maxreso<max_allowed) {
		TR << "Forcing res max to " << 1/max_allowed << std::endl;
		maxreso = max_allowed;
	}

	if (S2_bin) {
		minreso = minreso*minreso;
		maxreso = maxreso*maxreso;
	}

	// get bins
	Real step = (minreso-maxreso)/nbuckets;
	resbins.resize(nbuckets);
	counts.resize(nbuckets);
	for (Size i=1; i<=nbuckets; ++i) {
		counts[i] = 0;
		if (S2_bin) {
			resbins[i] = std::sqrt( maxreso + (i-0.5)*step );
		} else {
			resbins[i] =  maxreso + (i-0.5)*step;
		}
	}

	// get bin counts
	int H,K,L;
	for (int z=1; z<=(int)density.u3(); ++z) {
		H = (z < (int)density.u3()/2) ? z-1 : z-density.u3() - 1;
		for (int y=1; y<=(int)density.u2(); ++y) {
			K = (y < (int)density.u2()/2) ? y-1 : y-density.u2()-1;
			for (int x=1; x<=(int)density.u1(); ++x) {
				L = (x < (int)density.u1()/2) ? x-1 : x-density.u1()-1;
				Real s_i = (S2(H,K,L));
				if (!S2_bin) s_i=sqrt(s_i);
				int bucket_i = 1+(int)std::floor( (s_i-maxreso) / step );
				if ( bucket_i > 0 && bucket_i <= (int)nbuckets ) {
					counts[bucket_i]++;
				}
			}
		}
	}
}

/// @brief Compute intensities from model
utility::vector1< core::Real >
ElectronDensity::getIntensities( core::Size nbuckets, core::Real maxreso, core::Real minreso, bool S2_bin/*=false*/ ) {
	if (Fdensity.u1() == 0) numeric::fourier::fft3(density, Fdensity);

	Real min_allowed = sqrt(S2( density.u1()/2, density.u2()/2, density.u3()/2 ));
	Real max_allowed = std::min( sqrt(S2( 1,0,0 )), sqrt(S2( 0,1,0 )) );
	max_allowed = std::min( max_allowed, sqrt(S2( 0,0,1 )) );
	if (minreso>min_allowed) {
		TR << "Forcing res min to " << 1/min_allowed << std::endl;
		minreso = min_allowed;
	}
	if (maxreso<max_allowed) {
		TR << "Forcing res max to " << 1/max_allowed << std::endl;
		maxreso = max_allowed;
	}

	if (S2_bin) {
		minreso = minreso*minreso;
		maxreso = maxreso*maxreso;
	}
	Real step = (minreso-maxreso)/nbuckets;

	utility::vector1< core::Real > sum_I2(nbuckets, 0.0);
	utility::vector1< core::Size > counts(nbuckets, 0);

	int H,K,L;
	for (int z=1; z<=(int)density.u3(); ++z) {
		H = (z < (int)density.u3()/2) ? z-1 : z-density.u3() - 1;
		for (int y=1; y<=(int)density.u2(); ++y) {
			K = (y < (int)density.u2()/2) ? y-1 : y-density.u2()-1;
			for (int x=1; x<=(int)density.u1(); ++x) {
				L = (x < (int)density.u1()/2) ? x-1 : x-density.u1()-1;
				Real s_i = (S2(H,K,L));
				if (!S2_bin) s_i=sqrt(s_i);
				int bucket_i = 1+(int)std::floor( (s_i-maxreso) / step );
				if ( bucket_i > 0 && bucket_i <= (int)nbuckets ) {
					sum_I2[bucket_i] += std::real( Fdensity(x,y,z)*std::conj(Fdensity(x,y,z)) );
					counts[bucket_i]++;
				}
			}
		}
	}
	for (Size i=1; i<=nbuckets; ++i) {
		sum_I2[i] /= counts[i];
	}

	smooth_intensities(sum_I2);
	return sum_I2;
}

/// @brief Compute intensities from model, applying a mask first
utility::vector1< core::Real >
ElectronDensity::getIntensitiesMasked( poseCoords const &pose, core::Size nbuckets, core::Real maxreso, core::Real minreso, bool S2_bin/*=false*/ ) {
	//////////////
	// 1 compute mask
	core::Real radius = ATOM_MASK;
	ObjexxFCL::FArray3D< float > mask;
	mask.dimension(density.u1() , density.u2() , density.u3());
	mask=1.0;
	for (int i=1 ; i<=(int)pose.size(); ++i) {
		numeric::xyzVector< core::Real> cartX = pose[i].x_ - getTransform();
		numeric::xyzVector< core::Real> fracX = c2f*cartX;
		numeric::xyzVector< core::Real> atm_j, del_ij, atm_idx;
		atm_idx[0] = pos_mod (fracX[0]*grid[0] - origin[0] + 1 , (double)grid[0]);
		atm_idx[1] = pos_mod (fracX[1]*grid[1] - origin[1] + 1 , (double)grid[1]);
		atm_idx[2] = pos_mod (fracX[2]*grid[2] - origin[2] + 1 , (double)grid[2]);
		for (int z=1; z<=density.u3(); ++z) {
			atm_j[2] = z;
			del_ij[2] = (atm_idx[2] - atm_j[2]) / grid[2];
			if (del_ij[2] > 0.5) del_ij[2]-=1.0; if (del_ij[2] < -0.5) del_ij[2]+=1.0;
			del_ij[0] = del_ij[1] = 0.0;
			if ((f2c*del_ij).length_squared() > (radius+ATOM_MASK_PADDING)*(radius+ATOM_MASK_PADDING)) continue;
			for (int y=1; y<=density.u2(); ++y) {
				atm_j[1] = y;
				del_ij[1] = (atm_idx[1] - atm_j[1]) / grid[1] ;
				if (del_ij[1] > 0.5) del_ij[1]-=1.0; if (del_ij[1] < -0.5) del_ij[1]+=1.0;
				del_ij[0] = 0.0;
				if ((f2c*del_ij).length_squared() > (radius+ATOM_MASK_PADDING)*(radius+ATOM_MASK_PADDING)) continue;
				for (int x=1; x<=density.u1(); ++x) {
					atm_j[0] = x;
					del_ij[0] = (atm_idx[0] - atm_j[0]) / grid[0];
					if (del_ij[0] > 0.5) del_ij[0]-=1.0; if (del_ij[0] < -0.5) del_ij[0]+=1.0;
					numeric::xyzVector< core::Real > cart_del_ij = (f2c*del_ij);  // cartesian offset from (x,y,z) to atom_i
					core::Real d2 = (cart_del_ij).length_squared();
					if (d2 <= (radius+ATOM_MASK_PADDING)*(radius+ATOM_MASK_PADDING)) {
						core::Real sigmoid_msk = exp( d2 - (ATOM_MASK)*(ATOM_MASK)  );
						core::Real inv_msk = 1/(1+sigmoid_msk);
						mask(x,y,z) *= (1 - inv_msk);
					}
				}
			}
		}
	}
	//////////////

	//////////////
	// 2 apply to map
	ObjexxFCL::FArray3D< float > densityMask = density;
	for (int i=0; i<density.u1()*density.u2()*density.u3(); ++i)
		densityMask[i] *= (1-mask[i]);

	numeric::fourier::fft3(densityMask, Fdensity);
	//////////////

	//////////////
	// 3 get intensities
	Real min_allowed = sqrt(S2( density.u1()/2, density.u2()/2, density.u3()/2 ));
	Real max_allowed = std::min( sqrt(S2( 1,0,0 )), sqrt(S2( 0,1,0 )) );
	max_allowed = std::min( max_allowed, sqrt(S2( 0,0,1 )) );
	if (minreso>min_allowed) {
		TR << "Forcing res min to " << 1/min_allowed << std::endl;
		minreso = min_allowed;
	}
	if (maxreso<max_allowed) {
		TR << "Forcing res max to " << 1/max_allowed << std::endl;
		maxreso = max_allowed;
	}

	if (S2_bin) {
		minreso = minreso*minreso;
		maxreso = maxreso*maxreso;
	}
	Real step = (minreso-maxreso)/nbuckets;

	utility::vector1< core::Real > sum_I2(nbuckets, 0.0);
	utility::vector1< core::Size > counts(nbuckets, 0);

	int H,K,L;
	for (int z=1; z<=(int)density.u3(); ++z) {
		H = (z < (int)density.u3()/2) ? z-1 : z-density.u3() - 1;
		for (int y=1; y<=(int)density.u2(); ++y) {
			K = (y < (int)density.u2()/2) ? y-1 : y-density.u2()-1;
			for (int x=1; x<=(int)density.u1(); ++x) {
				L = (x < (int)density.u1()/2) ? x-1 : x-density.u1()-1;
				Real s_i = (S2(H,K,L));
				if (!S2_bin) s_i=sqrt(s_i);
				int bucket_i = 1+(int)std::floor( (s_i-maxreso) / step );
				if ( bucket_i > 0 && bucket_i <= (int)nbuckets ) {
					sum_I2[bucket_i] += std::real( Fdensity(x,y,z)*std::conj(Fdensity(x,y,z)) );
					counts[bucket_i]++;
				}
			}
		}
	}
	for (Size i=1; i<=nbuckets; ++i) {
		sum_I2[i] /= counts[i];
	}

	// clear Fdensity since it's masked
	Fdensity.clear();

	smooth_intensities(sum_I2);
	return sum_I2;
}



/// @brief Compute intensities from model
void
ElectronDensity::getIntensities(
	poseCoords const &pose, core::Size nbuckets, core::Real maxreso, core::Real minreso,
	utility::vector1< core::Real > &Imodel, bool S2_bin/*=false*/) {

	// rhoc & FrhoC
	calcRhoC( pose );

	// fft
	if (Fdensity.u1() == 0) numeric::fourier::fft3(density, Fdensity);

	Real min_allowed = sqrt(S2( density.u1()/2, density.u2()/2, density.u3()/2 ));
	Real max_allowed = std::min( sqrt(S2( 1,0,0 )), sqrt(S2( 0,1,0 )) );
	max_allowed = std::min( max_allowed, sqrt(S2( 0,0,1 )) );
	if (minreso>min_allowed) {
		TR << "Forcing res min to " << 1/min_allowed << std::endl;
		minreso = min_allowed;
	}
	if (maxreso<max_allowed) {
		TR << "Forcing res max to " << 1/max_allowed << std::endl;
		maxreso = max_allowed;
	}

	if (S2_bin) {
		minreso = minreso*minreso;
		maxreso = maxreso*maxreso;
	}
	Real step = (minreso-maxreso)/nbuckets;

	Imodel.clear(); Imodel.resize( nbuckets, 0.0 );
	utility::vector1< core::Size > counts(nbuckets, 0);

	int H,K,L;
	for (int z=1; z<=(int)density.u3(); ++z) {
		H = (z < (int)density.u3()/2) ? z-1 : z-density.u3() - 1;
		for (int y=1; y<=(int)density.u2(); ++y) {
			K = (y < (int)density.u2()/2) ? y-1 : y-density.u2()-1;
			for (int x=1; x<=(int)density.u1(); ++x) {
				L = (x < (int)density.u1()/2) ? x-1 : x-density.u1()-1;
				Real s_i = (S2(H,K,L));
				if (!S2_bin) s_i=sqrt(s_i);
				int bucket_i = 1+(int)std::floor( (s_i-maxreso) / step );
				if ( bucket_i > 0 && bucket_i <= (int)nbuckets ) {
					Imodel[bucket_i] += std::real( Frho_calc(x,y,z)*std::conj(Frho_calc(x,y,z)) );
					counts[bucket_i]++;
				}
			}
		}
	}
	for (Size i=1; i<=nbuckets; ++i) {
		Imodel[i] /= counts[i];
	}

	smooth_intensities(Imodel);
}

void
ElectronDensity::smooth_intensities(utility::vector1< core::Real > &Is) const {
	utility::vector1< core::Real > Is_in = Is;

	Is[1] = 0.65*Is_in[1]+0.23*Is_in[2]+0.12*Is_in[3];
	Is[2] = 0.3*Is_in[2]+0.35*Is_in[1]+0.23*Is_in[3]+0.12*Is_in[4];
	for (int i=2; i<=Is.size()-2; ++i) {
		Is[i] = 0.3*Is_in[i]+0.23*Is_in[i-1]+0.23*Is_in[i+1]+0.12*Is_in[i-2]+0.12*Is_in[i+2];
	}
	Is[Is.size()-1] = 0.3*Is_in[Is.size()-1]+0.23*Is_in[Is.size()]+0.23*Is_in[Is.size()-2]+0.12*Is_in[Is.size()-3];
	Is[Is.size()] = 0.3*Is_in[Is.size()]+0.23*Is_in[Is.size()-1]+0.12*Is_in[Is.size()-2];
}


void
ElectronDensity::scaleIntensities( utility::vector1< core::Real > scale_i, core::Real maxreso, core::Real minreso, bool S2_bin/*=false*/ ) {
	if (Fdensity.u1() == 0) numeric::fourier::fft3(density, Fdensity);
	Size nbuckets = scale_i.size();

	Real min_allowed = sqrt(S2( density.u1()/2, density.u2()/2, density.u3()/2 ));
	Real max_allowed = std::min( sqrt(S2( 1,0,0 )), sqrt(S2( 0,1,0 )) );
	max_allowed = std::min( max_allowed, sqrt(S2( 0,0,1 )) );
	if (minreso>min_allowed) {
		TR << "Forcing res min to " << 1/min_allowed << std::endl;
		minreso = min_allowed;
	}
	if (maxreso<max_allowed) {
		TR << "Forcing res max to " << 1/max_allowed << std::endl;
		maxreso = max_allowed;
	}

	if (S2_bin) {
		minreso = minreso*minreso;
		maxreso = maxreso*maxreso;
	}
	Real step = (minreso-maxreso)/nbuckets;

	int H,K,L;
	for (int z=1; z<=(int)density.u3(); ++z) {
		H = (z < (int)density.u3()/2) ? z-1 : z-density.u3() - 1;
		for (int y=1; y<=(int)density.u2(); ++y) {
			K = (y < (int)density.u2()/2) ? y-1 : y-density.u2()-1;
			for (int x=1; x<=(int)density.u1(); ++x) {
				L = (x < (int)density.u1()/2) ? x-1 : x-density.u1()-1;

				Real s_i = (S2(H,K,L));
				if (!S2_bin) s_i=sqrt(s_i);

				//fpd smooth interpolate between buckets
				Real bucket = 0.5+((s_i-maxreso) / step);
				int bucket_i = (int)std::floor(bucket);
				Real bucket_offset0 = bucket-bucket_i;
				Real bucket_offset1 = 1.0-bucket_offset0;

				if ( bucket_i > (int)nbuckets ) {
					Fdensity(x,y,z) = 0.0;
				} else if ( bucket_i == (int)nbuckets ) {
					// linear decay in last bin
					Fdensity(x,y,z) *= bucket_offset1*scale_i[nbuckets];
				} else if ( bucket_i < 0 ) {
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
ElectronDensity::maskDensityMap( poseCoords const &pose, core::Real radius ) {
	// get rho_c
	if (radius == 0) { radius = ATOM_MASK; }

	ObjexxFCL::FArray3D< float > mask;
	mask.dimension(density.u1() , density.u2() , density.u3());
	mask=1.0;

	for (int i=1 ; i<=(int)pose.size(); ++i) {
		numeric::xyzVector< core::Real> cartX = pose[i].x_ - getTransform();
		numeric::xyzVector< core::Real> fracX = c2f*cartX;
		numeric::xyzVector< core::Real> atm_j, del_ij, atm_idx;
		atm_idx[0] = pos_mod (fracX[0]*grid[0] - origin[0] + 1 , (double)grid[0]);
		atm_idx[1] = pos_mod (fracX[1]*grid[1] - origin[1] + 1 , (double)grid[1]);
		atm_idx[2] = pos_mod (fracX[2]*grid[2] - origin[2] + 1 , (double)grid[2]);
		for (int z=1; z<=density.u3(); ++z) {
			atm_j[2] = z;
			del_ij[2] = (atm_idx[2] - atm_j[2]) / grid[2];
			if (del_ij[2] > 0.5) del_ij[2]-=1.0; if (del_ij[2] < -0.5) del_ij[2]+=1.0;
			del_ij[0] = del_ij[1] = 0.0;
			if ((f2c*del_ij).length_squared() > (radius+ATOM_MASK_PADDING)*(radius+ATOM_MASK_PADDING)) continue;
			for (int y=1; y<=density.u2(); ++y) {
				atm_j[1] = y;
				del_ij[1] = (atm_idx[1] - atm_j[1]) / grid[1] ;
				if (del_ij[1] > 0.5) del_ij[1]-=1.0; if (del_ij[1] < -0.5) del_ij[1]+=1.0;
				del_ij[0] = 0.0;
				if ((f2c*del_ij).length_squared() > (radius+ATOM_MASK_PADDING)*(radius+ATOM_MASK_PADDING)) continue;
				for (int x=1; x<=density.u1(); ++x) {
					atm_j[0] = x;
					del_ij[0] = (atm_idx[0] - atm_j[0]) / grid[0];
					if (del_ij[0] > 0.5) del_ij[0]-=1.0; if (del_ij[0] < -0.5) del_ij[0]+=1.0;
					numeric::xyzVector< core::Real > cart_del_ij = (f2c*del_ij);  // cartesian offset from (x,y,z) to atom_i
					core::Real d2 = (cart_del_ij).length_squared();
					if (d2 <= (radius+ATOM_MASK_PADDING)*(radius+ATOM_MASK_PADDING)) {
						core::Real sigmoid_msk = exp( d2 - (ATOM_MASK)*(ATOM_MASK)  );
						core::Real inv_msk = 1/(1+sigmoid_msk);
						mask(x,y,z) *= (1 - inv_msk);
					}
				}
			}
		}
	}

	for (int i=0; i<density.u1()*density.u2()*density.u3(); ++i) {
		density[i] *= (1-mask[i]);
	}

	if (basic::options::option[ basic::options::OptionKeys::edensity::debug ]()) {
		ElectronDensity(density,1.0, numeric::xyzVector< core::Real >(0,0,0), false ).writeMRC( "rho_obs_masked.mrc");
		ElectronDensity(mask,1.0, numeric::xyzVector< core::Real >(0,0,0), false ).writeMRC( "rho_mask_inv.mrc");
	}

	// clear derived data
	density_change_trigger();
}


/// @brief Compute model-map FSC, masked by the pose
void
ElectronDensity::getFSC(
		poseCoords const &pose,
		core::Size nbuckets,
		core::Real maxreso,
		core::Real minreso,
		utility::vector1< core::Real > &modelmapFSC,
		utility::vector1< core::Real > &modelmapSigma,
		bool masked/*=false*/,
		bool S2_bin/*=false*/,
		Real mask_radius/*=0*/ ) {

	//////////////
	// 1 compute mask
	ObjexxFCL::FArray3D< float > *densityMasked=NULL;
	ObjexxFCL::FArray3D< float > densityCopy;

	if (masked) {
		core::Real radius = mask_radius>0 ? mask_radius : ATOM_MASK;
		ObjexxFCL::FArray3D< float > mask;
		mask.dimension(density.u1() , density.u2() , density.u3());
		mask=1.0;
		densityCopy = density;

		for (int i=1 ; i<=(int)pose.size(); ++i) {
			numeric::xyzVector< core::Real> cartX = pose[i].x_ - getTransform();
			numeric::xyzVector< core::Real> fracX = c2f*cartX;
			numeric::xyzVector< core::Real> atm_j, del_ij, atm_idx;
			atm_idx[0] = pos_mod (fracX[0]*grid[0] - origin[0] + 1 , (double)grid[0]);
			atm_idx[1] = pos_mod (fracX[1]*grid[1] - origin[1] + 1 , (double)grid[1]);
			atm_idx[2] = pos_mod (fracX[2]*grid[2] - origin[2] + 1 , (double)grid[2]);
			for (int z=1; z<=density.u3(); ++z) {
				atm_j[2] = z;
				del_ij[2] = (atm_idx[2] - atm_j[2]) / grid[2];
				if (del_ij[2] > 0.5) del_ij[2]-=1.0; if (del_ij[2] < -0.5) del_ij[2]+=1.0;
				del_ij[0] = del_ij[1] = 0.0;
				if ((f2c*del_ij).length_squared() > (radius+ATOM_MASK_PADDING)*(radius+ATOM_MASK_PADDING)) continue;
				for (int y=1; y<=density.u2(); ++y) {
					atm_j[1] = y;
					del_ij[1] = (atm_idx[1] - atm_j[1]) / grid[1] ;
					if (del_ij[1] > 0.5) del_ij[1]-=1.0; if (del_ij[1] < -0.5) del_ij[1]+=1.0;
					del_ij[0] = 0.0;
					if ((f2c*del_ij).length_squared() > (radius+ATOM_MASK_PADDING)*(radius+ATOM_MASK_PADDING)) continue;
					for (int x=1; x<=density.u1(); ++x) {
						atm_j[0] = x;
						del_ij[0] = (atm_idx[0] - atm_j[0]) / grid[0];
						if (del_ij[0] > 0.5) del_ij[0]-=1.0; if (del_ij[0] < -0.5) del_ij[0]+=1.0;
						numeric::xyzVector< core::Real > cart_del_ij = (f2c*del_ij);  // cartesian offset from (x,y,z) to atom_i
						core::Real d2 = (cart_del_ij).length_squared();
						if (d2 <= (radius+ATOM_MASK_PADDING)*(radius+ATOM_MASK_PADDING)) {
							core::Real sigmoid_msk = exp( d2 - (radius)*(radius)  );
							core::Real inv_msk = 1/(1+sigmoid_msk);
							mask(x,y,z) *= (1 - inv_msk);
						}
					}
				}
			}
		}

		//////////////
		// 2 apply to map
		for (int i=0; i<density.u1()*density.u2()*density.u3(); ++i)
			densityCopy[i] *= (1-mask[i]);

		densityMasked = &densityCopy;
	} else {
		densityMasked = &density;
	}

	numeric::fourier::fft3(*densityMasked, Fdensity);

	//////////////
	// 3 get FSC
	Real min_allowed = sqrt(S2( density.u1()/2, density.u2()/2, density.u3()/2 ));
	Real max_allowed = std::min( sqrt(S2( 1,0,0 )), sqrt(S2( 0,1,0 )) );
	max_allowed = std::min( max_allowed, sqrt(S2( 0,0,1 )) );
	if (minreso>min_allowed) {
		TR << "Forcing res min to " << 1/min_allowed << std::endl;
		minreso = min_allowed;
	}
	if (maxreso<max_allowed) {
		TR << "Forcing res max to " << 1/max_allowed << std::endl;
		maxreso = max_allowed;
	}

	if (S2_bin) {
		minreso = minreso*minreso;
		maxreso = maxreso*maxreso;
	}
	Real step = (minreso-maxreso)/nbuckets;

	calcRhoC( pose );
	numeric::xyzVector< core::Real > del_ij;
	utility::vector1<core::Real> num(nbuckets, 0.0), denom1(nbuckets, 0.0), denom2(nbuckets, 0.0), counter(nbuckets, 0.0);

	modelmapFSC.resize(nbuckets);
	modelmapSigma.resize(nbuckets);
	for (Size i=1; i<=nbuckets; ++i) modelmapSigma[i] = 0;

	int H,K,L;
	for (int z=1; z<=(int)density.u3(); ++z) {
		H = (z < (int)density.u3()/2) ? z-1 : z-density.u3() - 1;
		for (int y=1; y<=(int)density.u2(); ++y) {
			K = (y < (int)density.u2()/2) ? y-1 : y-density.u2()-1;
			for (int x=1; x<=(int)density.u1(); ++x) {
				L = (x < (int)density.u1()/2) ? x-1 : x-density.u1()-1;
				Real s_i = (S2(H,K,L));
				if (!S2_bin) s_i=sqrt(s_i);
				int bucket_i = 1+(int)std::floor( (s_i-maxreso) / step );
				if ( bucket_i > 0 && bucket_i <= (int)nbuckets ) {
					Real num_i = std::real( Frho_calc(x,y,z) * std::conj( Fdensity(x,y,z) ) );

					num[bucket_i] += num_i;
					denom1[bucket_i] += std::abs( Frho_calc(x,y,z) ) * std::abs( Frho_calc(x,y,z) );
					denom2[bucket_i] += std::abs( Fdensity(x,y,z) ) * std::abs( Fdensity(x,y,z) );

					Real denom_i = std::abs(Fdensity(x,y,z)) * std::abs( Frho_calc(x,y,z) );
					Real err = 0;
					if (denom_i > 0) {
						Real ratio_i = std::max( std::min( num_i/denom_i,1.0) , 0.0 );
						err =  acos( ratio_i );
					}

					modelmapSigma[bucket_i] += err;
					counter[bucket_i] += 1;
				}
			}
		}
	}

	// clear Fdensity since it's masked
	Fdensity.clear();

	for (Size i=1; i<=nbuckets; ++i) {
		denom1[i] = sqrt(denom1[i]);
		denom2[i] = sqrt(denom2[i]);
		modelmapFSC[i] = num[i] / (denom1[i]*denom2[i]);
		modelmapSigma[i] = sqrt( modelmapSigma[i] / counter[i] );
	}
}

/// @brief Compute model-map FSC, masked by the pose
void
ElectronDensity::getMLE(
		poseCoords const &pose,
		core::Size nbuckets,
		core::Real maxreso, core::Real minreso,
		utility::vector1< core::Real > const &mapmapSigma,
		utility::vector1< core::Real > &modelmapSigma,
		Real &errS2,
		bool masked/*=false*/, bool S2_bin/*=false*/,
		Real mask_radius/*=0*/ ) {
	//////////////
	// 1 compute mask
	ObjexxFCL::FArray3D< float > *densityMasked=NULL;
	ObjexxFCL::FArray3D< float > densityCopy;

	if (masked) {
		core::Real radius = mask_radius>0 ? mask_radius : ATOM_MASK;
		ObjexxFCL::FArray3D< float > mask;
		mask.dimension(density.u1() , density.u2() , density.u3());
		mask=1.0;
		densityCopy = density;

		for (int i=1 ; i<=(int)pose.size(); ++i) {
			numeric::xyzVector< core::Real> cartX = pose[i].x_ - getTransform();
			numeric::xyzVector< core::Real> fracX = c2f*cartX;
			numeric::xyzVector< core::Real> atm_j, del_ij, atm_idx;
			atm_idx[0] = pos_mod (fracX[0]*grid[0] - origin[0] + 1 , (double)grid[0]);
			atm_idx[1] = pos_mod (fracX[1]*grid[1] - origin[1] + 1 , (double)grid[1]);
			atm_idx[2] = pos_mod (fracX[2]*grid[2] - origin[2] + 1 , (double)grid[2]);
			for (int z=1; z<=density.u3(); ++z) {
				atm_j[2] = z;
				del_ij[2] = (atm_idx[2] - atm_j[2]) / grid[2];
				if (del_ij[2] > 0.5) del_ij[2]-=1.0; if (del_ij[2] < -0.5) del_ij[2]+=1.0;
				del_ij[0] = del_ij[1] = 0.0;
				if ((f2c*del_ij).length_squared() > (radius+ATOM_MASK_PADDING)*(radius+ATOM_MASK_PADDING)) continue;
				for (int y=1; y<=density.u2(); ++y) {
					atm_j[1] = y;
					del_ij[1] = (atm_idx[1] - atm_j[1]) / grid[1] ;
					if (del_ij[1] > 0.5) del_ij[1]-=1.0; if (del_ij[1] < -0.5) del_ij[1]+=1.0;
					del_ij[0] = 0.0;
					if ((f2c*del_ij).length_squared() > (radius+ATOM_MASK_PADDING)*(radius+ATOM_MASK_PADDING)) continue;
					for (int x=1; x<=density.u1(); ++x) {
						atm_j[0] = x;
						del_ij[0] = (atm_idx[0] - atm_j[0]) / grid[0];
						if (del_ij[0] > 0.5) del_ij[0]-=1.0; if (del_ij[0] < -0.5) del_ij[0]+=1.0;
						numeric::xyzVector< core::Real > cart_del_ij = (f2c*del_ij);  // cartesian offset from (x,y,z) to atom_i
						core::Real d2 = (cart_del_ij).length_squared();
						if (d2 <= (radius+ATOM_MASK_PADDING)*(radius+ATOM_MASK_PADDING)) {
							core::Real sigmoid_msk = exp( d2 - (radius)*(radius)  );
							core::Real inv_msk = 1/(1+sigmoid_msk);
							mask(x,y,z) *= (1 - inv_msk);
						}
					}
				}
			}
		}

		//////////////
		// 2 apply to map
		for (int i=0; i<density.u1()*density.u2()*density.u3(); ++i)
			densityCopy[i] *= (1-mask[i]);

		densityMasked = &densityCopy;
	} else {
		densityMasked = &density;
	}

	numeric::fourier::fft3(*densityMasked, Fdensity);

	//////////////
	// 3 get FSC
	Real min_allowed = sqrt(S2( density.u1()/2, density.u2()/2, density.u3()/2 ));
	Real max_allowed = std::min( sqrt(S2( 1,0,0 )), sqrt(S2( 0,1,0 )) );
	max_allowed = std::min( max_allowed, sqrt(S2( 0,0,1 )) );
	if (minreso>min_allowed) {
		TR << "Forcing res min to " << 1/min_allowed << std::endl;
		minreso = min_allowed;
	}
	if (maxreso<max_allowed) {
		TR << "Forcing res max to " << 1/max_allowed << std::endl;
		maxreso = max_allowed;
	}

	if (S2_bin) {
		minreso = minreso*minreso;
		maxreso = maxreso*maxreso;
	}
	Real step = (minreso-maxreso)/nbuckets;

	calcRhoC( pose );
	numeric::xyzVector< core::Real > del_ij;
	utility::vector1<core::Real> num(nbuckets, 0.0), counter(nbuckets, 0.0);

	modelmapSigma.resize(nbuckets);
	for (Size i=1; i<=nbuckets; ++i) modelmapSigma[i] = 0;

	int H,K,L;
	errS2 = 0.0;
	Real wt1=0, wtS=0, wtS2=0;
	for (int z=1; z<=(int)density.u3(); ++z) {
		H = (z < (int)density.u3()/2) ? z-1 : z-density.u3() - 1;
		for (int y=1; y<=(int)density.u2(); ++y) {
			K = (y < (int)density.u2()/2) ? y-1 : y-density.u2()-1;
			for (int x=1; x<=(int)density.u1(); ++x) {
				L = (x < (int)density.u1()/2) ? x-1 : x-density.u1()-1;
				Real s_i = (S2(H,K,L));
				if (!S2_bin) s_i=sqrt(s_i);
				int bucket_i = 1+(int)std::floor( (s_i-maxreso) / step );
				if ( bucket_i > 0 && bucket_i <= (int)nbuckets ) {
					Real num_i = std::real( Frho_calc(x,y,z) * std::conj( Fdensity(x,y,z) ) );

					num[bucket_i] += num_i;

					Real denom_i = std::abs(Fdensity(x,y,z)) * std::abs( Frho_calc(x,y,z) );
					Real err_j = 0;
					if (denom_i > 0) {
						Real ratio_i = std::max( std::min( num_i/denom_i,1.0) , 0.0 );
						err_j =  acos( ratio_i );
					}

					Real err_k = mapmapSigma[bucket_i];
					std::complex<Real> ml_error = exp( err_j*err_j / (std::complex<Real>(0,1)-2*err_k*err_k) );
					ml_error /= sqrt( 1.0+2.0*std::complex<Real>(0,1)*err_k*err_k );

					Real err = -std::arg( ml_error );  // E (phase_error)^2

					modelmapSigma[bucket_i] += err;

					errS2 += (1.0/(s_i*s_i))*err;
					wtS2+=(1.0/(s_i*s_i));

					counter[bucket_i] += 1;
				}
			}
		}
	}

	errS2 = sqrt( errS2 / wtS2 );

	// clear Fdensity since it's masked
	Fdensity.clear();

	for (Size i=1; i<=nbuckets; ++i) {
		modelmapSigma[i] = sqrt( modelmapSigma[i] / counter[i] );
	}
}



/// @brief Compute the FSC in the specified resolution range
void
ElectronDensity::getFSC(
		ObjexxFCL::FArray3D< float > const &density2,
		core::Size nbuckets,
		core::Real maxreso,
		core::Real minreso,
		utility::vector1< core::Real > &mapmapFSC,
		utility::vector1< core::Real > &mapmapSigma,
		bool S2_bin/*=false*/) {

	Real min_allowed = sqrt(S2( density.u1()/2, density.u2()/2, density.u3()/2 ));
	Real max_allowed = std::min( sqrt(S2( 1,0,0 )), sqrt(S2( 0,1,0 )) );
	max_allowed = std::min( max_allowed, sqrt(S2( 0,0,1 )) );
	if (minreso>min_allowed) {
		TR << "Forcing res min to " << 1/min_allowed << std::endl;
		minreso = min_allowed;
	}
	if (maxreso<max_allowed) {
		TR << "Forcing res max to " << 1/max_allowed << std::endl;
		maxreso = max_allowed;
	}

	if (S2_bin) {
		minreso = minreso*minreso;
		maxreso = maxreso*maxreso;
	}
	Real step = (minreso-maxreso)/nbuckets;

	runtime_assert( density.u1()==density2.u1() && density.u2()==density2.u2() && density.u3()==density2.u3() );

	// fft
	ObjexxFCL::FArray3D< std::complex<double> > Fdensity2;
	if (Fdensity.u1() == 0)
		numeric::fourier::fft3(density, Fdensity);
	numeric::fourier::fft3(density2, Fdensity2);

	// correl
	numeric::xyzVector< core::Real > del_ij;
	utility::vector1<core::Real> num(nbuckets, 0.0), denom1(nbuckets, 0.0), denom2(nbuckets, 0.0), counter(nbuckets, 0.0);
	mapmapFSC.resize(nbuckets);
	mapmapSigma.resize(nbuckets);

	for (Size i=1; i<=nbuckets; ++i) mapmapSigma[i] = 0;


	int H,K,L;
	for (int z=1; z<=(int)density.u3(); ++z) {
		H = (z < (int)density.u3()/2) ? z-1 : z-density.u3() - 1;
		for (int y=1; y<=(int)density.u2(); ++y) {
			K = (y < (int)density.u2()/2) ? y-1 : y-density.u2()-1;
			for (int x=1; x<=(int)density.u1(); ++x) {
				L = (x < (int)density.u1()/2) ? x-1 : x-density.u1()-1;
				Real s_i = (S2(H,K,L));
				if (!S2_bin) s_i=sqrt(s_i);
				int bucket_i = 1+(int)std::floor( (s_i-maxreso) / step );
				if ( bucket_i > 0 && bucket_i <= (int)nbuckets ) {
					Real num_i = std::real( Fdensity2(x,y,z) * std::conj( Fdensity(x,y,z) ) );

					num[bucket_i] += num_i;    // |E_1|*|E_2|*cos(a_1-a_2)
					denom1[bucket_i] += std::abs( Fdensity2(x,y,z) ) * std::abs( Fdensity2(x,y,z) );
					denom2[bucket_i] += std::abs( Fdensity(x,y,z) ) * std::abs( Fdensity(x,y,z) );

					Real denom_i = std::abs(Fdensity(x,y,z)) * std::abs( Fdensity2(x,y,z) );
					Real err = 0;
					if (denom_i > 0) {
						Real ratio_i = std::max( std::min( num_i/denom_i,1.0) , -1.0 );
						err =  acos( ratio_i );
					}

					mapmapSigma[bucket_i] += err*err;
					counter[bucket_i] += 1;
				}
			}
		}
	}

	for (Size i=1; i<=nbuckets; ++i) {
		denom1[i] = sqrt(denom1[i]);
		denom2[i] = sqrt(denom2[i]);
		mapmapFSC[i] = num[i] / (denom1[i]*denom2[i]);
		mapmapSigma[i] = sqrt( mapmapSigma[i] / counter[i] );
	}

	return;
}

core::Real
ElectronDensity::getRSCC( poseCoords const &pose ) {
	// rhoc & FrhoC
	calcRhoC( pose );

	core::Real sumC_i=0, sumO_i=0, sumCO_i=0, vol_i=0, CC_i=0;
 	core::Real sumO2_i=0.0, sumC2_i=0.0, varC_i=0, varO_i=0;
	core::Real clc_x, obs_x, eps_x;
	for (int x=0; x<density.u1()*density.u2()*density.u3(); ++x) {
		clc_x = rho_calc[x];
		obs_x = density[x];
		eps_x = 1-inv_rho_mask[x];

		sumCO_i += eps_x*clc_x*obs_x;
		sumO_i  += eps_x*obs_x;
		sumO2_i += eps_x*obs_x*obs_x;
		sumC_i  += eps_x*clc_x;
		sumC2_i += eps_x*clc_x*clc_x;
		vol_i   += eps_x;
	}
	varC_i = (sumC2_i - sumC_i*sumC_i / vol_i );
	varO_i = (sumO2_i - sumO_i*sumO_i / vol_i ) ;
	if (varC_i == 0 || varO_i == 0)
		CC_i = 0;
	else
		CC_i = (sumCO_i - sumC_i*sumO_i/ vol_i) / sqrt( varC_i * varO_i );

	return CC_i;
}

core::Real
ElectronDensity::getRSCC( ObjexxFCL::FArray3D< float > const &density2 ) {
	runtime_assert( density.u1()==density2.u1() && density.u2()==density2.u2() && density.u3()==density2.u3() );

	core::Real sumC_i=0, sumO_i=0, sumCO_i=0, vol_i=0, CC_i=0;
 	core::Real sumO2_i=0.0, sumC2_i=0.0, varC_i=0, varO_i=0;
	core::Real clc_x, obs_x, eps_x;
	for (int x=0; x<density.u1()*density.u2()*density.u3(); ++x) {
		clc_x = density2[x];
		obs_x = density[x];

		sumCO_i += clc_x*obs_x; sumO_i  += obs_x;
		sumO2_i += obs_x*obs_x; sumC_i  += clc_x;
		sumC2_i += clc_x*clc_x;
		vol_i   += 1.0;
	}
	varC_i = (sumC2_i - sumC_i*sumC_i / vol_i );
	varO_i = (sumO2_i - sumO_i*sumO_i / vol_i ) ;
	if (varC_i == 0 || varO_i == 0)
		CC_i = 0;
	else
		CC_i = (sumCO_i - sumC_i*sumO_i/ vol_i) / sqrt( varC_i * varO_i );

	return CC_i;
}

core::Real
ElectronDensity::maxNominalRes() {
	Real S = (1/sqrt(3.)) * sqrt(S2( density.u1()/2, density.u2()/2, density.u3()/2 ));
	return 1.0/S;
}


//
void
ElectronDensity::calcRhoC( poseCoords const &pose ) {
	// get rho_c
	rho_calc.dimension(density.u1() , density.u2() , density.u3());
	inv_rho_mask.dimension(density.u1() , density.u2() , density.u3());
	rho_calc=0.0;
	inv_rho_mask = 1.0;

	bool use_Bs = pose_has_nonzero_Bs( pose );  // warn on using this?
	if (!use_Bs) {
		TR << "Input pose has no nonzero B factors ... setting to default." << std::endl;
	}

	for (int i=1 ; i<=(int)pose.size(); ++i) {
		std::string elt_i = pose[i].elt_;
		OneGaussianScattering sig_j = get_A( elt_i );
		core::Real k = sig_j.k( pose[i].B_ );
		if (use_Bs) {
			k = std::min ( k, 4*M_PI*M_PI/minimumB );
		} else {
			k = std::min ( k, 4*M_PI*M_PI/effectiveB );
		}
		core::Real C = sig_j.C( k );
		if ( C < 1e-6 ) continue;

		numeric::xyzVector< core::Real> cartX = pose[i].x_ - getTransform();
		numeric::xyzVector< core::Real> fracX = c2f*cartX;
		numeric::xyzVector< core::Real> atm_j, del_ij, atm_idx;
		atm_idx[0] = pos_mod (fracX[0]*grid[0] - origin[0] + 1 , (double)grid[0]);
		atm_idx[1] = pos_mod (fracX[1]*grid[1] - origin[1] + 1 , (double)grid[1]);
		atm_idx[2] = pos_mod (fracX[2]*grid[2] - origin[2] + 1 , (double)grid[2]);
		for (int z=1; z<=density.u3(); ++z) {
			atm_j[2] = z;
			del_ij[2] = (atm_idx[2] - atm_j[2]) / grid[2];
			if (del_ij[2] > 0.5) del_ij[2]-=1.0; if (del_ij[2] < -0.5) del_ij[2]+=1.0;
			del_ij[0] = del_ij[1] = 0.0;
			if ((f2c*del_ij).length_squared() > (ATOM_MASK+ATOM_MASK_PADDING)*(ATOM_MASK+ATOM_MASK_PADDING)) continue;
			for (int y=1; y<=density.u2(); ++y) {
				atm_j[1] = y;
				del_ij[1] = (atm_idx[1] - atm_j[1]) / grid[1] ;
				if (del_ij[1] > 0.5) del_ij[1]-=1.0; if (del_ij[1] < -0.5) del_ij[1]+=1.0;
				del_ij[0] = 0.0;
				if ((f2c*del_ij).length_squared() > (ATOM_MASK+ATOM_MASK_PADDING)*(ATOM_MASK+ATOM_MASK_PADDING)) continue;
				for (int x=1; x<=density.u1(); ++x) {
					atm_j[0] = x;
					del_ij[0] = (atm_idx[0] - atm_j[0]) / grid[0];
					if (del_ij[0] > 0.5) del_ij[0]-=1.0; if (del_ij[0] < -0.5) del_ij[0]+=1.0;
					numeric::xyzVector< core::Real > cart_del_ij = (f2c*del_ij);  // cartesian offset from (x,y,z) to atom_i
					core::Real d2 = (cart_del_ij).length_squared();

					if (d2 <= (ATOM_MASK+ATOM_MASK_PADDING)*(ATOM_MASK+ATOM_MASK_PADDING)) {
						core::Real atm = C*exp(-k*d2);
						core::Real sigmoid_msk = exp( d2 - (ATOM_MASK)*(ATOM_MASK)  );
						core::Real inv_msk = 1/(1+sigmoid_msk);
						inv_rho_mask(x,y,z) *= (1 - inv_msk);
						rho_calc(x,y,z) += atm;
					}
				}
			}
		}
	}

	// ffts
	numeric::fourier::fft3(rho_calc, Frho_calc);
}


/*
//fpd
//   flat solvent model + bsol/ksol optimization
void
ElectronDensity::calcRhoCandSolvent( poseCoords const &pose ) {
	// get rho_c
	rho_calc.dimension(density.u1() , density.u2() , density.u3());
	for (int i=0; i<density.u1()*density.u2()*density.u3(); ++i) rho_calc[i]=0.0;

	rho_solv.dimension(density.u1() , density.u2() , density.u3());
	for (int i=0; i<density.u1()*density.u2()*density.u3(); ++i) rho_solv[i]=1.0;

	bool use_Bs = pose_has_nonzero_Bs( pose );
	if (!use_Bs) {
		TR << "Input pose has no nonzero B factors ... setting to default." << std::endl;
	}

	for (int i=1 ; i<=(int)pose.size(); ++i) {
		std::string elt_i = pose[i].elt_;
		OneGaussianScattering sig_j = get_A( elt_i );
		core::Real k = sig_j.k( pose[i].B_ );
		if (use_Bs) {
			k = std::min ( k, 4*M_PI*M_PI/minimumB );
		} else {
			k = std::min ( k, 4*M_PI*M_PI/effectiveB );
		}
		core::Real C = sig_j.C( k );
		if ( C < 1e-6 ) continue;

		numeric::xyzVector< core::Real> cartX = pose[i].x_ - getTransform();
		numeric::xyzVector< core::Real> fracX = c2f*cartX;
		numeric::xyzVector< core::Real> atm_j, del_ij, atm_idx;
		atm_idx[0] = pos_mod (fracX[0]*grid[0] - origin[0] + 1 , (double)grid[0]);
		atm_idx[1] = pos_mod (fracX[1]*grid[1] - origin[1] + 1 , (double)grid[1]);
		atm_idx[2] = pos_mod (fracX[2]*grid[2] - origin[2] + 1 , (double)grid[2]);
		for (int z=1; z<=density.u3(); ++z) {
			atm_j[2] = z;
			del_ij[2] = (atm_idx[2] - atm_j[2]) / grid[2];
			if (del_ij[2] > 0.5) del_ij[2]-=1.0; if (del_ij[2] < -0.5) del_ij[2]+=1.0;
			del_ij[0] = del_ij[1] = 0.0;
			if ((f2c*del_ij).length_squared() > (ATOM_MASK+ATOM_MASK_PADDING)*(ATOM_MASK+ATOM_MASK_PADDING)) continue;
			for (int y=1; y<=density.u2(); ++y) {
				atm_j[1] = y;
				del_ij[1] = (atm_idx[1] - atm_j[1]) / grid[1] ;
				if (del_ij[1] > 0.5) del_ij[1]-=1.0; if (del_ij[1] < -0.5) del_ij[1]+=1.0;
				del_ij[0] = 0.0;
				if ((f2c*del_ij).length_squared() > (ATOM_MASK+ATOM_MASK_PADDING)*(ATOM_MASK+ATOM_MASK_PADDING)) continue;
				for (int x=1; x<=density.u1(); ++x) {
					atm_j[0] = x;
					del_ij[0] = (atm_idx[0] - atm_j[0]) / grid[0];
					if (del_ij[0] > 0.5) del_ij[0]-=1.0; if (del_ij[0] < -0.5) del_ij[0]+=1.0;
					numeric::xyzVector< core::Real > cart_del_ij = (f2c*del_ij);  // cartesian offset from (x,y,z) to atom_i
					core::Real d2 = (cart_del_ij).length_squared();

					if (d2 <= (ATOM_MASK+ATOM_MASK_PADDING)*(ATOM_MASK+ATOM_MASK_PADDING)) {
						core::Real atm = C*exp(-k*d2);
						rho_calc(x,y,z) += atm;
					}
					rho_solv(x,y,z) *= (d2<ATOM_MASK*ATOM_MASK ? 0:1); // should use a probe to determine solvent accessability
				}
			}
		}
	}

	// ffts
	numeric::fourier::fft3(rho_calc, Frho_calc);
	numeric::fourier::fft3(rho_solv, Frho_solv);

	// add solvent contribution to Fmodel
	Real ksol = 0.35, bsol = 46.0;  // TO DO! fit these
	int H,K,L;
	for (int z=1; z<=(int)density.u3(); ++z) {
		H = (z < (int)density.u3()/2) ? z-1 : z-density.u3() - 1;
		for (int y=1; y<=(int)density.u2(); ++y) {
			K = (y < (int)density.u2()/2) ? y-1 : y-density.u2()-1;
			for (int x=1; x<=(int)density.u1(); ++x) {
				L = (x < (int)density.u1()/2) ? x-1 : x-density.u1()-1;
				Real s2_i = S2(H,K,L);
				Frho_calc(x,y,z) += ksol*exp(-bsol*s2_i/4)*Frho_solv(x,y,z);
			}
		}
	}
}
*/

numeric::xyzVector<core::Real> ElectronDensity::delt_cart(
					 numeric::xyzVector<core::Real> const & cartX1,
					 numeric::xyzVector<core::Real> const & cartX2) {
	numeric::xyzVector<core::Real> fracX1, fracX2, delt_cart;

	numeric::xyzVector<core::Real> delt_fracX(c2f*(cartX1-cartX2));
	for (Size i=0; i<3; ++i) {
		while (delt_fracX[i] > 0.5)  delt_fracX[i]-=1.0;
		while (delt_fracX[i] < -0.5) delt_fracX[i]+=1.0;
	}
	delt_cart = f2c*delt_fracX;

	return delt_cart;
}

numeric::xyzVector<core::Real> ElectronDensity::get_nearest_UC(
															   numeric::xyzVector<core::Real> const & cartX_in,
															   numeric::xyzVector<core::Real> const & cartX_ref) {
	// find a copy of cartX_in in any unit cell that is closest to cartX_ref
	core::Real distance_squared = cartX_in.distance_squared(cartX_ref);
	numeric::xyzVector<core::Real> cartX_out(cartX_in);

	bool searching(true);
	while (searching) {
		searching = false;
		numeric::xyzVector<core::Real> fracX(c2f*cartX_out);

		numeric::xyzVector<core::Real> shift(-1,-1,-1);
		for (shift[0] = -1; shift[0] < 1.1; shift[0] += 1) {
			for (shift[1] = -1; shift[1] < 1.1; shift[1] += 1) {
				for (shift[2] = -1; shift[2] < 1.1; shift[2] += 1) {
					if (shift.length_squared() < 0.1) continue;

					numeric::xyzVector<core::Real> fracX_copy(fracX+shift);
					numeric::xyzVector<core::Real> cartX_copy = f2c*fracX_copy;

					if (cartX_copy.distance_squared(cartX_ref) < distance_squared) {
						searching = true;
						cartX_out = cartX_copy;
						distance_squared = cartX_copy.distance_squared(cartX_ref);
					}
				}
			}
		}
	}
	return cartX_out;
}


numeric::xyzVector<core::Real> ElectronDensity::get_cart_unitCell(
														  numeric::xyzVector<core::Real> const & cartX) {
	numeric::xyzVector<core::Real> fracX(c2f*cartX);
	for (Size i=0; i<3; ++i) {
		while (fracX[i] > 0.5)  fracX[i]-=1.0;
		while (fracX[i] < -0.5) fracX[i]+=1.0;
	}
	numeric::xyzVector<core::Real> cartX_uc = f2c*fracX;
	return cartX_uc;
}


/////////////////////////////////////
/// get Fdrho_d(xyz)
///      compute if not already computed
///      returns pointers to FArray
utility::vector1< ObjexxFCL::FArray3D< std::complex<double> > * > ElectronDensity::getFdrhoc( OneGaussianScattering S ) {
	int atmid = S.a();
	std::map< int,ObjexxFCL::FArray3D< std::complex<double> > >::const_iterator itx = Fdrhoc_dx.find( atmid );

	if (itx == Fdrhoc_dx.end()) {
		ObjexxFCL::FArray3D<double> drhoc_dx,drhoc_dy,drhoc_dz;
		numeric::xyzVector< core::Real > del_ij;

		core::Real k = S.k( effectiveB );
		core::Real C = S.C( k );

		drhoc_dx.dimension(p_grid[0], p_grid[1], p_grid[2]);
		drhoc_dy.dimension(p_grid[0], p_grid[1], p_grid[2]);
		drhoc_dz.dimension(p_grid[0], p_grid[1], p_grid[2]);
		for (int z=1; z<=(int)p_grid[2]; ++z) {
			if (z < (int)p_grid[2]/2)
				del_ij[2] = ((core::Real)z - 1.0) / grid[2];
			else
				del_ij[2] = (((core::Real)z - p_grid[2] - 1.0)) / grid[2];
			for (int y=1; y<=(int)p_grid[1]; ++y) {
				if (y < (int)p_grid[1]/2)
					del_ij[1] = ((core::Real)y - 1.0) / grid[1] ;
				else
					del_ij[1] = (((core::Real)y - p_grid[1] - 1.0)) / grid[1];
				for (int x=1; x<=(int)p_grid[0]; ++x) {
					if (x < (int)p_grid[0]/2)
						del_ij[0] = ((core::Real)x - 1.0) / grid[0];
					else
						del_ij[0] = (((core::Real)x - p_grid[0] - 1.0)) / grid[0];

					numeric::xyzVector< core::Real > cart_del_ij = (f2c*del_ij);
					core::Real d2 = (cart_del_ij).length_squared();

					drhoc_dx(x,y,z) = cart_del_ij[0] * 4*C*k*exp( -k*d2 );
					drhoc_dy(x,y,z) = cart_del_ij[1] * 4*C*k*exp( -k*d2 );
					drhoc_dz(x,y,z) = cart_del_ij[2] * 4*C*k*exp( -k*d2 );
				}
			}
		}

		// ffts
		Fdrhoc_dx[atmid] = ObjexxFCL::FArray3D< std::complex<double> >();
		Fdrhoc_dy[atmid] = ObjexxFCL::FArray3D< std::complex<double> >();
		Fdrhoc_dz[atmid] = ObjexxFCL::FArray3D< std::complex<double> >();

		numeric::fourier::fft3(drhoc_dx, Fdrhoc_dx[atmid]);
		numeric::fourier::fft3(drhoc_dy, Fdrhoc_dy[atmid]);
		numeric::fourier::fft3(drhoc_dz, Fdrhoc_dz[atmid]);
	}

	utility::vector1< ObjexxFCL::FArray3D< std::complex<double> > * > retval;
	retval.push_back( &Fdrhoc_dx[atmid] );
	retval.push_back( &Fdrhoc_dy[atmid] );
	retval.push_back( &Fdrhoc_dz[atmid] );

	return retval;
}

/////////////////////////////////////
/// setup oversampled maps for fast density scoring
void ElectronDensity::setup_fastscoring_first_time(core::pose::Pose const &pose) {
	fastgrid = grid;
	fastorigin = origin;

	utility::vector1< core::Real > fastdens_params;
	if (basic::options::option[ basic::options::OptionKeys::edensity::fastdens_params ].user()) {
		fastdens_params = basic::options::option[ basic::options::OptionKeys::edensity::fastdens_params ]();
		runtime_assert( fastdens_params.size() == 2 );
	} else {
		fastdens_params.push_back( 0.4 );
		fastdens_params.push_back( 0.6 );
	}

	ObjexxFCL::FArray3D< double > rhoc, fastdens_score_i;
	ObjexxFCL::FArray3D< std::complex<double> > Frhoo, Frhoc;
	rhoc.dimension(fastgrid[0], fastgrid[1], fastgrid[2]);

	// atom count
	core::Size natms=0,nres=0;
	for (core::Size i=1; i<=pose.total_residue(); ++i) {
		if ( pose.residue(i).aa() == core::chemical::aa_vrt ) continue; nres++;
		core::conformation::Residue const& rsd_i = pose.residue(i);
		natms += rsd_i.nheavyatoms();
	}

	fastdens_score.dimension( density.u1(), density.u2(), density.u3(), nkbins_);

	core::Real max_val, min_val;

	// compute k limits (corresponding to B=0 and B=1000)
	OneGaussianScattering S_C = get_A( "C" );
	OneGaussianScattering S_O = get_A( "O" );
	OneGaussianScattering S_N = get_A( "N" );
	OneGaussianScattering S_S = get_A( "S" );

	kmin_ = S_C.k( 1000 ); // "infinite"
	kmax_ = std::max( std::max( std::max (S_S.k( 0 ), S_O.k( 0 )) , S_N.k( 0 )), S_C.k( 0 ) );
	kmax_ = std::min( kmax_, 4*M_PI*M_PI/(minimumB) );
	TR << "Setting [kmin_,kmax_] to [" << kmin_ << "," << kmax_ << "]" << std::endl;
	if (nkbins_ > 1) {
		kstep_ = (kmax_-kmin_) / (nkbins_-1);
	} else {
		kstep_ = 0.0;
	}

	// transform + standardize density map
	{
		ObjexxFCL::FArray3D< double > dens_std;
		dens_std.dimension( density.u1(), density.u2(), density.u3() );
		for (int i=0; i<density.u1()*density.u2()*density.u3(); ++i) {
			dens_std[i] = (density[i]-dens_mean)/dens_stdev;
		}
		numeric::fourier::fft3(dens_std, Frhoo);
	}

	for (int kbin=nkbins_; kbin>=1; --kbin) {
		rhoc = 0.0;
		numeric::xyzVector< core::Real > del_ij;
		core::Real k = (kbin-1)*kstep_ + kmin_;

		if (nkbins_ == 1) {
			OneGaussianScattering S = get_A( "C" );
			numeric::xyzVector< core::Real > del_ij;
			k = std::min ( S.k(0), 4*M_PI*M_PI/effectiveB );
		}

		if (k>0) {
			k = std::min ( k, 4*M_PI*M_PI/minimumB );

			core::Real ATOM_MASK_SQ = (18/k);  // very generous mask (rho<6e-6)
			core::Real FIXED_MASK_SQ = (3*3);
			core::Real C = pow(k, 1.5);
			core::Real VV = voxel_volume();

			core::Real rho_sum = 0.0, rho2_sum = 0.0;
			for (int z=1; z<=(int)fastgrid[2]; ++z) {
				if (z < (int)fastgrid[2]/2)
					del_ij[2] = ((core::Real)z - 1.0) / fastgrid[2];
				else
					del_ij[2] = (((core::Real)z - fastgrid[2] - 1.0)) / fastgrid[2];
				del_ij[0] = del_ij[1] = 0.0;
				if ((f2c*del_ij).length_squared() > ATOM_MASK_SQ) continue;  // early exit
				for (int y=1; y<=(int)fastgrid[1]; ++y) {
					if (y < (int)fastgrid[1]/2)
						del_ij[1] = ((core::Real)y - 1.0) / fastgrid[1] ;
					else
						del_ij[1] = (((core::Real)y - fastgrid[1] - 1.0)) / fastgrid[1];
					del_ij[0] = 0.0;
					if ((f2c*del_ij).length_squared() > ATOM_MASK_SQ) continue;  // early exit
					for (int x=1; x<=(int)fastgrid[0]; ++x) {
						if (x < (int)fastgrid[0]/2)
							del_ij[0] = ((core::Real)x - 1.0) / fastgrid[0];
						else
							del_ij[0] = (((core::Real)x - fastgrid[0] - 1.0)) / fastgrid[0];
							numeric::xyzVector< core::Real > cart_del_ij = (f2c*del_ij);
						core::Real d2 = (cart_del_ij).length_squared();
						if (d2 <= ATOM_MASK_SQ) {
							rhoc(x,y,z) = C*VV*exp(-k*d2);
						}
					}
				}
			}

			// standardize density & convert to an approximate CC
			// this has been minimally tuned but might be improved
			for (int i=0; i<fastgrid[0]*fastgrid[1]*fastgrid[2]; ++i) {
				rhoc[i] = rhoc[i] * (std::pow( k, -fastdens_params[1] ) - fastdens_params[2]);
			}
		}

		// ffts
		numeric::fourier::fft3(rhoc, Frhoc);
		Frhoc(1,1,1) = Frhoo(1,1,1) = 0.0;

		TR << "Bin " << kbin << ":  B(C/N/O/S)=" << S_C.B(k) << " / " << S_N.B(k) << " / " << S_O.B(k) << " / " << S_S.B(k) << "  sum=" << Frhoc(1,1,1) << std::endl;

		// convolute
		Frhoo(1,1,1) = 0.0;
		for (int i=1; i<=density.u1(); i++)
		for (int j=1; j<=density.u2(); j++)
		for (int k=1; k<=density.u3(); k++) {
			Frhoc(i,j,k) *= ( Frhoo(i,j,k) );
		}
		numeric::fourier::ifft3( Frhoc , fastdens_score_i );

		// copy to big array
		for (int i=1; i<=density.u1(); i++)
		for (int j=1; j<=density.u2(); j++)
		for (int k=1; k<=density.u3(); k++) {
 			fastdens_score(i,j,k,kbin) = fastdens_score_i(i,j,k);
		}

		// normalization
		//    [0.5% -- 99.5%] in min temp bin
		if (kbin == nkbins_) {
			std::vector<core::Real> scores_i( density.u1()*density.u2()*density.u3(),0.0 );
			for (int i=0; i<density.u1()*density.u2()*density.u3(); ++i)
				scores_i[i] = fastdens_score_i[i];

			core::Size cutat = legacy_? 1:5;
			std::nth_element (scores_i.begin(), scores_i.begin()+cutat-1, scores_i.end());
			std::nth_element (scores_i.begin(), scores_i.end()-cutat, scores_i.end());
			min_val = scores_i[cutat-1];
			max_val = scores_i[scores_i.size()-cutat];
		}

		if (basic::options::option[ basic::options::OptionKeys::edensity::debug ]()) {
			core::Real scalefactor = (natms/nres);
			core::Real mu=0.5*(max_val+min_val), sigma=0.5*scalefactor*(max_val-min_val);
			for (int i=1; i<=density.u1(); i++)
			for (int j=1; j<=density.u2(); j++)
			for (int k=1; k<=density.u3(); k++) {
				fastdens_score_i(i,j,k) = (fastdens_score(i,j,k,kbin)-mu)/sigma;
			}
			std::ostringstream oss; oss << "fastdens" << kbin << ".mrc";
			ElectronDensity(fastdens_score_i,1.0, numeric::xyzVector< core::Real >(0,0,0), false ).writeMRC( oss.str() );
		}
	}

	core::Real scalefactor = (natms/nres);
	core::Real mu=0.0, sigma=0.5*scalefactor*(max_val-min_val);
	if (legacy_) mu = 0.5*(max_val-min_val);
	for (int i=0; i<nkbins_*density.u1()*density.u2()*density.u3(); ++i) {
		fastdens_score[i] = (fastdens_score[i]-mu)/sigma;
	}

	// spline coeffs
	ObjexxFCL::FArray4D< double > temp_coeffs;
	spline_coeffs( fastdens_score, temp_coeffs);
	fastdens_score = temp_coeffs;
}


// use a pose to rescale the fastscoring bins so they correspond to RSCC
void ElectronDensity::rescale_fastscoring_temp_bins(core::pose::Pose const &pose, bool initBs) {
	if ( fastdens_score.u1()*fastdens_score.u2()*fastdens_score.u3()*fastdens_score.u4() == 0 )
		setup_fastscoring_first_time(pose);

	if (nkbins_ == 1) return;

	// pose->poseCoords
	poseCoords litePose;
	utility::vector1< core::Real> allBs;
	core::conformation::symmetry::SymmetryInfoCOP symm_info;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		core::conformation::symmetry::SymmetricConformation const & SymmConf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation const &> ( pose.conformation()) );
		symm_info = SymmConf.Symmetry_Info();
	}
	for (int i=1; i<=pose.total_residue(); ++i) {
		if (!remap_symm_ && symm_info && !symm_info->bb_is_independent( i ) ) continue;
		core::conformation::Residue const & rsd_i ( pose.residue(i) );
		if ( rsd_i.aa() == core::chemical::aa_vrt ) continue;

		core::Size natoms = rsd_i.nheavyatoms();
		for (int j=1; j<=natoms; ++j) {
			core::conformation::Atom const &atom_j( rsd_i.atom(j) );
			core::chemical::AtomTypeSet const & atom_type_set( rsd_i.atom_type_set() );

			poseCoord coord_j;
			coord_j.x_ = rsd_i.xyz( j );
			coord_j.B_ = 0;
			coord_j.elt_ = atom_type_set[ rsd_i.atom_type_index( j ) ].element();
			litePose.push_back( coord_j );

			core::Real B = 0.0;
			if (!initBs) B = pose.pdb_info()->temperature( i, j );
			allBs.push_back( B );
		}
	}

	// get avg K
	core::Real meanK=0.0;
	if (!initBs) {
		for (int i=1; i<=litePose.size(); ++i ) {
			std::string &elt_i = litePose[i].elt_;
			OneGaussianScattering sig_j = get_A( elt_i );
			meanK += sig_j.k( allBs[i] );
		}
		meanK /= litePose.size();
	}

	for (int RPT=1; RPT<=5; ++RPT) {
		core::Size ref=1;
		core::Real RSCC_ref=-999;
		utility::vector1< core::Real > RSCC(nkbins_), overlap(nkbins_);
		numeric::xyzVector< core::Real > fracX, idxX;
		for (int kbin=nkbins_; kbin>=1; --kbin) {
			core::Real k = (kbin-1)*kstep_ + kmin_;
			if (k==0) continue;

			overlap[kbin] = 0.0;
			for (int i=1; i<=litePose.size(); ++i ) {
				std::string &elt_i = litePose[i].elt_;
				OneGaussianScattering sig_j = get_A( elt_i );

				core::Real k_tgt = k;
				if (!initBs) k_tgt += sig_j.k( allBs[i] ) - meanK;
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
			RSCC[kbin] = getRSCC(litePose);
			if (RSCC[kbin] > RSCC_ref) {
				RSCC_ref = RSCC[kbin];
				ref = kbin;
			}
		}

		// rescale
		for (int kbin=nkbins_; kbin>=1; --kbin) {
			core::Real k = (kbin-1)*kstep_ + kmin_;
			if (k==0) continue;

			core::Real scale_i = RSCC[kbin]*overlap[ref] / (RSCC[ref]*overlap[kbin]);
			if (RPT == 5) TR << "BIN" << kbin << " : [" << k << "/" << meanK << "] " << RSCC[kbin] << "  " << overlap[kbin] << "  " << scale_i << std::endl;

			for (int i=1; i<=density.u1(); i++)
			for (int j=1; j<=density.u2(); j++)
			for (int k=1; k<=density.u3(); k++) {
				fastdens_score(i,j,k,kbin) *= scale_i;
			}
		}
	}
}



/////////////////////////////////////
/// Precomputes a bunch of stuff for patterson map scoring
void ElectronDensity::setup_patterson_first_time(core::pose::Pose const &pose) {
	bool is_set=false;
	int nres = pose.total_residue(); //reses.size();
	const core::Real FLUFF = 10.0;
	numeric::xyzVector< core::Real > CoM(0,0,0);
	core::Real nCoM=0;

	for (int i=1 ; i<=nres; ++i) {
		conformation::Residue const &rsd_i (pose.residue(i));
		if ( (rsd_i.aa() == core::chemical::aa_vrt) || (scoring_mask_.find(i) != scoring_mask_.end()) ) continue;
		int nheavyatoms = rsd_i.nheavyatoms();
		for (int j=1 ; j<=nheavyatoms; ++j) {
			numeric::xyzVector< core::Real > const &xyz_ij = rsd_i.atom(j).xyz();
			if (is_missing_density( xyz_ij )) continue;
			if (!is_set) {
				d_min = d_max = xyz_ij;
				is_set = true;
			}
			d_min[0] = std::min(d_min[0],xyz_ij[0]);
			d_min[1] = std::min(d_min[1],xyz_ij[1]);
			d_min[2] = std::min(d_min[2],xyz_ij[2]);
			d_max[0] = std::max(d_max[0],xyz_ij[0]);
			d_max[1] = std::max(d_max[1],xyz_ij[1]);
			d_max[2] = std::max(d_max[2],xyz_ij[2]);
			CoM+=xyz_ij;
			nCoM+=1;
		}
	}
	CoM /= nCoM;
	d_max -= CoM;
	d_min -= CoM;

	// extent <- range + minR
	p_extent = d_max - d_min + PattersonMaxR + 2*FLUFF;

	// use current grid spacing to find grid  in each dim
	//    grid has no prime factors >5 so that fft is fast
	//    grid is even so we can use real->complex fft shortcut (TO DO)
	//    mutiples depend on input spacegroup
	numeric::xyzVector< core::Real > f_extent = c2f*p_extent;
	p_grid[0] = findSampling5(f_extent[0]*grid[0] , MINMULT[0]);
	p_grid[1] = findSampling5(f_extent[1]*grid[1] , MINMULT[1]);
	p_grid[2] = findSampling5(f_extent[2]*grid[2] , MINMULT[2]);

	// center the molecule in our new unit cell
	f_extent = c2f*(d_min-FLUFF-PattersonMaxR/2.0);
	p_origin[0] = std::floor(f_extent[0]*grid[0]);
	p_origin[1] = std::floor(f_extent[1]*grid[1]);
	p_origin[2] = std::floor(f_extent[2]*grid[2]);

	// add fluff
	d_min = d_min-FLUFF;
	d_max = d_max+FLUFF;

	// DEBUGGING -- run p_o and p_c on same grid
	p_grid = grid;
	p_origin = origin;
	TR << "Setting up patterson map grid " << p_grid[0] << "x" << p_grid[1] << "x" << p_grid[2] << std::endl;
	TR << "                       origin " << p_origin[0] << "x" << p_origin[1] << "x" << p_origin[2] << std::endl;
	TR << "                       d_min  " << d_min[0] << "x" << d_min[1] << "x" << d_min[2] << std::endl;
	TR << "                       d_max  " << d_max[0] << "x" << d_max[1] << "x" << d_max[2] << std::endl;

	// setup precomputed symmetric matrices
	// for each matrix symmptr_i = symmptr[i]
	//    symmptr_i[x] = idx of symm_i(x)
	core::Size nsymmrot = symmRotOps.size();
	numeric::xyzVector< core::Real > frac_xyz;
	if (nsymmrot > 1 && !basic::options::option[ basic::options::OptionKeys::patterson::dont_use_symm_in_pcalc ]() ) {
		TR << "... precomputing " << nsymmrot << " rotations" << std::endl;
		symm_ptrs.resize(nsymmrot);

		// init arrays
		for (int i=1; i<= (int)nsymmrot; ++i) {
			symm_ptrs[i].dimension(p_grid[0], p_grid[1], p_grid[2]);
		}

		// build matrices
		for (int z=1; z<=(int)p_grid[2]; ++z) {
			if (z < (int)p_grid[2]/2) frac_xyz[2] = ((core::Real)z - 1.0) / grid[2];
			else                 frac_xyz[2] = (((core::Real)z - p_grid[2] - 1.0)) / grid[2];
			for (int y=1; y<=(int)p_grid[1]; ++y) {
				if (y < (int)p_grid[1]/2) frac_xyz[1] = ((core::Real)y - 1.0) / grid[1] ;
				else                 frac_xyz[1] = (((core::Real)y - p_grid[1] - 1.0)) / grid[1];
				for (int x=1; x<=(int)p_grid[0]; ++x) {
					if (x < (int)p_grid[0]/2) frac_xyz[0] = ((core::Real)x - 1.0) / grid[0];
					else                 frac_xyz[0] = (((core::Real)x - p_grid[0] - 1.0)) / grid[0];
					for (int i=1; i<= (int)nsymmrot; ++i) {
						// frac coords of symm copy
						numeric::xyzVector< core::Real > symm_i_xyz = symmRotOps[i]*frac_xyz;
						// idx of symm copy
						numeric::xyzVector< core::Size > symm_i_idx (
							pos_mod( (int)std::floor(symm_i_xyz[0]*grid[0]+0.5), p_grid[0] ),
							pos_mod( (int)std::floor(symm_i_xyz[1]*grid[1]+0.5), p_grid[1] ),
							pos_mod( (int)std::floor(symm_i_xyz[2]*grid[2]+0.5), p_grid[2] ) );
						symm_ptrs[i](x,y,z) = symm_i_idx[2]*p_grid[1]*p_grid[0] + symm_i_idx[1]*p_grid[0] + symm_i_idx[0];
					}
				}
			}
		}
	}

	/// compute patterson-space mask, compute S^2 array
	numeric::xyzVector< core::Real > del_ij;
	PattersonEpsilon.dimension(p_grid[0], p_grid[1], p_grid[2]);
	F_s2.dimension(p_grid[0], p_grid[1], p_grid[2]);
	core::Real eps_sum = 0;
	core::Real minS=999.0, maxS=0.0;
	int H,K,L;
	for (int z=1; z<=(int)p_grid[2]; ++z) {
		H = (z < (int)p_grid[2]/2) ? z-1 : z-p_grid[2] - 1;
		if (z < (int)p_grid[2]/2)
			del_ij[2] = ((core::Real)z - 1.0) / grid[2];
		else
			del_ij[2] = (((core::Real)z - p_grid[2] - 1.0)) / grid[2];
		for (int y=1; y<=(int)p_grid[1]; ++y) {
			K = (y < (int)p_grid[1]/2) ? y-1 : y-p_grid[1]-1;
			if (y < (int)p_grid[1]/2)
				del_ij[1] = ((core::Real)y - 1.0) / grid[1] ;
			else
				del_ij[1] = (((core::Real)y - p_grid[1] - 1.0)) / grid[1];
			for (int x=1; x<=(int)p_grid[0]; ++x) {
				L = (x < (int)p_grid[0]/2) ? x-1 : x-p_grid[0]-1;
				if (x < (int)p_grid[0]/2)
					del_ij[0] = ((core::Real)x - 1.0) / grid[0];
				else
					del_ij[0] = (((core::Real)x - p_grid[0] - 1.0)) / grid[0];

				numeric::xyzVector< core::Real > cart_del_ij = (f2c*del_ij);
				core::Real d2 = (cart_del_ij).length_squared();

				if (PattersonMinR <= 0) {
					PattersonEpsilon(x,y,z) = 1/(1+exp( d2 - PattersonMaxR*PattersonMaxR));
				} else {
					PattersonEpsilon(x,y,z) = 1/(1+exp( d2 - PattersonMaxR*PattersonMaxR))
											* 1/(1+exp( PattersonMinR*PattersonMinR - d2));
				}
				eps_sum += PattersonEpsilon(x,y,z);
				F_s2(x,y,z) = sqrt(S2(H,K,L));
				minS = std::min( minS, F_s2(x,y,z));
				maxS = std::max( maxS, F_s2(x,y,z));
			}
		}
	}

	// p_o resampled on p_c grid
	//ObjexxFCL::FArray3D< double >  eps_po;
	po_bar = 0.0;
	p_o.dimension( p_grid[0], p_grid[1], p_grid[2] );
	for (int z=1; z<=(int)p_grid[2]; ++z) {
		int dens_z = z;
		if (z > (int)p_grid[2]/2) dens_z = (int)grid[2] - (int)p_grid[2] + z;
		if (dens_z<=0 || dens_z > grid[2]) {
			continue; // patterson mask is bigger than unit cell
		}
		for (int y=1; y<=(int)p_grid[1]; ++y) {
			int dens_y = y;
			if (y>(int)p_grid[1]/2) dens_y = (int)grid[1] - (int)p_grid[1] + y;
			if (dens_y<=0 || dens_y > grid[1]) {
				continue; // patterson mask is bigger than unit cell
			}
			for (int x=1; x<=(int)p_grid[0]; ++x) {
				int dens_x = x;
				if (x>(int)p_grid[0]/2) dens_x = (int)grid[0] - (int)p_grid[0] + x;
				if (dens_x<=0 || dens_x > grid[0]) {
					continue; // patterson mask is bigger than unit cell
				}
				p_o(x,y,z) = density(dens_x,dens_y,dens_z);
			}
		}
	}

	////////////////////////
	////////////////////////
	// apply Ecalc to input patterson
	// bin reflections, normalize so that sum(F^2)=1 in each bin
	// truncate at -patterson::resolution_cutoffs
	// bin after truncating
	// apply sigmaA
	lowres_cut = 999.0;
	hires_cut = 0.01;
	if (basic::options::option[ basic::options::OptionKeys::patterson::resolution_cutoffs ].user() ) {
		utility::vector1<core::Real> rescut_in = basic::options::option[ basic::options::OptionKeys::patterson::resolution_cutoffs ]();
		if (rescut_in.size() == 1)
			hires_cut = rescut_in[1];
		if (rescut_in.size() >= 2) {
			hires_cut = std::min(rescut_in[1],rescut_in[2]);
			lowres_cut = std::max(rescut_in[1],rescut_in[2]);
		}
	}

	// find resolution bin cutoffs
	core::Size nbuckets = basic::options::option[ basic::options::OptionKeys::patterson::nshells ]();
	core::Real max_S3 = std::pow( std::min( 1/hires_cut , maxS) , 3.0);
	core::Real min_S3 = std::pow( 1/lowres_cut , 3.0);

	// now categorize each point
	// save these for normalizing Pcalc's later
	bucket_id.dimension( p_grid[0], p_grid[1], p_grid[2] );
	for (int i=0; i< (int)(p_grid[0]*p_grid[1]*p_grid[2]) ; ++i) {
		// buckets = 1+floor(nb*(S2.^(3/2)-minS^3)/(maxS^3-minS^3));
		int bucket_i = 1+(int)std::floor( nbuckets*( std::pow( F_s2[i], 3.0 ) - min_S3 ) / (max_S3 - min_S3) );

		if ( bucket_i <= 0 || bucket_i > (int)nbuckets )
			bucket_id[i]=(core::Size)0;
		else
			bucket_id[i]=(core::Size)bucket_i;
	}

	if ( !basic::options::option[ basic::options::OptionKeys::patterson::no_ecalc ]() ) {
		if (basic::options::option[ basic::options::OptionKeys::edensity::debug ]()) {
			ElectronDensity(p_o,1.0,numeric::xyzVector< core::Real >(0,0,0), true).writeMRC( "p_obs_unnormalized.mrc");
		}

		// per-bucket normalization
		ObjexxFCL::FArray3D< std::complex<double> > Fp_o;
		numeric::fourier::fft3(p_o, Fp_o);
		bucket_counts.resize(nbuckets,0);
		utility::vector1<core::Real> F2(nbuckets, 0);
		for (int i=0; i<(int)(p_grid[0]*p_grid[1]*p_grid[2]) ; ++i) {
			if (bucket_id[i] > 0) {
				F2[ bucket_id[i] ] += std::abs(Fp_o[i]);  // should be norm?
				bucket_counts[ bucket_id[i] ]++;
			}
		}
		//std::cerr << "HIRES CUT " << hires_cut << "  LOWRES CUT " << lowres_cut << std::endl;
		//std::cerr << "MAX_S " << std::pow( max_S3 , 1.0/3.0 ) << "  MIN_S " << std::pow( min_S3 , 1.0/3.0 ) << std::endl;
		for (int i=1; i<=(int)nbuckets; ++i) {
			F2[i] /= bucket_counts[ i ];
		}

		// multiply F^2 by sigmaA^2
		core::Real rmsd = basic::options::option[ basic::options::OptionKeys::patterson::rmsd ]();
		core::Real fsol = basic::options::option[ basic::options::OptionKeys::patterson::Fsol ]();
		core::Real bsol = basic::options::option[ basic::options::OptionKeys::patterson::Bsol ]();
		for (int i=0; i<(int)(p_grid[0]*p_grid[1]*p_grid[2]) ; ++i) {
			if (bucket_id[i] > 0) {
				core::Real S2_i = F_s2[i]*F_s2[i];
				core::Real sigmaa_i = sqrt( (1-fsol*std::exp(-bsol*S2_i))*std::exp(-8*(M_PI*M_PI/3)*rmsd*rmsd*S2_i) );
				Fp_o[i] = sigmaa_i * (std::abs(Fp_o[i]) / F2[ bucket_id[i] ]);
			} else {
				Fp_o[i] = 0;
			}
		}

		// back to patterson space
		numeric::fourier::ifft3(Fp_o, p_o);
	}

	// standardize
	// core::Real delta,mean=0,M2=0;
	// for (int i=0; i< (int)(p_grid[0]*p_grid[1]*p_grid[2]) ; ++i) {
	// 	delta = p_o[i] - mean;
	// 	mean = mean + delta/(i+1);
	// 	M2 = M2 + delta*(p_o[i] - mean);
	// }
	// M2 = sqrt(M2/((int)p_grid[0]*p_grid[1]*p_grid[2]));
	// for (int i=0; i< (int)(p_grid[0]*p_grid[1]*p_grid[2]) ; ++i) {
	// 	p_o[i] = (p_o[i]-mean)/M2;
	// 	po_bar += p_o[i] * PattersonEpsilon[i];
	// }
	// po_bar /= eps_sum;

	// normalization in reciprocal space now
	po_bar = 0.0;

	if (basic::options::option[ basic::options::OptionKeys::edensity::debug ]()) {
		ElectronDensity(p_o,1.0,numeric::xyzVector< core::Real >(0,0,0), true).writeMRC( "p_obs.mrc" );
	}

	// collect scatterers
	for (int i=1 ; i<=nres; ++i) {
		conformation::Residue const &rsd_i (pose.residue(i));
		if ( (rsd_i.aa() == core::chemical::aa_vrt) || (scoring_mask_.find(i) != scoring_mask_.end()) ) continue;
		for (int j=1 ; j<=(int)rsd_i.nheavyatoms(); ++j) {
			// ensure that this element exists in the precomputed scatterers
			getFdrhoc( get_A( rsd_i.atom_type_set()[ rsd_i.atom_type_index( j ) ].element() ) );
		}
	}
}


/////////////////////////////////////
/// Match pose to patterson map
core::Real ElectronDensity::matchPoseToPatterson(
		core::pose::Pose const &pose,
        bool cacheCCs /*=false */ ) {
	using namespace numeric::statistics;

	// make sure map is loaded
	if (!isLoaded) {
		TR << "[ ERROR ]  ElectronDensity::matchPose called but no map is loaded!\n";
		return 0.0;
	}

	int nres = pose.total_residue(); //reses.size();
	numeric::xyzVector< core::Real > cartX, fracX;
	numeric::xyzVector< core::Real > atm_i, atm_j, del_ij, atm_idx_ij;

	core::Real SC_scaling = basic::options::option[ basic::options::OptionKeys::patterson::sc_scaling ]();

	///////////////////////////
	//// first time this function is called, we need to set up a bunch of stuff
	//// if the pose has moved a lot we need to set up a new grid
	// quick check of bounding box, CoM
	// bool is_set = false;
	// core::Real nCoM=0;
	numeric::xyzVector< core::Real > newd_min, newd_max, CoM(0,0,0);
	// for (int i=1 ; i<=nres; ++i) {
	// 	conformation::Residue const &rsd_i (pose.residue(i));
	// 	if ( (rsd_i.aa() == core::chemical::aa_vrt) || (scoring_mask_.find(i) != scoring_mask_.end()) ) continue;
	// 	int nheavyatoms = rsd_i.nheavyatoms();
	// 	for (int j=1 ; j<=nheavyatoms; ++j) {
	// 		numeric::xyzVector< core::Real > const &xyz_ij = rsd_i.atom(j).xyz();
	// 		if (is_missing_density( xyz_ij )) continue;
	// 		if (!is_set) {
	// 			newd_min = newd_max = xyz_ij;
	// 			is_set = true;
	// 		}
	// 		newd_min[0] = std::min(newd_min[0],xyz_ij[0]);
	// 		newd_min[1] = std::min(newd_min[1],xyz_ij[1]);
	// 		newd_min[2] = std::min(newd_min[2],xyz_ij[2]);
	// 		newd_max[0] = std::max(newd_max[0],xyz_ij[0]);
	// 		newd_max[1] = std::max(newd_max[1],xyz_ij[1]);
	// 		newd_max[2] = std::max(newd_max[2],xyz_ij[2]);
	// 		CoM+=xyz_ij;
	// 		nCoM+=1;
	// 	}
	// }
	// CoM /= nCoM;
	// newd_max -= CoM;
	// newd_min -= CoM;

	// save CoM for repacking
	//p_CoM = CoM;

	bool needToSetup = (p_extent.length_squared() == 0);
	p_CoM = numeric::xyzVector< core::Real >(0,0,0);
	//for (int j=0; j<3; ++j) {
	//	needToSetup |= newd_min[j] < d_min[j];
	//	needToSetup |= newd_max[j] > d_max[j];
	//}
	// do not let the grid change while minimizing
	if (needToSetup && ! pose.energies().use_nblist() ) {
		setup_patterson_first_time(pose);
	}

	// set up arrays sized by p_grid
	rho_calc.dimension(p_grid[0], p_grid[1], p_grid[2]);

	for (int i=0; i<(int)(p_grid[0]*p_grid[1]*p_grid[2]); ++i) {
		rho_calc[i] = 0.0;
	}

	///////////////////////////
	/// COMPUTE RHO_C, MASK
	///    remember atoms that were used to calculate rho_c
	rho_calc_atms.resize(nres);
	rho_calc_as.resize(nres);
	rho_calc_sum = 0.0;
	for (int i=1 ; i<=nres; ++i) {
		conformation::Residue const &rsd_i (pose.residue(i));

		// skip vrts & masked reses
		if ( (rsd_i.aa() == core::chemical::aa_vrt) || (scoring_mask_.find(i) != scoring_mask_.end()) ) continue;

		int nheavyatoms = rsd_i.nheavyatoms();
		rho_calc_atms[i].resize(nheavyatoms);
		rho_calc_as[i].resize(nheavyatoms);

		for (int j=1 ; j<=nheavyatoms; ++j) {
			conformation::Atom const &atm_i( rsd_i.atom(j) ); //(pose.residue(i).atom("CA"));
			chemical::AtomTypeSet const & atom_type_set( rsd_i.atom_type_set() );
			std::string elt_i = atom_type_set[ rsd_i.atom_type_index( j ) ].element();
			OneGaussianScattering sig_j = get_A( elt_i );
			core::Real k = sig_j.k( effectiveB );
			core::Real C = sig_j.C( k );

			if ( (Size) j > rsd_i.last_backbone_atom()+1)  // count CB as bb atom
				C *= SC_scaling;

			rho_calc_atms[i][j] = atm_i.xyz()-CoM;
			rho_calc_as[i][j] = sig_j;

			// skip randomized residues / residues with no mass (e.g. centroids)
			if ( is_missing_density( atm_i.xyz() ) || C < 1e-6) {
				rho_calc_as[i][j] = OneGaussianScattering(); // weight == 0
				continue;
			}

			cartX = atm_i.xyz()-CoM;
			fracX = c2f*cartX;
			// atm_idx_ij[0] = fracX[0]*grid[0] - p_origin[0] + 1;
			// atm_idx_ij[1] = fracX[1]*grid[1] - p_origin[1] + 1;
			// atm_idx_ij[2] = fracX[2]*grid[2] - p_origin[2] + 1;
			atm_idx_ij[0] = pos_mod (fracX[0]*p_grid[0] - p_origin[0] + 1 , (double)p_grid[0]);
			atm_idx_ij[1] = pos_mod (fracX[1]*p_grid[1] - p_origin[1] + 1 , (double)p_grid[1]);
			atm_idx_ij[2] = pos_mod (fracX[2]*p_grid[2] - p_origin[2] + 1 , (double)p_grid[2]);

			for (int z=1; z<=(int)p_grid[2]; ++z) {
				del_ij[2] = (atm_idx_ij[2] - z) / grid[2];
				if (del_ij[2] > 0.5) del_ij[2]-=1.0;
				if (del_ij[2] < -0.5) del_ij[2]+=1.0;
				del_ij[0] = del_ij[1] = 0.0;
				if ((f2c*del_ij).length_squared() > square(ATOM_MASK+1)) continue;
				for (int y=1; y<=(int)p_grid[1]; ++y) {
					del_ij[1] = (atm_idx_ij[1] - y) / grid[1] ;
					if (del_ij[1] > 0.5) del_ij[1]-=1.0;
					if (del_ij[1] < -0.5) del_ij[1]+=1.0;
					del_ij[0] = 0.0;
					if ((f2c*del_ij).length_squared() > square(ATOM_MASK+1)) continue;
					for (int x=1; x<=(int)p_grid[0]; ++x) {
						del_ij[0] = (atm_idx_ij[0] - x) / grid[0];
						if (del_ij[0] > 0.5) del_ij[0]-=1.0;
						if (del_ij[0] < -0.5) del_ij[0]+=1.0;

						numeric::xyzVector< core::Real > cart_del_ij = (f2c*del_ij);  // cartesian offset from (x,y,z) to atom_i
						core::Real d2 = (cart_del_ij).length_squared();
						if (d2 <= (ATOM_MASK+1)*(ATOM_MASK+1)) {
							core::Real atm = C*exp(-k*d2);
							rho_calc(x,y,z) += atm;
							rho_calc_sum+= atm;
						}
					}
				}
			}
		}
	}

	// rho_c -> patterson
	for (int i=0; i< (int)(p_grid[0]*p_grid[1]*p_grid[2]) ; ++i) {
		rho_calc[i] /= rho_calc_sum;
	}
	numeric::fourier::fft3(rho_calc, Frho_calc);

	ObjexxFCL::FArray3D< std::complex<double> > Fcalc2;
	Fcalc2.dimension( p_grid[0],p_grid[1],p_grid[2] );

	if ( !basic::options::option[ basic::options::OptionKeys::patterson::no_ecalc ]() ) {
		if (basic::options::option[ basic::options::OptionKeys::edensity::debug ]()) {
			for (int i=0; i<(int)(p_grid[0]*p_grid[1]*p_grid[2]) ; ++i)
				if (bucket_id[i] > 0) {
					Fcalc2[i] = std::norm(Frho_calc[i]);
				} else {
					Fcalc2[i] = 0;
				}
			numeric::fourier::ifft3(Fcalc2, Pcalc);
			ElectronDensity(Pcalc,1.0,numeric::xyzVector< core::Real >(0,0,0), true).writeMRC( "p_calc_unnormalized.mrc");
		}

		// normalize by resolution shell
		// not while minimizing --> assume normalization fixed during minimizaion
		if (! pose.energies().use_nblist()) {
			//utility::vector1<core::Real> F2(bucket_counts.size(), 0);
			F2 = utility::vector1<core::Real>(bucket_counts.size(), 0);
			for (int i=0; i<(int)(p_grid[0]*p_grid[1]*p_grid[2]) ; ++i) {
				if (bucket_id[i] > 0)
					F2[ bucket_id[i] ] += std::norm(Frho_calc[i]);
			}
			for (int i=1; i<=(int)bucket_counts.size(); ++i) {
				F2[i] /= bucket_counts[ i ];
			}
		}

		for (int i=0; i<(int)(p_grid[0]*p_grid[1]*p_grid[2]) ; ++i) {
			if (bucket_id[i] > 0) {
				Fcalc2[i] = std::norm(Frho_calc[i]) / F2[ bucket_id[i] ];
				Frho_calc[i] = Frho_calc[i] / F2[ bucket_id[i] ];
			} else {
				Fcalc2[i] = 0;
				Frho_calc[i] = 0;
			}
		}
	} else {
		for (int i=0; i<(int)(p_grid[0]*p_grid[1]*p_grid[2]) ; ++i) {
			if (bucket_id[i] > 0) {
				Fcalc2[i] = std::norm(Frho_calc[i]);
			} else {
				Fcalc2[i] = 0;
			}
		}
	}

	// back to realspace
	numeric::fourier::ifft3(Fcalc2, Pcalc);

	// sum over symmetric orientations
 	if (symm_ptrs.size() > 1 && !basic::options::option[ basic::options::OptionKeys::patterson::dont_use_symm_in_pcalc ]() ) {
 		ObjexxFCL::FArray3D< double > Pcalcmonomer = Pcalc;
 		for (int i=0; i< (int)(p_grid[0]*p_grid[1]*p_grid[2]) ; ++i) {
 			for (int j=2; j<=(int)symm_ptrs.size(); ++j) {  // j==1 is identity
 				Pcalc[i] += Pcalcmonomer[ symm_ptrs[j][i] ];
 			}
 		}
 	}

	// correlate
	core::Real sumC=0.0, sumC2=0.0, sumO=0.0, sumO2=0.0, sumCO=0.0, vol=0.0;

	sumC=0.0;
	sumC2=0.0;
	for (int z=1; z<=(int)p_grid[2]; ++z) {
		for (int y=1; y<=(int)p_grid[1]; ++y) {
			for (int x=1; x<=(int)p_grid[0]; ++x) {
				core::Real eps_x = PattersonEpsilon(x,y,z);
				if (eps_x < 1e-6) continue;

				core::Real clc_x = Pcalc(x,y,z);

				// find corresponding point in dens_grid
				core::Real obs_x = p_o(x,y,z);

				// SMOOTHED
				sumCO += eps_x*clc_x*obs_x;
				sumO  += eps_x*obs_x;
				sumO2 += eps_x*obs_x*obs_x;
				sumC  += eps_x*clc_x;
				sumC2 += eps_x*clc_x*clc_x;
				vol   += eps_x;
			}
		}
	}
	core::Real varC = (sumC2 - sumC*sumC / vol );
	core::Real varO = (sumO2 - sumO*sumO / vol ) ;

	if (basic::options::option[ basic::options::OptionKeys::edensity::debug ]()) {
			ElectronDensity(Pcalc,1.0,numeric::xyzVector< core::Real >(0,0,0), true).writeMRC( "p_calc.mrc");
			ElectronDensity(PattersonEpsilon,1.0,numeric::xyzVector< core::Real >(0,0,0), true).writeMRC( "p_eps.mrc");
	}

	///////////////////////////////
	// precompute derivatives
	numeric::xyzVector< core::Real > sum_dcc(0,0,0);
	core::Size natoms(0);
	if (cacheCCs) {
		p_sumC = sumC; p_sumC2 = sumC2;
		p_sumO = sumO; p_sumO2 = sumO2;
		p_sumCO = sumCO; p_vol = vol;

		ObjexxFCL::FArray3D< double > dpcc_dx, dpcc_dy, dpcc_dz;

		for (int i=1 ; i<=nres; ++i) {
			conformation::Residue const &rsd_i (pose.residue(i));
			dCCdxs_pat[i].resize( rsd_i.nheavyatoms() );
			std::fill (dCCdxs_pat[i].begin(), dCCdxs_pat[i].end(), numeric::xyzVector< core::Real >(0.0,0.0,0.0));
		}

		core::Real f = (p_sumCO - p_sumC*p_sumO/p_vol);
		core::Real g = (varC * varO);
		core::Real pc_bar = p_sumC/p_vol;

		// for each scatterer
		for (std::map< int,ObjexxFCL::FArray3D< std::complex<double> > >::const_iterator it =  Fdrhoc_dx.begin();
		     it !=  Fdrhoc_dx.end(); ++it) {
			ObjexxFCL::FArray3D< std::complex< double > > &Fdrhoc_dx_i = Fdrhoc_dx[ it->first ];
			ObjexxFCL::FArray3D< std::complex< double > > &Fdrhoc_dy_i = Fdrhoc_dy[ it->first ];
			ObjexxFCL::FArray3D< std::complex< double > > &Fdrhoc_dz_i = Fdrhoc_dz[ it->first ];

			//OLD conv_p_co = real( ifft( fft(epsilon.*p_o)   .* fft(rho_c) .* conj(fft(drhox)) ) );
			//OLD conv_p_c  = real( ifft( fft(epsilon)        .* fft(rho_c) .* conj(fft(drhox)) ) );
			//OLD conv_p_c2 = real( ifft( fft(2*epsilon.*p_c) .*  fft(rho_c) .* conj(fft(drhox)) ) );
			// f = ( sumCO - sumC*sumO/sumV  );
			// g = ( varC*varO );
			// Fconv_p_pcc = fft2( p_mask .* ( (p_o-po_bar)/sqrt(g) - f * varO*(p_c-pc_bar) / (g^1.5) ) ) ...
			//    .* frho_c .* conj(fft2(drhox));
			{
				ObjexxFCL::FArray3D< std::complex<double> > Fdpcc_dx, Fdpcc_dy, Fdpcc_dz;
				dpcc_dx.dimension(p_grid[0],p_grid[1],p_grid[2]);
				for (int i=0; i<(int)(p_grid[2]*p_grid[1]*p_grid[0]); ++i) {
					dpcc_dx[i] = PattersonEpsilon[i] *
					             ( (p_o[i]-po_bar)/sqrt(g) - f * varO * (Pcalc[i]-pc_bar) / std::pow(g,1.5) );
				}

				numeric::fourier::fft3(dpcc_dx, Fdpcc_dx);
				Fdpcc_dy = Fdpcc_dx;
				Fdpcc_dz = Fdpcc_dx;

				for (int i=0; i<(int)(p_grid[2]*p_grid[1]*p_grid[0]); ++i) {
					Fdpcc_dx[i] = 1.0/rho_calc_sum * Fdpcc_dx[i] * Frho_calc[i] * std::conj(Fdrhoc_dx_i[i]);
					Fdpcc_dy[i] = 1.0/rho_calc_sum * Fdpcc_dy[i] * Frho_calc[i] * std::conj(Fdrhoc_dy_i[i]);
					Fdpcc_dz[i] = 1.0/rho_calc_sum * Fdpcc_dz[i] * Frho_calc[i] * std::conj(Fdrhoc_dz_i[i]);
				}

				// inverse FFT
				//   if accuracy is poor may need to do an fft_resize here
				//   to resample on a finer grid
				numeric::fourier::ifft3(Fdpcc_dx, dpcc_dx);
				numeric::fourier::ifft3(Fdpcc_dy, dpcc_dy);
				numeric::fourier::ifft3(Fdpcc_dz, dpcc_dz);

				if (basic::options::option[ basic::options::OptionKeys::edensity::debug ]()) {
					ElectronDensity(dpcc_dx,1.0,numeric::xyzVector< core::Real >(0,0,0), true).writeMRC( "dpcc_dx.mrc");
					ElectronDensity(dpcc_dy,1.0,numeric::xyzVector< core::Real >(0,0,0), true).writeMRC( "dpcc_dy.mrc");
					ElectronDensity(dpcc_dz,1.0,numeric::xyzVector< core::Real >(0,0,0), true).writeMRC( "dpcc_dz.mrc");
				}
			}

			// spline coeffs
			ObjexxFCL::FArray3D< double > dpcc_dx_coeffs, dpcc_dy_coeffs, dpcc_dz_coeffs;
			if ( basic::options::option[ basic::options::OptionKeys::patterson::use_spline_interpolation ]()) {
				spline_coeffs( dpcc_dx, dpcc_dx_coeffs );
				spline_coeffs( dpcc_dy, dpcc_dy_coeffs );
				spline_coeffs( dpcc_dz, dpcc_dz_coeffs );
			}

			// for each atom
			for (int i=1 ; i<=nres; ++i) {
				conformation::Residue const &rsd_i (pose.residue(i));
				if ( rho_calc_as[i].size() == 0 ) continue;

				for (int j=1 ; j<=(int)rsd_i.nheavyatoms(); ++j) {
					if (rho_calc_as[i][j].a() != it->first) continue;

					cartX = rsd_i.atom(j).xyz()-CoM;
					fracX = c2f*cartX;

					atm_idx_ij[0] = fracX[0]*grid[0] - p_origin[0] + 1;
					atm_idx_ij[1] = fracX[1]*grid[1] - p_origin[1] + 1;
					atm_idx_ij[2] = fracX[2]*grid[2] - p_origin[2] + 1;

					// may want to do spline interpolation here
					numeric::xyzVector< core::Real > dpcc;
					if ( basic::options::option[ basic::options::OptionKeys::patterson::use_spline_interpolation ]()) {
						dpcc[0] =  interp_spline( dpcc_dx_coeffs , atm_idx_ij );
						dpcc[1] =  interp_spline( dpcc_dy_coeffs , atm_idx_ij );
						dpcc[2] =  interp_spline( dpcc_dz_coeffs , atm_idx_ij );
					} else {
						dpcc[0] =  interp_linear( dpcc_dx , atm_idx_ij );
						dpcc[1] =  interp_linear( dpcc_dy , atm_idx_ij );
						dpcc[2] =  interp_linear( dpcc_dz , atm_idx_ij );
					}
					dCCdxs_pat[i][j] = dpcc;

					// sum over symmetries
					// i _think_ if we don't care about self interactions we just need to multiply by # of symm
					if (symm_ptrs.size() > 1 && !basic::options::option[ basic::options::OptionKeys::patterson::dont_use_symm_in_pcalc ]()) {
						dCCdxs_pat[i][j] *= symm_ptrs.size();
					}
					//std::cerr << i << "   " << j << "   " << atm_idx_ij << " ==> " << dpcc << std::endl;
					sum_dcc += dCCdxs_pat[i][j];
					natoms++;
				}
			}
		}
		sum_dcc /= natoms;
		TR.Debug << "sum(grad) = " << sum_dcc << std::endl;

		// renormalize so grads sum to 0
		for (int i=1 ; i<=nres; ++i) {
			conformation::Residue const &rsd_i (pose.residue(i));
			if ( rho_calc_as[i].size() == 0 ) continue;
			for (int j=1 ; j<=(int)rsd_i.nheavyatoms(); ++j) {
				dCCdxs_pat[i][j] -= sum_dcc;
			}
		}
	}

	///////////////////////////////
	/// DUMP MAPS FOR DEBUGGING
	if (basic::options::option[ basic::options::OptionKeys::edensity::debug ]()) {
		//TR << "Patterson correl = " << ( (p_sumCO - p_sumC*p_sumO/ p_vol) / sqrt( varC * varO ) ) << std::endl;
		//TR << "          sum_co =  " << p_sumCO << std::endl;
		//TR << "          sum_c  =  " << p_sumC << std::endl;
		//TR << "          sum_o  =  " << p_sumO << std::endl;
		//TR << "          sum_c2 =  " << p_sumC2 << std::endl;
		//TR << "          sum_o2 =  " << p_sumO2 << std::endl;

		//if (useMask) ElectronDensity(Pmask,1.0,numeric::xyzVector< core::Real >(0,0,0), true).writeMRC( "p_mask.mrc" );
		ElectronDensity(rho_calc,1.0, numeric::xyzVector< core::Real >(0,0,0), false).writeMRC( "rho_calc.mrc" );
	}
	return ( (sumCO - sumC*sumO/vol) / sqrt( varC * varO ) );
}



///  Rematch the pose to a patterson map, using previous rho_calc with only rsd changed
///   do not change rho_calc or update cache
///  DOES NOT WORK WITH MASK!
core::Real ElectronDensity::rematchResToPatterson( core::conformation::Residue const &rsd ) const {
	int resid = rsd.seqpos();

	// constants
	numeric::xyzVector< core::Real > cartX, fracX;
	numeric::xyzVector< core::Real > atm_i, atm_j, del_ij, atm_idx_ij;

	// copy old data
	ObjexxFCL::FArray3D< double > rho_calc_copy = rho_calc;

	// subtract old residue's density
	int nheavyatoms_old = rho_calc_atms[resid].size();
	for (int j=1 ; j<=nheavyatoms_old; ++j) {
		OneGaussianScattering const &sig_j = rho_calc_as[resid][j];
		core::Real k = sig_j.k( effectiveB );
		core::Real C = sig_j.C( k );

		// if this atom's weight is 0 continue
		if ( C < 1e-6 ) continue;

		cartX = rho_calc_atms[resid][j];
		fracX = c2f*cartX;
		atm_idx_ij[0] = fracX[0]*grid[0] - p_origin[0] + 1;
		atm_idx_ij[1] = fracX[1]*grid[1] - p_origin[1] + 1;
		atm_idx_ij[2] = fracX[2]*grid[2] - p_origin[2] + 1;

		for (int z=1; z<=(int)p_grid[2]; ++z) {
			del_ij[2] = (atm_idx_ij[2] - z) / grid[2];
			del_ij[0] = del_ij[1] = 0.0;
			if ((f2c*del_ij).length_squared() > square(ATOM_MASK+1)) continue;
			for (int y=1; y<=(int)p_grid[1]; ++y) {
				del_ij[1] = (atm_idx_ij[1] - y) / grid[1] ;
				del_ij[0] = 0.0;
				if ((f2c*del_ij).length_squared() > square(ATOM_MASK+1)) continue;
				for (int x=1; x<=(int)p_grid[0]; ++x) {
					del_ij[0] = (atm_idx_ij[0] - x) / grid[0];

					numeric::xyzVector< core::Real > cart_del_ij = (f2c*del_ij);  // cartesian offset from (x,y,z) to atom_i
					core::Real d2 = (cart_del_ij).length_squared();
					if (d2 <= (ATOM_MASK+1)*(ATOM_MASK+1)) {
						core::Real atm = C*exp(-k*d2);
						rho_calc_copy(x,y,z) -= atm/rho_calc_sum;  // reuse old rho_calc_sum
					}
				}
			}
		}
	}


	// add density to new, updating cache
	int nheavyatoms_new = rsd.nheavyatoms();

	// for (int j=1 ; j<=rsd.natoms(); ++j) {
	//  	cartX = rsd.atom(j).xyz() - p_CoM;
	//  	fracX = rho_calc_atms[resid][j];
	//  	std::cout << j << "  " << cartX[0] << "," << cartX[1] << "," << cartX[2] << "  " << fracX[0] << "," << fracX[1] << "," << fracX[2];
	//  	if (j<= nheavyatoms_new) std::cout <<  " *";
	//  	std::cout << std::endl;
	// }

	for (int j=1 ; j<=nheavyatoms_new; ++j) {
		conformation::Atom const &atm_i( rsd.atom(j) );
		chemical::AtomTypeSet const & atom_type_set( rsd.atom_type_set() );
		std::string elt_i = atom_type_set[ rsd.atom_type_index( j ) ].element();

		OneGaussianScattering sig_j = get_A( elt_i );
		core::Real k = sig_j.k( effectiveB );
		core::Real C = sig_j.C( k );

		// skip randomized residues / residues with no mass (e.g. centroids)
		if ( is_missing_density( atm_i.xyz() ) || C < 1e-6) continue;

		cartX = atm_i.xyz() - p_CoM;
		fracX = c2f*cartX;
		atm_idx_ij[0] = fracX[0]*grid[0] - p_origin[0] + 1;
		atm_idx_ij[1] = fracX[1]*grid[1] - p_origin[1] + 1;
		atm_idx_ij[2] = fracX[2]*grid[2] - p_origin[2] + 1;

		for (int z=1; z<=(int)p_grid[2]; ++z) {
			del_ij[2] = (atm_idx_ij[2] - z) / grid[2];
			del_ij[0] = del_ij[1] = 0.0;
			if ((f2c*del_ij).length_squared() > square(ATOM_MASK+1)) continue;
			for (int y=1; y<=(int)p_grid[1]; ++y) {
				del_ij[1] = (atm_idx_ij[1] - y) / grid[1] ;
				del_ij[0] = 0.0;
				if ((f2c*del_ij).length_squared() > square(ATOM_MASK+1)) continue;
				for (int x=1; x<=(int)p_grid[0]; ++x) {
					del_ij[0] = (atm_idx_ij[0] - x) / grid[0];

					numeric::xyzVector< core::Real > cart_del_ij = (f2c*del_ij);  // cartesian offset from (x,y,z) to atom_i
					core::Real d2 = (cart_del_ij).length_squared();
					if (d2 <= (ATOM_MASK+1)*(ATOM_MASK+1)) {
						core::Real atm = C*exp(-k*d2);
						rho_calc_copy(x,y,z) += atm/rho_calc_sum;  // reuse old rho_calc_sum
					}
				}
			}
		}
	}

	// update all saved stuff derived from this
	ObjexxFCL::FArray3D< std::complex<double> > Fcalc2;
	ObjexxFCL::FArray3D< double > Pcalc_copy;

	numeric::fourier::fft3(rho_calc_copy, Fcalc2);
	// normalize by resolution shell
	utility::vector1<core::Real> F2(bucket_counts.size(), 0);
	for (int i=0; i<(int)(p_grid[0]*p_grid[1]*p_grid[2]) ; ++i) {
		if (bucket_id[i] > 0)
			F2[ bucket_id[i] ] += std::norm(Fcalc2[i]);
	}
	for (int i=1; i<=(int)bucket_counts.size(); ++i) {
		F2[i] /= bucket_counts[ i ];
	}

	Fcalc2.dimension( p_grid[0],p_grid[1],p_grid[2] );
	for (int i=0; i<(int)(p_grid[0]*p_grid[1]*p_grid[2]) ; ++i) {
		if (bucket_id[i] > 0) {
			Fcalc2[i] = std::norm(Fcalc2[i]) / F2[ bucket_id[i] ];
		} else {
			Fcalc2[i] = 0;
		}
	}

	numeric::fourier::ifft3(Fcalc2, Pcalc_copy);

	// sum over symmetric orientations
 	if (symm_ptrs.size() > 1 && !basic::options::option[ basic::options::OptionKeys::patterson::dont_use_symm_in_pcalc ]() ) {
 		ObjexxFCL::FArray3D< double > Pcalcmonomer = Pcalc_copy;
 		for (int i=0; i< (int)(p_grid[0]*p_grid[1]*p_grid[2]) ; ++i) {
 			for (int j=2; j<=(int)symm_ptrs.size(); ++j) {  // j==1 is identity
 				Pcalc_copy[i] += Pcalcmonomer[ symm_ptrs[j][i] ];
 			}
 		}
 	}

	// correlate
	core::Real sumC=0.0, sumC2=0.0, sumO=0.0, sumO2=0.0, sumCO=0.0, vol=0.0;

	sumC=0.0;
	sumC2=0.0;
	for (int z=1; z<=(int)p_grid[2]; ++z) {
		for (int y=1; y<=(int)p_grid[1]; ++y) {
			for (int x=1; x<=(int)p_grid[0]; ++x) {
				core::Real eps_x = PattersonEpsilon(x,y,z);
				if (eps_x < 1e-6) continue;

				core::Real clc_x = Pcalc_copy(x,y,z);

				// find corresponding point in dens_grid
				core::Real obs_x = p_o(x,y,z);

				// SMOOTHED
				sumCO += eps_x*clc_x*obs_x;
				sumO  += eps_x*obs_x;
				sumO2 += eps_x*obs_x*obs_x;
				sumC  += eps_x*clc_x;
				sumC2 += eps_x*clc_x*clc_x;
				vol   += eps_x;
			}
		}
	}
	core::Real varC = (sumC2 - sumC*sumC / vol );
	core::Real varO = (sumO2 - sumO*sumO / vol ) ;
	return ( (sumCO - sumC*sumO/vol) / sqrt( varC * varO ) );
}

////////////////////////////////
////////////////////////////////

void ElectronDensity::updateCachedDensity( core::conformation::Residue const &rsd ) {
	int resid = rsd.seqpos();

	// constants
	numeric::xyzVector< core::Real > cartX, fracX;
	numeric::xyzVector< core::Real > atm_i, atm_j, del_ij, atm_idx_ij;

	// subtract old residue's density
	int nheavyatoms_old = rho_calc_atms[resid].size();
	for (int j=1 ; j<=nheavyatoms_old; ++j) {
		OneGaussianScattering const &sig_j = rho_calc_as[resid][j];
		core::Real k = sig_j.k( effectiveB );
		core::Real C = sig_j.C( k );

		// if this atom's weight is 0 continue
		if ( C < 1e-6 ) continue;

		cartX = rho_calc_atms[resid][j];
		fracX = c2f*cartX;
		atm_idx_ij[0] = fracX[0]*grid[0] - p_origin[0] + 1;
		atm_idx_ij[1] = fracX[1]*grid[1] - p_origin[1] + 1;
		atm_idx_ij[2] = fracX[2]*grid[2] - p_origin[2] + 1;

		for (int z=1; z<=(int)p_grid[2]; ++z) {
			del_ij[2] = (atm_idx_ij[2] - z) / grid[2];
			del_ij[0] = del_ij[1] = 0.0;
			if ((f2c*del_ij).length_squared() > square(ATOM_MASK+1)) continue;
			for (int y=1; y<=(int)p_grid[1]; ++y) {
				del_ij[1] = (atm_idx_ij[1] - y) / grid[1] ;
				del_ij[0] = 0.0;
				if ((f2c*del_ij).length_squared() > square(ATOM_MASK+1)) continue;
				for (int x=1; x<=(int)p_grid[0]; ++x) {
					del_ij[0] = (atm_idx_ij[0] - x) / grid[0];

					numeric::xyzVector< core::Real > cart_del_ij = (f2c*del_ij);  // cartesian offset from (x,y,z) to atom_i
					core::Real d2 = (cart_del_ij).length_squared();
					if (d2 <= (ATOM_MASK+1)*(ATOM_MASK+1)) {
						core::Real atm = C*exp(-k*d2);
						rho_calc(x,y,z) -= atm;
					}
				}
			}
		}
	}

	// add density to new, updating cache
	int nheavyatoms_new = rsd.nheavyatoms();
	rho_calc_atms[resid].resize(nheavyatoms_new);
	rho_calc_as[resid].resize(nheavyatoms_new);
	for (int j=1 ; j<=nheavyatoms_new; ++j) {
		conformation::Atom const &atm_i( rsd.atom(j) );
		chemical::AtomTypeSet const & atom_type_set( rsd.atom_type_set() );
		std::string elt_i = atom_type_set[ rsd.atom_type_index( j ) ].element();

		OneGaussianScattering sig_j = get_A( elt_i );
		core::Real k = sig_j.k( effectiveB );
		core::Real C = sig_j.C( k );

		rho_calc_atms[resid][j] = atm_i.xyz()-p_CoM;
		rho_calc_as[resid][j] = sig_j;

		// skip randomized residues / residues with no mass (e.g. centroids)
		if ( is_missing_density( atm_i.xyz() ) || C < 1e-6) {
			rho_calc_as[resid][j] = OneGaussianScattering(); // weight == 0
			continue;
		}

		cartX = atm_i.xyz();
		fracX = c2f*cartX;
		atm_idx_ij[0] = fracX[0]*grid[0] - p_origin[0] + 1;
		atm_idx_ij[1] = fracX[1]*grid[1] - p_origin[1] + 1;
		atm_idx_ij[2] = fracX[2]*grid[2] - p_origin[2] + 1;

		for (int z=1; z<=(int)p_grid[2]; ++z) {
			del_ij[2] = (atm_idx_ij[2] - z) / grid[2];
			del_ij[0] = del_ij[1] = 0.0;
			if ((f2c*del_ij).length_squared() > square(ATOM_MASK+1)) continue;
			for (int y=1; y<=(int)p_grid[1]; ++y) {
				del_ij[1] = (atm_idx_ij[1] - y) / grid[1] ;
				del_ij[0] = 0.0;
				if ((f2c*del_ij).length_squared() > square(ATOM_MASK+1)) continue;
				for (int x=1; x<=(int)p_grid[0]; ++x) {
					del_ij[0] = (atm_idx_ij[0] - x) / grid[0];

					numeric::xyzVector< core::Real > cart_del_ij = (f2c*del_ij);  // cartesian offset from (x,y,z) to atom_i
					core::Real d2 = (cart_del_ij).length_squared();
					if (d2 <= (ATOM_MASK+1)*(ATOM_MASK+1)) {
						core::Real atm = C*exp(-k*d2);
						rho_calc(x,y,z) += atm;
					}
				}
			}
		}
	}
}




/////////////////////////////////////
/// Computes the symmatric rotation matrices.  Stores mapping in 'symmap'.
void ElectronDensity::compute_symm_rotations(
	core::pose::Pose const &pose,
	core::conformation::symmetry::SymmetryInfoCOP symmInfo /*=NULL*/
) {
	// symmetry
	bool isSymm = (symmInfo.get() != NULL);
	bool remapSymm = remap_symm_;

	if (!isSymm || !remapSymm) return;

	int nsubunits = symmInfo->subunits();
	int nres_per = symmInfo->num_independent_residues();

	// mapping at the (non-virtual) leaves
	// MUST BE DONE EVERY TIME THE POSE CHANGES
	for (int subunit_i=1 ;  subunit_i<=nsubunits; ++subunit_i) {
		// symm
		int startRes = (subunit_i-1)*nres_per+1;
		int source_subunit = subunit_i;
		numeric::xyzMatrix< core::Real > R_i = numeric::xyzMatrix<core::Real>::rows(1,0,0, 0,1,0, 0,0,1);

		if (!symmInfo->bb_is_independent(startRes)) {
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
		for (int i=1 ; i<=nvrts; ++i) {
			int resid = i+vrtStart;
			if (vrts_mapped[i]) {
				continue;
			}
			utility::vector1< core::kinematics::Edge > edges_i = pose.fold_tree().get_outgoing_edges(resid);
			int nchildren = edges_i.size();

			if (nchildren == 0) {
				//TR.Debug << "[ WARNING ] VRT (" << resid << ") at leaf ... ignoring" << std::endl;
				symmap[ resid ] = make_pair( utility::vector1<int>(nsubunits,0), numeric::xyzMatrix<core::Real>::rows(1,0,0, 0,1,0, 0,0,1) );
				vrts_mapped[i] = true;
			} else if ( nchildren == 1) {
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
				for (int j=1; j<=nchildren; ++j) {
					int downstream = edges_i[j].stop();
					if (downstream <= vrtStart) {
						TR.Error << "[ Error ] VRT (" << resid << ") contains multiple jumps to subunits!  Exiting." << std::endl;
						exit(1);
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
					for (int j=2; j<=nchildren; ++j) {
						utility::vector1<int> const &mapping_j = symmap[ edges_i[j].stop() ].first;
						for (int k=1; k<=nsubunits; ++k) {
							if (mapping_j[k] == 0) continue;   // subunit k is not under child j

							// find the subunit to which R_i maps k
							numeric::xyzMatrix< core::Real > R_ik = (symmap[ -k ].second) * R_i;  // R_i * R_k
							core::Real bestFit = 999999;
							int bestL=-1;
							for (int l=1; l<=nsubunits; ++l) {
								numeric::xyzMatrix< core::Real > R_diff = R_ik - (symmap[ -l ].second);
								core::Real thisErr = R_diff.col_x().length_squared() +  R_diff.col_y().length_squared() +  R_diff.col_z().length_squared();
								if (thisErr < bestFit) {
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
	for (int i=1, iend=pose.total_residue(); i<=iend; ++i) {
		if ( pose.residue(i).aa() == core::chemical::aa_vrt ) continue;
		core::Size nAtms = pose.residue(i).nheavyatoms();
		dCCdxs_res[i].resize( nAtms );
		std::fill (dCCdxs_res[i].begin(), dCCdxs_res[i].end(), numeric::xyzVector< core::Real >(0.0,0.0,0.0));
	}
}

/////////////////////////////////////
/// compute rho_c based on a list of atoms
void ElectronDensity::compute_rho(core::pose::Pose const & pose,
								  utility::vector1<core::id::AtomID> const & atom_ids,
								  ObjexxFCL::FArray3D< double > & calculated_density,
								  ObjexxFCL::FArray3D< double > & inv_rho_mask) {
	core::Real SC_scaling = basic::options::option[ basic::options::OptionKeys::edensity::sc_scaling ]();
	core::Real ATOM_MASK_SQ = (ATOM_MASK+ATOM_MASK_PADDING)*(ATOM_MASK+ATOM_MASK_PADDING);
	TR.Warning << "ATOM_MASK: " << ATOM_MASK << std::endl;

	for (int x=0; x<(int)calculated_density.size(); ++x) {
		calculated_density[x] = 0.0;
		inv_rho_mask[x] = 1.0;
	}

	for (core::Size i_id=1; i_id<=atom_ids.size(); ++i_id) {
		core::Size ires = atom_ids[i_id].rsd();
		core::Size iatom = atom_ids[i_id].atomno();

		conformation::Atom const & atm_i( pose.residue(ires).atom(iatom) );
		chemical::AtomTypeSet const & atom_type_set( pose.residue(ires).atom_type_set() );
		std::string elt_i = atom_type_set[ pose.residue(ires).atom_type_index( iatom ) ].element();

		OneGaussianScattering sig_j = get_A( elt_i );
		core::Real k = sig_j.k( effectiveB );
		core::Real C = sig_j.C( k );

		if ( (Size) iatom > pose.residue(ires).last_backbone_atom())
			C *= SC_scaling;
		if ( is_missing_density( atm_i.xyz() ) ) continue;
		if ( C < 1e-6 ) continue;

		numeric::xyzVector< core::Real > cartX, fracX;
		numeric::xyzVector< core::Real > atm_j, del_ij, atm_idx_ij;

		cartX = atm_i.xyz() - getTransform();
		fracX = c2f*cartX;
		atm_idx_ij[0] = pos_mod (fracX[0]*grid[0] - origin[0] + 1 , (double)grid[0]);
		atm_idx_ij[1] = pos_mod (fracX[1]*grid[1] - origin[1] + 1 , (double)grid[1]);
		atm_idx_ij[2] = pos_mod (fracX[2]*grid[2] - origin[2] + 1 , (double)grid[2]);

		using namespace ObjexxFCL::format;
		for (int z=1; z<=calculated_density.u3(); ++z) {
			atm_j[2] = z;
			del_ij[2] = (atm_idx_ij[2] - atm_j[2]) / grid[2];
			if (del_ij[2] > 0.5) del_ij[2]-=1.0;
			if (del_ij[2] < -0.5) del_ij[2]+=1.0;

			del_ij[0] = del_ij[1] = 0.0;
			if ((f2c*del_ij).length_squared() > ATOM_MASK_SQ) continue;  // early exit

			for (int y=1; y<=calculated_density.u2(); ++y) {
				atm_j[1] = y;

				del_ij[1] = (atm_idx_ij[1] - atm_j[1]) / grid[1] ;
				if (del_ij[1] > 0.5) del_ij[1]-=1.0;
				if (del_ij[1] < -0.5) del_ij[1]+=1.0;
				del_ij[0] = 0.0;
				if ((f2c*del_ij).length_squared() > ATOM_MASK_SQ) continue;  // early exit

				for (int x=1; x<=calculated_density.u1(); ++x) {
					atm_j[0] = x;
					del_ij[0] = (atm_idx_ij[0] - atm_j[0]) / grid[0];
					if (del_ij[0] > 0.5) del_ij[0]-=1.0;
					if (del_ij[0] < -0.5) del_ij[0]+=1.0;

					numeric::xyzVector< core::Real > cart_del_ij = (f2c*del_ij);  // cartesian offset from (x,y,z) to atom_i
					core::Real d2 = (cart_del_ij).length_squared();
					if (d2 <= ATOM_MASK_SQ) {
						core::Real atm = C*exp(-k*d2);
						calculated_density(x,y,z) += atm;

						core::Real sigmoid_msk = exp( d2 - (ATOM_MASK)*(ATOM_MASK)  );
						core::Real inv_msk = 1/(1+sigmoid_msk);
						inv_rho_mask(x,y,z) *= (1 - inv_msk);

					}
				}
			}
		}

	}
}

/*
exapansion of eqn 4 in Cowtan, Acta Cryst. (1998) D54. 750-6
Note: conjugation is flipped in this implementation because we use the correlation function between \rho_f(y) and \rho(y+x), instead of \rho_f(y) and \rho(y-x)
 LaTeX equation:
 \overline{\rho _f} &= \frac{1}{\sum_y\epsilon_f(y)} \sum_y\epsilon_f(y)\rho _f(y) \indent (mean \ \rho_f)
 \\ \overline{\rho}(x) &= \frac{1}{\sum_y\epsilon_f(y)} \sum_y\epsilon_f(y)\rho(y-x) \indent (convolute\ mask\ and\ map)
 \\ \sigma^2_{\rho _f} &= \frac{1}{\sum_y\epsilon_f(y)} \sum_y\epsilon_f(y)[\rho_f(y)-\overline{\rho _f}]^2
 \\ &= \begin{bmatrix}\frac{1}{\sum_y\epsilon_f(y)}\sum_y\epsilon_f(y)\rho^2_y(y)\end{bmatrix} - [\overline{\rho _f}]^2 \indent (mean \ \rho^2_f)
 \\ \sigma^2_\rho(x) &= \frac{1}{\sum_y\epsilon_f(y)} \sum_y\epsilon_f(y)[\rho(y-x)-\overline{\rho}(x)]^2
 \\ &= \begin{bmatrix}\frac{1}{\sum_y\epsilon_f(y)} \sum_y\epsilon_f(y) \rho ^2(y-x) \end{bmatrix} - [\overline{\rho}(x)]^2 \indent (convolute\ mask\ and\ map\ squared)
 \\
 \\ t(x)&=\sum_y\epsilon_f(y)\{[\rho_f(y)-\overline{\rho _f}]-\frac{\sigma_{\rho _f}}{\sigma_\rho(x)}[\rho(y-x)-\overline{\rho}(x)]\}^2
 \\ \sigma^2_\rho(x) \cdot t(x) &= \sigma^2_\rho(x) \cdot \sum_y\epsilon_f(y)[\rho_f(y)-\overline{\rho _f}]^2
 \\ & \indent + \sigma^2_{\rho _f} \cdot \sum_y\epsilon_f(y)[\rho(y-x)-\overline{\rho}(x)]^2
 \\ & \indent - 2 \cdot \sigma_{\rho _f} \sigma_\rho(x) \cdot \sum_y \epsilon_f(y)\{[\rho_f(y)-\overline{\rho _f}]\times[\rho(y-x)-\overline{\rho}(x)]\}
 \\ &= 2 \cdot \sigma^2_{\rho _f} \sigma^2_\rho(x) \cdot \sum_y \epsilon_f(y)
 \\ & \indent - 2 \cdot \sigma_{\rho _f} \sigma_\rho(x) \cdot \{
 \\ & \indent \indent \indent \indent \indent \indent \sum_y \epsilon_f(y) \rho_f(y) \rho(y-x)
 \\ & \indent \indent \indent \indent \indent \indent + \overline{\rho _f} \cdot \bar{\rho}(x) \cdot \sum_y \epsilon_f(y)
 \\ & \indent \indent \indent \indent \indent \indent - \overline{\rho _f} \cdot \sum_y \epsilon_f(y)  \rho(y-x)
 \\ & \indent \indent \indent \indent \indent \indent - \bar{\rho}(x) \cdot \sum_y \epsilon_f(y) \rho_f(y)
 \\ & \indent \indent \indent \indent \indent \indent \}
 \\ &= 2 \cdot \sigma^2_{\rho _f} \sigma^2_\rho(x) \cdot \sum_y \epsilon_f(y)
 \\ & \indent - 2 \cdot \sigma_{\rho _f} \sigma_\rho(x) \cdot \{ \sum_y \epsilon_f(y) \rho_f(y) \rho(y-x) - \overline{\rho _f} \cdot \bar{\rho}(x) \cdot \sum_y \epsilon_f(y) \}
 */

numeric::xyzVector< double > ElectronDensity::match_fragment(
							  ObjexxFCL::FArray3D< double > const & rho_calc,
							  ObjexxFCL::FArray3D< double > const & mask,
							  ObjexxFCL::FArray3D< double > const & rho_obs,
							  core::Real radius /* = 6.0 */) {
	using namespace ObjexxFCL::format;

	numeric::xyzVector<core::Real> box_cartX(radius, radius, radius);
	numeric::xyzVector<core::Real> box_fracX(c2f*box_cartX);
	numeric::xyzVector<int> box_grid(0,0,0);
	if (radius > 1e-3) {
		for (Size i=0;i<3;++i) box_grid[i] = (int) ceil(box_fracX[i]*grid[i]);
	}

	core::Real sum_mask = 0.;//·[epsilon_f]
	for (Size x=0; x < mask.size(); ++x) {
		sum_mask += mask[x];  // mask sum
	}
	//TR << "sum_mask " << F(8,3,sum_mask) << std::endl;

	core::Real sum_mask_rhoc = 0.;//·[epsilon_f*rho_f]
	for (Size x=0; x < mask.size(); ++x) {
		sum_mask_rhoc += mask[x]*rho_calc[x];
	}
	//TR << "sum_mask_rho " << F(8,3,sum_mask_rhoc) << std::endl;
	core::Real mean_rhoc = sum_mask_rhoc / sum_mask;//·[epsilon_f*rho_f] / ·(epsilon)
	//TR << "mean_rhoc " << F(8,3,mean_rhoc) << std::endl;

	core::Real sum_mask_rhoc_sq(0.0);//·[epsilon_f*rho_f*rho_f]
	for (Size x=0; x < rho_obs.size(); ++x) {
		sum_mask_rhoc_sq += mask[x]*rho_calc[x]*rho_calc[x];
	}

	ObjexxFCL::FArray3D< double > rhoo_squared(rho_obs.u1(),rho_obs.u2(),rho_obs.u3());
	for (Size x=0; x < rho_obs.size(); ++x) {
		rhoo_squared[x] = rho_obs[x]*rho_obs[x];
	}

	ObjexxFCL::FArray3D< double > rho_calc_mask(rho_obs.u1(),rho_obs.u2(),rho_obs.u3());
	for (Size x=0; x < mask.size(); ++x) {
		rho_calc_mask[x] = mask[x]*rho_calc[x];
	}

	// FFT of rho, rho squared, and mask
	ObjexxFCL::FArray3D< std::complex<double> > Fmask;
	numeric::fourier::fft3(mask, Fmask);

	ObjexxFCL::FArray3D< std::complex<double> > Frhoo;
	numeric::fourier::fft3(rho_obs, Frhoo);

	ObjexxFCL::FArray3D< std::complex<double> > Frhoo_squared;
	numeric::fourier::fft3(rhoo_squared, Frhoo_squared);

	ObjexxFCL::FArray3D< std::complex<double> > Frho_calc_mask;
	numeric::fourier::fft3(rho_calc_mask, Frho_calc_mask);

	// convolution: conjugate map times mask
	ObjexxFCL::FArray3D< std::complex<double> > Fconv_mask__rhoo;
	conj_map_times(Fconv_mask__rhoo, Fmask, Frhoo);
	ObjexxFCL::FArray3D< double > conv_mask__rhoo;
	numeric::fourier::ifft3(Fconv_mask__rhoo, conv_mask__rhoo);
	ObjexxFCL::FArray3D< double > mean_conv_rhoo(rho_obs.u1(),rho_obs.u2(),rho_obs.u3()); // \bar{\rho}(x) in the paper
	for (Size x=0; x < mask.size(); ++x) {
		mean_conv_rhoo[x] = conv_mask__rhoo[x] / sum_mask;
		//TR << "mean_conv_rhoo " << I(10,x) << F(8,3,mean_conv_rhoo[x]) << std::endl;
	}

	// convolution: conjugate map squared times mask
	ObjexxFCL::FArray3D< std::complex<double> > Fconv_mask__rhoo_sq;
	conj_map_times(Fconv_mask__rhoo_sq, Fmask, Frhoo_squared);
	ObjexxFCL::FArray3D< double > conv_mask__rhoo_sq;
	numeric::fourier::ifft3(Fconv_mask__rhoo_sq, conv_mask__rhoo_sq);

	// convolute fragment map with observed map
	ObjexxFCL::FArray3D< std::complex<double> > Fconv_mask_rhoc__rhoo;
	conj_map_times(Fconv_mask_rhoc__rhoo, Frho_calc_mask, Frhoo);
	ObjexxFCL::FArray3D< double > conv_mask_rhoc__rhoo;
	numeric::fourier::ifft3(Fconv_mask_rhoc__rhoo, conv_mask_rhoc__rhoo);

	// standard deviations
	core::Real sigma_rho_frag_sq(sum_mask_rhoc_sq/sum_mask - mean_rhoc*mean_rhoc);
	//TR << "sigma_rho_frag_sq " << F(8,3,sigma_rho_frag_sq) << std::endl;
	ObjexxFCL::FArray3D< double > sigma_rhoo_sq(rho_obs.u1(),rho_obs.u2(),rho_obs.u3());
	for (Size x=0; x < mask.size(); ++x) {
		sigma_rhoo_sq[x] = conv_mask__rhoo_sq[x]/sum_mask - mean_conv_rhoo[x]* mean_conv_rhoo[x];
		//TR << "sigma_rhoo_sq " << I(10,x) << F(8,3,sigma_rhoo_sq[x]) << std::endl;
	}

	numeric::xyzVector< int > grid_idx (0,0,0);
	numeric::xyzVector< int > max_score_idx (0,0,0);
	double max_score(0.);
	double min_score(0.);
	ObjexxFCL::FArray3D< double > match_score_density(rho_obs.u1(),rho_obs.u2(),rho_obs.u3());

	for (grid_idx[2] = 1; grid_idx[2] <= density.u3(); grid_idx[2]++) {
		if (box_grid[2] != 0) { if (grid_idx[2] - 1 > box_grid[2] && (density.u3()-grid_idx[2]) > box_grid[2]) continue; }
		for (grid_idx[1] = 1; grid_idx[1] <= density.u2(); grid_idx[1]++) {
			if (box_grid[1] != 0) { if (grid_idx[1] - 1 > box_grid[1] && (density.u2()-grid_idx[1]) > box_grid[1]) continue; }
			for (grid_idx[0] = 1; grid_idx[0] <= density.u1(); grid_idx[0]++) {
				if (box_grid[0] != 0) { if (grid_idx[0] - 1 > box_grid[0] && (density.u1()-grid_idx[0]) > box_grid[0]) continue; }

				core::Real sigma_sq = sigma_rho_frag_sq * sigma_rhoo_sq(grid_idx[0],grid_idx[1],grid_idx[2]);
				if ( sigma_sq < 1e-30) continue;
                core::Real sigma = sqrt(sigma_sq);
                core::Real term1 =  sigma * mean_rhoc * mean_conv_rhoo(grid_idx[0],grid_idx[1],grid_idx[2]);
                core::Real term2 = -sigma * conv_mask_rhoc__rhoo(grid_idx[0],grid_idx[1],grid_idx[2]) / sum_mask;

				core::Real match_score = -2. * (sigma_sq + term1 + term2) / sigma_rhoo_sq(grid_idx[0],grid_idx[1],grid_idx[2]) ;
				match_score_density(grid_idx[0],grid_idx[1],grid_idx[2]) = match_score;

				//TR << "Match " << I(4,grid_idx[0]) << I(4,grid_idx[1]) << I(4,grid_idx[2]) << F(8,3,sigma_sq) << F(8,3,sigma) << F(8,3, mean_conv_rhoo(grid_idx[0],grid_idx[1],grid_idx[2])) << F(8,3, conv_mask_rhoc__rhoo(grid_idx[0],grid_idx[1],grid_idx[2])) << F(8,3,term1) << F(8,3,term2) << F(8,3,sigma_sq + term1 + term2) << F(10,3,match_score) << std::endl;

				if (max_score_idx == numeric::xyzVector< int > (0,0,0)
					|| max_score < match_score) {
					max_score_idx = grid_idx;
					max_score = match_score;
				}
				if (max_score_idx == numeric::xyzVector< int > (0,0,0)
					|| min_score < match_score) {
					min_score = match_score;
				}

			}
		}
	}

	/*
	for (grid_idx[2] = 1; grid_idx[2] <= density.u3(); grid_idx[2]++) {
		for (grid_idx[1] = 1; grid_idx[1] <= density.u2(); grid_idx[1]++) {
			for (grid_idx[0] = 1; grid_idx[0] <= density.u1(); grid_idx[0]++) {
				core::Real sigma_sq = sigma_rho_frag_sq * sigma_rhoo_sq(grid_idx[0],grid_idx[1],grid_idx[2]);
				if ( sigma_sq < 1e-30) {
					match_score_density(grid_idx[0],grid_idx[1],grid_idx[2]) = min_score;
				}
			}
		}
	}
	core::scoring::electron_density::ElectronDensity temp_edensity_map(*this);
	temp_edensity_map.set_data(match_score_density);
	temp_edensity_map.writeMRC("match_score.mrc");
	 */

	numeric::xyzVector<core::Real> fracX(
										 ( (core::Real)max_score_idx[0] - 1 ) / grid[0],
										 ( (core::Real)max_score_idx[1] - 1 ) / grid[1],
										 ( (core::Real)max_score_idx[2] - 1 ) / grid[2] );

	for (core::Size i=0; i<3; ++i) {
		if ( fracX[i] > 0.5 ) {
			fracX[i] -= 1.0;
		}
	}
	numeric::xyzVector_double cartX = f2c*fracX;

	//TR << "Match idx " << I(8,max_score_idx[0]) << I(8,max_score_idx[1]) << I(8,max_score_idx[2]) << F(10,3,max_score) << std::endl;
	//TR << "Match frac" << F(8,3,fracX[0]) << F(8,3,fracX[1]) << F(8,3,fracX[2]) << F(10,3,max_score) << std::endl;
	//TR << "Match cart" << F(8,3,cartX[0]) << F(8,3,cartX[1]) << F(8,3,cartX[2]) << F(10,3,max_score) << std::endl;

	return cartX;
}

/////////////////////////////////////
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
	if (!isLoaded) {
		TR << "[ ERROR ]  ElectronDensity::matchRes called but no map is loaded!\n";
		return 0.0;
	}

	if ( scoring_mask_.find(resid) != scoring_mask_.end() ) return 0.0;

	// symmetry
	bool isSymm = (symmInfo.get() != NULL);
	bool remapSymm = remap_symm_;

	// symm
	if (isSymm && !symmInfo->bb_is_independent(resid) && !remapSymm) return 0.0; // only score monomer

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
	for (int i=HALFWINDOW-1; i>=0; --i) {
		if ( (resid-i)>1 && pose.fold_tree().is_cutpoint( resid-i-1 ) ) win_start = resid-i;
		if ( (resid+i)<=nres && pose.fold_tree().is_cutpoint( resid+i ) ) win_stop = resid+i;
	}
	win_start = std::max(1, (int)win_start);
	win_stop = std::min((int)win_stop, nres);

	// 1: grab context atom ids
	//    only extend left/right to the next chainbreak
	for (Size i=win_start; i<=win_stop; ++i) {
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
				if (j<win_start || j>win_stop) {
					neighborResids.insert( j );
				}
			}
		}

		if (i==(Size)resid) continue;  // already added these atoms
		if (!rsd_i.is_polymer()) continue;

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
	for (std::set< core::Size >::iterator it=neighborResids.begin(); it!=neighborResids.end(); it++) {
		core::conformation::Residue const &rsd_i( pose.residue(*it) );
		if (!rsd_i.is_polymer()) continue;
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

		if ( is_missing_density( atom.xyz() ) ) continue;  // check if coords were randomized

		cartX = atom.xyz() - getTransform();
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
		idxX = numeric::xyzVector< core::Real >( fracX[0]*grid[0] - origin[0] + 1  ,
		                                         fracX[1]*grid[1] - origin[1] + 1  ,
		                                         fracX[2]*grid[2] - origin[2] + 1  );
		idxX_high[0] = std::max(idxX[0],idxX_high[0]);
		idxX_high[1] = std::max(idxX[1],idxX_high[1]);
		idxX_high[2] = std::max(idxX[2],idxX_high[2]);
	}
	
	int firstMaskedAtom = 1; 
	int lastMaskedAtom = atmList.size();
	if (basic::options::option[ basic::options::OptionKeys::edensity::unmask_bb ]())
		firstMaskedAtom = rsd.last_backbone_atom()+1;

	// add in neighbor atoms
	// don't let them modify size of bounding box
	for ( int j=1; j<=nNeighborAtoms; ++j ) {
		conformation::Atom const &atom (neighborAtoms[j]);
		numeric::xyzVector< core::Real > cartX, idxX, fracX;
		if ( is_missing_density( atom.xyz() ) ) continue;  // check if coords were randomized
		cartX = atom.xyz() - getTransform();
		fracX = c2f*cartX;
		idxX = numeric::xyzVector< core::Real >( fracX[0]*grid[0] - origin[0] + 1 ,
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

	for (int x=0; x<bbox_dims[0]*bbox_dims[1]*bbox_dims[2]; ++x) {
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

	for (int i=1; i<=(int)atmList.size(); ++i) {
		numeric::xyzVector< core::Real > & atm_i = atmList[i];
		atm_i[0] -= bbox_min[0];
		atm_i[1] -= bbox_min[1];
		atm_i[2] -= bbox_min[2];

		numeric::xyzVector< core::Real > atm_j, del_ij;

		core::Real k,C;
		if ( i <= nResAtms) {
			std::string elt_i = atom_type_set[ rsd.atom_type_index( i ) ].element();
			OneGaussianScattering sig_j = get_A( elt_i );
			k = sig_j.k( effectiveB );
			C = sig_j.C( k );
		} else if (i<= lastMaskedAtom) {
			OneGaussianScattering const &sig_j = contextAtomAs[i-nResAtms];
			k = sig_j.k( effectiveB );
			C = sig_j.C( k );
		} else {
			OneGaussianScattering const &sig_j = neighborAtomAs[i-lastMaskedAtom];
			k = sig_j.k( effectiveB );
			C = sig_j.C( k );
		}

		for (int z=1; z<=bbox_dims[2]; ++z) {
			atm_j[2] = z;
			del_ij[2] = (atm_i[2] - atm_j[2]) / grid[2];
			// wrap-around??
			if (del_ij[2] > 0.5) del_ij[2]-=1.0;
			if (del_ij[2] < -0.5) del_ij[2]+=1.0;

			del_ij[0] = del_ij[1] = 0.0;
			if ((f2c*del_ij).length_squared() > (ATOM_MASK+ATOM_MASK_PADDING)*(ATOM_MASK+ATOM_MASK_PADDING)) continue;

			mapZ = (z+bbox_min[2]) % grid[2];
			if (mapZ <= 0) mapZ += grid[2];
			if (mapZ > density.u3()) continue;

			for (int y=1; y<=bbox_dims[1]; ++y) {
				atm_j[1] = y;

				// early exit?
				del_ij[1] = (atm_i[1] - atm_j[1]) / grid[1] ;
				// wrap-around??
				if (del_ij[1] > 0.5) del_ij[1]-=1.0;
				if (del_ij[1] < -0.5) del_ij[1]+=1.0;
				del_ij[0] = 0.0;
				if ((f2c*del_ij).length_squared() > (ATOM_MASK+ATOM_MASK_PADDING)*(ATOM_MASK+ATOM_MASK_PADDING)) continue;

				mapY = (y+bbox_min[1]) % grid[1];
				if (mapY <= 0) mapY += grid[1];
				if (mapY > density.u2()) continue;

				for (int x=1; x<=bbox_dims[0]; ++x) {
					atm_j[0] = x;

					// early exit?
					del_ij[0] = (atm_i[0] - atm_j[0]) / grid[0];
					// wrap-around??
					if (del_ij[0] > 0.5) del_ij[0]-=1.0;
					if (del_ij[0] < -0.5) del_ij[0]+=1.0;

					mapX = (x+bbox_min[0]) % grid[0];
					if (mapX <= 0) mapX += grid[0];
					if (mapX > density.u1()) continue;

					numeric::xyzVector< core::Real > cart_del_ij = (f2c*del_ij);  // cartesian offset from atom_i to (x,y,z)
					core::Real d2 = (cart_del_ij).length_squared();
					if (d2 <= (ATOM_MASK+ATOM_MASK_PADDING)*(ATOM_MASK+ATOM_MASK_PADDING)) {
						core::Real atm = C*exp(-k*d2);
						core::Real sigmoid_msk = exp( d2 - (ATOM_MASK)*(ATOM_MASK)  );
						core::Real inv_msk = 1/(1+sigmoid_msk);

						rho_obs(x,y,z) = density(mapX,mapY,mapZ);
						if (i>=firstMaskedAtom && i<=lastMaskedAtom) {
							rho_calc_fg(x,y,z) += atm;
							inv_rho_mask(x,y,z) *= (1 - inv_msk);
							if (cacheCCs) {
								int idx = (z-1)*rho_calc_fg.u2()*rho_calc_fg.u1() + (y-1)*rho_calc_fg.u1() + x-1;

								core::Real eps_i = (1-inv_msk), inv_eps_i;
								if (eps_i == 0) // divide-by-zero
									inv_eps_i = sigmoid_msk;
								else
									inv_eps_i = 1/eps_i;

								rho_dx_pt[i].push_back  ( idx );
								rho_dx_atm[i].push_back ( (-2*k*atm)*cart_del_ij );
								rho_dx_mask[i].push_back( (-2*sigmoid_msk*inv_msk*inv_msk*inv_eps_i)*cart_del_ij );
							}
						} else {
							rho_calc_bg(x,y,z) += atm;
						}
					}
				}
			}
		}
	}

	
	if (basic::options::option[ basic::options::OptionKeys::edensity::debug ]() && resid == 1) {
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

	for (int x=0; x<bbox_dims[0]*bbox_dims[1]*bbox_dims[2]; ++x) {
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
	if (varC_i == 0 || varO_i == 0)
		CC_i = 0;
	else
		CC_i = (sumCO_i - sumC_i*sumO_i/ vol_i) / sqrt( varC_i * varO_i );

	if (cacheCCs) {
		CCs[resid] = CC_i;
	}

	///////////////////////////
	/// 4  CALCULATE PER-ATOM DERIVATIVES
	if (cacheCCs) {
		for (int j=1 ; j<=(int)lastMaskedAtom; ++j) {
			utility::vector1< int > const &rho_dx_pt_ij   = rho_dx_pt[j];
			utility::vector1< numeric::xyzVector<core::Real> > const &rho_dx_mask_ij = rho_dx_mask[j];
			utility::vector1< numeric::xyzVector<core::Real> > const &rho_dx_atm_ij  = rho_dx_atm[j];

			numeric::xyzVector< core::Real > dVdx_ij(0,0,0), dOdx_ij(0,0,0), dO2dx_ij(0,0,0), dCOdx_ij(0,0,0), dC2dx_ij(0,0,0), dCdx_ij(0,0,0);

			int npoints = rho_dx_pt_ij.size();
			for (int n=1; n<=npoints; ++n) {
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
	bool ignoreBs /*=false*/
) {
	// make sure map is loaded
	if (!isLoaded) {
		TR << "[ ERROR ]  ElectronDensity::matchResFast called but no map is loaded!\n";
		return 0.0;
	}

	if ( fastdens_score.u1()*fastdens_score.u2()*fastdens_score.u3()*fastdens_score.u4() == 0 )
		setup_fastscoring_first_time(pose);

	if ( scoring_mask_.find(resid) != scoring_mask_.end() ) return 0.0;

	// symmetry
	bool isSymm = (symmInfo.get() != NULL);
	bool remapSymm = remap_symm_;

	// symm
	if (isSymm && !symmInfo->bb_is_independent(resid) && !remapSymm) return 0.0; // only score monomer

	core::Real score = 0;
	numeric::xyzVector< core::Real > fracX, idxX;
	for (Size i=1; i<=rsd.nheavyatoms(); ++i) {
		chemical::AtomTypeSet const & atom_type_set( rsd.atom_type_set() );
		std::string elt_i = atom_type_set[ rsd.atom_type_index( i ) ].element();
		OneGaussianScattering sig_j = get_A( elt_i );
		core::Real B = pose.pdb_info() ? pose.pdb_info()->temperature( rsd.seqpos(), i ) : effectiveB;
		core::Real k = sig_j.k( B );

		if (ignoreBs) k = 4*M_PI*M_PI/effectiveB;

		core::Real kbin = 1;
		if (nkbins_>1)
			kbin = (k - kmin_)/kstep_ + 1;
		if (kbin<1) kbin=1;
		if (kbin>nkbins_) kbin=nkbins_;

		fracX = c2f*rsd.atom(i).xyz();

		idxX[0] = pos_mod (fracX[0]*fastgrid[0] - fastorigin[0] + 1 , (double)fastgrid[0]);
		idxX[1] = pos_mod (fracX[1]*fastgrid[1] - fastorigin[1] + 1 , (double)fastgrid[1]);
		idxX[2] = pos_mod (fracX[2]*fastgrid[2] - fastorigin[2] + 1 , (double)fastgrid[2]);
		idxX = numeric::xyzVector<core::Real>( fracX[0]*fastgrid[0] - fastorigin[0] + 1,
		                                       fracX[1]*fastgrid[1] - fastorigin[1] + 1,
		                                       fracX[2]*fastgrid[2] - fastorigin[2] + 1);
		core::Real score_i = interp_spline( fastdens_score , kbin, idxX );
		core::Real W = sig_j.a(  ) / 6.0;
		score += W*score_i;
	}

	return score;
}


//  Compute the gradient (fast density score)
void ElectronDensity::dCCdx_fastRes(
				int atmid, int resid,
				numeric::xyzVector<core::Real> const &X,
        core::conformation::Residue const &rsd,
        core::pose::Pose const &pose,
				numeric::xyzVector<core::Real> &dCCdX ){

	// make sure map is loaded
	if (!isLoaded) {
		TR << "[ ERROR ]  ElectronDensity::dCCdx_fastRes called but no map is loaded!\n";
		return;
	}

	if ( fastdens_score.u1()*fastdens_score.u2()*fastdens_score.u3()*fastdens_score.u4() == 0 )
		setup_fastscoring_first_time(pose);

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
	if (nkbins_>1)
		kbin = (k - kmin_)/kstep_ + 1;
	if (kbin<1) kbin=1;
	if (kbin>nkbins_) kbin=nkbins_;

	numeric::xyzVector<core::Real> dCCdX_grid;
	core::Real dkbin;

	core::Real W = sig_j.a(  ) / 6.0;
	interp_dspline( fastdens_score , idxX, kbin,  dCCdX_grid, dkbin );
	dCCdX[0] = W*dCCdX_grid[0]*c2f(1,1)*fastgrid[0] + W*dCCdX_grid[1]*c2f(2,1)*fastgrid[1] + W*dCCdX_grid[2]*c2f(3,1)*fastgrid[2];
	dCCdX[1] = W*dCCdX_grid[0]*c2f(1,2)*fastgrid[0] + W*dCCdX_grid[1]*c2f(2,2)*fastgrid[1] + W*dCCdX_grid[2]*c2f(3,2)*fastgrid[2];
	dCCdX[2] = W*dCCdX_grid[0]*c2f(1,3)*fastgrid[0] + W*dCCdX_grid[1]*c2f(2,3)*fastgrid[1] + W*dCCdX_grid[2]*c2f(3,3)*fastgrid[2];

	if (ExactDerivatives) {
		numeric::xyzVector<core::Real> dCCdX1 = dCCdX;

		core::conformation::Residue rsd_copy = rsd;

		rsd_copy.atom( atmid ).xyz( numeric::xyzVector<core::Real>( X[0]+NUM_DERIV_H_CEN,X[1],X[2] ) );
		core::Real CC_px = getDensityMap().matchResFast( resid, rsd_copy, pose, NULL );
		rsd_copy.atom( atmid ).xyz( numeric::xyzVector<core::Real>( X[0]-NUM_DERIV_H_CEN,X[1],X[2] ) );
		core::Real CC_mx = getDensityMap().matchResFast( resid, rsd_copy, pose, NULL );
		rsd_copy.atom( atmid ).xyz( numeric::xyzVector<core::Real>( X[0],X[1]+NUM_DERIV_H_CEN,X[2] ) );
		core::Real CC_py = getDensityMap().matchResFast( resid, rsd_copy, pose, NULL );
		rsd_copy.atom( atmid ).xyz( numeric::xyzVector<core::Real>( X[0],X[1]-NUM_DERIV_H_CEN,X[2] ) );
		core::Real CC_my = getDensityMap().matchResFast( resid, rsd_copy, pose, NULL );
		rsd_copy.atom( atmid ).xyz( numeric::xyzVector<core::Real>( X[0],X[1],X[2]+NUM_DERIV_H_CEN ) );
		core::Real CC_pz = getDensityMap().matchResFast( resid, rsd_copy, pose, NULL );
		rsd_copy.atom( atmid ).xyz( numeric::xyzVector<core::Real>( X[0],X[1],X[2]-NUM_DERIV_H_CEN ) );
		core::Real CC_mz = getDensityMap().matchResFast( resid, rsd_copy, pose, NULL );

		// rescore with orig pose
		getDensityMap().matchRes( resid, rsd, pose, NULL, false );

		dCCdX[0] = (CC_px-CC_mx)/(2*NUM_DERIV_H_CEN); // * dCCdxs_res_multiplier[resid][atmid];
		dCCdX[1] = (CC_py-CC_my)/(2*NUM_DERIV_H_CEN); // * dCCdxs_res_multiplier[resid][atmid];
		dCCdX[2] = (CC_pz-CC_mz)/(2*NUM_DERIV_H_CEN); // * dCCdxs_res_multiplier[resid][atmid];

		TR << "   " <<  dCCdX<< "  ;  " <<  dCCdX1 << std::endl;
	}
}

//  Compute the gradient (fast density score) w.r.t B factors
Real
ElectronDensity::dCCdB_fastRes(
	int atmid, int resid,
	core::conformation::Residue const &rsd,
	core::pose::Pose const &pose
)
{
	// make sure map is loaded
	if (!isLoaded) {
		TR << "[ ERROR ]  ElectronDensity::dCCdB_fast called but no map is loaded!\n";
		return 0;
	}

	if ( fastdens_score.u1()*fastdens_score.u2()*fastdens_score.u3()*fastdens_score.u4() == 0 )
		setup_fastscoring_first_time(pose);

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
	if (nkbins_ == 1)
		return 0;
	core::Real kbin = (k - kmin_)/kstep_ + 1;
	if (kbin<1 || kbin>nkbins_)
		return 0;

	numeric::xyzVector<core::Real> dCCdX_grid;
	core::Real dscore_dkbin=0.0, dscore_dk;
	interp_dspline( fastdens_score , idxX, kbin,  dCCdX_grid, dscore_dkbin );
	dscore_dk = dscore_dkbin / kstep_;

	core::Real W = sig_j.a(  ) / 6.0;
	Real dCCdb = W*dscore_dk*sig_j.dk( B );

	return dCCdb;
}

//  Compute the gradient (fast density score) w.r.t B factors
void
ElectronDensity::dCCdBs(
	core::pose::Pose const &pose,
	utility::vector1< core::Real>  & dE_dvars
) {
	// 1 get rhoc
	core::scoring::electron_density::poseCoords litePose;

	core::conformation::symmetry::SymmetryInfoCOP symm_info;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		core::conformation::symmetry::SymmetricConformation const & SymmConf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation const &> ( pose.conformation()) );
		symm_info = SymmConf.Symmetry_Info();
	}

	for (int i=1; i<=pose.total_residue(); ++i) {
		if (symm_info && !symm_info->bb_is_independent( i ) ) continue;
		core::conformation::Residue const & rsd_i ( pose.residue(i) );
		if ( rsd_i.aa() == core::chemical::aa_vrt ) continue;

		core::Size natoms = rsd_i.nheavyatoms();
		for (int j=1; j<=natoms; ++j) {
			core::conformation::Atom const &atom_j( rsd_i.atom(j) );
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

	utility::vector1< utility::vector1< int > > rho_dx_pt(natoms);
	utility::vector1< utility::vector1< core::Real > > rho_dx_bs(natoms);

	for (int i=1 ; i<=(int)natoms; ++i) {
		std::string elt_i = litePose[i].elt_;
		OneGaussianScattering sig_j = get_A( elt_i );
		core::Real B_i = litePose[i].B_;
		core::Real k = sig_j.k( B_i );
		core::Real dk = sig_j.dk( B_i );
		k = std::min ( k, 4*M_PI*M_PI/minimumB );
		core::Real C = sig_j.C( k ); if ( C < 1e-6 ) continue;

		numeric::xyzVector< core::Real> cartX = litePose[i].x_ - getTransform();
		numeric::xyzVector< core::Real> fracX = c2f*cartX;
		numeric::xyzVector< core::Real> atm_j, del_ij, atm_idx;
		atm_idx[0] = pos_mod (fracX[0]*grid[0] - origin[0] + 1 , (double)grid[0]);
		atm_idx[1] = pos_mod (fracX[1]*grid[1] - origin[1] + 1 , (double)grid[1]);
		atm_idx[2] = pos_mod (fracX[2]*grid[2] - origin[2] + 1 , (double)grid[2]);
		for (int z=1; z<=density.u3(); ++z) {
			atm_j[2] = z;
			del_ij[2] = (atm_idx[2] - atm_j[2]) / grid[2];
			if (del_ij[2] > 0.5) del_ij[2]-=1.0; if (del_ij[2] < -0.5) del_ij[2]+=1.0;
			del_ij[0] = del_ij[1] = 0.0;
			if ((f2c*del_ij).length_squared() > (ATOM_MASK+ATOM_MASK_PADDING)*(ATOM_MASK+ATOM_MASK_PADDING)) continue;
			for (int y=1; y<=density.u2(); ++y) {
				atm_j[1] = y;
				del_ij[1] = (atm_idx[1] - atm_j[1]) / grid[1] ;
				if (del_ij[1] > 0.5) del_ij[1]-=1.0; if (del_ij[1] < -0.5) del_ij[1]+=1.0;
				del_ij[0] = 0.0;
				if ((f2c*del_ij).length_squared() > (ATOM_MASK+ATOM_MASK_PADDING)*(ATOM_MASK+ATOM_MASK_PADDING)) continue;
				for (int x=1; x<=density.u1(); ++x) {
					atm_j[0] = x;
					del_ij[0] = (atm_idx[0] - atm_j[0]) / grid[0];
					if (del_ij[0] > 0.5) del_ij[0]-=1.0; if (del_ij[0] < -0.5) del_ij[0]+=1.0;
					numeric::xyzVector< core::Real > cart_del_ij = (f2c*del_ij);  // cartesian offset from (x,y,z) to atom_i
					core::Real d2 = (cart_del_ij).length_squared();

					if (d2 <= (ATOM_MASK+ATOM_MASK_PADDING)*(ATOM_MASK+ATOM_MASK_PADDING)) {
						core::Real atm = C*exp(-k*d2);
						core::Real sigmoid_msk = exp( d2 - (ATOM_MASK)*(ATOM_MASK)  );
						core::Real inv_msk = 1/(1+sigmoid_msk);
						inv_rho_mask(x,y,z) *= (1 - inv_msk);
						rho_calc(x,y,z) += atm;

						int idx = (z-1)*density.u2()*density.u1() + (y-1)*density.u1() + x-1;
						rho_dx_pt[i].push_back ( idx );
						rho_dx_bs[i].push_back ( -atm*(2*d2*k - 3)*dk / k );
					}
				}
			}
		}
	}

	// CCs
	core::Real sumC_i=0, sumO_i=0, sumCO_i=0, vol_i=0, CC_i=0, sumO2_i=0.0, sumC2_i=0.0, varC_i=0, varO_i=0;
	for (int x=0; x<density.u1()*density.u2()*density.u3(); ++x) {
		core::Real clc_x = rho_calc[x];
		core::Real obs_x = density[x];
		core::Real eps_x = 1-inv_rho_mask[x];

		sumCO_i += eps_x*clc_x*obs_x;
		sumO_i  += eps_x*obs_x;
		sumO2_i += eps_x*obs_x*obs_x;
		sumC_i  += eps_x*clc_x;
		sumC2_i += eps_x*clc_x*clc_x;
		vol_i   += eps_x;
	}
	varC_i = (sumC2_i - sumC_i*sumC_i / vol_i );
	varO_i = (sumO2_i - sumO_i*sumO_i / vol_i ) ;
	if (varC_i > 0 && varO_i > 0)
		CC_i = (sumCO_i - sumC_i*sumO_i/ vol_i) / sqrt( varC_i * varO_i );

	//std::cerr << CC_i << std::endl;

	// dCCdbs
	for (int i=1 ; i<=natoms; ++i) {
		utility::vector1< int > const &rho_dx_pt_i = rho_dx_pt[i];
		utility::vector1< core::Real > const &rho_dx_bs_i = rho_dx_bs[i];
		int npoints = rho_dx_pt_i.size();

		core::Real dCOdx_i=0, dC2dx_i=0;
		for (int n=1; n<=npoints; ++n) {
			const int x(rho_dx_pt_i[n]);
			core::Real clc_x = rho_calc[x];
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
	numeric::xyzVector< core::Real > const & X,
	core::conformation::Residue const &rsd,
	core::pose::Pose const & pose,
	numeric::xyzVector<core::Real> &dCCdX
)
{
	// make sure map is loaded
	if (!isLoaded) {
		TR << "[ ERROR ]  ElectronDensity::dCCdx_res called but no map is loaded!\n";
		dCCdX = numeric::xyzVector<core::Real>(0.0,0.0,0.0);
		exit(1);
	}

	static bool warned=false;
	if (!DensScoreInMinimizer) {
		if (!warned)
			TR << "[ WARNING ] dCCdx_res called but DensityScoreInMinimizer = false ... returning 0" << std::endl;
		//warned = true;
		dCCdX = numeric::xyzVector<core::Real>(0.0,0.0,0.0);
		return;
	}

	// if we didn't score this residue
	//      (because it was missing density or masked or a symm copy)
	// then don't compute its derivative
	if ( (int)dCCdxs_res[resid].size() < atmid ) {
		//std::cerr << "no " << resid << "." << atmid << std::endl;
 		dCCdX = numeric::xyzVector<core::Real>(0.0,0.0,0.0);
		return;
	}

	if (! ExactDerivatives ) {
		dCCdX = dCCdxs_res[resid][atmid];
	} else {
		////////////////////////////////
		//////  debug : compare d_true to d_exact
		//
		core::conformation::Residue rsd_copy = rsd;

		rsd_copy.atom( atmid ).xyz( numeric::xyzVector<core::Real>( X[0]+NUM_DERIV_H_CEN,X[1],X[2] ) );
		core::Real CC_px = getDensityMap().matchRes( resid, rsd_copy, pose, NULL, false );

		rsd_copy.atom( atmid ).xyz( numeric::xyzVector<core::Real>( X[0]-NUM_DERIV_H_CEN,X[1],X[2] ) );
		core::Real CC_mx = getDensityMap().matchRes( resid, rsd_copy, pose, NULL, false );

		rsd_copy.atom( atmid ).xyz( numeric::xyzVector<core::Real>( X[0],X[1]+NUM_DERIV_H_CEN,X[2] ) );
		core::Real CC_py = getDensityMap().matchRes( resid, rsd_copy, pose, NULL, false );

		rsd_copy.atom( atmid ).xyz( numeric::xyzVector<core::Real>( X[0],X[1]-NUM_DERIV_H_CEN,X[2] ) );
		core::Real CC_my = getDensityMap().matchRes( resid, rsd_copy, pose, NULL, false );

		rsd_copy.atom( atmid ).xyz( numeric::xyzVector<core::Real>( X[0],X[1],X[2]+NUM_DERIV_H_CEN ) );
		core::Real CC_pz = getDensityMap().matchRes( resid, rsd_copy, pose, NULL, false );

		rsd_copy.atom( atmid ).xyz( numeric::xyzVector<core::Real>( X[0],X[1],X[2]-NUM_DERIV_H_CEN ) );
		core::Real CC_mz = getDensityMap().matchRes( resid, rsd_copy, pose, NULL, false );

		// rescore with orig pose
		getDensityMap().matchRes( resid, rsd, pose, NULL, false );

		dCCdX[0] = (CC_px-CC_mx)/(2*NUM_DERIV_H_CEN);
		dCCdX[1] = (CC_py-CC_my)/(2*NUM_DERIV_H_CEN);
		dCCdX[2] = (CC_pz-CC_mz)/(2*NUM_DERIV_H_CEN);

		//////
		//////  debug : compare d_true to d_exact
		numeric::xyzVector<core::Real> dCCdX1 = dCCdxs_res[resid][atmid];
		TR << "   " <<  dCCdX<< "  ;  " <<  dCCdX1 << std::endl;
	}
}


/////////////////////////////////////
//  Compute the gradient (whole-structure CA density score)
void ElectronDensity::dCCdx_cen( int resid,
                                 numeric::xyzVector<core::Real> const &X,
                                 core::pose::Pose const &pose,
                                 numeric::xyzVector<core::Real> &dCCdX ) {
	// make sure map is loaded
	if (!isLoaded) {
		TR << "[ ERROR ]  ElectronDensity::dCCdx_cen called but no map is loaded!\n";
		dCCdX = numeric::xyzVector<core::Real>(0.0,0.0,0.0);
		exit(1);
	}
	// make sure map is loaded
	if (!isLoaded) {
		TR << "[ ERROR ]  ElectronDensity::dCCdx_cen called but no map is loaded!\n";
		dCCdX = numeric::xyzVector<core::Real>(0.0,0.0,0.0);
		exit(1);
	}

	static bool warned=false;
	if (!DensScoreInMinimizer) {
		if (!warned)
			TR << "[ WARNING ] dCCdx_cen called but DensityScoreInMinimizer = false ... returning 0" << std::endl;
		//warned = true;
		dCCdX = numeric::xyzVector<core::Real>(0.0,0.0,0.0);
		return;
	}

	if (! ExactDerivatives ) {
		dCCdX = dCCdxs_cen[resid];
	} else {
		//
		core::pose::Pose pose_copy = pose;
		id::AtomID id(2,resid);

		pose_copy.set_xyz(id, numeric::xyzVector<core::Real>( X[0]+NUM_DERIV_H_CEN,X[1],X[2] ) );
		core::Real CC_px = getDensityMap().matchPose( pose_copy , NULL, true );

		pose_copy.set_xyz(id, numeric::xyzVector<core::Real>( X[0]-NUM_DERIV_H_CEN,X[1],X[2] ) );
		core::Real CC_mx = getDensityMap().matchPose( pose_copy , NULL, true );

		pose_copy.set_xyz(id, numeric::xyzVector<core::Real>( X[0],X[1]+NUM_DERIV_H_CEN,X[2] ) );
		core::Real CC_py = getDensityMap().matchPose( pose_copy , NULL, true );

		pose_copy.set_xyz(id, numeric::xyzVector<core::Real>( X[0],X[1]-NUM_DERIV_H_CEN,X[2] ) );
		core::Real CC_my = getDensityMap().matchPose( pose_copy , NULL, true );

		pose_copy.set_xyz(id, numeric::xyzVector<core::Real>( X[0],X[1],X[2]+NUM_DERIV_H_CEN ) );
		core::Real CC_pz = getDensityMap().matchPose( pose_copy , NULL, true );

		pose_copy.set_xyz(id, numeric::xyzVector<core::Real>( X[0],X[1],X[2]-NUM_DERIV_H_CEN ) );
		core::Real CC_mz = getDensityMap().matchPose( pose_copy , NULL, true );

		// rescore with orig pose
		getDensityMap().matchPose( pose , NULL, true );

		dCCdX[0] = (CC_px-CC_mx)/(2*NUM_DERIV_H_CEN);
		dCCdX[1] = (CC_py-CC_my)/(2*NUM_DERIV_H_CEN);
		dCCdX[2] = (CC_pz-CC_mz)/(2*NUM_DERIV_H_CEN);

		//////  debug : compare d_true to d_exact
		numeric::xyzVector<core::Real> dCCdX1 = dCCdxs_cen[resid];
		TR << "   " <<  dCCdX<< "  ;  " <<  dCCdX1 << std::endl;
	}
}


/////////////////////////////////////
//  Compute the gradient (whole-structure allatom density score)
void ElectronDensity::dCCdx_aacen( int atmid, int resid,
                             numeric::xyzVector<core::Real> const &X,
                             core::pose::Pose const &pose,
                             numeric::xyzVector<core::Real> &dCCdX ) {
	// make sure map is loaded
	if (!isLoaded) {
		TR << "[ ERROR ]  ElectronDensity::dCCdx_aacen called but no map is loaded!\n";
		dCCdX = numeric::xyzVector<core::Real>(0.0,0.0,0.0);
		exit(1);
	}

	static bool warned=false;
	if (!DensScoreInMinimizer) {
		if (!warned)
			TR << "[ WARNING ] dCCdx_aacen called but DensityScoreInMinimizer = false ... returning 0" << std::endl;
		//warned = true;
		dCCdX = numeric::xyzVector<core::Real>(0.0,0.0,0.0);
		return;
	}

	// if we didn't score this residue
	//      (because it was missing density or masked or a symm copy)
	// then don't compute its derivative
	if ( (int)dCCdxs_aacen[resid].size() < atmid ) {
 		dCCdX = numeric::xyzVector<core::Real>(0.0,0.0,0.0);
		return;
	}

	if (! ExactDerivatives ) {
		dCCdX = dCCdxs_aacen[resid][atmid];
	} else {
		//
		core::pose::Pose pose_copy = pose;
		id::AtomID id(atmid,resid);

		pose_copy.set_xyz(id, numeric::xyzVector<core::Real>( X[0]+NUM_DERIV_H_CEN,X[1],X[2] ) );
		core::Real CC_px = getDensityMap().matchPose( pose_copy , NULL, true );

		pose_copy.set_xyz(id, numeric::xyzVector<core::Real>( X[0]-NUM_DERIV_H_CEN,X[1],X[2] ) );
		core::Real CC_mx = getDensityMap().matchPose( pose_copy , NULL, true );

		pose_copy.set_xyz(id, numeric::xyzVector<core::Real>( X[0],X[1]+NUM_DERIV_H_CEN,X[2] ) );
		core::Real CC_py = getDensityMap().matchPose( pose_copy , NULL, true );

		pose_copy.set_xyz(id, numeric::xyzVector<core::Real>( X[0],X[1]-NUM_DERIV_H_CEN,X[2] ) );
		core::Real CC_my = getDensityMap().matchPose( pose_copy , NULL, true );

		pose_copy.set_xyz(id, numeric::xyzVector<core::Real>( X[0],X[1],X[2]+NUM_DERIV_H_CEN ) );
		core::Real CC_pz = getDensityMap().matchPose( pose_copy , NULL, true );

		pose_copy.set_xyz(id, numeric::xyzVector<core::Real>( X[0],X[1],X[2]-NUM_DERIV_H_CEN ) );
		core::Real CC_mz = getDensityMap().matchPose( pose_copy , NULL, true );

		// rescore with orig pose
		getDensityMap().matchPose( pose , NULL, true );

		dCCdX[0] = (CC_px-CC_mx)/(2*NUM_DERIV_H_CEN);
		dCCdX[1] = (CC_py-CC_my)/(2*NUM_DERIV_H_CEN);
		dCCdX[2] = (CC_pz-CC_mz)/(2*NUM_DERIV_H_CEN);

		//////
		//////  debug : compare d_true to d_exact
		numeric::xyzVector<core::Real> dCCdX1 = dCCdxs_aacen[resid][atmid];
		TR << "   " <<  dCCdX<< "  ;  " <<  dCCdX1 << std::endl;
	}
}

/////////////////////////////////////
//  Compute the gradient w.r.t. patterson correlation
void ElectronDensity::dCCdx_pat( int atmid, int resid,
                                  numeric::xyzVector<core::Real> const &X,
                                  core::pose::Pose const &pose,
                                  numeric::xyzVector<core::Real> &dCCdX ) {
	// make sure map is loaded
	if (!isLoaded) {
		TR << "[ ERROR ]  ElectronDensity::dCCdx_aacen called but no map is loaded!\n";
		exit(1);
	}

	static bool warned=false;
	if (!DensScoreInMinimizer) {
		if (!warned)
			TR << "[ WARNING ] dCCdx_pat called but DensityScoreInMinimizer = false ... returning 0" << std::endl;
		//warned = true;
		dCCdX = numeric::xyzVector<core::Real>(0.0,0.0,0.0);
		return;
	}

	// if we didn't score this residue
	//      (because it was missing density or masked or a symm copy)
	// then don't compute its derivative
	if ( (int)dCCdxs_pat[resid].size() < atmid ) {
 		dCCdX = numeric::xyzVector<core::Real>(0.0,0.0,0.0);
		return;
	}

	if (! ExactDerivatives ) {
		dCCdX = dCCdxs_pat[resid][atmid];
	} else {
		//
		core::pose::Pose pose_copy = pose;
		id::AtomID id(atmid,resid);

		pose_copy.set_xyz(id, numeric::xyzVector<core::Real>( X[0]+NUM_DERIV_H_CEN,X[1],X[2] ) );
		core::Real CC_px = getDensityMap().matchPoseToPatterson( pose_copy , false );

		pose_copy.set_xyz(id, numeric::xyzVector<core::Real>( X[0]-NUM_DERIV_H_CEN,X[1],X[2] ) );
		core::Real CC_mx = getDensityMap().matchPoseToPatterson( pose_copy , false );

		pose_copy.set_xyz(id, numeric::xyzVector<core::Real>( X[0],X[1]+NUM_DERIV_H_CEN,X[2] ) );
		core::Real CC_py = getDensityMap().matchPoseToPatterson( pose_copy , false );

		pose_copy.set_xyz(id, numeric::xyzVector<core::Real>( X[0],X[1]-NUM_DERIV_H_CEN,X[2] ) );
		core::Real CC_my = getDensityMap().matchPoseToPatterson( pose_copy , false );

		pose_copy.set_xyz(id, numeric::xyzVector<core::Real>( X[0],X[1],X[2]+NUM_DERIV_H_CEN ) );
		core::Real CC_pz = getDensityMap().matchPoseToPatterson( pose_copy , false );

		pose_copy.set_xyz(id, numeric::xyzVector<core::Real>( X[0],X[1],X[2]-NUM_DERIV_H_CEN ) );
		core::Real CC_mz = getDensityMap().matchPoseToPatterson( pose_copy , false );

		// rescore with orig pose
		getDensityMap().matchPoseToPatterson( pose , false );

		dCCdX[0] = (CC_px-CC_mx)/(2*NUM_DERIV_H_CEN);
		dCCdX[1] = (CC_py-CC_my)/(2*NUM_DERIV_H_CEN);
		dCCdX[2] = (CC_pz-CC_mz)/(2*NUM_DERIV_H_CEN);

		//////
		//////  debug : compare d_true to d_exact
		numeric::xyzVector<core::Real> dCCdX1 = dCCdxs_pat[resid][atmid];
		TR << resid << "   " << atmid << "  ::  " <<  dCCdX << "  ::  " <<  dCCdX1 << std::endl;
	}
}

/////////////////////////////////////
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

/////////////////////////////////////
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

	if (!mapin) {
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
	if ( !mapin.read(mapString, 3) || (std::string(mapString) != "MAP")) {
		TR << "[ ERROR ]  'MAP' string missing, not a valid MRC map.  Not loading map." << std::endl;
		return false;
	}
	// Check the file endianness
	if (mode != 2) {
		swap4_aligned(&mode, 1);
		if (mode != 2) {
			TR << "[ ERROR ]   Non-real (32-bit float) data types are unsupported.  Not loading map." << std::endl;
			return false;
		} else {
			swap = true; // enable byte swapping
		}
	}

	// Swap all the information obtained from the header
	if (swap) {
		swap4_aligned(extent, 3);
		swap4_aligned(&origin[0], 3);
		swap4_aligned(&altorigin[0], 3);
		swap4_aligned(&grid[0], 3);
		swap4_aligned(&cellDimensions[0], 3);
		swap4_aligned(&cellAngles[0], 3);
		swap4_aligned(crs2xyz, 3);
		swap4_aligned(&symBytes, 1);
	}

	if (reso != 0) TR << " Setting resolution to " << reso << "A" << std::endl;
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
	dataOffset = filesize - 4L*((long long)extent[0]*(long long)extent[1]*(long long)extent[2]);
	if (dataOffset != (CCP4HDSIZE + symBytes)) {
		if (dataOffset == CCP4HDSIZE) {
			// Bogus symmetry record information
			TR << "[ WARNING ] File contains bogus symmetry record.  Continuing." << std::endl;
			symBytes = 0;
		} else if (dataOffset < CCP4HDSIZE) {
			TR << "[ ERROR ] File appears truncated and doesn't match header.  Not loading map." << std::endl;
			return false;
		} else if ((dataOffset > CCP4HDSIZE) && (dataOffset < (1024*1024))) {
				// Fix for loading SPIDER files which are larger than usual
		 		// In this specific case, we must absolutely trust the symBytes record
				dataOffset = CCP4HDSIZE + symBytes;
				TR << "[ WARNING ]  File is larger than expected and doesn't match header.  Reading anyway." << std::endl;
		} else {
			TR << "[ ERROR ] File is MUCH larger than expected and doesn't match header.  Not loading map." << std::endl;
			TR << dataOffset  << std::endl;
			return false;
		}
	}

	// Read symmetry records -- organized as 80-byte lines of text.
	utility::vector1< std::string > symList;
	symData[80]='\0';
	if (symBytes != 0) {
		TR << "Symmetry records found:" << std::endl;
		mapin.seekg(CCP4HDSIZE, std::ios::beg);
		for (int i = 0; i < symBytes/80; i++) {
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
	if (grid[0] == 0 && extent[0] > 0) {
		grid[0] = extent[0] - 1;
		TR << "[ WARNING ] Fixed X interval count.  Continuing." << std::endl;
  	}
	if (grid[1] == 0 && extent[1] > 0) {
		grid[1] = extent[1] - 1;
		TR << "[ WARNING ]  Fixed Y interval count.  Continuing." << std::endl;
  	}
	if (grid[2] == 0 && extent[2] > 0) {
		grid[2] = extent[2] - 1;
		TR << "[ WARNING ]  Fixed Z interval count.  Continuing." << std::endl;
	}

	// Mapping between CCP4 column, row, section and Cartesian x, y, z.
	if (crs2xyz[0] == 0 && crs2xyz[1] == 0 && crs2xyz[2] == 0) {
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

	for (coord[2] = 1; coord[2] <= extent[2]; coord[2]++) {
		for (coord[1] = 1; coord[1] <= extent[1]; coord[1]++) {
			// Read an entire row of data from the file, then write it into the
			// datablock with the correct slice ordering.
			if ( mapin.eof() ) {
				TR << "[ ERROR ] Unexpected end-of-file. Not loading map." << std::endl;
				return false;
			}
			if ( mapin.fail() ) {
				TR << "[ ERROR ] Problem reading the file. Not loading map." << std::endl;
				return false;
			}
			if ( !mapin.read( reinterpret_cast< char* >(rowdata), sizeof(float)*extent[0]) ) {
				TR << "[ ERROR ] Error reading data row. Not loading map." << std::endl;
				return false;
			}

			for (coord[0] = 1; coord[0] <= extent[0]; coord[0]++) {
				density( coord[xyz2crs[0]], coord[xyz2crs[1]], coord[xyz2crs[2]]) = rowdata[coord[0]-1];
			}
		}
	}

	if (swap == 1)
		swap4_aligned( &density[0], vol_xySize * vol_zsize);
	delete [] rowdata;

	this->origin[0] = origin[xyz2crs[0]];
	this->origin[1] = origin[xyz2crs[1]];
	this->origin[2] = origin[xyz2crs[2]];

	// grid doesnt seemed to get remapped in ccp4 maps
	this->grid[0] = grid[0];
	this->grid[1] = grid[1];
	this->grid[2] = grid[2];

	// advanced: force the apix value different than what is provided
	numeric::xyzVector< core::Real > ori_scale(0,0,0);
	if (force_apix_ > 0) {
		ori_scale[0] = (force_apix_ * this->grid[0]) / cellDimensions[0];
		ori_scale[1] = (force_apix_ * this->grid[1]) / cellDimensions[1];
		ori_scale[2] = (force_apix_ * this->grid[2]) / cellDimensions[2];
		cellDimensions[0] = force_apix_ * this->grid[0];
		cellDimensions[1] = force_apix_ * this->grid[1];
		cellDimensions[2] = force_apix_ * this->grid[2];
		TR << "Forcing apix to " << force_apix_ << std::endl;
	}

	///////////////////////////////////
	/// POST PROCESSING
	// expand to unit cell
	this->computeCrystParams();
	TR << " voxel vol.: " << voxel_volume() << std::endl;

	// mrc format maps occasionally specify a real-valued origin in a different spot in the header
 	if (  altorigin[0]!=0 &&  altorigin[0]!=0 &&  altorigin[0]!=0 &&
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
 	}	else {
		use_altorigin = false;
	}

	// if we force the apix, adjust the origin accordingly
	if (force_apix_ > 0) {
		this->origin = numeric::xyzVector<core::Real>(
			this->origin[0]*ori_scale[0] ,
			this->origin[1]*ori_scale[1] ,
			this->origin[2]*ori_scale[2] );

 		TR << "Force apix repositioning the origin:\n";
 		TR << "     origin =" << this->origin[0] << " x " << this->origin[1] << " x " << this->origin[2] << std::endl;
	}

	this->efforigin = this->origin;
	this->expandToUnitCell();

	// resample the map
	if (gridSpacing > 0) this->resize( gridSpacing );

	// >>> figure out effective B <<<
	core::Real max_del_grid = std::max( cellDimensions[0]/((double)this->grid[0]) , cellDimensions[1]/((double)this->grid[1]) );
	max_del_grid = std::max( max_del_grid , cellDimensions[2]/((double)this->grid[2]) );

	// if input resolution is given, use that as long as it is greater than grid spacing
	TR << "Minimum resolution = " << 2*max_del_grid << std::endl;

	// finally, set effective B factor
	minimumB = 16*max_del_grid*max_del_grid;
	TR << "Minimum B factor = " << minimumB << std::endl;

	// if "auto" input resolution is given (and there is no b factor fitting), assume oversampling
	if (reso == 0 && !legacy_)
		max_del_grid *= 1.5;
	else
		max_del_grid = std::max( max_del_grid, reso/2 );

	// if input resolution is given, use that as long as it is greater than grid spacing
	TR << "Effective resolution = " << 2*max_del_grid << std::endl;

	// finally, set effective B factor
	if (effectiveB == 0)  // may be directly set via a flag
		effectiveB = 16*max_del_grid*max_del_grid;
	else
		effectiveB = std::max( effectiveB, minimumB );
	TR << "Effective B factor = " << effectiveB << std::endl;

	// adjust mask widths based on resolution
	OneGaussianScattering cscat = get_A( "C" );

	// make sure mask extends >= 3 carbon STDEVS
	core::Real mask_min = 3.0 * sqrt( effectiveB / (2*M_PI*M_PI) );
	//core::Real mask_min = 2; 				// ptc - use this threshold for electron density rotamer recovery
	if (ATOM_MASK < mask_min) {
		TR << "Override ATOM_MASK (was " << ATOM_MASK << ", now " << mask_min << ")" << std::endl;
		ATOM_MASK = mask_min;
	}
	mask_min = 3.0 * sqrt(2.0) * std::max( 2.4+1.6*reso , reso ) / M_PI;
	if (CA_MASK < mask_min) {
		TR << "Override CA_MASK (was " << CA_MASK << ", now " << mask_min << ")" << std::endl;
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


/////////////////////////////////////
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
	for (int i=1; i<=(int)symList.size(); ++i) {
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
		numeric::xyzMatrix< core::Real > rot(0);
		numeric::xyzVector< core::Real > trans(0,0,0);
		int k;
		for (int j=1; j<=3; ++j) {
			k = rows[j].find('/');
			if (k != (int)std::string::npos) {
				// numerator
				int startNum = rows[j].find_last_not_of("0123456789",k-1) + 1;
				int startDenom = k+1;
				float denom = std::atof( &rows[j][startDenom]);

				// make sure this shift corresponds to a point in the map
				core::Size oldMinMult = MINMULT[j-1];
				while ( std::fmod(MINMULT[j-1] , denom) > 1e-6 ) MINMULT[j-1] += oldMinMult;

				trans[j-1] = std::atof( &rows[j][startNum]) / denom;
			} else {
				trans[j-1] = 0;
			}

			if (rows[j].find("-X") != std::string::npos)
				rot(j,1) = -1;
			else if (rows[j].find("X") != std::string::npos)
				rot(j,1) = 1;

			if (rows[j].find("-Y") != std::string::npos)
				rot(j,2) = -1;
			else if (rows[j].find("Y") != std::string::npos)
				rot(j,2) = 1;

			if (rows[j].find("-Z") != std::string::npos)
				rot(j,3) = -1;
			else if (rows[j].find("Z") != std::string::npos)
				rot(j,3) = 1;
		}

		symmOps.push_back( RT( rot, trans ) );

		// remember rotations only for patterson case
		//   only put one of (-aX,-bY,-cZ) and ( aX,bY,cZ)
		// ignore translations
		bool dup = false;
		for (int i=1; i<=(int)symmRotOps.size(); ++i) {
			numeric::xyzMatrix<core::Real> R_i = symmRotOps[i];
			if ( (   rot(1,1) == R_i(1,1) && rot(1,2) == R_i(1,2) && rot(1,3) == R_i(1,3)
			      && rot(2,1) == R_i(2,1) && rot(2,2) == R_i(2,2) && rot(2,3) == R_i(2,3)
			      && rot(3,1) == R_i(3,1) && rot(3,2) == R_i(3,2) && rot(3,3) == R_i(3,3) ) ||
			   (     rot(1,1) == -R_i(1,1) && rot(1,2) == -R_i(1,2) && rot(1,3) == -R_i(1,3)
			      && rot(2,1) == -R_i(2,1) && rot(2,2) == -R_i(2,2) && rot(2,3) == -R_i(2,3)
			      && rot(3,1) == -R_i(3,1) && rot(3,2) == -R_i(3,2) && rot(3,3) == -R_i(3,3) ) ) {
				dup = true;
				break;
			}
		}
		if (!dup) {
			//std::cerr << "   ADD ROTATION " << std::endl
			//          << rot(1,1) << "," << rot(1,2) << "," << rot(1,3) << std::endl
			//          << rot(2,1) << "," << rot(2,2) << "," << rot(2,3) << std::endl
			//          << rot(3,1) << "," << rot(3,2) << "," << rot(3,3) << std::endl;
			symmRotOps.push_back( rot );
		}
	}

	return;
}

/////////////////////////////////////
// expand density to cover unit cell
// maintain origin
void ElectronDensity::expandToUnitCell() {
	numeric::xyzVector< int > extent( density.u1(), density.u2(), density.u3() );

	// if it already covers unit cell do nothing
	if ( grid[0] == extent[0] && grid[1] == extent[1] && grid[2] == extent[2] )
		return;

	ObjexxFCL::FArray3D< float > newDensity( grid[0],grid[1],grid[2], 0.0 );

	// copy the block
	int limX=std::min(extent[0],grid[0]),
	    limY=std::min(extent[1],grid[1]),
	    limZ=std::min(extent[2],grid[2]);
	for (int x=1; x<=limX; ++x)
	for (int y=1; y<=limY; ++y)
	for (int z=1; z<=limZ; ++z) {
		newDensity( x,y,z ) = density( x,y,z );
	}

	// apply symmetry
	// why backwards? it is a mystery
	for (int x=grid[0]; x>=1; --x)
	for (int y=grid[1]; y>=1; --y)
	for (int z=grid[2]; z>=1; --z) {
		if (x <= limX && y <= limY && z <= limZ)
			continue;

		numeric::xyzVector<core::Real> fracX(
			( (core::Real)x + origin[0] - 1 ) / grid[0],
			( (core::Real)y + origin[1] - 1 ) / grid[1],
			( (core::Real)z + origin[2] - 1 ) / grid[2] );

		for (int symm_idx=1; symm_idx<=(int)symmOps.size(); symm_idx++) {
			numeric::xyzVector<core::Real> SfracX =
				symmOps[symm_idx].get_rotation() * fracX +  symmOps[symm_idx].get_translation();

			// indices of symm copy
			int Sx = pos_mod((int)floor(SfracX[0]*grid[0]+0.5 - origin[0]) , grid[0]) + 1;
			int Sy = pos_mod((int)floor(SfracX[1]*grid[1]+0.5 - origin[1]) , grid[1]) + 1 ;
			int Sz = pos_mod((int)floor(SfracX[2]*grid[2]+0.5 - origin[2]) , grid[2]) + 1 ;

			if (Sx <= limX && Sy <= limY && Sz <= limZ) {
				newDensity( x,y,z ) = density(Sx,Sy,Sz);
			}
		}
	}

	if (basic::options::option[ basic::options::OptionKeys::edensity::debug ]()) {
		ElectronDensity(density,1.0, numeric::xyzVector< core::Real >(0,0,0), false ).writeMRC( "rho_before_expand.mrc");
		ElectronDensity(newDensity,1.0, numeric::xyzVector< core::Real >(0,0,0), false ).writeMRC( "rho_after_expand.mrc");
	}

	// new map!
	density = newDensity;
}


void ElectronDensity::showCachedScores( utility::vector1< int > const &reses ) {
	for (int i=1; i<=(int)reses.size(); ++i ) {
		TR << "   " << reses[i] << ": " << CCs[reses[i]] << std::endl;
	}
}

/////////////////////////////////////
//  DEBUGGING: write _this_ map in MATLAB v5 format
bool ElectronDensity::writeMAT(std::string filestem) {
	using namespace std;

	string outfile = filestem + ".mat";
	numeric::xyzVector<int> dims( density.u1(), density.u2(), density.u3() );

	ofstream out;
	out.open(outfile.c_str(), ios::out);

	// write file header
	short buff_s;
	long buff_l;
	out.write(
		"MATLAB 5.0 MAT-file, Platform: GLNX86   "
		"                                        "
		"                                        "
		"    ", 124);
	buff_s=0x0100;            // endian-ness test
	out.write((char*)&buff_s, 2);
	buff_s=0x4D49;            // code(?)
	out.write((char*)&buff_s, 2);

	// write data block header
	buff_l = 0x0000000E;      // data type (array)
	out.write((char*)&buff_l, 4);
	buff_l = 4*dims[0]*dims[1]*dims[2] + 64;      // number of bytes
	out.write((char*)&buff_l, 4);

	// write array header
	buff_l = 0x00000006;      // array data type (single)
	out.write((char*)&buff_l, 4);
	buff_l = 0x00000008;      // # bytes in array data block
	out.write((char*)&buff_l, 4);
	buff_l = 0x00000007;      // global/logical flags
	out.write((char*)&buff_l, 4);
	buff_l = 0x00000000;      // padding
	out.write((char*)&buff_l, 4);

	// write dims array
	buff_l = 0x00000005;      // dims data type (int32)
	out.write((char*)&buff_l, 4);
	buff_l = 0x0000000C;      // # bytes in dims data block
	out.write((char*)&buff_l, 4);
	out.write((char*)&dims[0], 12);
	buff_l = 0x00000000;      // padding
	out.write((char*)&buff_l, 4);

	// write name array
	buff_l=0x0001;            // name data type (int8)
	out.write((char*)&buff_l, 4);
	buff_l=0x0005;            // # bytes in name
	out.write((char*)&buff_l, 4);
	out.write("match   ", 8); // pad with spaces so next write is word-aligned

	// write array
	buff_l = 0x00000007;      // array data type (single)
	out.write((char*)&buff_l, 4);
	buff_l = 4*dims[0]*dims[1]*dims[2];
	out.write((char*)&buff_l, 4);

	for (int k=1; k<=dims[2]; k++) {
		for (int j=1; j<=dims[1]; j++) {
			for (int i=1; i<=dims[0]; i++) {
				float buff_f = (float) density(i,j,k);
				out.write((char*)&buff_f, 4);
			}
		}
	}

	out.close();

	return true;
}

/////////////////////////////////////
//  DEBUGGING: write _this_ map in MRC format
bool ElectronDensity::writeMRC(std::string mapfilename, bool writeRhoCalc, bool writeRhoMask) {
	std::fstream outx( (mapfilename).c_str() , std::ios::binary | std::ios::out );

	float buff_f;
	int buff_i;
	float buff_vf[3];
	int buff_vi[3];
	int symBytes = 0;

	if (!outx ) {
		TR.Error << "[ ERROR ]  Error opening MRC map for writing." << std::endl;
		return false;
	}

	runtime_assert( !( writeRhoCalc && writeRhoMask ) ); // for now

	if (writeRhoCalc)
		runtime_assert( density.u1() == rho_calc.u1() && density.u2() == rho_calc.u2() && density.u3() == rho_calc.u3() );
	if (writeRhoMask)
		runtime_assert( density.u1() == inv_rho_mask.u1() && density.u2() == inv_rho_mask.u2() && density.u3() == inv_rho_mask.u3() );

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
	for (int i=0; i<25; ++i) {
		outx.write(reinterpret_cast <char*>(&buff_i), sizeof(int));
	}

	// alt origin (MRC)
	float ori_float[3]={0,0,0};
	numeric::xyzVector<core::Real> origin_realspace(0,0,0);
	if (use_altorigin) {
		idx2cart ( numeric::xyzVector<core::Real>(1,1,1), origin_realspace );
		ori_float[0] = (float)( origin_realspace[0] );
		ori_float[1] = (float)( origin_realspace[1] );
		ori_float[2] = (float)( origin_realspace[2] );
	}
 	outx.write(reinterpret_cast <char*>(ori_float), sizeof(float)*3);

	// Write "MAP" at byte 208, indicating a CCP4 file.
	char buff_s[80]; strcpy(buff_s, "MAP DD");
	outx.write(reinterpret_cast <char*>(buff_s), 8);

	// fill remainder of head with junk
	int nJunkWords = (CCP4HDSIZE - 216) /4;
	buff_i = 0;
	for (int i=0; i<nJunkWords; ++i) {
		outx.write(reinterpret_cast <char*>(&buff_i), sizeof(int));
	}

	// data
	int coord[3];
	for (coord[2] = 1; coord[2] <= density.u3(); coord[2]++) {
		for (coord[1] = 1; coord[1] <= density.u2(); coord[1]++) {
			for (coord[0] = 1; coord[0] <= density.u1(); coord[0]++) {
				buff_f = (float) density(coord[0],coord[1],coord[2]);
				if (writeRhoCalc) buff_f = (float) rho_calc(coord[0],coord[1],coord[2]);
				if (writeRhoMask) buff_f = 1.0 - (float) inv_rho_mask(coord[0],coord[1],coord[2]);
				outx.write(reinterpret_cast <char*>(&buff_f), sizeof(float));
			}
		}
	}

	return true;
}


void ElectronDensity::density_change_trigger() {
	fastdens_score.clear();
	fastgrid = numeric::xyzVector< core::Real >(0,0,0);

	PattersonEpsilon.clear(); p_o.clear();
	p_extent = numeric::xyzVector< core::Real >(0,0,0);

	Fdensity.clear();

	// not strictly necessary
	rho_calc.clear();
	rho_solv.clear();
	inv_rho_mask.clear();
	Frho_calc.clear();
	Frho_solv.clear();

	computeStats();
	computeGradients();
}


/////////////////////////////////////
// resize a map (using FFT-interpolation)
void ElectronDensity::resize( core::Real approxGridSpacing ) {
	// potentially expand map to cover entire unit cell
	if ( grid[0] != density.u1() || grid[1] != density.u2() || grid[2] != density.u3() ){
		TR << "[ ERROR ] resize() not supported for maps not covering the entire unit cell."<< std::endl;
		TR << "   " << grid[0] << " != " << density.u1()
		   << " || " << grid[1] << " != " << density.u2()
		   << " || " << grid[2] << " != " << density.u3() << std::endl;
		exit(1);
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
	for (int i=0; i< newDims[0]*newDims[1]*newDims[2] ; ++i)
		density[i] = (float)newDensity[i];
	this->grid = newGrid;
	this->efforigin = origin = newOri;

	TR << " new extent: " << density.u1() << " x " << density.u2() << " x " << density.u3() << std::endl;
	TR << " new origin: " << origin[0] << " x " << origin[1] << " x " << origin[2] << std::endl;
	TR << "   new grid: " << grid[0] << " x " << grid[1] << " x " << grid[2] << std::endl;
	TR << "     new VV: " << voxel_volume() << std::endl;
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
	core::Real sum=0, sum2=0, sumCoM=0;
	dens_max = -std::numeric_limits< core::Real >::max();
	dens_min = std::numeric_limits< core::Real >::max();
	centerOfMass = numeric::xyzVector<core::Real>(0,0,0);

	int N = density.u1()*density.u2()*density.u3();
	//core::Real N2 = square((core::Real)N);

	for (int i=1; i<=density.u1(); i++)
	for (int j=1; j<=density.u2(); j++)
	for (int k=1; k<=density.u3(); k++) {
		sum += density(i,j,k);
		sum2 += density(i,j,k)*density(i,j,k);
		dens_max = std::max(dens_max, (core::Real)density(i,j,k));
		dens_min = std::min(dens_min, (core::Real)density(i,j,k));

 		centerOfMass[0] += i*density(i,j,k);
 		centerOfMass[1] += j*density(i,j,k);
 		centerOfMass[2] += k*density(i,j,k);
 		sumCoM += density(i,j,k);
	}
	centerOfMass /= sumCoM;
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
	if (D == 0) {
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
