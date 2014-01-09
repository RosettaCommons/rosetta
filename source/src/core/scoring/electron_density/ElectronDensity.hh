// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/electron_density/ElectronDensity.hh
/// @brief  Scoring a structure against an electron density map
/// @author Frank DiMaio

#ifndef INCLUDED_core_scoring_electron_density_ElectronDensity_hh
#define INCLUDED_core_scoring_electron_density_ElectronDensity_hh

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/Atom.fwd.hh>
#include <core/scoring/electron_density/util.hh>
#include <core/scoring/electron_density/xray_scattering.hh>
#include <core/scoring/electron_density/ElectronDensity.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

// Utility headers
#include <utility/exit.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/FArray4D.hh>
#include <numeric/xyzVector.hh>

// C++ headers
#include <string>
#include <map>
#include <complex>

#include <core/kinematics/RT.hh>
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>


namespace core {
namespace scoring {
namespace electron_density {

const core::Real MAX_FLT = 1e37;

float pos_mod(float x,float y);

class ElectronDensity : public utility::pointer::ReferenceCount {
public:
	///@brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~ElectronDensity();
	/// @brief constructor
	ElectronDensity();

	/// @brief calulated density from a vector of poses
	ElectronDensity( utility::vector1< core::pose::PoseOP > poses, core::Real reso, core::Real apix );

	/// @brief constructor from an FArray3D (debugging only)
	template<class T>
	ElectronDensity( ObjexxFCL::FArray3D< T > const &map,
	                 core::Real apix = 1.0,
	                 numeric::xyzVector< core::Real > new_origin=numeric::xyzVector< core::Real >(0,0,0),
	                 bool fftshift=false) {
		init();

		isLoaded = true;
		efforigin = origin = new_origin;
		grid = numeric::xyzVector< int >(map.u1(),map.u2(),map.u3());
		cellDimensions = numeric::xyzVector< float >(apix*map.u1(),apix*map.u2(),apix*map.u3());
		cellAngles = numeric::xyzVector< float >(90,90,90);
		density.dimension( map.u1(),map.u2(),map.u3() );
		if (fftshift) origin -= grid/2;
		for (int i=1; i<=(int)map.u1(); ++i) {
			int fi = (int)(fftshift ? pos_mod( i-(map.u1()/2)-1 , map.u1())+1 : i);
			for (int j=1; j<=(int)map.u2(); ++j) {
				int fj = (int)(fftshift ? pos_mod( j-(map.u2()/2)-1 , map.u2())+1 : j);
				for (int k=1; k<=(int)map.u3(); ++k) {
					int fk = (int)(fftshift ? pos_mod( k-(map.u3()/2)-1 , map.u3())+1 : k);
					density(fi,fj,fk) = (float)map(i,j,k);
				}
			}
		}
	}

	/// @brief Initialize map from cmd line options
	void
	init();

	/// @brief Load an MRC (="new-CCP4") density map
	bool
	readMRCandResize(
		std::string mapfile,
		core::Real reso=5.0,
		core::Real gridSpacing=0.0);

	/// @brief Load an MRC (="new-CCP4") density map
	bool readMRCandResize(
		std::istream & mapin,
		std::string mapfile, // just used as information about the map
		core::Real reso=5.0,
		core::Real gridSpacing=0.0);

	/// @brief (debugging) Write MRC mapfile
	bool writeMRC(std::string mapfilestem, bool writeRhoCalc=false, bool writeRhoMask=false);

	/// @brief (debugging) Write MATLAB v5 mapfile
	bool writeMAT(std::string mapfilestem);

	/// @brief Align a pose about a 2D rotation axis
	numeric::xyzMatrix< core::Real > rotAlign2DPose( core::pose::Pose const &pose, std::string axis );

	/// @brief Quickly matches a centroid pose into a low-resolution density map
	///   by placing a single Gaussian at each CA
	core::Real
	matchCentroidPose(
		core::pose::Pose const &pose,
		core::conformation::symmetry::SymmetryInfoCOP symmInfo = NULL,
		bool cacheCCs=false
	);

	/// @brief Match a pose into a medium-resolution density map
	///   by placing a single Gaussian at each atom
	core::Real
	matchPose(
		core::pose::Pose const &pose,
		core::conformation::symmetry::SymmetryInfoCOP symmInfo=NULL,
		bool cacheCCs=false
	);


	/////  "Density tools"
	/////    -- all use "light pose" representation
	/////
	/// @brief  mask the density map (like chimera 'zone' but sigmoidal mask)
	void
	maskDensityMap( poseCoords const &pose, core::Real radius );

	/// @brief get resolution bins (informational)
	void
	getResolutionBins(
			core::Size nbuckets, core::Real maxreso, core::Real minreso,
 			utility::vector1< core::Real > &,
		 	utility::vector1< core::Size > &,
			bool S2_bin=false );

	/// @brief Compute map intensities
	utility::vector1< core::Real >
	getIntensities( core::Size nbuckets, core::Real maxreso, core::Real minreso, bool S2_bin=false );

	/// @brief Compute map intensities, masked by a pose
	utility::vector1< core::Real >
	getIntensitiesMasked( poseCoords const &pose, core::Size nbuckets, core::Real maxreso, core::Real minreso, bool S2_bin=false );

	/// @brief Compute intensities from a pose
	void
	getIntensities( poseCoords const &pose, core::Size nbuckets, core::Real maxreso, core::Real minreso, utility::vector1< core::Real > &Imodel, bool S2_bin=false);

	/// @brief Compute map-map FSC & errors
	void
	getMapMapError(
			ObjexxFCL::FArray3D< float > const &map2,
			core::Size nbuckets, core::Real maxreso, core::Real minreso,
			utility::vector1< core::Real >& FSC,
			utility::vector1< core::Real >& phaseError,
			utility::vector1< core::Real >& complexPlaneError,
			bool S2_bin=false);

	/// @brief Compute model-map FSC & errors
	void
	getModelMapError(
			poseCoords const &pose,
			core::Size nbuckets, core::Real maxreso, core::Real minreso,
			utility::vector1< core::Real > const &mapmapPhaseError,
			utility::vector1< core::Real > const &mapmapComplexPlaneError,
			utility::vector1< core::Real >& FSC,
			utility::vector1< core::Real >& phaseError,
			core::Real &phaseErrorSum,
			core::Real &complexPlaneErrorK,
			core::Real &complexPlaneErrorProb,
			bool masked=false, bool S2_bin=false, core::Real mask_radius=0);

	/// @brief Compute model-map FSC & errors
	void
	getFSC(
			poseCoords const &pose,
			core::Size nbuckets, core::Real maxreso, core::Real minreso,
			utility::vector1< core::Real >& FSC,
			bool masked=false, bool S2_bin=false, core::Real mask_radius=0) {
		utility::vector1< core::Real > mapmapPhaseError(nbuckets,0.1);
		utility::vector1< core::Real > mapmapComplexPlaneError(nbuckets,0.1);
		utility::vector1< core::Real > phaseError;
		core::Real phaseErrorSum, complexPlaneErrorK, complexPlaneErrorProb;
		getModelMapError(
				pose, nbuckets, maxreso, minreso,
				mapmapPhaseError, mapmapComplexPlaneError,
				FSC, phaseError, phaseErrorSum, complexPlaneErrorK, complexPlaneErrorProb, masked, S2_bin, mask_radius );
	}

/*  OLD VERSIONS
	void
	getFSC(
		poseCoords const &pose,
		core::Size nbuckets,
		core::Real maxreso, core::Real minreso,
		utility::vector1< core::Real >&, utility::vector1< core::Real > &,
		bool masked=false, bool S2_bin=false, core::Real mask_radius=0.0, bool verbose=false );

	/// @brief Compute model-map expected phase error
	void
	getMLE(
		poseCoords const &pose,
		core::Size nbuckets,
		core::Real maxreso, core::Real minreso,
		utility::vector1< core::Real > const &,
		utility::vector1< core::Real >&, Real &errS2,
		bool masked=false, bool S2_bin=false, bool phase_only=false, core::Real mask_radius=0.0 );

	/// @brief Compute map-map FSC
	void
	getFSC(
			ObjexxFCL::FArray3D< float > const &map2,
			core::Size nbuckets,
			core::Real maxreso, core::Real minreso,
			utility::vector1< core::Real >&, utility::vector1< core::Real > &,
			bool S2_bin=false, bool phase_only=false, bool verbose=false);
*/

	/// @brief Compute model-map RSCC
	core::Real
	getRSCC( poseCoords const &pose );

	/// @brief Compute map-map RSCC
	core::Real
	getRSCC( ObjexxFCL::FArray3D< float > const &map2 );

	/// @brief Compute intensities, update density
	void
	scaleIntensities( utility::vector1< core::Real > I_tgt, core::Real maxreso, core::Real minreso, bool S2_bin=false );

	void
	calcRhoC( poseCoords const &pose );

	// for now use masking to handle solvent
	/*
	void
	calcRhoCandSolvent( poseCoords const &pose );
	*/

	core::Real
	maxNominalRes();


	/////
	/////  Scorefunction stuff
	/////
	/// @brief Match a pose to a patterson map
	core::Real matchPoseToPatterson( core::pose::Pose const &pose, bool cacheCCs=false );

	/// @brief Rematch the pose to a patterson map, using previous rho_calc with only rsd changed
	core::Real rematchResToPatterson( core::conformation::Residue const &rsd ) const;

	/// @brief Update cached rho_calc by changing residue 'rsd'
	void updateCachedDensity( core::conformation::Residue const &rsd );

	/// @brief Match a residue's conformation to the density map.
	///   Backbone atoms from adjacent residues are also used for scoring.
	///   Returns the correlation coefficient between map and pose
	///   Internally stores per-res CCs, per-atom dCC/dxs
	core::Real
	matchRes(
		int resid,
		core::conformation::Residue const &rsd,
		core::pose::Pose const &pose,
		core::conformation::symmetry::SymmetryInfoCOP symmInfo=NULL,
		bool cacheCCs=false
	);

	/// @brief Match a residue's conformation to the density map.
	///    Same as matchRes, but using a fast approximation to the match function
	core::Real
	matchResFast( int resid,
		core::conformation::Residue const &rsd,
		core::pose::Pose const &pose,
		core::conformation::symmetry::SymmetryInfoCOP symmInfo=NULL,
		bool ignoreBs=false
	);

	/// @brief Computes the symmatric rotation matrices
	void
	compute_symm_rotations(
		core::pose::Pose const &pose,
		core::conformation::symmetry::SymmetryInfoCOP symmInfo=NULL
	);

	/// @brief Return the gradient of CC w.r.t. atom X's movement
	/// Uses information stored from the previous call to matchRes with this resid
	void  dCCdx_res(
		int atmid,
		int resid,
		numeric::xyzVector<core::Real> const &X,
		core::conformation::Residue const &rsd,
		core::pose::Pose const &pose,
		numeric::xyzVector<core::Real> &gradX
	);

	/// @brief Return the gradient of "fast CC" w.r.t. atom X's movement
	/// Uses information stored from the previous call to matchRes with this resid
	void
	dCCdx_fastRes(
		int atmid,
		int resid,
		numeric::xyzVector<core::Real> const &X,
		core::conformation::Residue const &rsd,
		core::pose::Pose const &pose,
		numeric::xyzVector<core::Real> &gradX
	);

	/// @brief Gradient of CC w.r.t B factors
	Real
	dCCdB_fastRes(
		int atmid, int resid,
		core::conformation::Residue const &rsd,
		core::pose::Pose const &pose
	);

	/// @brief calculate all gradients of CC w.r.t B factors (slow more precise version)
	void
	dCCdBs(
		core::pose::Pose const &pose,
		utility::vector1< core::Real>  & dE_dvars
	);

	//@brief Called to initialize scorefunction scaling for B-factor refinement
	void rescale_fastscoring_temp_bins(core::pose::Pose const &pose, bool init=false);

	/// @brief Return the gradient of CC w.r.t. res X's CA's movement
	/// Centroid-mode analogue of dCCdx
	void
	dCCdx_cen(
		int resid,
		numeric::xyzVector<core::Real> const &X,
		core::pose::Pose const &pose,
		numeric::xyzVector<core::Real> &gradX
	);

	/// @brief Return the gradient of whole-structure-CC w.r.t. atom X's movement
	/// non-sliding-window analogue of dCCdx
	void
	dCCdx_aacen(
		int atmid,
		int resid,
		numeric::xyzVector<core::Real> const &X,
		core::pose::Pose const &pose,
		numeric::xyzVector<core::Real> &gradX
	);

	/// @brief Return the gradient of patterson-CC w.r.t. atom X's movement
	void
	dCCdx_pat(
		int atmid,
		int resid,
		numeric::xyzVector<core::Real> const &X,
		core::pose::Pose const &pose,
		numeric::xyzVector<core::Real> &gradX
	);

	/// @brief Resets the counters used for derivative computation in
	///   sliding-window/fast scoring
	void clear_dCCdx_res_cache( core::pose::Pose const &pose );

	void compute_rho(
		core::pose::Pose const & pose,
		utility::vector1<core::id::AtomID> const & atom_ids,
		ObjexxFCL::FArray3D< double > & rho_calc,
		ObjexxFCL::FArray3D< double > & inv_rho_mask);

	numeric::xyzVector< double > match_fragment(
		ObjexxFCL::FArray3D< double > const & rho_calc,
		ObjexxFCL::FArray3D< double > const & mask,
		ObjexxFCL::FArray3D< double > const & rho_obs,
		core::Real radius = 6.0 );

	/// @brief Get the transformation from indices to Cartesian coords using 'real' origin
	numeric::xyzVector<core::Real> getTransform() {
		numeric::xyzVector<core::Real> idxX(  grid[0]-origin[0]+1 , grid[1]-origin[1]+1 , grid[2]-origin[2]+1 ) , cartX;
		idx2cart( idxX , cartX );
		return cartX;
	}

	/// @brief set # of residues
	void set_nres(int nres) {
		if ( (int) CCs.size() != nres ) {
			CCs.resize(nres, 0.0);
			dCCdxs_res.resize(nres);
			dCCdxs_pat.resize(nres);
			dCCdxs_cen.resize(nres);
			dCCdxs_aacen.resize(nres);
			symmap.clear();   // reset if #residues changes
		}
	}

	/// @brief Print cached CCs
	void showCachedScores( utility::vector1< int > const &reses );
	inline core::Real getCachedScore( core::Size resid ) {
		runtime_assert( resid <= CCs.size() );
		return CCs[resid];
	}

	//////////////////////////////////
	// property getters and setters
	inline void setUseDensityInMinimizer( bool newVal ) { DensScoreInMinimizer = newVal; }
	inline bool getUseDensityInMinimizer() const { return DensScoreInMinimizer; }

	inline void setUseExactDerivatives( bool newVal ) { ExactDerivatives = newVal; }
	inline bool getUseExactDerivatives() const { return ExactDerivatives; }

	inline void setSCscaling( core::Real SC_scalingin) { SC_scaling_=SC_scalingin; }
	inline core::Real getSCscaling() { return SC_scaling_; }

	inline void setWindow( core::Size window_in ) { WINDOW_ = window_in; }
	inline core::Size getWindow( ) { return WINDOW_; }

	inline void setScoreWindowContext( bool newVal ) { score_window_context_ = newVal; }
	inline bool getScoreWindowContext() const { return score_window_context_; }

	//////////////////////////////////
	//  map properties
	inline bool isMapLoaded() const { return this->isLoaded; };
	inline core::Real getNumDerivH() const { return NUM_DERIV_H; }
	inline core::Real getMean() const { return dens_mean; }
	inline core::Real getMin()  const { return dens_min;  }
	inline core::Real getMax()  const { return dens_max;  }
	inline core::Real getStdev() const { return dens_stdev; }
	inline core::Real getResolution( ) const { return this->reso; }
	inline numeric::xyzVector<core::Real> getCoM() const { return centerOfMass; }
	inline numeric::xyzVector<core::Real> getOrigin() const { return origin; }
	inline numeric::xyzVector<core::Real> getEffOrigin() const { return efforigin; }
	inline utility::vector1< core::kinematics::RT > getsymmOps() const { return symmOps; }

	inline core::Real getAtomMask( ) const { return ATOM_MASK; }

	void maskResidues( int scoring_mask ) {
		scoring_mask_[ scoring_mask ] = 1;
	}
	void maskResidues( utility::vector1< int > const & scoring_mask ) {
		for (core::Size i=1; i<= scoring_mask.size(); ++i)
			scoring_mask_[ scoring_mask[i] ] = 1;
	}
	void clearMask( ) {
		scoring_mask_.clear();
	}


	//////////////////////////////////
	//////////////////////////////////
	// raw data pointer
	inline ObjexxFCL::FArray3D< float > const & data() const { return density; };
	inline void set_data(ObjexxFCL::FArray3D< double > const & density_in) {
		assert(density.u1() == density_in.u1());
		assert(density.u2() == density_in.u2());
		assert(density.u3() == density_in.u3());
		for (Size i=0;i<density_in.size(); ++i) density[i] = density_in[i];
	};


	//////////////////////////////////
	//////////////////////////////////
	// helper functions to convert between indices and cartesian coords
	inline void cart2idx( numeric::xyzVector<core::Real> const & cartX , numeric::xyzVector<core::Real> &idxX ) const {
		numeric::xyzVector<core::Real> fracX = c2f*cartX;
		idxX = numeric::xyzVector<core::Real>( fracX[0]*grid[0] - efforigin[0] + 1,
		                                       fracX[1]*grid[1] - efforigin[1] + 1,
		                                       fracX[2]*grid[2] - efforigin[2] + 1);
	}

	template<class Q>
	void idx2cart( numeric::xyzVector<Q> const & idxX , numeric::xyzVector<core::Real> &cartX ) const {
		numeric::xyzVector<core::Real> fracX( (idxX[0]  + efforigin[0] -1 ) / grid[0],
		                                      (idxX[1]  + efforigin[1] -1 ) / grid[1],
		                                      (idxX[2]  + efforigin[2] -1 ) / grid[2] );
		cartX = f2c*fracX;
	}

	template<class Q>
	void idxoffset2cart( numeric::xyzVector<Q> const & idxX , numeric::xyzVector<core::Real> &cartX ) const {
		numeric::xyzVector<core::Real> fracX( ( (core::Real) idxX[0] ) / grid[0],
		                                      ( (core::Real) idxX[1] ) / grid[1],
		                                      ( (core::Real) idxX[2] ) / grid[2] );
		cartX = f2c*fracX;
	}

	//////////////////////////////////
	//////////////////////////////////
	// helper functions to convert between fractional and cartesian coords
	inline void cart2frac( numeric::xyzVector<core::Real> const & cartX , numeric::xyzVector<core::Real> & fracX ) const {
		fracX = c2f*(cartX);
	}
	inline void frac2cart( numeric::xyzVector<core::Real> const & fracX , numeric::xyzVector<core::Real> &cartX ) const {
		cartX = f2c*fracX;
	}

	numeric::xyzVector<core::Real> delt_cart(numeric::xyzVector<core::Real> const & cartX1, numeric::xyzVector<core::Real> const & cartX2);
	numeric::xyzVector<core::Real> get_cart_unitCell(numeric::xyzVector<core::Real> const & cartX);
	numeric::xyzVector<core::Real> get_nearest_UC(numeric::xyzVector<core::Real> const & cartX_in, numeric::xyzVector<core::Real> const & cartX_ref);


	numeric::xyzVector<core::Real> dens_grad ( numeric::xyzVector<core::Real> const & idxX ) const;

	/// resize the map via FFT
	void resize( core::Real approxGridSpacing );

	//// access cached data from last scored pose
	void get_symmMap(int vrtid, utility::vector1<int> &X_map, numeric::xyzMatrix<core::Real> &R) {
		runtime_assert( symmap.find( vrtid ) != symmap.end() );
		X_map = symmap[ vrtid ].first;
		R = symmap[ vrtid ].second;
	}

	// gets rotation vactor for subunit 'subunit' in last-scored pose (Rosetta symmetry)
	void get_R(int subunit, numeric::xyzMatrix<core::Real> &R) {
		runtime_assert( symmap.find( -subunit ) != symmap.end() );
		R = symmap[ -subunit ].second;
	}

///////////
// PRIVATE MEMBER FUNCTIONS
///////////
private:
	/// @brief The function is called everytime the density changes
	void
	density_change_trigger();

	// smooth an intensity spectrum
	void
	smooth_intensities(utility::vector1< core::Real > &) const;

	// helper functions for map statistics
	void computeGradients();
	void computeStats();
	int suggestRadius();

	// helper functions for symmetry
	void initializeSymmOps( utility::vector1< std::string > const & symList );
	void computeCrystParams();
	void expandToUnitCell();

	// setup patterson map scoring data
	void setup_patterson_first_time(core::pose::Pose const &pose);

	// setup fast density scoring data
	void setup_fastscoring_first_time(core::pose::Pose const &pose);

	// get Fdrho_d(xyz)
	// compute if not already computed
	utility::vector1< ObjexxFCL::FArray3D< std::complex<double> > * > getFdrhoc( OneGaussianScattering S );

	// get S2 (reciprocal space dist^2)
	double S2(int h, int k, int l) {
		return ( h*h*RcellDimensions[0]*RcellDimensions[0]
		           + k*k*RcellDimensions[1]*RcellDimensions[1]
		           + l*l*RcellDimensions[2]*RcellDimensions[2]
		           + 2*h*k*RcellDimensions[0]*RcellDimensions[1]*cosRcellAngles[2]
		           + 2*h*l*RcellDimensions[0]*RcellDimensions[2]*cosRcellAngles[1]
		           + 2*k*l*RcellDimensions[1]*RcellDimensions[2]*cosRcellAngles[0] );
	}

	// volume of 1 voxel
	double voxel_volume( ) {
		return V / (grid[0]*grid[1]*grid[2]);
	}

///////////
// DATA
///////////
private:
	// do we have a map loaded?
	bool isLoaded;

	// the density data array
	ObjexxFCL::FArray3D< float > density;

	// fft of density -- only used in FSC calc
	ObjexxFCL::FArray3D< std::complex<double> > Fdensity;

	// Controllable parameters
	std::map< core::Size, bool > scoring_mask_;
	core::Real reso, ATOM_MASK, CA_MASK, force_apix_, SC_scaling_;
	core::Real ATOM_MASK_PADDING;
	core::Size WINDOW_;
	bool score_window_context_, remap_symm_;

	// turn on legacy options
	bool legacy_;

	bool DensScoreInMinimizer, ExactDerivatives;
	core::Real NUM_DERIV_H, NUM_DERIV_H_CEN, PattersonB, PattersonMinR, PattersonMaxR;
	ObjexxFCL::FArray3D< float > PattersonEpsilon;

	// (patterson only) map resamped on p_calc grid
	ObjexxFCL::FArray3D< double >  p_o;
	double po_bar;

	// (fast scoring) precomputed rhocrhoo, d_rhocrhoo
	ObjexxFCL::FArray4D< double > fastdens_score;
	core::Size nkbins_;
	core::Real kmin_,kmax_,kstep_;
	numeric::xyzVector< int > fastgrid;           // grid & origin
	numeric::xyzVector< core::Real > fastorigin;  // for resampled maps

	///////////////////
	/// TONS OF CACHED STUFF <<< this should live in its own class!
	///////////////////
	// previously scored computed density map, fft(rho_calc), and patterson map
	ObjexxFCL::FArray3D< double > rho_calc, rho_solv, inv_rho_mask;

	core::Real rho_calc_sum;
	ObjexxFCL::FArray3D< std::complex<double> > Frho_calc, Frho_solv;
	utility::vector1<core::Size> bucket_counts;
	ObjexxFCL::FArray3D< core::Size > bucket_id;
	core::Real lowres_cut, hires_cut;
	ObjexxFCL::FArray3D< double > F_s2;
	ObjexxFCL::FArray3D< double > Pcalc;
	utility::vector1<core::Real> F2;

	// symm pointer matrices
	utility::vector1< ObjexxFCL::FArray3D< core::Size > > symm_ptrs;

	// (precomputed)
	// FFT of gradient of rho_calc with respec to an atom at the origin's movement in x/y/z
	// computed for each scatterer
	std::map< int , ObjexxFCL::FArray3D< std::complex<double> > > Fdrhoc_dx;
	std::map< int , ObjexxFCL::FArray3D< std::complex<double> > > Fdrhoc_dy;
	std::map< int , ObjexxFCL::FArray3D< std::complex<double> > > Fdrhoc_dz;

	// atoms, scattering used to calculate density map rho_calc
	utility::vector1< utility::vector1< numeric::xyzVector< core::Real > > > rho_calc_atms;
	utility::vector1< utility::vector1< OneGaussianScattering > > rho_calc_as;

	// cached patterson map statistics
	core::Real p_sumC, p_sumC2, p_sumO, p_sumO2, p_sumCO, p_vol;

	// patterson map is calculated in P1, in an alternate (padded) grid to avoid self peaks
	numeric::xyzVector< core::Real > d_min,d_max;  // (remember bounding coords)
	numeric::xyzVector< core::Real > p_extent, p_origin;
	numeric::xyzVector< core::Size > p_grid;
	numeric::xyzVector< core::Real > p_CoM;

	// map info
	core::Real minimumB;    // minimum B factor allowed by map
	core::Real effectiveB;  // B factor blurring based on map resolution
	numeric::xyzVector< int > grid;
	numeric::xyzVector< core::Real > origin, efforigin;
	bool use_altorigin;   // which field to write origin to ... only affects map outputting

	// cache scoring-related statistics
	utility::vector1<core::Real>  CCs;
	core::Real CC_cen, CC_aacen, CC_pat;
	utility::vector1< numeric::xyzVector< core::Real > > dCCdxs_cen;
	utility::vector1< utility::vector1< numeric::xyzVector< core::Real > > > dCCdxs_aacen;
	utility::vector1< utility::vector1< numeric::xyzVector< core::Real > > > dCCdxs_res;
	utility::vector1< utility::vector1< numeric::xyzVector< core::Real > > > dCCdxs_pat;

	///////////////////
	/// SYMMETRY (Rosetta's symmetry, not necessarily crystal symmetry)
	///////////////////
	// map vrtid -> subunit mapping, rotation
	// if (vrtid < 0) then it refers to the mapping from a non-vrt in subunit# -vrtid
	std::map< int , std::pair< utility::vector1<int> , numeric::xyzMatrix<core::Real> > > symmap;

	///////////////////
	/// VISUALIZATION-SPECIFIC
	///////////////////
	// gradients, used for displaying isocontoured surface
	//     ... mutable for the viewer to access
	// in non-graphics builds this never gets initialized
	mutable ObjexxFCL::FArray3D< double > coeff_grad_x, coeff_grad_y, coeff_grad_z;

	///////////////////
	/// CRYSTAL INFO
	///////////////////
	// TO DO --- put all this in a self-contained class
	// converting fractional to cartesian coords
	numeric::xyzMatrix<core::Real> f2c, c2f;

	// unit cell, reciprocal unit cell parameters, volume
	numeric::xyzVector<float> cellDimensions, cellAngles;
	numeric::xyzVector<float> RcellDimensions, cosRcellAngles;
	core::Real V, RV;

	// symmetric transforms (in frac. coords)
	// this is only used to expand the density data outside the ASU
	//    and is unrelated Rosetta's symmetric modelling
	utility::vector1< core::kinematics::RT > symmOps;
	utility::vector1< numeric::xyzMatrix<core::Real> > symmRotOps;

	// min multiples in each dim
	numeric::xyzVector<core::Size> MINMULT;

	// map statistics
	numeric::xyzVector<core::Real> centerOfMass;
	core::Real dens_mean, dens_min, dens_max, dens_stdev;
};

/// @brief The EDM instance
ElectronDensity& getDensityMap(std::string filename ="", bool force_reload =false);

/// @brief The EDM instance
ElectronDensity& getDensityMap_legacy(std::string filename ="", bool force_reload =false);


// x mod y, returns z in [0,y-1]
inline int pos_mod(int x,int y) {
	int r=x%y; if (r<0) r+=y;
	return r;
}
inline float pos_mod(float x,float y) {
	float r=std::fmod(x,y); if (r<0) r+=y;
	return r;
}
inline double pos_mod(double x,double y) {
	double r=std::fmod(x,y); if (r<0) r+=y;
	return r;
}

}
}
}


#endif

