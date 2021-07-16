// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief Utils for docking into density
/// @details
/// @author Frank DiMaio
/// @author Danny Farrell

#include <protocols/electron_density/DockIntoDensityUtils.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>

#include <core/scoring/electron_density/ElectronDensity.hh>
#include <core/scoring/electron_density/xray_scattering.hh>

#include <basic/Tracer.hh>


static basic::Tracer TR( "protocols.electron_density.DockIntoDensityUtils" );

namespace protocols {
namespace electron_density {


void
pose_spherical_samples(
	core::pose::Pose const &pose,
	ObjexxFCL::FArray3D< core::Real > & sigR,
	ObjexxFCL::FArray3D< core::Real > & epsR,
	protocols::electron_density::PoseSphericalSamplesOptions const & params
) {
	using namespace core;

	scoring::electron_density::ElectronDensity &density = core::scoring::electron_density::getDensityMap();

	numeric::xyzVector< Real > reference_atm;
	utility::vector1< numeric::xyzVector< Real > > atmList;
	utility::vector1< Real > all_K, all_C;

	numeric::xyzVector< Real > massSum(0.0,0.0,0.0), centerCA(0.0,0.0,0.0);

	// atom mask ... 3sigma from carbon
	Real const ATOM_MASK = 3.0 * sqrt( density.getEffectiveBfactor() / (2*M_PI*M_PI) );

	for ( Size i=1; i<= pose.size(); ++i ) {
		conformation::Residue const & rsd( pose.residue(i) );
		if ( rsd.aa() == core::chemical::aa_vrt ) continue;
		for ( Size j=1; j<= rsd.nheavyatoms(); ++j ) {
			conformation::Atom const & atom( rsd.atom(j) );
			atmList.push_back(atom.xyz());
			massSum += atom.xyz();

			if ( i==(pose.size()+1)/2 && j==2 ) {
				centerCA = atom.xyz();
			}

			chemical::AtomTypeSet const & atom_type_set( rsd.atom_type_set() );
			std::string const elt_i = atom_type_set[ rsd.atom_type_index( j ) ].element();
			scoring::electron_density::OneGaussianScattering sig_j = core::scoring::electron_density::get_A( elt_i );

			core::Real const K_i = sig_j.k( density.getEffectiveBfactor() );
			all_K.push_back( K_i );
			all_C.push_back( sig_j.C( K_i ) );
		}
	}
	int nAtms=atmList.size();
	massSum /= nAtms; // center_of_mass = mass_sum / nAtms

	if ( params.center_on_middle_ca_ ) {
		massSum = centerCA;
	}

	// precompute sines & cosines
	utility::vector1<core::Real> cT,cG, sT,sG;
	cT.resize(2*params.B_); cG.resize(2*params.B_); sT.resize(2*params.B_); sG.resize(2*params.B_);
	for ( core::Size t=1; t<=2*params.B_; ++t ) {
		core::Real const theta = (2.0*t-1.0)*M_PI/(4*params.B_);
		sT[t] = sin(theta);
		cT[t] = cos(theta);
	}
	for ( core::Size p=1; p<=2*params.B_; ++p ) {
		core::Real const phi = (2.0*p-2.0)*M_PI/(2*params.B_);
		sG[p] = sin(phi);
		cG[p] = cos(phi);
	}

	//////////////////
	// pose -> spherical-sampled density
	// 1. one models each atom with a Gaussian sphere of density
	// 2. interpolate this calculated density in cencentric spherical shells
	// (extending out to D Ang in 1 Ang steps)
	//////////////////
	sigR.dimension( 2*params.B_, 2*params.B_, params.nRsteps_ );
	sigR = 0.0;

	epsR.dimension( 2*params.B_, 2*params.B_, params.nRsteps_ );
	epsR = 0.0;

	// for each atom
	for ( int i=1; i<=nAtms; ++i ) {
		core::Real k=all_K[i];
		core::Real C=all_C[i];

		atmList[i] -= massSum;

		core::Real const atomR = atmList[i].length();
		if ( atomR < 1e-5 ) {
			// uniform contribution to inner shells
			for ( Size ridx=1; ridx<=params.nRsteps_; ++ridx ) {
				core::Real const atomD = ridx * params.delRsteps_;
				if ( atomD < ATOM_MASK ) {
					core::Real const atomH = C * exp(-k*atomD*atomD); // <-- this is the place to calculate density
					for ( core::Size t=1; t<=2*params.B_; ++t ) {
						for ( core::Size p=1; p<=2*params.B_; ++p ) {
							sigR(p,t,ridx) += atomH;
							epsR(p,t,ridx) = 1.0;
						}
					}
				}
			}
			continue;
		}

		core::Real const beta = acos( atmList[i][2] / atomR );
		core::Real const gamma = atan2( atmList[i][0] , atmList[i][1] );   // x and y switched from usual convention

		core::Real const st1 = sin(beta);
		core::Real const sg1 = sin(gamma);
		core::Real const ct1 = cos(beta);
		core::Real const cg1 = cos(gamma);

		if ( params.laplacian_offset_ != 0 ) {
			TR << "Applying laplacian filter with offset of: " << params.laplacian_offset_ << " A" << std::endl;
			for ( core::Size ridx=1; ridx<=params.nRsteps_; ++ridx ) {
				core::Real const shellR = ridx * params.delRsteps_;
				for ( core::Size t=1; t<=2*params.B_; ++t ) {
					core::Real const minAtomD = atomR*atomR + shellR*shellR - 2*atomR*shellR*(st1*sT[t]+ct1*cT[t]);
					if ( minAtomD>ATOM_MASK*ATOM_MASK ) continue; // this just loops back to 2xB so we still get to sigR
					for ( core::Size p=1; p<=2*params.B_; ++p ) {
						core::Real const atomD = atomR*atomR + shellR*shellR - 2*atomR*shellR*(st1*sT[t]*(sg1*sG[p]+cg1*cG[p])+ct1*cT[t]);
						if ( atomD < ATOM_MASK*ATOM_MASK ) {
							core::Real const atomH = C * exp(-k*atomD);
							sigR(p,t,ridx) += (-6 * atomH);
							epsR(p,t,ridx) = 1.0;
						}
					}
				}
			}
			// compute laplacian for surrounding coordinates
			for ( int xyz = 0; xyz < 3; ++xyz ) {
				for ( int lapl = 0; lapl < 2; ++lapl ) {
					reference_atm = atmList[i];
					atmList[i][xyz] = atmList[i][xyz] + ( ( (lapl==0) ? 1.0 : -1.0 ) * params.laplacian_offset_ );
					core::Real const lap_atomR = atmList[i].length();
					core::Real const beta_2 = acos( atmList[i][2] / lap_atomR );
					core::Real const gamma_2 = atan2( atmList[i][0] , atmList[i][1] );   // x and y switched from usual convention
					core::Real const st1_2 = sin(beta_2);
					core::Real const sg1_2 = sin(gamma_2);
					core::Real const ct1_2 = cos(beta_2);
					core::Real const cg1_2 = cos(gamma_2);
					// residue index
					for ( core::Size ridx=1; ridx<=params.nRsteps_; ++ridx ) {
						core::Real shellR = ridx * params.delRsteps_;
						for ( core::Size t=1; t<=2*params.B_; ++t ) {
							core::Real minAtomD =  lap_atomR*lap_atomR + shellR*shellR - 2*lap_atomR*shellR*(st1_2*sT[t]+ct1_2*cT[t]);
							if ( minAtomD>ATOM_MASK*ATOM_MASK ) continue;
							for ( core::Size p=1; p<=2*params.B_; ++p ) {
								core::Real const atomD = lap_atomR*lap_atomR + shellR*shellR - 2*lap_atomR*shellR*(st1_2*sT[t]*(sg1_2*sG[p]+cg1_2*cG[p])+ct1_2*cT[t]);
								if ( atomD < ATOM_MASK*ATOM_MASK ) {
									core::Real const atomH = C * exp(-k*atomD);
									sigR(p,t,ridx) += atomH;
									epsR(p,t,ridx) = 1.0;
								}
							}
						}
					}
					// set atm back to original value
					atmList[i] = reference_atm;
				}
			}
		} else {
			for ( core::Size ridx=1; ridx<=params.nRsteps_; ++ridx ) {
				core::Real shellR = ridx * params.delRsteps_;
				for ( core::Size t=1; t<=2*params.B_; ++t ) {
					core::Real minAtomD =  atomR*atomR + shellR*shellR - 2*atomR*shellR*(st1*sT[t]+ct1*cT[t]);
					if ( minAtomD>ATOM_MASK*ATOM_MASK ) continue;
					for ( core::Size p=1; p<=2*params.B_; ++p ) {
						core::Real atomD = atomR*atomR + shellR*shellR - 2*atomR*shellR*(st1*sT[t]*(sg1*sG[p]+cg1*cG[p])+ct1*cT[t]);
						if ( atomD < ATOM_MASK*ATOM_MASK ) {
							core::Real atomH = C * exp(-k*atomD);
							sigR(p,t,ridx) += atomH;
							epsR(p,t,ridx) = 1.0;
						}
					}
				}
			}
		}
	} // loop through each atom
}

} // electron_density
} // protocols
