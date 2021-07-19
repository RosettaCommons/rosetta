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

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>

#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>

#include <core/conformation/Residue.hh>
#include <core/pose/util.hh>

#include <core/scoring/electron_density/ElectronDensity.hh>
#include <core/scoring/electron_density/xray_scattering.hh>
#include <core/scoring/ScoreFunction.hh>

#include <core/id/AtomID.hh>

#include <core/io/pdb/pdb_writer.hh>
#include <core/io/StructFileRep.hh>
#include <core/io/StructFileRepOptions.hh>

#include <core/types.hh>

#include <utility/io/ozstream.hh>
#include <utility/vector1.hh>

#include <numeric/fourier/SHT.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzMatrix.hh>

#include <basic/Tracer.hh>


static basic::Tracer TR( "protocols.electron_density.DockIntoDensityUtils" );

namespace protocols {
namespace electron_density {

// non-superposed RMS
core::Real
get_rms(core::pose::PoseOP const r1, core::pose::PoseOP const r2, DensitySymmInfo const &d) {
	runtime_assert( r1->size() == r2->size() );
	core::Size nres = r1->size();
	core::Real rms=0.0;
	core::Size N=0;
	for ( int i=1; i<=(int)nres; ++i ) {
		if ( !r1->residue(i).is_protein() ) continue;
		rms += d.min_symm_dist2( r1->residue(i).xyz(2), r2->residue(i).xyz(2) );
		N++;
	}

	return (std::sqrt(rms/N));
}

// non-superposed RMS
core::Real
get_rms(RefinementResult const & r1, RefinementResult const & r2, DensitySymmInfo const & d ) {
	return get_rms(r1.pose_, r2.pose_, d);
}


/// CLASSES
numeric::xyzVector< core::Real >
RefinementResult::center() {
	core::Size const centerres = (pose_->size()+1)/2;
	return (pose_->residue(centerres).xyz(2));
}


/// FUNCTIONS

// move the pose given rotation/translation matrix
void
apply_transform(
	core::pose::Pose & pose,
	RBfitResult const& transform
) {
	// then each atom x can be transformed to this optimal configuration by:
	utility::vector1< core::id::AtomID > ids;
	utility::vector1< numeric::xyzVector< core::Real > > positions;
	numeric::xyzVector< core::Real > xyz_rot;

	for ( core::Size irsd = 1; irsd <= pose.size(); ++irsd ) {
		for ( core::Size iatom = 1; iatom <= pose.residue_type( irsd ).natoms(); ++iatom ) {
			numeric::xyzVector< core::Real > const atom_xyz = pose.xyz( core::id::AtomID( iatom, irsd ) );
			ids.push_back( core::id::AtomID( iatom, irsd ) );
			xyz_rot = transform.rotation_*( atom_xyz + transform.pre_trans_ ) + transform.post_trans_;
			positions.push_back( xyz_rot );
		}
	}
	pose.batch_set_xyz( ids, positions );
}


core::Real get_rot_angle( numeric::xyzMatrix<core::Real> const & R ) {
	core::Real const trace = R.xx() + R.yy() + R.zz();
	return acos( std::min( std::max(0.5* (trace-1),-1.0), 1.0) );
}


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


core::Real
get_radius(
	core::pose::Pose const & pose,
	numeric::xyzVector< core::Real > &com,
	bool const center_on_middle_ca ){
	numeric::xyzVector< core::Real > centerCA(0.0,0.0,0.0);
	com = numeric::xyzVector< core::Real >(0.0,0.0,0.0);
	core::Size nAtms=0;
	for ( core::Size i=1; i<= pose.size(); ++i ) {
		core::conformation::Residue const & rsd( pose.residue(i) );
		if ( rsd.aa() == core::chemical::aa_vrt ) continue;
		for ( core::Size j=1; j<= rsd.nheavyatoms(); ++j ) {
			core::conformation::Atom const & atom( rsd.atom(j) );
			com += atom.xyz();
			if ( i==(pose.size()+1)/2 && j==2 ) {
				centerCA = atom.xyz();
			}
			nAtms++;
		}
	}
	com /= nAtms;

	if ( center_on_middle_ca ) {
		com = centerCA;
	}

	core::Real maxSize = 0;
	for ( core::Size i=1; i<= pose.size(); ++i ) {
		core::conformation::Residue const & rsd( pose.residue(i) );
		if ( rsd.aa() == core::chemical::aa_vrt ) continue;
		for ( core::Size j=1; j<= rsd.nheavyatoms(); ++j ) {
			core::conformation::Atom const & atom( rsd.atom(j) );
			maxSize = std::max( maxSize, (com-atom.xyz()).length_squared() );
		}
	}
	maxSize = std::sqrt( maxSize ); // this the radius
	return maxSize;
}


core::Real
get_spectrum(
	core::pose::Pose const& pose,
	utility::vector1< core::Real > &pose_1dspec,
	core::Real const delR,
	core::Real const fragDens,
	bool const convolute_single_residue,
	bool const center_on_middle_ca
) {
	numeric::xyzVector< core::Real > com(0.0,0.0,0.0);
	core::Real const extent = get_radius( pose, com, center_on_middle_ca );

	// grid spacing == delR
	auto const ngrid = (core::Size) std::ceil( extent / delR + 2);
	pose_1dspec.clear();
	pose_1dspec.resize(ngrid, 0.0);
	utility::vector1< core::Real > pose_1dspec_for_nrsteps = pose_1dspec;

	core::Real massSum=0.0;
	// for fine point selection... generate pose_1dspec based on the middle residue of the pose.
	// The smoothing of the map is just based on whatever the extent of the middle residue.
	// It doesn't seem to matter but seems strange that it is dependent on that.
	if ( convolute_single_residue == true ) {
		core::Size const midres = [&]{
			core::Size midres = (pose.size()+1)/2;
			while ( pose.residue(midres).is_virtual_residue() ) {
				++midres;
			}
			return midres;
		}();
		core::conformation::Residue const & rsd( pose.residue(midres) );
		core::conformation::Atom const & residue_CA( rsd.atom(2) );
		for ( core::Size j=1; j<= rsd.nheavyatoms(); ++j ) {
			core::conformation::Atom const & atom( rsd.atom(j) );
			core::Real const binnum = ( atom.xyz()-residue_CA.xyz() ).length() / delR + 1.0;
			core::Real const fpart = binnum - std::floor(binnum);
			auto const binint = (core::Size) std::floor(binnum);
			pose_1dspec[binint] += (1-fpart);
			pose_1dspec[binint+1] += (fpart);
		}
	} else {
		// for coarse point selection... generate pose_1dspec based on whole pose
		for ( core::Size i=1; i<=pose.size(); ++i ) {
			core::conformation::Residue const & rsd( pose.residue(i) );
			if ( rsd.aa() == core::chemical::aa_vrt ) continue;
			for ( core::Size j=1; j<= rsd.nheavyatoms(); ++j ) {
				core::conformation::Atom const & atom( rsd.atom(j) );
				core::Real const binnum = ( atom.xyz()-com ).length() / delR + 1.0;
				core::Real const fpart = binnum - std::floor(binnum);
				auto const binint = (core::Size) std::floor(binnum);
				pose_1dspec[binint] += (1-fpart);
				pose_1dspec[binint+1] += (fpart);
			}
		}
	}

	// this is to calculate massSum via full pose for nRsteps_
	for ( core::Size i=1; i<=pose.size(); ++i ) {
		core::conformation::Residue const & rsd( pose.residue(i) );
		if ( rsd.aa() == core::chemical::aa_vrt ) continue;
		for ( core::Size j=1; j<= rsd.nheavyatoms(); ++j ) {
			core::conformation::Atom const & atom( rsd.atom(j) );
			core::Real const binnum = ( atom.xyz()-com ).length() / delR + 1.0;
			core::Real const fpart = binnum - std::floor(binnum);
			auto const binint = (core::Size) std::floor(binnum);
			pose_1dspec_for_nrsteps[binint] += (1-fpart);
			pose_1dspec_for_nrsteps[binint+1] += (fpart);
			massSum += 1;
		}
	}

	// now setting nRsteps_
	core::Real running_total=0.0;
	core::Size const nRsteps = [&] {
		core::Size nRsteps = ngrid;
		for ( int i=1; i<=(int)ngrid; ++i ) {
			running_total += pose_1dspec_for_nrsteps[i] / massSum;
			if ( running_total > fragDens ) {
				nRsteps = i;
				break;
			}
		}
		return nRsteps;
	}();

	running_total=0.0;
	for ( int i=1; i<=(int)ngrid; ++i ) {
		running_total += pose_1dspec[i] / massSum;
		pose_1dspec[i] /= (i*i); // normalize
		TR << "spectrum " << i << ": " << pose_1dspec[i] << " " << pose_1dspec[i]*i*i/massSum << " " << running_total << std::endl;
	}
	return nRsteps;
}



void
map_from_spectrum(
	utility::vector1< core::Real > const& pose_1dspec,
	ObjexxFCL::FArray3D< core::Real > &rot,
	core::Real const delR ) {
	// ugly
	for ( int z=1; z<=(int)rot.u3(); ++z ) {
		for ( int y=1; y<=(int)rot.u2(); ++y ) {
			for ( int x=1; x<=(int)rot.u1(); ++x ) {
				numeric::xyzVector< core::Real > idxX, cartX;
				idxX[0] = x<=rot.u1()/2 ? x-1 : x-1-rot.u1();
				idxX[1] = y<=rot.u2()/2 ? y-1 : y-1-rot.u2();
				idxX[2] = z<=rot.u3()/2 ? z-1 : z-1-rot.u3();

				core::scoring::electron_density::getDensityMap().idxoffset2cart( idxX, cartX );
				core::Real const d = cartX.length() / delR;
				core::Real const fpart = d - std::floor(d);
				core::Size const dint = (core::Size) std::floor(d) + 1;
				if ( dint<pose_1dspec.size() ) { // last entry is always 0 so this check is valid
					rot(x,y,z) = (1-fpart)*pose_1dspec[dint] + (fpart)*pose_1dspec[dint+1];
				} else {
					rot(x,y,z) = 0.0;
				}
			}
		}
	}
}




utility::vector1< std::pair< numeric::xyzVector< core::Real >, core::Real > >
create_and_sort_point_score_pairs(
	ObjexxFCL::FArray3D< float > const & densdata,
	ObjexxFCL::FArray3D< core::Real > const & rot,
	core::Size const grid_step) {
	// sort points
	utility::vector1< std::pair< numeric::xyzVector< core::Real >, core::Real > > point_score_pairs;
	for ( int z=1; z<=(int)densdata.u3(); z+=grid_step ) {
		for ( int y=1; y<=(int)densdata.u2(); y+=grid_step ) {
			for ( int x=1; x<=(int)densdata.u1(); x+=grid_step ) {
				numeric::xyzVector< core::Real > const x_idx(x,y,z);
				core::Real const dens_value = -rot(x,y,z);
				std::pair< numeric::xyzVector< core::Real >, core::Real > point_score_pair = std::make_pair(x_idx, dens_value);
				point_score_pairs.push_back(point_score_pair);
			}
		}
	}
	std::sort(point_score_pairs.begin(), point_score_pairs.end(), PointScoreComparator());
	return point_score_pairs;
}


void
dump_points_to_search_to_pdb(
	utility::vector1< numeric::xyzVector< core::Real > > const & points_to_search,
	std::string const & filename ) {
	//dump the points
	core::io::StructFileRep sfr;
	auto & chains = sfr.chains();
	chains.push_back(utility::vector0<core::io::AtomInformation>());
	core::io::StructFileRepOptionsCOP options = utility::pointer::make_shared< core::io::StructFileRepOptions >();
	core::io::AtomInformation ai;
	for ( core::Size i=1; i<=points_to_search.size(); i++ ) {
		ai.serial = i;
		ai.name = "MG  ";
		ai.resName = " MG";
		ai.chainID = 'A';
		numeric::xyzVector< core::Real > x_cart;
		numeric::xyzVector< core::Real > x_idx = points_to_search[i];
		core::scoring::electron_density::getDensityMap().idx2cart( x_idx, x_cart );
		ai.x = x_cart[0];
		ai.y = x_cart[1];
		ai.z = x_cart[2];
		ai.occupancy = 1.0;
		ai.element = "MG";
		chains[0].push_back(ai);
	}

	std::string const pdb_contents(core::io::pdb::create_pdb_contents_from_sfr(sfr, options));
	utility::io::ozstream file(std::string(filename).c_str(), std::ios::out | std::ios::binary);
	file.write( pdb_contents.c_str(), pdb_contents.size() );
	file.close();
}

utility::vector1< numeric::xyzVector< core::Real > >
select_density_points( core::pose::Pose const & pose,
	SelectDensityPointsOptions const & params,
	core::Size & nRsteps ) {

	utility::vector1< numeric::xyzVector< core::Real > > points_to_search;

	ObjexxFCL::FArray3D< float > const & densdata = core::scoring::electron_density::getDensityMap().get_data();
	ObjexxFCL::FArray3D< std::complex<core::Real> > Fdens, Frot;
	numeric::fourier::fft3(densdata, Fdens);

	// make rotationally averaged pose map
	utility::vector1< core::Real > pose_1dspec;

	nRsteps = get_spectrum( pose, pose_1dspec, params.delRsteps_, params.fragDens_, params.convolute_single_residue_, params.center_on_middle_ca_ );

	ObjexxFCL::FArray3D< core::Real > rot;
	rot.dimension( densdata.u1(), densdata.u2(), densdata.u3() );

	map_from_spectrum( pose_1dspec, rot, params.delRsteps_);
	if ( basic::options::option[ basic::options::OptionKeys::edensity::debug ]() ) {
		core::scoring::electron_density::ElectronDensity(rot,1.0, numeric::xyzVector< core::Real >(0,0,0), false).writeMRC( "spectrum.mrc" );
	}

	numeric::fourier::fft3(rot, Frot);

	// FFT convolution
	core::Size const npts = densdata.u1()*densdata.u2()*densdata.u3();
	for ( core::Size i=0; i< npts; ++i ) {
		Frot[i] = Fdens[i] * std::conj( Frot[i] );  // reuse
	}
	numeric::fourier::ifft3(Frot, rot);

	// mask asu
	if ( params.symminfo_.enabled() ) {
		params.symminfo_.mask_asu( rot , core::scoring::electron_density::getDensityMap(), 0);
	}

	auto const point_score_pairs = create_and_sort_point_score_pairs(densdata, rot, params.gridStep_);

	core::Real minDistNative = 1e4;
	for ( core::Size i=1; i<=point_score_pairs.size(); i++ ) {
		bool hasneighbor = false;
		numeric::xyzVector< core::Real > x_idx = point_score_pairs[i].first;
		numeric::xyzVector< core::Real > x_cart(x_idx[0],x_idx[1],x_idx[2]);
		core::scoring::electron_density::getDensityMap().idx2cart( x_idx, x_cart );
		for ( core::Size j=1; j<=points_to_search.size(); j++ ) {
			numeric::xyzVector< core::Real > x_idx_stored = points_to_search[j];
			numeric::xyzVector< core::Real > x_cart_stored(x_idx_stored[0],x_idx_stored[1],x_idx_stored[2]);
			core::scoring::electron_density::getDensityMap().idx2cart( x_idx_stored, x_cart_stored );
			core::Real distance = (x_cart - x_cart_stored).length();
			if ( distance < params.point_radius_ ) hasneighbor = true;
		}
		if ( !hasneighbor ) {
			points_to_search.push_back(point_score_pairs[i].first);
			if ( params.native_ ) {
				numeric::xyzVector< core::Real > x_cart2;
				core::scoring::electron_density::getDensityMap().idx2cart( x_idx, x_cart2 );
				core::Real const distNative = params.symminfo_.min_symm_dist2(x_cart2, params.native_com_);
				minDistNative = std::min( minDistNative, distNative );
			}
		}
		if ( points_to_search.size() >= params.topNtrans_ ) break;
	}


	if ( basic::options::option[ basic::options::OptionKeys::edensity::debug ]() ) {
		core::scoring::electron_density::ElectronDensity mapdmp = core::scoring::electron_density::getDensityMap();
		mapdmp.set_data(rot);
		mapdmp.writeMRC( "filter.mrc" );

		dump_points_to_search_to_pdb( points_to_search, "selectedpoints.pdb");
	}

	if ( params.native_ ) TR << "Closest point to native: " << std::sqrt(minDistNative) << std::endl;
	return points_to_search;
}


// read results, normalize scores, remove redundant (CA rms)
void
do_filter( RefinementResultDB & results, DensitySymmInfo const & symminfo, core::Real const cluster_radius) {
	// dumb -- copy to vector
	utility::vector1< RefinementResult > results_sort;
	while ( results.size() > 0 )
			results_sort.push_back( results.pop() );
	utility::vector1< bool > selector(results_sort.size(), false);

	for ( int i=results_sort.size(); i>=1; --i ) {
		selector[i] = true;
		for ( int j=i+1; j<=(int)results_sort.size() && selector[i]; ++j ) {
			if ( selector[j] && get_rms(results_sort[i], results_sort[j], symminfo)<cluster_radius ) {
				selector[i] = false;
			}
		}
		if ( selector[i] ) {
			results.add_element( results_sort[i] );
		}
	}

	return;
}


// filter RB results
void
do_filter(
	utility::vector1< core::pose::PoseOP > const &poses,
	RBfitResultDB & results,
	bool const rescore,
	DensitySymmInfo const & symminfo,
	core::Real const cluster_radius
) {
	// dumb -- copy to vector
	utility::vector1< RBfitResult > results_sort;
	utility::vector1< core::pose::PoseOP > selected;
	while ( results.size() > 0 )
			results_sort.push_back( results.pop() );
	utility::vector1< bool > selector(results_sort.size(), false);

	core::scoring::ScoreFunctionOP scorefxn_dens( new core::scoring::ScoreFunction() );
	scorefxn_dens->set_weight( core::scoring::elec_dens_whole_structure_allatom, 1.0); //?

	for ( int i=results_sort.size(); i>=1; --i ) {
		selector[i] = true;
		// make (& rescore) the pose:
		core::pose::PoseOP const posecopy( utility::pointer::make_shared< core::pose::Pose >( *(poses[ results_sort[i].pose_idx_ ]) ) );
		apply_transform( *posecopy, results_sort[i] );

		if ( rescore ) {
			TR << "[" << i << "/" << results_sort.size() << "]" << "\r" << std::flush;

			core::pose::addVirtualResAsRoot( *posecopy );
			results_sort[i].score_ = (*scorefxn_dens)(*posecopy);
		}
		for ( int j=1; j<=(int)selected.size() && selector[i]; ++j ) {
			if ( get_rms(posecopy, selected[j], symminfo) <cluster_radius ) {
				selector[i] = false;
			}
		}

		if ( selector[i] ) {
			results.add_element( results_sort[i] );
			selected.push_back( posecopy );
		}
	}

	return;
}

// fast (rot-only) filter
void
do_filter(
	RBfitResultDB & results,
	core::Size const delR,
	core::Size const nRsteps,
	core::Real const cluster_radius) {
	// dumb -- copy to vector
	utility::vector1< RBfitResult > results_sort;
	while ( results.size() > 0 )
			results_sort.push_back( results.pop() );
	utility::vector1< bool > selector(results_sort.size(), false);

	core::Real nsel=0;
	for ( int i=results_sort.size(); i>=1; --i ) {
		selector[i] = true;

		for ( int j=i+1; j<=(int)selector.size() && selector[i]; ++j ) {
			core::Real est_rms = 0.6*delR*nRsteps * get_rot_angle( results_sort[i].rotation_ * numeric::inverse(results_sort[j].rotation_) );

			if ( selector[j] && est_rms < cluster_radius ) {
				selector[i] = false;
			}
		}

		if ( selector[i] ) {
			results.add_element( results_sort[i] );
			nsel++;
		}
	}

	return;
}


// do the main search over the map
void
density_grid_search(
	core::Size pose_idx,
	core::pose::Pose const & pose,
	RBfitResultDB & results,
	utility::vector1< numeric::xyzVector< core::Real > > const & points_to_search,
	DensityGridSearchOptions const & params
) {
	core::scoring::electron_density::ElectronDensity &density = core::scoring::electron_density::getDensityMap();

	// allocate space for SHT
	numeric::fourier::SHT SOFT(params.B_, params.nRsteps_);

	// get com of pose
	numeric::xyzMatrix<core::Real> rot;
	numeric::xyzVector<core::Real> pretrans, posttrans;
	get_radius( pose, pretrans, params.center_on_middle_ca_ );

	pretrans=-1.0*pretrans;

	// get pose SPHARM
	ObjexxFCL::FArray3D< core::Real > poseSig, poseCoefR, poseCoefI;
	ObjexxFCL::FArray3D< core::Real > poseEps, poseEpsCoefR, poseEpsCoefI;

	/// step 1: map pose to spherically sampled density + mask
	pose_spherical_samples( pose, poseSig, poseEps,
		PoseSphericalSamplesOptions{params.B_, params.nRsteps_, params.delRSteps_, params.laplacian_offset_, params.center_on_middle_ca_});

	SOFT.sharm_transform( poseSig, poseCoefR, poseCoefI );
	//poseEps = 1.0; // uncomment to turn off masking
	SOFT.sharm_transform( poseEps, poseEpsCoefR, poseEpsCoefI );

	if ( basic::options::option[ basic::options::OptionKeys::edensity::debug ]() ) {
		core::scoring::electron_density::ElectronDensity(poseSig, 1.0, numeric::xyzVector< core::Real >(0,0,0), false).writeMRC( "Pose_sigR.mrc" );
		core::scoring::electron_density::ElectronDensity(poseEps, 1.0, numeric::xyzVector< core::Real >(0,0,0), false).writeMRC( "Pose_epsR.mrc" );
	}

	SOFT.sharm_invTransform( poseSig, poseCoefR, poseCoefI );        // generate bandlimited signal
	SOFT.sharm_invTransform( poseEps, poseEpsCoefR, poseEpsCoefI );  // generate bandlimited signal

	// *** OPTIONAL 1 ***
	// normalize poseEps to [0,1]
	/*
	core::Real minE = 1e6, maxE = -1e6;
	for (core::Size j=0; j<4*B_*B_*nRsteps_; ++j) {
	minE = std::min( minE, poseEps[j] );
	maxE = std::max( maxE, poseEps[j] );
	}
	core::Real scale = 1.0/(maxE-minE);
	for (core::Size j=0; j<4*B_*B_*nRsteps_; ++j) {
	poseEps[j] = (poseEps[j]-minE) * scale;
	}
	for (core::Size j=0; j<nRsteps_; ++j) {
	poseCoefR[B_*B_*j] = poseCoefR[B_*B_*j] - minE;
	}
	for (core::Size j=0; j<B_*B_*nRsteps_; ++j) {
	poseCoefR[j] = (poseCoefR[j]) * scale;
	poseCoefI[j] = (poseCoefI[j]) * scale;
	}
	*/

	double sumRWt = 0.0, sumTwt = 0.0;
	for ( core::Size r=0; r<params.nRsteps_; ++r ) {
		sumRWt += 4.0 * numeric::square((r+1.0));
	}
	for ( core::Size t=0; t<2*params.B_; ++t ) {
		sumTwt += std::sin( (t+0.5) * numeric::constants::d::pi / (2*params.B_) );
	}

	core::Real sumEps=0.0, sumSig=0.0, sumSig2=0.0;
	for ( core::Size r=0; r<params.nRsteps_; ++r ) {
		double const thisRWt = 4.0 * numeric::square((r+1.0)) / sumRWt;
		for ( core::Size t=0; t<2*params.B_; ++t ) {
			double const thisTWt = std::sin( (t+0.5) * numeric::constants::d::pi / (2*params.B_) ) / sumTwt;
			for ( core::Size p=0; p<2*params.B_; ++p ) {
				double const thisWt = thisRWt*thisTWt/(2*params.B_);
				sumEps += thisWt*poseEps(p+1,t+1,r+1);
				sumSig += thisWt*poseEps(p+1,t+1,r+1)*poseSig(p+1,t+1,r+1);
				sumSig2 += thisWt*poseEps(p+1,t+1,r+1)*poseSig(p+1,t+1,r+1)*poseSig(p+1,t+1,r+1);
			}
		}
	}

	if ( basic::options::option[ basic::options::OptionKeys::edensity::debug ]() ) {
		core::scoring::electron_density::ElectronDensity(poseSig, 1.0, numeric::xyzVector< core::Real >(0,0,0), false).writeMRC( "Pose_sigRbl.mrc" );
		core::scoring::electron_density::ElectronDensity(poseEps, 1.0, numeric::xyzVector< core::Real >(0,0,0), false).writeMRC( "Pose_epsRbl.mrc" );
	}

	ObjexxFCL::FArray3D< core::Real > mapSig, mapCoefR, mapCoefI, map2Sig, map2CoefR, map2CoefI;
	ObjexxFCL::FArray3D< core::Real > so3_correl, mask_correl, mask2_correl;
	for ( core::Size i=1; i<=points_to_search.size(); ++i ) {
		// get cartesian coords of ths point
		density.idx2cart( points_to_search[i], posttrans );

		density.mapSphericalSamples( mapSig, params.nRsteps_, params.delRSteps_, params.B_, points_to_search[i], params.laplacian_offset_ );
		map2Sig = mapSig;
		for ( core::Size j=0; j<4*params.B_*params.B_*params.nRsteps_; ++j ) {
			map2Sig[j] *= mapSig[j];
		}
		if ( basic::options::option[ basic::options::OptionKeys::edensity::debug ]() ) {
			core::scoring::electron_density::ElectronDensity(
				mapSig, 1.0, numeric::xyzVector< core::Real >(0,0,0), false).writeMRC("Map"+utility::to_string(i)+"_sigR.mrc" );
		}

		// new version
		SOFT.sharm_transform( mapSig, mapCoefR, mapCoefI );
		SOFT.sharm_transform( map2Sig, map2CoefR, map2CoefI );

		SOFT.so3_correlate(so3_correl, mapCoefR,mapCoefI,  poseCoefR,poseCoefI);
		SOFT.so3_correlate(mask_correl, mapCoefR,mapCoefI,  poseEpsCoefR,poseEpsCoefI);
		SOFT.so3_correlate(mask2_correl, map2CoefR,map2CoefI,  poseEpsCoefR,poseEpsCoefI);

		// we initially oversample since the set is so clustered
		// no matter how many we want, don't take more than 1/8 of everything (which still might be a lot)
		core::Size const nperRot = std::min( 100*params.max_rot_per_trans_ , params.B_*params.B_*params.B_);
		core::Size const nperRotCl = params.max_rot_per_trans_;

		RBfitResultDB local_results( nperRot );

		for ( core::Size j=0; j<8*params.B_*params.B_*params.B_; ++j ) {
			core::Real const sumSigMap = so3_correl[j] / (4*numeric::constants::d::pi);
			core::Real const sumMap = mask_correl[j] / (4*numeric::constants::d::pi);
			core::Real const sumMap2 = mask2_correl[j] / (4*numeric::constants::d::pi);

			double const CC = (
				(sumEps*sumSigMap - sumSig*sumMap) / (
				std::sqrt( sumEps*sumSig2 - sumSig*sumSig )
				* std::sqrt( sumEps*sumMap2 - sumMap*sumMap )
				) );

			if ( local_results.to_add_element( CC ) ) {
				SOFT.idx_to_rot(j , rot);
				local_results.add_element( RBfitResult( pose_idx, CC, rot, pretrans, posttrans ) );
			}
		}

		do_filter( local_results, params.delRSteps_, params.nRsteps_, params.cluster_radius_ );
		while ( local_results.size() > nperRotCl ) local_results.pop();

		core::Size const nclust = local_results.size();

		core::Real bestrms=1e6,bestscore=0;
		while ( local_results.size() > 0 ) {
			RBfitResult sol_i = local_results.pop();

			if ( params.native_ ) {
				core::pose::PoseOP const posecopy( utility::pointer::make_shared< core::pose::Pose >( pose ) );
				apply_transform( *posecopy, sol_i );
				core::pose::addVirtualResAsRoot( *posecopy );
				core::Real const rms_i = get_rms(params.native_, posecopy, params.symminfo_);
				if ( rms_i < bestrms ) {
					bestrms = rms_i;
					bestscore = sol_i.score_;
				}
			}

			results.add_element( sol_i );
		}

		core::Real minDistNative=1e6;
		if ( params.native_ ) {
			core::Real const distNative = (posttrans-params.native_com_).length_squared();
			minDistNative = std::min( minDistNative, distNative );
		}

		if ( std::sqrt(minDistNative) <5.0 ) {
			TR << "[" << i << "/" << points_to_search.size() << "]" << " nmdls " << nclust << " pointrms " << std::sqrt(minDistNative)
				<< " rms " << bestrms << " score " << bestscore << std::endl;
		}
		if ( i%100 == 0 ) {
			TR << "[" << i << "/" << points_to_search.size() << "] " << results.top().score_ << " / " << results.size() << std::endl;
		}
	}

	//TR << "[" << points_to_search_.size() << "/" << points_to_search_.size() << "]" << std::endl;
}


} // electron_density
} // protocols
