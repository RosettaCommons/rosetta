// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief protocols for folding into density
/// @details
/// @author Frank DiMaio

#include <protocols/electron_density/DockIntoDensityMover.hh>

#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>

#include <core/chemical/AA.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/io/AtomInformation.hh>
#include <core/io/StructFileRep.hh>
#include <core/io/StructFileRepOptions.fwd.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/util.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/types.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/fourier/SHT.hh>


#include <protocols/electron_density/DockIntoDensityMover.hh>
#include <protocols/jd2/util.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

#include <utility/file/file_sys_util.hh>
#include <utility/file/FileName.hh>
#include <utility/io/ozstream.hh>
#include <utility/string_util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


#include <basic/Tracer.hh>

#include <core/chemical/AtomTypeSet.hh> // AUTO IWYU For AtomTypeSet
#include <core/io/StructFileRepOptions.hh> // AUTO IWYU For StructFileRepOptions
#include <core/scoring/electron_density/xray_scattering.hh> // AUTO IWYU For OneGaussianScattering, get_A




namespace protocols {
namespace electron_density {

using core::Size;

static basic::Tracer TR( "protocols.electron_density.DockIntoDensityMover" );

struct
	ScoredPoint {
	core::Real score;
	numeric::xyzVector< core::Real > pos;
};

struct
	ScoredPointCmptr {
	bool operator()(ScoredPoint a, ScoredPoint b) {
		return a.score > b.score;
	}
};
// non-superposed RMS
core::Real
get_rms(core::pose::PoseOP r1, core::pose::PoseOP r2, DensitySymmInfo const &d) {
	runtime_assert( r1->size() == r2->size() );
	core::Size nres = r1->size();
	core::Real rms=0.0;
	core::Size N=0;
	for ( int i=1; i<=(int)nres; ++i ) {
		if ( !r1->residue(i).is_protein() ) continue;
		//rms += (r1->residue(i).xyz(2) - r2->residue(i).xyz(2)).length_squared();
		rms += d.min_symm_dist2( r1->residue(i).xyz(2), r2->residue(i).xyz(2) );
		N++;
	}

	return (std::sqrt(rms/N));
}

// non-superposed RMS
core::Real
get_rms(RefinementResult r1, RefinementResult r2, DensitySymmInfo const &d) {
	return get_rms(r1.pose_, r2.pose_, d);
}

core::Real get_rot_angle( numeric::xyzMatrix<core::Real> R ) {
	core::Real trace = R.xx() + R.yy() + R.zz();
	return acos( std::min( std::max(0.5* (trace-1),-1.0), 1.0) );
}


void
DockIntoDensityMover::print_best_rms( core::pose::Pose const &pose, RBfitResultDB const &results ) {
	core::Real bestrms=1e4;
	core::Size bestrank=0;
	RBfitResultDB resultscopy = results;
	while ( resultscopy.size() > 0 ) {
		RBfitResult sol_i = resultscopy.pop();
		core::pose::PoseOP posecopy ( new core::pose::Pose(pose) );
		apply_transform( *posecopy, sol_i );
		core::pose::addVirtualResAsRoot( *posecopy );
		core::Real rms_i = get_rms(native_, posecopy, symminfo_);
		if ( rms_i < bestrms ) {
			bestrms = rms_i; bestrank=resultscopy.size()+1;
		}
	}
	TR << "Best RMS = " << bestrms << " at rank " << bestrank << std::endl;
}



void
DockIntoDensityMover::setNative( core::pose::PoseOP native ) {
	native_ = native;
	core::pose::addVirtualResAsRoot( *native_ );

	native_com_ = numeric::xyzVector< core::Real >(0,0,0);
	core::Size N=0;
	for ( int i=1; i<=(int)native_->size(); ++i ) {
		if ( !native_->residue(i).is_protein() ) continue;
		native_com_ += native_->residue(i).xyz(2);
		N++;
	}
	native_com_ /= N;
}


// get_radius: calculate the extent of a pose
core::Real
DockIntoDensityMover::get_radius( core::pose::Pose const & pose, numeric::xyzVector< core::Real > &com ){
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

	if ( center_on_middle_ca_ ) {
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


// move the pose given rotation/translation matrix
void
DockIntoDensityMover::apply_transform (
	core::pose::Pose & pose,
	RBfitResult const& transform
) {
	// then each atom x can be transformed to this optimal configuration by:
	utility::vector1< core::id::AtomID > ids;
	utility::vector1< numeric::xyzVector< core::Real > > positions;
	numeric::xyzVector< core::Real > xyz_rot;

	for ( core::Size irsd = 1; irsd <= pose.size(); ++irsd ) {
		for ( core::Size iatom = 1; iatom <= pose.residue_type( irsd ).natoms(); ++iatom ) {
			numeric::xyzVector< core::Real > atom_xyz = pose.xyz( core::id::AtomID( iatom, irsd ) );
			ids.push_back( core::id::AtomID( iatom, irsd ) );
			xyz_rot = transform.rotation_*( atom_xyz + transform.pre_trans_ ) + transform.post_trans_;
			positions.push_back( xyz_rot );
		}
	}
	pose.batch_set_xyz( ids, positions );
}


void
DockIntoDensityMover::get_spectrum( core::pose::Pose const& pose, utility::vector1< core::Real > &pose_1dspec ) {
	numeric::xyzVector< core::Real > com(0.0,0.0,0.0);
	core::Real extent = get_radius( pose, com );

	// grid spacing == delR
	auto ngrid = (core::Size) std::ceil( extent / delR_ + 2);
	pose_1dspec.clear();
	pose_1dspec.resize(ngrid, 0.0);
	utility::vector1< core::Real > pose_1dspec_for_nrsteps = pose_1dspec;

	core::Real massSum=0.0;
	// for fine point selection... generate pose_1dspec based on the middle residue of the pose.
	if ( convolute_single_residue_ == true ) {
		core::Size midres = (pose.size()+1)/2;
		while ( pose.residue(midres).is_virtual_residue() ) {
			++midres;
		}
		core::conformation::Residue const & rsd( pose.residue(midres) );
		core::conformation::Atom const & residue_CA( rsd.atom(2) );
		for ( core::Size j=1; j<= rsd.nheavyatoms(); ++j ) {
			core::conformation::Atom const & atom( rsd.atom(j) );
			core::Real binnum = ( atom.xyz()-residue_CA.xyz() ).length() / delR_ + 1.0;
			core::Real fpart = binnum - std::floor(binnum);
			auto binint = (core::Size) std::floor(binnum);
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
				core::Real binnum = ( atom.xyz()-com ).length() / delR_ + 1.0;
				core::Real fpart = binnum - std::floor(binnum);
				auto binint = (core::Size) std::floor(binnum);
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
			core::Real binnum = ( atom.xyz()-com ).length() / delR_ + 1.0;
			core::Real fpart = binnum - std::floor(binnum);
			auto binint = (core::Size) std::floor(binnum);
			pose_1dspec_for_nrsteps[binint] += (1-fpart);
			pose_1dspec_for_nrsteps[binint+1] += (fpart);
			massSum += 1;
		}
	}

	// now setting nRsteps_
	core::Real const fracDens=fragDens_; // choose radius covering this fraction of density mass
	core::Real running_total=0.0;
	nRsteps_ = ngrid;
	for ( int i=1; i<=(int)ngrid; ++i ) {
		running_total += pose_1dspec_for_nrsteps[i] / massSum;
		if ( running_total > fracDens ) {
			nRsteps_ = i;
			break;
		}
	}

	running_total=0.0;
	for ( int i=1; i<=(int)ngrid; ++i ) {
		running_total += pose_1dspec[i] / massSum;
		pose_1dspec[i] /= (i*i); // normalize
		TR << "spectrum " << i << ": " << pose_1dspec[i] << " " << pose_1dspec[i]*i*i/massSum << " " << running_total << std::endl;
	}
}


void
DockIntoDensityMover::map_from_spectrum( utility::vector1< core::Real > const& pose_1dspec, ObjexxFCL::FArray3D< core::Real > &rot ) {
	// ugly
	for ( int z=1; z<=(int)rot.u3(); ++z ) {
		for ( int y=1; y<=(int)rot.u2(); ++y ) {
			for ( int x=1; x<=(int)rot.u1(); ++x ) {
				numeric::xyzVector< core::Real > idxX, cartX;
				idxX[0] = x<=rot.u1()/2 ? x-1 : x-1-rot.u1();
				idxX[1] = y<=rot.u2()/2 ? y-1 : y-1-rot.u2();
				idxX[2] = z<=rot.u3()/2 ? z-1 : z-1-rot.u3();

				core::scoring::electron_density::getDensityMap().idxoffset2cart( idxX, cartX );
				core::Real d = cartX.length() / delR_;
				core::Real fpart = d - std::floor(d);
				core::Size dint = (core::Size) std::floor(d) + 1;
				if ( dint<pose_1dspec.size() ) { // last entry is always 0 so this check is valid
					rot(x,y,z) = (1-fpart)*pose_1dspec[dint] + (fpart)*pose_1dspec[dint+1];
				} else {
					rot(x,y,z) = 0.0;
				}
			}
		}
	}

	if ( basic::options::option[ basic::options::OptionKeys::edensity::debug ]() ) {
		core::scoring::electron_density::ElectronDensity(rot,1.0, numeric::xyzVector< core::Real >(0,0,0), false).writeMRC( "spectrum.mrc" );
	}
}


void
DockIntoDensityMover::predefine_search( utility::vector1< numeric::xyzVector<core::Real> > &pts_in ) {
	points_defined_ = true;

	points_to_search_.clear();

	for ( int i=1; i<=(int)pts_in.size(); ++i ) {
		numeric::xyzVector< core::Real > x_idx;
		core::scoring::electron_density::getDensityMap().cart2idx( pts_in[i] , x_idx );
		points_to_search_.push_back( x_idx );
	}
}


void
DockIntoDensityMover::select_points( core::pose::Pose & pose ) {
	ObjexxFCL::FArray3D< float > const & densdata = core::scoring::electron_density::getDensityMap().get_data();
	ObjexxFCL::FArray3D< std::complex<core::Real> > Fdens, Frot;
	numeric::fourier::fft3(densdata, Fdens);

	// make rotationally averaged pose map
	utility::vector1< core::Real > pose_1dspec;
	get_spectrum( pose, pose_1dspec );

	if ( points_defined_ ) return; // points were predefined, don't change
	points_to_search_.clear();

	ObjexxFCL::FArray3D< core::Real > rot;
	rot.dimension( densdata.u1(), densdata.u2(), densdata.u3() );
	map_from_spectrum( pose_1dspec, rot );
	numeric::fourier::fft3(rot, Frot);

	// FFT convolution
	core::Size npts = densdata.u1()*densdata.u2()*densdata.u3();
	for ( core::Size i=0; i< npts; ++i ) {
		Frot[i] = Fdens[i] * std::conj( Frot[i] );  // reuse
	}
	numeric::fourier::ifft3(Frot, rot);

	// mask asu
	if ( symminfo_.enabled() ) {
		symminfo_.mask_asu( rot , core::scoring::electron_density::getDensityMap(), 0);
	}

	// sort points
	utility::vector1< std::pair< numeric::xyzVector< core::Real >, core::Real > > point_score_pairs;
	for ( int z=1; z<=(int)densdata.u3(); z+=gridStep_ ) {
		for ( int y=1; y<=(int)densdata.u2(); y+=gridStep_ ) {
			for ( int x=1; x<=(int)densdata.u1(); x+=gridStep_ ) {
				numeric::xyzVector< core::Real > x_idx(x,y,z);
				core::Real dens_value = -rot(x,y,z);
				std::pair< numeric::xyzVector< core::Real >, core::Real > point_score_pair = std::make_pair(x_idx, dens_value);
				point_score_pairs.push_back(point_score_pair);
			}
		}
	}
	core::Real minDistNative = 1e4;
	std::sort(point_score_pairs.begin(), point_score_pairs.end(), PointScoreComparator());
	for ( core::Size i=1; i<=point_score_pairs.size(); i++ ) {
		bool hasneighbor = false;
		numeric::xyzVector< core::Real > x_idx = point_score_pairs[i].first;
		numeric::xyzVector< core::Real > x_cart(x_idx[0],x_idx[1],x_idx[2]);
		core::scoring::electron_density::getDensityMap().idx2cart( x_idx, x_cart );
		for ( core::Size j=1; j<=points_to_search_.size(); j++ ) {
			numeric::xyzVector< core::Real > x_idx_stored = points_to_search_[j];
			numeric::xyzVector< core::Real > x_cart_stored(x_idx_stored[0],x_idx_stored[1],x_idx_stored[2]);
			core::scoring::electron_density::getDensityMap().idx2cart( x_idx_stored, x_cart_stored );
			core::Real distance = (x_cart - x_cart_stored).length();
			if ( distance < point_radius_ ) hasneighbor = true;
		}
		if ( !hasneighbor ) {
			points_to_search_.push_back(point_score_pairs[i].first);
			if ( native_ ) {
				numeric::xyzVector< core::Real > x_cart2;
				core::scoring::electron_density::getDensityMap().idx2cart( x_idx, x_cart2 );
				core::Real distNative = symminfo_.min_symm_dist2(x_cart2, native_com_);
				minDistNative = std::min( minDistNative, distNative );
			}
		}
		if ( points_to_search_.size() >= topNtrans_ ) break;
	}


	if ( basic::options::option[ basic::options::OptionKeys::edensity::debug ]() ) {
		core::scoring::electron_density::ElectronDensity mapdmp = core::scoring::electron_density::getDensityMap();
		mapdmp.set_data(rot);
		mapdmp.writeMRC( "filter.mrc" );

		//dump the points
		core::io::StructFileRep sfr;
		auto & chains = sfr.chains();
		chains.push_back(utility::vector0<core::io::AtomInformation>());
		core::io::StructFileRepOptionsCOP options = utility::pointer::make_shared< core::io::StructFileRepOptions >();
		core::io::AtomInformation ai;
		for ( core::Size i=1; i<=points_to_search_.size(); i++ ) {
			ai.serial = i;
			ai.name = "MG  ";
			ai.resName = " MG";
			ai.chainID = 'A';
			numeric::xyzVector< core::Real > x_cart;
			numeric::xyzVector< core::Real > x_idx = points_to_search_[i];
			core::scoring::electron_density::getDensityMap().idx2cart( x_idx, x_cart );
			ai.x = x_cart[0];
			ai.y = x_cart[1];
			ai.z = x_cart[2];
			ai.occupancy = 1.0;
			ai.element = "MG";
			chains[0].push_back(ai);
		}

		std::string const pdb_contents(core::io::pdb::create_pdb_contents_from_sfr(sfr, options));
		utility::io::ozstream file(std::string("selectedpoints.pdb").c_str(), std::ios::out | std::ios::binary);
		file.write( pdb_contents.c_str(), pdb_contents.size() );
		file.close();
	}

	if ( native_ ) TR << "Closest point to native: " << std::sqrt(minDistNative) << std::endl;
}


void
DockIntoDensityMover::poseSphericalSamples(
	core::pose::Pose const &pose,
	ObjexxFCL::FArray3D< core::Real > & sigR,
	ObjexxFCL::FArray3D< core::Real > & epsR
) {
	using namespace core;

	core::scoring::electron_density::ElectronDensity &density = core::scoring::electron_density::getDensityMap();

	core::Size B=B_;
	core::Real delRsteps=delR_;
	core::Size nRsteps=nRsteps_;

	numeric::xyzVector< Real > reference_atm;
	utility::vector1< numeric::xyzVector< Real > > atmList;
	utility::vector1< Real > all_K, all_C;

	numeric::xyzVector< Real > massSum(0.0,0.0,0.0), centerCA(0.0,0.0,0.0);

	// atom mask ... 3sigma from carbon
	core::Real ATOM_MASK = 3.0 * sqrt( density.getEffectiveBfactor() / (2*M_PI*M_PI) );

	for ( core::Size i=1; i<= pose.size(); ++i ) {
		conformation::Residue const & rsd( pose.residue(i) );
		if ( rsd.aa() == core::chemical::aa_vrt ) continue;
		for ( core::Size j=1; j<= rsd.nheavyatoms(); ++j ) {
			conformation::Atom const & atom( rsd.atom(j) );
			atmList.push_back(atom.xyz());
			massSum += atom.xyz();

			if ( i==(pose.size()+1)/2 && j==2 ) {
				centerCA = atom.xyz();
			}

			chemical::AtomTypeSet const & atom_type_set( rsd.atom_type_set() );
			std::string elt_i = atom_type_set[ rsd.atom_type_index( j ) ].element();
			core::scoring::electron_density::OneGaussianScattering sig_j = core::scoring::electron_density::get_A( elt_i );

			core::Real K_i = sig_j.k( density.getEffectiveBfactor() );
			all_K.push_back( K_i );
			all_C.push_back( sig_j.C( K_i ) );
		}
	}
	int nAtms=atmList.size();
	massSum /= nAtms; // center_of_mass = mass_sum / nAtms

	if ( center_on_middle_ca_ ) {
		massSum = centerCA;
	}

	// precompute sines & cosines
	utility::vector1<core::Real> cT,cG, sT,sG;
	cT.resize(2*B); cG.resize(2*B); sT.resize(2*B); sG.resize(2*B);
	for ( core::Size t=1; t<=2*B; ++t ) {
		core::Real theta = (2.0*t-1.0)*M_PI/(4*B);
		sT[t] = sin(theta);
		cT[t] = cos(theta);
	}
	for ( core::Size p=1; p<=2*B; ++p ) {
		core::Real phi = (2.0*p-2.0)*M_PI/(2*B);
		sG[p] = sin(phi);
		cG[p] = cos(phi);
	}

	//////////////////
	// pose -> spherical-sampled density
	// 1. one models each atom with a Gaussian sphere of density
	// 2. interpolate this calculated density in cencentric spherical shells
	// (extending out to D Ang in 1 Ang steps)
	//////////////////
	if ( laplacian_offset_ != 0 ) {
		TR << "Applying laplacian filter with offset of: " << laplacian_offset_ << " A" << std::endl;
	}
	sigR.dimension( 2*B, 2*B, nRsteps );
	sigR = 0.0;

	epsR.dimension( 2*B, 2*B, nRsteps );
	epsR = 0.0;

	// for each atom
	for ( int i=1; i<=nAtms; ++i ) {
		core::Real k=all_K[i];
		core::Real C=all_C[i];

		atmList[i] -= massSum;

		core::Real atomR = atmList[i].length();
		if ( atomR < 1e-5 ) {
			// uniform contribution to inner shells
			for ( core::Size ridx=1; ridx<=nRsteps; ++ridx ) {
				core::Real atomD = ridx * delRsteps;
				if ( atomD < ATOM_MASK ) {
					core::Real atomH = C * exp(-k*atomD*atomD); // <-- this is the place to calculate density
					for ( core::Size t=1; t<=2*B; ++t ) {
						for ( core::Size p=1; p<=2*B; ++p ) {
							sigR(p,t,ridx) += atomH;
							epsR(p,t,ridx) = 1.0;
						}
					}
				}
			}
			continue;
		}


		core::Real beta = acos( atmList[i][2] / atomR );
		core::Real gamma = atan2( atmList[i][0] , atmList[i][1] );   // x and y switched from usual convention

		core::Real st1 = sin(beta);
		core::Real sg1 = sin(gamma);
		core::Real ct1 = cos(beta);
		core::Real cg1 = cos(gamma);

		if ( laplacian_offset_ != 0 ) {
			for ( core::Size ridx=1; ridx<=nRsteps; ++ridx ) {
				core::Real shellR = ridx * delRsteps;
				for ( core::Size t=1; t<=2*B; ++t ) {
					core::Real minAtomD =  atomR*atomR + shellR*shellR - 2*atomR*shellR*(st1*sT[t]+ct1*cT[t]);
					if ( minAtomD>ATOM_MASK*ATOM_MASK ) continue; // this just loops back to 2xB so we still get to sigR
					for ( core::Size p=1; p<=2*B; ++p ) {
						core::Real atomD = atomR*atomR + shellR*shellR - 2*atomR*shellR*(st1*sT[t]*(sg1*sG[p]+cg1*cG[p])+ct1*cT[t]);
						if ( atomD < ATOM_MASK*ATOM_MASK ) {
							core::Real atomH = C * exp(-k*atomD);
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
					atmList[i][xyz] = atmList[i][xyz] + ( ( (lapl==0) ? 1.0 : -1.0 ) * laplacian_offset_ );
					atomR = atmList[i].length();
					core::Real beta_2 = acos( atmList[i][2] / atomR );
					core::Real gamma_2 = atan2( atmList[i][0] , atmList[i][1] );   // x and y switched from usual convention
					core::Real st1_2 = sin(beta_2);
					core::Real sg1_2 = sin(gamma_2);
					core::Real ct1_2 = cos(beta_2);
					core::Real cg1_2 = cos(gamma_2);
					// residue index
					for ( core::Size ridx=1; ridx<=nRsteps; ++ridx ) {
						core::Real shellR = ridx * delRsteps;
						for ( core::Size t=1; t<=2*B; ++t ) {
							core::Real minAtomD =  atomR*atomR + shellR*shellR - 2*atomR*shellR*(st1_2*sT[t]+ct1_2*cT[t]);
							if ( minAtomD>ATOM_MASK*ATOM_MASK ) continue;
							for ( core::Size p=1; p<=2*B; ++p ) {
								core::Real atomD = atomR*atomR + shellR*shellR - 2*atomR*shellR*(st1_2*sT[t]*(sg1_2*sG[p]+cg1_2*cG[p])+ct1_2*cT[t]);
								if ( atomD < ATOM_MASK*ATOM_MASK ) {
									core::Real atomH = C * exp(-k*atomD);
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
			for ( core::Size ridx=1; ridx<=nRsteps; ++ridx ) {
				core::Real shellR = ridx * delRsteps;
				for ( core::Size t=1; t<=2*B; ++t ) {
					core::Real minAtomD =  atomR*atomR + shellR*shellR - 2*atomR*shellR*(st1*sT[t]+ct1*cT[t]);
					if ( minAtomD>ATOM_MASK*ATOM_MASK ) continue;
					for ( core::Size p=1; p<=2*B; ++p ) {
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


// do the main search over the map
void
DockIntoDensityMover::density_grid_search (
	core::Size pose_idx,
	core::pose::Pose & pose,
	RBfitResultDB & results
) {
	core::scoring::electron_density::ElectronDensity &density = core::scoring::electron_density::getDensityMap();
	//numeric::xyzVector<int> grid = density.getGrid();

	// allocate space for SHT
	numeric::fourier::SHT SOFT(B_, nRsteps_);

	// get com of pose
	numeric::xyzMatrix<core::Real> rot;
	numeric::xyzVector<core::Real> pretrans, posttrans;
	get_radius( pose, pretrans );

	pretrans=-1.0*pretrans;

	runtime_assert( points_to_search_.size() >= 1 ); // sanity check

	// get pose SPHARM
	ObjexxFCL::FArray3D< core::Real > poseSig, poseCoefR, poseCoefI;
	ObjexxFCL::FArray3D< core::Real > poseEps, poseEpsCoefR, poseEpsCoefI;
	poseSphericalSamples( pose, poseSig, poseEps );

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
	for ( core::Size r=0; r<nRsteps_; ++r ) {
		sumRWt += 4.0 * numeric::square((r+1.0));
	}
	for ( core::Size t=0; t<2*B_; ++t ) {
		sumTwt += std::sin( (t+0.5) * numeric::constants::d::pi / (2*B_) );
	}

	core::Real sumEps=0.0, sumSig=0.0, sumSig2=0.0;
	for ( core::Size r=0; r<nRsteps_; ++r ) {
		double thisRWt = 4.0 * numeric::square((r+1.0)) / sumRWt;
		for ( core::Size t=0; t<2*B_; ++t ) {
			double thisTWt = std::sin( (t+0.5) * numeric::constants::d::pi / (2*B_) ) / sumTwt;
			for ( core::Size p=0; p<2*B_; ++p ) {
				double thisWt = thisRWt*thisTWt/(2*B_);
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
	for ( core::Size i=1; i<=points_to_search_.size(); ++i ) {
		// get cartesian coords of ths point
		density.idx2cart( points_to_search_[i], posttrans );

		density.mapSphericalSamples( mapSig, nRsteps_, delR_, B_, points_to_search_[i], laplacian_offset_ );
		map2Sig = mapSig;
		for ( core::Size j=0; j<4*B_*B_*nRsteps_; ++j ) {
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
		core::Size nperRot = std::min( 100*max_rot_per_trans_ , B_*B_*B_);
		core::Size nperRotCl = max_rot_per_trans_;

		RBfitResultDB local_results( nperRot );

		for ( core::Size j=0; j<8*B_*B_*B_; ++j ) {
			core::Real sumSigMap = so3_correl[j] / (4*numeric::constants::d::pi);
			core::Real sumMap = mask_correl[j] / (4*numeric::constants::d::pi);
			core::Real sumMap2 = mask2_correl[j] / (4*numeric::constants::d::pi);

			double CC = (
				(sumEps*sumSigMap - sumSig*sumMap) / (
				std::sqrt( sumEps*sumSig2 - sumSig*sumSig )
				* std::sqrt( sumEps*sumMap2 - sumMap*sumMap )
				) );

			if ( local_results.to_add_element( CC ) ) {
				SOFT.idx_to_rot(j , rot);
				local_results.add_element( RBfitResult( pose_idx, CC, rot, pretrans, posttrans ) );
			}
		}


		do_filter( local_results );
		while ( local_results.size() > nperRotCl ) local_results.pop();

		core::Size nclust = local_results.size();

		core::Real bestrms=1e6,bestscore=0;
		while ( local_results.size() > 0 ) {
			RBfitResult sol_i = local_results.pop();

			if ( native_ ) {
				core::pose::PoseOP posecopy ( new core::pose::Pose(pose) );
				apply_transform( *posecopy, sol_i );
				core::pose::addVirtualResAsRoot( *posecopy );
				core::Real rms_i = get_rms(native_, posecopy, symminfo_);
				if ( rms_i < bestrms ) {
					bestrms = rms_i; bestscore = sol_i.score_;
				}
			}

			results.add_element( sol_i );
		}

		core::Real minDistNative=1e6;
		if ( native_ ) {
			core::Real distNative = (posttrans-native_com_).length_squared();
			minDistNative = std::min( minDistNative, distNative );
		}

		if ( std::sqrt(minDistNative) <5.0 ) {
			TR << "[" << i << "/" << points_to_search_.size() << "]" << " nmdls " << nclust << " pointrms " << std::sqrt(minDistNative)
				<< " rms " << bestrms << " score " << bestscore << std::endl;
		}
		if ( i%100 == 0 ) {
			TR << "[" << i << "/" << points_to_search_.size() << "] " << results.top().score_ << " / " << results.size() << std::endl;
		}
	}

	//TR << "[" << points_to_search_.size() << "/" << points_to_search_.size() << "]" << std::endl;
}


void
DockIntoDensityMover::do_refinement (
	utility::vector1< core::pose::PoseOP > const &poses,
	RBfitResultDB & results_in,
	RefinementResultDB & results_out
) {
	using namespace core::pack;
	using namespace core::pack::task;
	using namespace core::pack::task::operation;

	core::scoring::ScoreFunctionOP scorefxn_dens( new core::scoring::ScoreFunction() );
	scorefxn_dens->set_weight( core::scoring::elec_dens_whole_structure_allatom, 1.0);

	core::scoring::ScoreFunctionOP scorefxn_refine = core::scoring::get_score_function();
	scorefxn_refine->set_weight( core::scoring::elec_dens_fast, dens_wt_);
	core::scoring::ScoreFunctionOP scorefxn_refine_rb( new core::scoring::ScoreFunction() );
	scorefxn_refine_rb->set_weight( core::scoring::elec_dens_fast, dens_wt_);

	core::kinematics::MoveMapOP bbmm( new core::kinematics::MoveMap );
	bbmm->set_bb( true ); bbmm->set_chi( true ); bbmm->set_jump( true );
	protocols::minimization_packing::MinMoverOP bbmin( new protocols::minimization_packing::MinMover(
		bbmm, scorefxn_refine, "lbfgs_armijo_nonmonotone", 1e-5, true ) );
	bbmin->max_iter(200); // make a parameter?

	// packer
	TaskFactoryOP tf( new TaskFactory() );
	tf->push_back( utility::pointer::make_shared< InitializeFromCommandline >()); // get extra rotamer flags from command line
	tf->push_back( utility::pointer::make_shared< operation::IncludeCurrent >()); // include current rotamer by default
	tf->push_back( utility::pointer::make_shared< RestrictToRepacking >()); // do not design
	protocols::minimization_packing::PackRotamersMoverOP packer( new protocols::minimization_packing::PackRotamersMover() );
	packer->task_factory( tf );
	packer->score_function( scorefxn_refine );

	protocols::simple_moves::SwitchResidueTypeSetMover to_all_atom( core::chemical::FA_STANDARD );

	// do refinement
	core::Size ntotal=results_in.size();
	while ( results_in.size() > 0 ) {
		RBfitResult sol_i = results_in.pop();
		//TR << "[" << ntotal-results_in.size() << "/" << ntotal << "]" << "\r" << std::flush;

		core::pose::PoseOP posecopy ( new core::pose::Pose( *(poses[ sol_i.pose_idx_ ]) ) );
		apply_transform( *posecopy, sol_i );
		core::pose::addVirtualResAsRoot( *posecopy );

		// Setup rigid-body movemap now! (we need to know root jump number)
		core::kinematics::MoveMapOP rbmm = utility::pointer::make_shared< core::kinematics::MoveMap >();
		rbmm->set_bb( false ); rbmm->set_chi( false ); rbmm->set_jump( false );
		int root = posecopy->fold_tree().root();
		utility::vector1< core::kinematics::Edge > root_edges = posecopy->fold_tree().get_outgoing_edges(root);
		for ( core::Size i=1; i<=root_edges.size(); ++i ) rbmm->set_jump ( root_edges[i].label() , true );

		// Setup rigid-body min now!
		protocols::minimization_packing::MinMoverOP rbmin(
			new protocols::minimization_packing::MinMover( rbmm, scorefxn_refine_rb, "lbfgs_armijo_nonmonotone", 1e-5, true ) );
		rbmin->max_iter(200); // make a parameter?

		core::Real scoreb = (*scorefxn_dens)(*posecopy), scorei=0;

		// rbmin/pack/fullmin
		if ( do_refine_ ) {
			rbmin->apply( *posecopy );

			if ( posecopy->is_centroid() ) {
				to_all_atom.apply( *posecopy );
			}

			for ( core::Size i=1; i<=ncyc_; ++i ) {
				packer->apply( *posecopy );
				scorei = (*scorefxn_dens)(*posecopy);

				if ( min_backbone_ ) {
					//scorefxn_refine->show( *posecopy );
					bbmin->apply( *posecopy );
					//scorefxn_refine->show( *posecopy );
				} else {
					rbmin->apply( *posecopy );
				}
			}
		}

		// rescore (expensive sf)
		core::Real scoref = (*scorefxn_dens)(*posecopy);

		core::Real rms=0.0;
		if ( native_ ) {
			rms = get_rms( RefinementResult( 0.0,0.0,0.0, posecopy ), RefinementResult( 0.0,0.0,0.0, native_ ), symminfo_ );
		}

		TR << "[" << ntotal-results_in.size() << "/" << ntotal << "] " << scoreb << " : " << scorei << " : " << scoref << "  rms=" << rms << std::endl;

		// store
		results_out.add_element( RefinementResult( -scoref,-scoreb,sol_i.score_, posecopy ) );
	}
	TR << std::endl;
}


// read results, normalize scores, remove redundant (CA rms)
void
DockIntoDensityMover::do_filter( RefinementResultDB & results ) {
	// dumb -- copy to vector
	utility::vector1< RefinementResult > results_sort;
	while ( results.size() > 0 )
			results_sort.push_back( results.pop() );
	utility::vector1< bool > selector(results_sort.size(), false);

	for ( int i=results_sort.size(); i>=1; --i ) {
		selector[i] = true;
		for ( int j=i+1; j<=(int)results_sort.size() && selector[i]; ++j ) {
			if ( selector[j] && get_rms(results_sort[i], results_sort[j], symminfo_)<cluster_radius_ ) {
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
DockIntoDensityMover::do_filter( utility::vector1< core::pose::PoseOP > const &poses, RBfitResultDB & results, bool rescore ) {
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
		core::pose::PoseOP posecopy ( new core::pose::Pose( *(poses[ results_sort[i].pose_idx_ ]) ) );
		apply_transform( *posecopy, results_sort[i] );

		if ( rescore ) {
			TR << "[" << i << "/" << results_sort.size() << "]" << "\r" << std::flush;

			core::pose::addVirtualResAsRoot( *posecopy );
			results_sort[i].score_ = (*scorefxn_dens)(*posecopy);
		}
		for ( int j=1; j<=(int)selected.size() && selector[i]; ++j ) {
			if ( get_rms(posecopy, selected[j], symminfo_) <cluster_radius_ ) {
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
DockIntoDensityMover::do_filter( RBfitResultDB & results ) {
	// dumb -- copy to vector
	utility::vector1< RBfitResult > results_sort;
	while ( results.size() > 0 )
			results_sort.push_back( results.pop() );
	utility::vector1< bool > selector(results_sort.size(), false);

	core::Real nsel=0;
	for ( int i=results_sort.size(); i>=1; --i ) {
		selector[i] = true;

		for ( int j=i+1; j<=(int)selector.size() && selector[i]; ++j ) {
			core::Real est_rms = 0.6*delR_*nRsteps_ * get_rot_angle( results_sort[i].rotation_ * numeric::inverse(results_sort[j].rotation_) );

			if ( selector[j] && est_rms < cluster_radius_ ) {
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


void
DockIntoDensityMover::apply( core::pose::Pose & pose) {
	// call multipose mover
	utility::vector1< core::pose::PoseOP > temp;
	temp.push_back( utility::pointer::make_shared< core::pose::Pose >(pose) );
	apply_multi( temp );
	pose = *(temp[1]);
}


void
DockIntoDensityMover::apply_multi( utility::vector1< core::pose::PoseOP > & poses) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pose;

	std::string base_name = tag_;
	if ( base_name.size() == 0 ) {
		// get output name ... assumes jd2 ... anything that doesn't use jd2 must call setTag
		base_name = protocols::jd2::current_input_tag();
		utility::vector1< std::string > temp_out_names= utility::split( base_name );
		utility::file::FileName out_name = utility::file::combine_names( temp_out_names );
		base_name = out_name.base();
	}

	RBfitResultDB results_filtered(cluster_oversample_*topNfilter_);    // oversample for clustering
	RefinementResultDB results_refine(cluster_oversample_*topNfinal_);  // oversample for clustering

	if ( passthrough_ ) {
		poses.clear();
		poses.push_back( native_ );

		RBfitResultDB results(poses.size());
		for ( int i=1; i<=(int)poses.size(); ++i ) {
			results.add_element(
				RBfitResult(1, -1000, numeric::xyzMatrix<core::Real>::identity(), numeric::xyzVector<core::Real>(0,0,0), numeric::xyzVector<core::Real>(0,0,0))
			);
			do_refinement ( poses, results, results_refine );
		}

		// force score to -1000 <---- HACK
		RefinementResult sol_i = results_refine.pop();
		sol_i.score_ = -1000;
		results_refine.add_element(sol_i);
	} else {
		// for each input fragment...
		for ( int i=1; i<=(int)poses.size(); ++i ) {
			core::pose::PoseOP pose_i = poses[i];

			// set nRsteps
			TR << " *** Input " << i << " of " << poses.size() << " ***" << std::endl;
			numeric::xyzVector< core::Real > com;

			// 1: get points to sample
			//  NOTE: this sets nRsteps_!
			select_points( *pose_i );
			TR << "Searching a max radius of " << nRsteps_*delR_ << " (" << nRsteps_ << " steps)" << std::endl;
			TR << "Selected " << points_to_search_.size() << " translations to search" << std::endl;

			// 2: do the grid search
			//    top N per position
			RBfitResultDB results(cluster_oversample_*topNfilter_);

			density_grid_search( i, *pose_i, results );
			TR << "Have " << results.size() << " solutions with score >= " << results.top().score_ << std::endl;
			if ( native_ ) print_best_rms( *pose_i, results );

			// 3: rescore
			if ( cluster_radius_>0.0 ) {
				do_filter( poses, results, false ); // rescore=true
			}

			while ( results.size() > topNfilter_ ) results.pop();

			TR << "Have " << results.size() << " solutions after clustering" << std::endl;
			if ( native_ ) print_best_rms( *pose_i, results );

			// 4: combine
			while ( results.size() > 0 ) {
				results_filtered.add_element( results.pop() );
			}
		}
	}

	TR << "Final output" << std::endl;

	// now cluster the aggregate pool (if necessary)
	if ( poses.size()>1 ) {
		TR << "Have " << results_filtered.size() << " total solutions before clustering" << std::endl;
		do_filter( poses, results_filtered, false ); // rescore=false
	}

	TR << "Have " << results_filtered.size() << " total solutions after clustering" << std::endl;

	// 5: refinement
	do_refinement ( poses, results_filtered, results_refine );
	TR << "Have " << results_refine.size() << " total solutions after refinement" << std::endl;

	// cluster a final time (in case structures minimized into the same energy well)
	// also apply local rescoring
	do_filter( results_refine );

	// only dump the top N
	while ( results_refine.size() > topNfinal_ ) { results_refine.pop(); }

	// dump the final structres
	core::io::silent::SilentFileOptions opts;
	core::io::silent::SilentFileData silent_file_data( silent_, false, false, "binary", opts ); //true to store argv in silent file

	while ( results_refine.size() > 0 ) {
		// retrieve data from results_list
		RefinementResult sol_i = results_refine.pop();
		core::pose::PoseOP posecopy = sol_i.pose_;

		// TO DO: fullatom refinement
		if ( silent_.size() > 0 ) {
			core::Real rms=0.0;
			if ( native_ ) {
				rms = get_rms( sol_i, RefinementResult( 0.0,0.0,0.0, native_ ), symminfo_ );
			}

			std::string silent_fn = base_name+"_"+utility::to_string( results_refine.size()+1 );
			core::io::silent::BinarySilentStruct silent_stream( opts, *posecopy, silent_fn );
			silent_stream.add_energy( "spharm_score", sol_i.spharm_score_ );
			silent_stream.add_energy( "init_score", sol_i.prerefine_score_ );
			silent_stream.add_energy( "dens_score", sol_i.score_ );
			silent_stream.add_energy( "rms", rms );
			silent_file_data.write_silent_struct( silent_stream, silent_ );
		} else {
			// tag
			core::io::RemarkInfo remark;
			std::ostringstream oss;

			if ( posecopy->pdb_info() ) {
				oss << "RANK = " << results_refine.size()+1;

				remark.num = 1; remark.value = oss.str();
				posecopy->pdb_info()->remarks().push_back( remark );

				oss.str(""); oss.clear();
				oss << "SCORE = " << sol_i.score_;
				remark.num = 1; remark.value = oss.str();
				posecopy->pdb_info()->remarks().push_back( remark );
			}

			std::string outname = base_name+"_"+ObjexxFCL::right_string_of( results_refine.size()+1, 6, '0' )+".pdb";
			posecopy->dump_pdb( outname );
		}
	}
}




}
}

