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
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>

#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>

#include <core/conformation/Residue.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>

#include <core/scoring/electron_density/ElectronDensity.hh>
#include <core/scoring/electron_density/xray_scattering.hh>
#include <core/scoring/ScoreFunction.hh>

#include <core/id/AtomID.hh>

#include <core/io/pdb/pdb_writer.hh>
#include <core/io/StructFileRep.hh>
#include <core/io/StructFileRepOptions.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinarySilentStruct.hh>

#include <core/types.hh>

#include <utility/io/ozstream.hh>
#include <utility/vector1.hh>
#include <utility/exit.hh>

#include <numeric/fourier/SHT.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzMatrix.hh>

#include <basic/Tracer.hh>


static basic::Tracer TR( "protocols.electron_density.DockIntoDensityUtils" );

namespace protocols {
namespace electron_density {


std::ostream& operator<< (std::ostream& out, const RBfitResult& result) {
	out << "RBfitResult( pose_idx: " << result.pose_idx_ << ", score: " << result.score_
		<< ", pre_trans: (" << result.pre_trans_.x()  << ", " << result.pre_trans_.y()  << ", " << result.pre_trans_.z()
		<< "), post_trans: (" << result.post_trans_.x()  << ", " << result.post_trans_.y() << ", " << result.post_trans_.z()
		<< "), rotation: ()" << std::endl;
	return out;
}
// for classes
//
template< typename T >
void
write_RBfitResultDB( RBfitResultDB fit_result_DB, T & outresults ) {
	while ( fit_result_DB.size() > 0 ) {
		core::Size const rank = fit_result_DB.size();
		RBfitResult const sol_i = fit_result_DB.pop();
		// write out RBfitResult, so we can combine them later
		outresults << sol_i.pose_idx_ << "\t" << rank << "\t" << sol_i.score_
			<< "\t" << sol_i.rotation_.xx() << "\t" << sol_i.rotation_.xy() << "\t" << sol_i.rotation_.xz()
			<< "\t" << sol_i.rotation_.yx() << "\t" << sol_i.rotation_.yy() << "\t" << sol_i.rotation_.yz()
			<< "\t" << sol_i.rotation_.zx() << "\t" << sol_i.rotation_.zy() << "\t" << sol_i.rotation_.zz()
			<< "\t" << sol_i.pre_trans_[0] << "\t" << sol_i.pre_trans_[1] << "\t" <<  sol_i.pre_trans_[2]
			<< "\t" << sol_i.post_trans_[0] << "\t" << sol_i.post_trans_[1] << "\t" << sol_i.post_trans_[2] << std::endl;
	}
}


template< typename T >
void
dump_RefinementDB_to_silent(
	T resultDB,
	std::string const & outfile,
	std::string const & tag_prefix,
	std::string const & final_chain,
	bool const centroid_output,
	bool const append_to_outfile,
	utility::vector1< core::pose::PoseCOP > const & natives,
	DensitySymmInfo const & symminfo,
	bool const legacy_rms
) {
	if ( !append_to_outfile ) {
		remove( outfile.c_str() );
	}

	core::io::silent::SilentFileOptions sfopts;
	core::io::silent::SilentFileData silent_file_data( outfile, false, false, "binary", sfopts );
	protocols::simple_moves::SwitchResidueTypeSetMover to_cen("centroid");
	protocols::simple_moves::SwitchResidueTypeSetMover to_fa("fa_standard");
	while ( resultDB.size() > 0 ) {
		core::Size const rank = resultDB.size();
		RefinementResult const sol_i = resultDB.pop();

		std::string const current_tag( tag_prefix + "_" + ObjexxFCL::right_string_of(rank, 6, '0') );

		core::pose::Pose const posecopy = [&]{
			core::pose::Pose tmp_pose(*sol_i.pose_);

			if ( centroid_output ) to_cen.apply(tmp_pose);
			else to_fa.apply(tmp_pose);

			if ( tmp_pose.pdb_info() == nullptr ) {
				tmp_pose.pdb_info(utility::pointer::make_shared< core::pose::PDBInfo >(tmp_pose.size()));
			}
			// Temporary until chains are strings i guess
			if ( final_chain.at(0) != '^' ) tmp_pose.pdb_info()->set_chains( final_chain.at(0) );

			return tmp_pose;
		}();

		core::io::silent::BinarySilentStruct silent_stream( sfopts, posecopy, current_tag );
		silent_stream.clear_energies();

		silent_stream.add_energy( "score", -sol_i.score_ ); // REQUIRED!! order is IMPORTANT
		silent_stream.add_energy( "spharm_score", sol_i.spharm_score_ );
		silent_stream.add_energy( "init_score", sol_i.prerefine_score_ );
		silent_stream.add_energy( "dens_score", sol_i.score_ );
		silent_stream.add_energy( "dens_rank", rank );
		core::Real best_rms = 9999.0, best_gdt = 0;
		core::Size best_rms_idx = 0;
		for ( core::Size i=1; i<=natives.size(); ++i ) {
			core::Real const rms = [&]{
				if ( legacy_rms ) return get_rms( *sol_i.pose_, *natives.at(i), symminfo );
				else { return get_rms( *sol_i.pose_, *natives.at(i), symminfo, true ); }
			}();
			// this isn't perfect but i guess it's ok
			if ( rms < best_rms ) {
				best_rms = rms;
				best_rms_idx = i;
				if ( !legacy_rms ) best_gdt = get_gdt( *sol_i.pose_, *natives[i], symminfo, true );
			}
		}

		silent_stream.add_energy( "rms", best_rms );
		if ( !legacy_rms ) silent_stream.add_energy( "gdt", best_gdt );
		silent_stream.add_energy( "native_idx", best_rms_idx );
		silent_file_data.write_silent_struct( silent_stream, outfile );
	}
}

template
void
dump_RefinementDB_to_silent<RevRefinementResultDB>(
	RevRefinementResultDB resultDB,
	std::string const & outfile,
	std::string const & tag_prefix,
	std::string const & final_chain,
	bool const centroid_output,
	bool const append_to_outfile,
	utility::vector1< core::pose::PoseCOP > const & natives,
	DensitySymmInfo const & symminfo,
	bool const legacy_rms
);


template
void
dump_RefinementDB_to_silent<RefinementResultDB>(
	RefinementResultDB resultDB,
	std::string const & outfile,
	std::string const & tag_prefix,
	std::string const & final_chain,
	bool const centroid_output,
	bool const append_to_outfile,
	utility::vector1< core::pose::PoseCOP > const & natives,
	DensitySymmInfo const & symminfo,
	bool const legacy_rms
);


void
compare_RBfitDB_to_native(
	RBfitResultDB resultDB,
	core::pose::Pose const & pose,
	core::pose::PoseCOPs const & natives,
	utility::vector1< numeric::xyzVector< core::Real > > const & native_coms,
	utility::vector1< numeric::xyzVector< core::Real > > const & native_middle_cas,
	DensitySymmInfo const & symminfo,
	bool const rot_middle_ca,
	core::Real const rms_cutoff ) {
	core::Size initial_DB_size = resultDB.size();
	core::Size num_results_w_rms_under_cutoff = 0;

	TR.Debug << "Checking the RMS of " << resultDB.size() << " results to the native" << std::endl;

	while ( resultDB.size() > 0 ) {
		core::Size const result_idx = resultDB.size();
		RBfitResult const sol_i = resultDB.pop();

		for ( core::Size i = 1; i <= natives.size(); ++i ) {
			// only check if translation within 10 A of native rotation centers
			if ( rot_middle_ca && sqrt( symminfo.min_symm_dist2( sol_i.post_trans_, native_middle_cas[i] ) ) > 10 ) continue;
			else if ( sqrt( symminfo.min_symm_dist2( sol_i.post_trans_, native_coms[i] ) ) > 10 ) continue;

			core::pose::Pose const posecopy = [&]{
				core::pose::Pose tmp_pose(pose);
				apply_transform( tmp_pose, sol_i );
				core::pose::addVirtualResAsRoot( tmp_pose );
				return tmp_pose;
			}();
			core::Real const rms = get_rms(posecopy, *natives[i], symminfo, true);
			core::Real const gdt_result = get_gdt(posecopy, *natives[i], symminfo, true);
			if ( rms <= rms_cutoff ) {
				TR << "Search kept result: rank: " << result_idx << "/" << initial_DB_size << " rms: " << rms << " GDT: " << gdt_result << " score: " << sol_i.score_
					<< " closest to com_idx " << i << std::endl;
				++num_results_w_rms_under_cutoff;
			}
		}
	}
	TR << "Completed native comparison between " << initial_DB_size << " results. Have " << num_results_w_rms_under_cutoff << " results <= " << rms_cutoff << " rms." << std::endl;
}


// non-superposed RMS
core::Real
get_rms(core::pose::Pose const & r1, core::pose::Pose const & r2, DensitySymmInfo const &d) {
	runtime_assert( r1.size() == r2.size() );
	core::Size const nres = r1.size();
	core::Real rms=0.0;
	core::Size N=0;
	for ( core::Size i=1; i<=nres; ++i ) {
		if ( !r1.residue(i).is_protein() ) continue;
		rms += d.min_symm_dist2( r1.residue(i).xyz(2), r2.residue(i).xyz(2) );
		N++;
	}

	return (std::sqrt(rms/N));
}


template< typename T >
void
dump_and_raise_bad_pose_alignment(
	core::pose::Pose const & r1,
	char const r1_chain,
	core::Size const r1_resnum,
	core::Size const r1_posenum,
	core::pose::Pose const & r2,
	core::Size const r2_posenum,
	T & out) {
	out << "r1 seq: " << r1.sequence() << std::endl;
	out << "r2 seq: " << r2.sequence() << std::endl;
	out << "r1 chain: " << r1_chain << " r1 pdbinfo resnum: " << r1_resnum << std::endl;
	out << "r2 chain: " << r2.pdb_info()->chain(r2_posenum) << " r2 pdbinfo resnum: " << r2.pdb_info()->number(r2_posenum) << std::endl;
	out << "r1 name: " << r1.residue(r1_posenum).name3() << " r2 name: " << r2.residue(r2_posenum).name3() << std::endl;
	out << "r1 pose num: " << r1_posenum << " r2 pose num: " << r2_posenum << std::endl;
	throw CREATE_EXCEPTION(utility::excn::BadInput, "Found that in get_rms or get_gdt, poses did not align in residue numbering or sequence");
}


template< typename T >
void
dump_and_raise_no_pose_alignment(
	core::pose::Pose const & r1,
	core::pose::Pose const & r2,
	T & out) {
	out << "r1 seq: " << r1.sequence() << std::endl;
	out << "r1 start " << r1.pdb_info()->number(1) << std::endl;
	out << "r1 chain: " << r1.pdb_info()->chain(1) << std::endl;
	out << "r2 seq: " << r2.sequence() << std::endl;
	out << "r2 start " << r2.pdb_info()->number(1) << std::endl;
	out << "r2 chain: " << r2.pdb_info()->chain(2) << std::endl;
	throw CREATE_EXCEPTION(
		utility::excn::BadInput,
		"0 residues aligned when getting rms or gdt. something must be wrong with your inputs");
}


core::Real
get_rms(core::pose::Pose const & r1, core::pose::Pose const & r2, DensitySymmInfo const &d, bool const native) {
	core::Real rms=0.0;
	core::Size N=0;
	for ( core::Size i = 1; i <= r1.size(); ++i ) {
		int const r1_resnum = r1.pdb_info()->number(i);
		char const r1_chain = r1.pdb_info()->chain(i);
		core::Size const r2_posenum = r2.pdb_info()->pdb2pose( r1_chain, r1_resnum );
		if ( r2_posenum == 0 ) continue;

		if ( r1.residue(i).name3() != r2.residue(r2_posenum).name3() ) {
			dump_and_raise_bad_pose_alignment(
				r1, r1_chain, r1_resnum, i,
				r2, r2_posenum,
				TR.Error
			);
		}

		core::Real const distance = d.min_symm_dist2( r1.residue(i).xyz(2), r2.residue(r2_posenum).xyz(2) );
		rms += distance;
		++N;
	}

	if ( N == 0 ) {
		if ( native ) {
			// Sometimes we have a native that is the same chain, but is out of frame of our input pose.
			// So if we're doing a native check, we can just return 1000, but if we're doing clustering
			// this should never happen and we fail
			TR.Warning << "pose did not align at all N==0, does your native pose contain your query seq?" << std::endl;
			return 1000.0;
		} else {
			dump_and_raise_no_pose_alignment(r1, r2, TR.Error);
			return 1000.0;  // tidy wants this
		}
	} else if ( N <= 5 ) {
		TR.Warning << "WARNING only aligned to native with " << N << " residues! that's too small to make sense! returning 1000" << std::endl;
		return 1000.0;
	} else return std::sqrt(rms/N);
}


// non-superposed RMS
core::Real
get_rms(RefinementResult const & r1, RefinementResult const & r2, DensitySymmInfo const & d ) {
	return get_rms(*r1.pose_, *r2.pose_, d);
}



core::Real
get_gdt(
	core::pose::Pose const & r1,
	core::pose::Pose const & r2,
	DensitySymmInfo const &d,
	bool native)
{
	core::Real one=0.0, two=0.0, four=0.0, eight=0.0;
	core::Size N=0;
	for ( core::Size i = 1; i <= r1.size(); ++i ) {
		int const r1_resnum = r1.pdb_info()->number(i);
		char const r1_chain = r1.pdb_info()->chain(i);
		core::Size const r2_posenum = r2.pdb_info()->pdb2pose( r1_chain, r1_resnum );
		if ( r2_posenum == 0 ) continue;

		if ( r1.residue(i).name3() != r2.residue(r2_posenum).name3() ) {
			dump_and_raise_bad_pose_alignment(
				r1, r1_chain, r1_resnum, i,
				r2, r2_posenum,
				TR.Error
			);
		}

		core::Real const distance = d.min_symm_dist2( r1.residue(i).xyz(2), r2.residue(r2_posenum).xyz(2) );
		if ( distance <= 1 ) ++one;
		if ( distance <= 2 ) ++two;
		if ( distance <= 4 ) ++four;
		if ( distance <= 8 ) ++eight;
		++N;
	}

	if ( N == 0 ) {
		if ( native ) {
			// Sometimes we have a native that is the same chain, but is out of frame of our input pose.
			// So if we're doing a native check, we can just return 0.0, but if we're doing clustering
			// this should never happen and we fail
			TR.Warning << "pose did not align at all N==0, does your native pose contain your query seq?" << std::endl;
			return 0.0;
		} else {
			dump_and_raise_no_pose_alignment(r1, r2, TR.Error);
			return 0.0;
		}
	} else if ( N <= 5 ) {
		TR.Warning << "WARNING only aligned to native with " << N << " residues! that's too small to make sense! returning 0.0" << std::endl;
		return 0.0;
	} else {
		core::Real const gdt = 100 * (one + two + four + eight) / (4 * N);
		return gdt;
	}
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
dump_points_to_search_to_pdb_or_txt(
	utility::vector1< numeric::xyzVector< core::Real > > const & points_to_search,
	std::string const & pdb_filename,
	std::string const & txt_filename
) {
	//dump the points
	if ( !txt_filename.empty() ) {
		//dump the points
		std::ofstream outpoints;
		outpoints.open(txt_filename.c_str());
		for ( core::Size i=1; i<=points_to_search.size(); i++ ) {
			outpoints << std::to_string(i) << "\t" << points_to_search[i][0] << "\t" << points_to_search[i][1] << "\t" << points_to_search[i][2] << std::endl;
		}
		outpoints.close();
	}

	if ( !pdb_filename.empty() ) {
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
		utility::io::ozstream file(std::string(pdb_filename).c_str(), std::ios::out | std::ios::binary);
		file.write( pdb_contents.c_str(), pdb_contents.size() );
		file.close();
	}
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


	utility::vector1< core::Real > distance_to_native_points( params.natives_.size(), 10000 );
	utility::vector1< core::Size > best_native_point_ranks( params.natives_.size(), 10000 );

	for ( core::Size i=1; i<=point_score_pairs.size(); i++ ) {
		bool hasneighbor = false;
		numeric::xyzVector< core::Real > const x_idx = point_score_pairs[i].first;
		numeric::xyzVector< core::Real > x_cart(x_idx[0],x_idx[1],x_idx[2]);
		core::scoring::electron_density::getDensityMap().idx2cart( x_idx, x_cart );
		for ( core::Size j=1; j<=points_to_search.size(); j++ ) {
			numeric::xyzVector< core::Real > const x_idx_stored = points_to_search[j];
			numeric::xyzVector< core::Real > x_cart_stored(x_idx_stored[0],x_idx_stored[1],x_idx_stored[2]);
			core::scoring::electron_density::getDensityMap().idx2cart( x_idx_stored, x_cart_stored );
			core::Real const distance = (x_cart - x_cart_stored).length();
			if ( distance < params.point_radius_ ) hasneighbor = true;
		}
		if ( !hasneighbor ) {
			points_to_search.push_back(point_score_pairs[i].first);

			for ( core::Size native_i = 1; native_i <= params.natives_.size(); ++native_i ) {

				core::Real const distNative = [&]{
					if ( params.center_on_middle_ca_ ) {
						return params.symminfo_.min_symm_dist2(x_cart, params.native_middle_cas_[native_i]);
					} else {
						return params.symminfo_.min_symm_dist2(x_cart, params.native_coms_[native_i]);
					}
				}();

				if ( distNative < distance_to_native_points[ native_i ] ) {
					best_native_point_ranks[ native_i ] = points_to_search.size();
					distance_to_native_points[ native_i ] = distNative;
				}
			}
		}
		if ( points_to_search.size() >= params.topNtrans_ ) break;
	}


	if ( basic::options::option[ basic::options::OptionKeys::edensity::debug ]() ) {
		core::scoring::electron_density::ElectronDensity mapdmp = core::scoring::electron_density::getDensityMap();
		mapdmp.set_data(rot);
		mapdmp.writeMRC( "filter.mrc" );

		dump_points_to_search_to_pdb_or_txt( points_to_search, "selectedpoints.pdb", "");
	}

	for ( core::Size i = 1; i <= params.natives_.size(); ++i ) {
		TR << "Closest point to ";
		if ( params.center_on_middle_ca_ ) {
			TR << "native middle CA";
		} else {
			TR << "native COM";
		}
		TR << " index: " << i << " is: " <<  std::sqrt(distance_to_native_points[i]) <<
			" at rank: " << best_native_point_ranks[i] << std::endl;
	}

	return points_to_search;
}


// read results, normalize scores, remove redundant (CA rms)
void
cluster_RefinementDB( RefinementResultDB & results, DensitySymmInfo const & symminfo, core::Real const cluster_radius, core::Size const target_size ) {
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
		if ( results.size() == target_size && target_size != 0 ) return;
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
		core::pose::PoseOP posecopy( poses[ results_sort[i].pose_idx_ ]->clone() );
		apply_transform( *posecopy, results_sort[i] );

		if ( rescore ) {
			TR << "[" << i << "/" << results_sort.size() << "]" << "\r" << std::flush;

			core::pose::addVirtualResAsRoot( *posecopy );
			results_sort[i].score_ = (*scorefxn_dens)(*posecopy);
		}
		for ( int j=1; j<=(int)selected.size() && selector[i]; ++j ) {
			if ( get_rms(*posecopy, *selected[j], symminfo) <cluster_radius ) {
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
cluster_RBfitResultDB_fast(
	RBfitResultDB & results,
	core::Size const delR,
	core::Size const nRsteps,
	core::Real const cluster_radius,
	core::Size const max_results,
	bool const include_distance,
	core::scoring::electron_density::ElectronDensity const & dens) {
	// dumb -- copy to vector
	utility::vector1< RBfitResult > results_sort;
	while ( results.size() > 0 )
			results_sort.push_back( results.pop() );
	utility::vector1< bool > selector(results_sort.size(), false);

	numeric::xyzVector< core::Real > i_posttrans_cart, j_posttrans_cart;
	core::Real nsel=0, distance=0.0;
	for ( core::Size i=results_sort.size(); i>=1; --i ) {
		if ( max_results != 0 && max_results == results.size() ) break;
		selector[i] = true;
		if ( include_distance ) {
			numeric::xyzVector< core::Real > const & i_posttrans_idx = results_sort[i].post_trans_;
			dens.idx2cart( i_posttrans_idx, i_posttrans_cart );
		}

		for ( core::Size j=i+1; j<=selector.size() && selector[i]; ++j ) {
			if ( include_distance ) {
				numeric::xyzVector< core::Real > const & j_posttrans_idx = results_sort[j].post_trans_;
				dens.idx2cart( j_posttrans_idx, j_posttrans_cart );
				distance = (i_posttrans_cart - j_posttrans_cart).length();
			}
			core::Real const est_rms = distance + 0.6*delR*nRsteps * get_rot_angle( results_sort[i].rotation_ * numeric::inverse(results_sort[j].rotation_) );

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
	core::Size const pose_idx,
	core::pose::Pose const & pose,
	RBfitResultDB & results,
	utility::vector1< numeric::xyzVector< core::Real > > const & points_to_search,
	DensityGridSearchOptions const & params
) {
	if ( params.nRsteps_ == 0  ) {
		throw CREATE_EXCEPTION(utility::excn::Exception, "Code error: params.nRsteps_ cannot be 0 when passed into"
			" densiy_grid_search.  Try using get_spectrum() to return nRsteps before calling this function");
	}
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

	std::stringstream outresults;
	outresults << "index" << "\t" << "rank" << "\t" << "score" << "\t" << "rot_xx" << "\t" <<  "rot_xy" << "\t" << "rot_xz"
		<< "\t" << "rot_yx" << "\t" << "rot_yy" << "\t" << "rot_yz"
		<< "\t" << "rot_zx" << "\t" << "rot_zy" << "\t" << "rot_zz"
		<< "\t" << "pre_x" << "\t" << "pre_y" << "\t" << "pre_z"
		<< "\t" << "post_x" << "\t" << "post_y" << "\t" << "post_z" << std::endl;

	TR.Debug << "Beginning point search from " << params.point_search_start_ << " to " << params.point_search_end_ << " now." << std::endl;

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

		RBfitResultDB local_results( nperRot );

		for ( core::Size j=0; j<8*params.B_*params.B_*params.B_; ++j ) {
			core::Real const sumSigMap = so3_correl[j] / (4*numeric::constants::d::pi);
			core::Real const sumMap = mask_correl[j] / (4*numeric::constants::d::pi);
			core::Real const sumMap2 = mask2_correl[j] / (4*numeric::constants::d::pi);

			core::Real const CC = (
				(sumEps*sumSigMap - sumSig*sumMap) / (
				std::sqrt( sumEps*sumSig2 - sumSig*sumSig )
				* std::sqrt( sumEps*sumMap2 - sumMap*sumMap )
				) );

			if ( local_results.to_add_element( CC ) ) {
				SOFT.idx_to_rot(j , rot);
				local_results.add_element( RBfitResult( pose_idx, CC, rot, pretrans, posttrans ) );
			}
		}

		cluster_RBfitResultDB_fast(
			local_results,
			params.delRSteps_,
			params.nRsteps_,
			params.cluster_radius_,
			params.max_rot_per_trans_,
			params.include_distance_during_fast_cluster_,
			density);

		while ( local_results.size() > params.max_rot_per_trans_ ) local_results.pop();

		if ( !params.output_fn_.empty() ) write_RBfitResultDB( local_results, outresults );

		while ( local_results.size() > 0 ) {
			RBfitResult const sol_i = local_results.pop();

			if ( !params.natives_.empty() ) {
				for ( core::Size native_i = 1; native_i <= params.natives_.size(); ++native_i ) {
					core::Real const distNative = ( posttrans - params.native_coms_[ native_i ] ).length();
					if ( distNative < params.rms_cutoff_ ) {
						TR << "[" << i << "/" << points_to_search.size() << "] COM dist to native " << distNative << " at com_idx " << native_i << std::endl;
						compare_RBfitDB_to_native(
							local_results,
							pose,
							params.natives_,
							params.native_coms_,
							params.native_middle_cas_,
							params.symminfo_,
							params.center_on_middle_ca_,
							params.rms_cutoff_
						);
					}
				}
			}

			results.add_element( sol_i );
		}
		if ( i%100 == 0 ) {
			TR << "[" << i << "/" << points_to_search.size() << "] " << results.top().score_ << " / " << results.size() << std::endl;
		}
	}

	if ( !params.output_fn_.empty() ) {
		std::ofstream out_file( params.output_fn_.c_str() );
		out_file << outresults.str() << std::endl;
	}
}


} // electron_density
} // protocols
