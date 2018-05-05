// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/io/CrystInfo.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>

#include <protocols/cryst/spacegroup.hh>
#include <protocols/cryst/refinable_lattice.hh>

#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/mpistream.hh>
#include <utility/string_util.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>

#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/format.hh>

#include <utility/excn/Exceptions.hh>

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/option_macros.hh>
#include <basic/basic.hh>
#include <basic/database/open.hh>

// option includes
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <fstream>
#include <iostream>
#include <cmath>

#include <sstream>
#include <string>
#include <queue>

#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.io.hh>
#include <numeric/random/random.hh>
#include <numeric/NumericTraits.hh>
#include <numeric/fourier/FFT.hh>
#include <numeric/constants.hh>
#include <numeric/model_quality/rms.hh>

namespace protocols {
namespace cryst {


core::Size
get_nres_asu( core::pose::Pose const & pose ) {
	core::Size nres_tgt = pose.size();
	core::conformation::symmetry::SymmetryInfoCOP symm_info;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		auto const & SymmConf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation const &> ( pose.conformation()) );
		symm_info = SymmConf.Symmetry_Info();
		nres_tgt = symm_info->num_independent_residues();
	}
	if ( pose.residue(nres_tgt).aa() == core::chemical::aa_vrt ) nres_tgt--;
	while ( !pose.residue(nres_tgt).is_protein() ) nres_tgt--;
	return nres_tgt;
}

core::Real
crystRMSfast (core::pose::Pose &pose_native, core::pose::Pose &pose_decoy) {
	core::Real sumdistfast=0;

	core::Size nres_native = get_nres_asu(pose_native);
	core::Size nres_decoy  = get_nres_asu(pose_decoy);

	runtime_assert( nres_native==nres_decoy );

	numeric::xyzVector< core::Real > com_native(0,0,0), com_decoy(0,0,0);
	core::Real radius_native=0, radius_decoy=0;
	for ( core::Size i=1; i<= nres_native; ++i ) com_native += pose_native.residue(i).xyz(2);
	com_native /= nres_native;
	for ( core::Size i=1; i<= nres_native; ++i ) radius_native = std::max( radius_native, (pose_native.residue(i).xyz(2)-com_native).length() );
	for ( core::Size i=1; i<= nres_decoy; ++i ) com_decoy += pose_decoy.residue(i).xyz(2);
	com_decoy /= nres_decoy;
	for ( core::Size i=1; i<= nres_decoy; ++i ) radius_decoy = std::max( radius_decoy, (pose_decoy.residue(i).xyz(2)-com_decoy).length() );

	// apply symmops
	core::io::CrystInfo ci_native = pose_native.pdb_info()->crystinfo();
	protocols::cryst::Spacegroup sg_native;
	sg_native.set_spacegroup(ci_native.spacegroup());
	sg_native.set_parameters(ci_native.A(),ci_native.B(),ci_native.C(), ci_native.alpha(), ci_native.beta(), ci_native.gamma());

	utility::vector1< numeric::xyzVector<core::Real> > coms_native;
	coms_native.push_back( com_native );

	for ( int s=1; s<=(int)sg_native.nsymmops(); ++s ) {
		numeric::xyzMatrix<core::Real> R_i = sg_native.symmop(s).get_rotation();
		numeric::xyzVector<core::Real> T_i = sg_native.symmop(s).get_translation();
		for ( int i=-(int)2; i<=(int)2; ++i ) {
			for ( int j=-(int)2; j<=(int)2; ++j ) {
				for ( int k=-(int)2; k<=(int)2; ++k ) {
					if ( s==1 && i==0 && j==0 && k==0 ) continue;
					numeric::xyzVector<core::Real> com_ijk = sg_native.f2c() * (R_i*(sg_native.c2f()*com_native) + T_i + numeric::xyzVector<core::Real>(i,j,k));
					core::Size distance = (com_ijk-com_native).length();
					if ( distance<radius_native + 18 ) {
						coms_native.push_back(com_ijk);
					}
				}
			}
		}
	}

	core::io::CrystInfo ci_decoy = pose_decoy.pdb_info()->crystinfo();
	protocols::cryst::Spacegroup sg_decoy;
	sg_decoy.set_spacegroup(ci_decoy.spacegroup());
	sg_decoy.set_parameters(ci_decoy.A(),ci_decoy.B(),ci_decoy.C(), ci_decoy.alpha(), ci_decoy.beta(), ci_decoy.gamma());

	utility::vector1< numeric::xyzVector<core::Real> > coms_decoy;
	coms_decoy.push_back( com_decoy );

	for ( int s=1; s<=(int)sg_decoy.nsymmops(); ++s ) {
		numeric::xyzMatrix<core::Real> R_i = sg_decoy.symmop(s).get_rotation();
		numeric::xyzVector<core::Real> T_i = sg_decoy.symmop(s).get_translation();
		for ( int i=-(int)2; i<=(int)2; ++i ) {
			for ( int j=-(int)2; j<=(int)2; ++j ) {
				for ( int k=-(int)2; k<=(int)2; ++k ) {
					if ( s==1 && i==0 && j==0 && k==0 ) continue;
					numeric::xyzVector<core::Real> com_ijk = sg_decoy.f2c() * (R_i*(sg_decoy.c2f()*com_decoy) + T_i + numeric::xyzVector<core::Real>(i,j,k));
					core::Size distance = (com_ijk-com_decoy).length();
					if ( distance<radius_decoy + 12 ) {
						coms_decoy.push_back(com_ijk);
					}
				}
			}
		}
	}

	// 1B - get chainA->chainA transformation
	numeric::xyzMatrix< core::Real > R;
	{
		ObjexxFCL::FArray2D< core::Real > init_coords( 3, nres_native ), final_coords( 3, nres_native );
		ObjexxFCL::FArray1D< core::Real > ww( nres_native, 1.0 );
		ObjexxFCL::FArray2D< core::Real > uu( 3, 3, 0.0 );
		core::Real ctx;
		for ( core::Size i=1; i<=nres_native; ++i ) {
			if ( pose_decoy.residue(i).is_protein() ) {
				numeric::xyzVector< core::Real > const &x_decoy = pose_decoy.residue(i).atom(2).xyz();
				numeric::xyzVector< core::Real > const &x_native = pose_native.residue(i).atom(2).xyz();
				init_coords(1,i) = x_decoy[0]; init_coords(2,i) = x_decoy[1]; init_coords(3,i) = x_decoy[2];
				final_coords(1,i) = x_native[0]; final_coords(2,i) = x_native[1]; final_coords(3,i) = x_native[2];
			}
		}
		numeric::model_quality::findUU( final_coords, init_coords, ww, nres_native, uu, ctx );
		R.xx( uu(1,1) ); R.xy( uu(2,1) ); R.xz( uu(3,1) );
		R.yx( uu(1,2) ); R.yy( uu(2,2) ); R.yz( uu(3,2) );
		R.zx( uu(1,3) ); R.zy( uu(2,3) ); R.zz( uu(3,3) );
	}

	// 2 get decoy->native correspondence
	for ( core::Size j=2; j<=coms_decoy.size(); ++j ) {
		core::Real mindist = 1e6;
		numeric::xyzVector< core::Real > decoy_j_off = coms_decoy[j]-coms_decoy[1];
		for ( core::Size k=2; k<=coms_native.size(); ++k ) {
			numeric::xyzVector< core::Real > native_k_off = coms_native[k]-coms_native[1];
			core::Real dist = (native_k_off - R*decoy_j_off).length();
			if ( dist < mindist ) { mindist=dist; }
		}
		sumdistfast += mindist*mindist;
	}
	sumdistfast = std::sqrt( sumdistfast / (coms_decoy.size()) );
	return sumdistfast;
}

core::Real
crystRMS (core::pose::Pose &pose_native, core::pose::Pose &pose_decoy, bool allatom) {
	using namespace core;
	using namespace core::scoring;
	using namespace core::conformation::symmetry;

	protocols::cryst::MakeLatticeMover setup_native, setup_decoy;

	///
	///
	if ( !core::pose::symmetry::is_symmetric( pose_native ) ) {
		setup_native.contact_dist(16);
		setup_native.apply( pose_native );
	}
	if ( !core::pose::symmetry::is_symmetric( pose_decoy ) ) {
		setup_decoy.contact_dist(10);
		setup_decoy.apply( pose_decoy );
	}

	SymmetryInfoOP symm_info_native = (dynamic_cast<SymmetricConformation &> ( pose_native.conformation())).Symmetry_Info();
	core::Size nres_asu = symm_info_native->num_independent_residues();
	core::Size nsubunits_native = symm_info_native->subunits();
	SymmetryInfoOP symm_info_decoy = (dynamic_cast<SymmetricConformation &> ( pose_decoy.conformation())).Symmetry_Info();
	core::Size nres_asu_decoy = symm_info_decoy->num_independent_residues();
	core::Size nsubunits_decoy = symm_info_decoy->subunits();
	runtime_assert( nres_asu == nres_asu_decoy);

	// 1 get CoM of all subunits
	utility::vector1< core::Vector > coms_native(nsubunits_native, core::Vector(0,0,0) );
	for ( core::Size j=1; j<=nsubunits_native; ++j ) {
		core::Size count=0;
		for ( core::Size i=1; i<=nres_asu; ++i ) {
			core::Size residx = (j-1)*nres_asu+i;
			core::conformation::Residue const &rsd = pose_native.residue(residx);
			if ( allatom ) {
				for ( core::Size k=1; k<=rsd.nheavyatoms(); ++k ) {
					coms_native[j] += rsd.atom(k).xyz();
					count++;
				}
			} else {
				if ( rsd.is_protein() ) {
					coms_native[j] += rsd.atom(2).xyz();
					count++;
				}
			}
		}
		coms_native[j] /= count;
	}

	utility::vector1< core::Vector > coms_decoy(nsubunits_decoy, core::Vector(0,0,0) );
	for ( core::Size j=1; j<=nsubunits_decoy; ++j ) {
		core::Size count=0;
		for ( core::Size i=1; i<=nres_asu; ++i ) {
			core::Size residx = (j-1)*nres_asu+i;
			core::conformation::Residue const &rsd = pose_decoy.residue(residx);
			if ( allatom ) {
				for ( core::Size k=1; k<=rsd.nheavyatoms(); ++k ) {
					coms_decoy[j] += rsd.atom(k).xyz();
					count++;
				}
			} else {
				if ( rsd.is_protein() ) {
					coms_decoy[j] += rsd.atom(2).xyz();
					count++;
				}
			}
		}
		coms_decoy[j] /= count;
	}

	// get natoms in asu
	core::Size nres_aln = 0;
	if ( allatom ) {
		for ( core::Size i=1; i<=nres_asu; ++i ) {
			core::conformation::Residue const &rsd = pose_native.residue(i);
			nres_aln += rsd.nheavyatoms();
		}
	} else {
		for ( core::Size i=1; i<=nres_asu; ++i ) {
			core::conformation::Residue const &rsd = pose_native.residue(i);
			if ( rsd.is_protein() ) nres_aln++;
		}
	}

	// 1B - get chainA->chainA transformation
	numeric::xyzMatrix< core::Real > R;
	{
		ObjexxFCL::FArray2D< core::Real > init_coords( 3, nres_aln ), final_coords( 3, nres_aln );
		ObjexxFCL::FArray1D< core::Real > ww( nres_aln, 1.0 );
		ObjexxFCL::FArray2D< core::Real > uu( 3, 3, 0.0 );
		numeric::Real ctx;
		core::Size idx=1;
		for ( core::Size i=1; i<=nres_asu; ++i ) {
			if ( allatom ) {
				if ( pose_decoy.residue(i).is_protein() ) {
					Vector const &x_decoy = pose_decoy.residue(i).atom(2).xyz();
					Vector const &x_native = pose_native.residue(i).atom(2).xyz();
					init_coords(1,idx) = x_decoy[0]; init_coords(2,idx) = x_decoy[1]; init_coords(3,idx) = x_decoy[2];
					final_coords(1,idx) = x_native[0]; final_coords(2,idx) = x_native[1]; final_coords(3,idx) = x_native[2];
					idx++;
				}
			} else {
				core::conformation::Residue const &rsd_n = pose_native.residue(i);
				core::conformation::Residue const &rsd_d = pose_decoy.residue(i);
				runtime_assert( rsd_n.nheavyatoms() == rsd_d.nheavyatoms() );
				for ( core::Size k=1; k<=rsd_n.nheavyatoms(); ++k ) {
					Vector const &x_decoy = rsd_d.atom(k).xyz();
					Vector const &x_native = rsd_n.atom(k).xyz();
					init_coords(1,idx) = x_decoy[0]; init_coords(2,idx) = x_decoy[1]; init_coords(3,idx) = x_decoy[2];
					final_coords(1,idx) = x_native[0]; final_coords(2,idx) = x_native[1]; final_coords(3,idx) = x_native[2];
					idx++;
				}
			}
		}
		numeric::model_quality::findUU( final_coords, init_coords, ww, nres_asu, uu, ctx );
		R.xx( uu(1,1) ); R.xy( uu(2,1) ); R.xz( uu(3,1) );
		R.yx( uu(1,2) ); R.yy( uu(2,2) ); R.yz( uu(3,2) );
		R.zx( uu(1,3) ); R.zy( uu(2,3) ); R.zz( uu(3,3) );
	}

	// 2 get decoy->native correspondence
	utility::vector1<core::Size> map_decoy2native(nsubunits_decoy,1); // 1->1 by construction
	core::Real sumdist = 0.0;
	for ( core::Size j=2; j<=nsubunits_decoy; ++j ) {
		Real mindist = 1e6;
		Size minidx = 0;
		Vector decoy_j_off = coms_decoy[j]-coms_decoy[1];
		for ( core::Size k=2; k<=nsubunits_native; ++k ) {
			Vector native_k_off = coms_native[k]-coms_native[1];
			Real dist = (native_k_off - R*decoy_j_off).length();
			if ( dist < mindist ) { mindist=dist; minidx=k; }
		}
		sumdist += mindist*mindist;
		map_decoy2native[j] = minidx;
	}

	// 3 construct corresponding CA arrays and compute RMS
	float rms;
	{
		ObjexxFCL::FArray2D< core::Real > init_coords( 3, nsubunits_decoy*nres_aln ), final_coords( 3, nsubunits_decoy*nres_aln );
		core::Size idx=1;
		for ( core::Size j=1; j<=nsubunits_decoy; ++j ) {
			core::Size resdecoy = (j-1)*nres_asu;
			core::Size resnative = (map_decoy2native[j]-1)*nres_asu;
			for ( core::Size i=1; i<=nres_asu; ++i ) {
				if ( allatom ) {
					core::conformation::Residue const &rsd_n = pose_native.residue(resnative+i);
					core::conformation::Residue const &rsd_d = pose_decoy.residue(resdecoy+i);
					for ( core::Size k=1; k<=rsd_n.nheavyatoms(); ++k ) {
						Vector const &x_decoy = rsd_d.atom(k).xyz();
						Vector const &x_native = rsd_n.atom(k).xyz();
						init_coords(1,idx) = x_decoy[0]; init_coords(2,idx) = x_decoy[1]; init_coords(3,idx) = x_decoy[2];
						final_coords(1,idx) = x_native[0]; final_coords(2,idx) = x_native[1]; final_coords(3,idx) = x_native[2];
						idx++;
					}
				} else {
					if ( pose_decoy.residue(resdecoy).is_protein() ) {
						Vector const &x_decoy = pose_decoy.residue(resdecoy+i).atom(2).xyz();
						Vector const &x_native = pose_native.residue(resnative+i).atom(2).xyz();
						init_coords(1,idx) = x_decoy[0]; init_coords(2,idx) = x_decoy[1]; init_coords(3,idx) = x_decoy[2];
						final_coords(1,idx) = x_native[0]; final_coords(2,idx) = x_native[1]; final_coords(3,idx) = x_native[2];
						idx++;
					}
				}
			}
		}

		rms = numeric::model_quality::rms_wrapper( nsubunits_decoy*nres_aln, final_coords, init_coords );
	}

	return rms;
}

}
}

