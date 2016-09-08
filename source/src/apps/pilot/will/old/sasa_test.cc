// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /src/apps/pilat/will/coiled_coil.cc
/// @brief samples coiled coils

#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/sasa.hh>
#include <numeric/random/random.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/prof.hh>
#include <core/chemical/ResidueType.hh>
#include <numeric/xyz.functions.hh>
#include <basic/options/option.hh>
#include <core/scoring/packing/surf_vol.hh>
#include <core/scoring/ScoreFunction.hh>


#include <time.h>

//Auto Headers
#include <core/import_pose/import_pose.hh>


using core::Size;
using core::Real;
using numeric::random::uniform;
using numeric::random::gaussian;
using core::pose::Pose;
typedef numeric::xyzVector<Real> Vec;
typedef utility::vector1<Vec>    Vecs;
typedef numeric::xyzMatrix<Real> Mat;

inline void rot_pose( core::pose::Pose & pose, Mat const & rot ) {
	for(Size ir = 1; ir <= pose.size(); ++ir) {
		for(Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
			core::id::AtomID const aid(core::id::AtomID(ia,ir));
			pose.set_xyz( aid, rot * pose.xyz(aid) );
		}
	}
}

// inline void rot_pose( core::pose::Pose & pose, Mat const & rot, Vec const & cen ) {
// 	trans_pose(pose,-cen);
// 	rot_pose(pose,rot);
// 	trans_pose(pose,cen);
// }

inline void rot_pose( core::pose::Pose & pose, Vec const & axis, Real const & ang ) {
	rot_pose(pose,rotation_matrix_degrees(axis,ang));
}

// inline void rot_pose( core::pose::Pose & pose, Vec const & axis, Real const & ang, Vec const & cen ) {
// 	rot_pose(pose,rotation_matrix_degrees(axis,ang),cen);
// }


Real dot_sasa(Pose & pose, clock_t & t) {
	PROF_START( basic::STAGE1 );
	t = clock();
	Real s = core::scoring::calc_total_sasa(pose,1.4);
	t = clock()-t;
	PROF_STOP( basic::STAGE1 );
	return s;
}

Real dab_sasa(Pose & pose, clock_t & t){
	core::scoring::ScoreFunction sf;
	sf.set_weight(core::scoring::dab_sasa,1.0);
	PROF_START( basic::STAGE2 );
	t = clock();
	// Real s = sf(pose); // use to profile score function
	Real s = core::scoring::packing::get_surf_tot(pose,1.4); // use to profile direct call
	t = clock()-t;
	PROF_STOP( basic::STAGE2 );
	return s;
}


int main (int argc, char *argv[])
{

	try {


	devel::init(argc,argv);

	Pose pose;
	core::import_pose::pose_from_file(pose,core::options::option[core::options::OptionKeys::in::file::s]()[1], core::import_pose::PDB_file);

	clock_t tot1=0,tot2=0;
	for(Size i = 1; i <= 10; ++i){
		rot_pose(pose,numeric::xyzVector<Real>(uniform(),uniform(),uniform()),uniform()*360);
		// for(Size j = 1; j <= pose.size(); ++j) {
		// 	pose.set_phi(j,pose.phi(j)+gaussian()/10.0);
		// 	pose.set_psi(j,pose.psi(j)+gaussian()/10.0);
		// }
		clock_t t1,t2;
		Real SASAdab = dab_sasa(pose,t1);
		Real SASAdot = dot_sasa(pose,t2);
		tot1 += t1; tot2 += t2;
		std::cerr << "DOTS vs DAB: dab/dot " << (Real)tot1/(Real)tot2 << " " << (Real)t1/(Real)t2 << " " << t1 << " " << t2 << " values: " << SASAdab << " " << SASAdot << std::endl;
	}
	std::cerr << "PROFILE info:" << std::endl;
	basic::prof_show();


	return 0;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}


