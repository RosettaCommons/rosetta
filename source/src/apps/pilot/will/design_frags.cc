// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief

#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <devel/init.hh>
#include <core/import_pose/import_pose.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/kinematics/Stub.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyz.functions.hh>
#include <utility/string_util.hh>

#include <utility/vector1.hh>


// core::kinematics::Stub getxform(core::conformation::Residue const & downstream_resi, core::conformation::Residue const & upstream_resi) {
// 	core::kinematics::Stub s;
// 	numeric::xyzVector<Real> A1 = downstream_resi.xyz("N" ) - downstream_resi.xyz("C");
// 	numeric::xyzVector<Real> B1 = downstream_resi.xyz("CA") - downstream_resi.xyz("C");
// 	numeric::xyzVector<Real> A2 =   upstream_resi.xyz("C" ) -   upstream_resi.xyz("N");
// 	numeric::xyzVector<Real> B2 =   upstream_resi.xyz("CA") -   upstream_resi.xyz("N");
// 	s.M = alignVectorSets( A1,B1, A2,B2 );
// 	s.v = upstream_resi.xyz("N")-s.M*downstream_resi.xyz("C");
// 	return s;
// }

using core::Real;
using core::Size;

struct Quat {
	float w,x,y,z;
	Quat(numeric::xyzMatrix<Real> M) {
		Real r = 2.0 * std::sqrt( 1 + M.xx() - M.yy() - M.zz() );
		w = (M.zy()-M.yz()) / r;
		x = r / 4.0;
		y = (M.xy()+M.yx()) / r;
		z = (M.zx()+M.xz()) / r;
	}
};

core::kinematics::Stub getxform(core::conformation::Residue const & move_resi, core::conformation::Residue const & fixd_resi) {
	core::kinematics::Stub s;
	s.M = alignVectorSets(move_resi.xyz(1)-move_resi.xyz(2),move_resi.xyz(3)-move_resi.xyz(2),fixd_resi.xyz(1)-fixd_resi.xyz(2),fixd_resi.xyz(3)-fixd_resi.xyz(2));
	s.v = fixd_resi.xyz(2)-s.M*move_resi.xyz(2);
	return s;
}



void pose2frags(core::pose::Pose const & pose, Size minlen, Size maxlen, std::string tag) {
	for(Size i = 1; i <= pose.n_residue()-minlen; ++i) {
		for(Size l = minlen; l <= maxlen; ++l) {
			Size j = i+l+1;
			if(j > pose.n_residue()) continue;
			core::kinematics::Stub s = getxform( pose.residue(i), pose.residue(j) );
			Quat q(s.M);
			if( fabs(std::sqrt(q.w*q.w+q.x*q.x+q.y*q.y+q.z*q.z) - 1.0) > 0.00001 ) utility_exit_with_message("quaternion is not lenght 1!");
			std::cout << "loop " << tag << " " << i << " " << j << "    " << s.v.x() << " " << s.v.y() << " " << s.v.z() << "    " << q.w << " " << q.x << " " << q.y << " " << q.z << std::endl;
			Real e1 = -asin(s.M.zx());
			Real e2 = atan2(s.M.zy(),s.M.zz());
			Real e3 = atan2(s.M.yx()/cos(e1),s.M.xx()/cos(e1));
			std::cout << "euler " << e1 << " " << e2 << " " << e3 << std::endl;
		}
	}
}

int
main (int argc, char *argv[]) {

	try {

	using namespace basic::options;

	devel::init( argc, argv );


	utility::vector1<std::string> files = option[OptionKeys::in::file::s]();
	for(Size i = 1; i <= files.size(); ++i) {
		core::pose::Pose pose;
		core::import_pose::pose_from_pdb(pose,files[i]);
		std::string tag = utility::file_basename(files[i]).substr(3,4);
		pose2frags(pose,1,20,tag);
	}


	return 0;



	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
