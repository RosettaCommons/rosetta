// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:

#include <numeric/geometry/hashing/xyzStripeHashWithMeta.hh>
#include <protocols/sic_dock/xyzStripeHashPoseWithMeta.hh>
//#include <apps/pilot/will/gpu/gpu_refold.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
//#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <devel/init.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.io.hh>
#include <core/scoring/sasa.hh>
#include <fstream>

#include <core/kinematics/Stub.hh>

#include <apps/pilot/will/SphereSample.hh>

using core::Real;
using core::Size;
using core::id::AtomID;
typedef numeric::xyzVector<core::Real> Vec;

OPT_KEY( Real, clash_dis )
void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	OPT( in::file::s );
	NEW_OPT( clash_dis, "clash dist", 3.5 );
}


utility::vector1<Vec>
get_surface_points(core::pose::Pose const & p, Real /*DIST*/) {
	SphereSample ss;
	ss.pdb_from_level(1,"test1.pdb");
	ss.pdb_from_level(2,"test2.pdb");
	ss.pdb_from_level(3,"test3.pdb");
	ss.pdb_from_level(4,"test4.pdb");
	ss.pdb_from_level(5,"test5.pdb");
	ss.pdb_from_level(6,"test6.pdb");
	ss.pdb_from_level(7,"test7.pdb");						

	utility::vector1<Vec>	points;
	for(Size ir = 1; ir <= p.n_residue(); ++ir) {
		for(Size ia = 1; ia <= p.residue(ir).nheavyatoms(); ++ia){
			points.push_back(p.xyz(AtomID(ia,ir)));
		}
	}
	
//	for(Size )
	return points;
}


int main(int argc, char *argv[]) {

	try {

	register_options();
	devel::init(argc,argv);
	using namespace basic::options;


	Real const DIST(option[OptionKeys::clash_dis]());
	//Real const DIST2(DIST*DIST);  // unused ~Labonte

	core::pose::Pose p;
	core::import_pose::pose_from_pdb(p,option[OptionKeys::in::file::s]()[1]);
	if(false) {
		for(Size ir = 1; ir <= p.n_residue(); ++ir) {
			if( p.residue(ir).is_lower_terminus() ) core::pose::remove_lower_terminus_type_from_pose_residue(p,ir);
			if( p.residue(ir).is_upper_terminus() ) core::pose::remove_upper_terminus_type_from_pose_residue(p,ir);		
		}
	}
	
	utility::vector1<Vec> points = get_surface_points(p,DIST);
	
	utility_exit_with_message("TEST GET SURF PTS");
	
	protocols::sic_dock::xyzStripeHashPoseWithMeta xyzhash(DIST,p,protocols::sic_dock::ALL);
	xyzhash.sanity_check();

	using namespace ObjexxFCL::format;
	std::cout << xyzhash.grid_size() << std::endl;
	std::cout << xyzhash.natom() << std::endl;
	std::cout << xyzhash.xdim() << " " << xyzhash.ydim() << " " << xyzhash.zdim() << std::endl;
	for(Size j = 0; j < (Size)xyzhash.natom(); ++j) {
		std::cout << F(12,7,xyzhash.grid_atoms()[j].x) << " ";
		std::cout << F(12,7,xyzhash.grid_atoms()[j].y) << " ";
		std::cout << F(12,7,xyzhash.grid_atoms()[j].z) << std::endl;
	}
	for(Size j = 0; j < Size(xyzhash.xdim()*xyzhash.ydim()*xyzhash.zdim()); ++j) {
		std::cout << xyzhash.grid_stripe()[j].x << " ";
		std::cout << xyzhash.grid_stripe()[j].y << std::endl;
	}
	return 0;


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

}
