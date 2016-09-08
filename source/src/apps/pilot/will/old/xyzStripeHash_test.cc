// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#include <numeric/geometry/hashing/xyzStripeHash.hh>
#include <protocols/sic_dock/types.hh>
#include <protocols/sic_dock/xyzStripeHashPose.hh>
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
#include <numeric/random/random.hh>
#include <numeric/xyz.io.hh>
#include <fstream>
#include <protocols/scoring/ImplicitFastClashCheck.hh>

#include <core/kinematics/Stub.hh>

#include <core/conformation/AtomGraphData.hh>
#include <core/conformation/AtomGraph.hh>
#include <core/conformation/find_neighbors.hh>

#include <set>

using core::Size;
using core::id::AtomID;
typedef numeric::xyzVector<float> Vec;

#include <time.h>
#include <sys/time.h>
#include <stdio.h>
#ifdef __MACH__
#include <mach/mach_time.h>
#endif
double time_highres() {
#ifdef __MACH__
  mach_timebase_info_data_t info;
  mach_timebase_info(&info);
  return mach_absolute_time() / 1000000000.0;
  //uint64_t duration = mach_absolute_time();
  //duration *= info.numer;
  //duration /= info.denom;
#else
  timespec tp;
  clock_gettime(CLOCK_REALTIME, &tp);
  return tp.tv_sec + tp.tv_nsec/1000000000.0;
#endif
  return 0;
}

OPT_KEY( Boolean, dump_hash )
void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	OPT( in::file::s );
	NEW_OPT( dump_hash,"dump hash data and quit", false );
}

int slow_nbcount( protocols::sic_dock::xyzStripeHashPose const & xyzhash, float DIST, Vec const & rv ){
	Vec tmp = rv+xyzhash.translation();
	int safe_nbcount = 0;
	for(int j = 0; j < (int)xyzhash.natom(); ++j) {
		float const & hx = xyzhash.grid_atoms()[j].x();
		float const & hy = xyzhash.grid_atoms()[j].y();
		float const & hz = xyzhash.grid_atoms()[j].z();
		if( (hx-tmp.x())*(hx-tmp.x()) + (hy-tmp.y())*(hy-tmp.y()) + (hz-tmp.z())*(hz-tmp.z()) <= DIST*DIST ) safe_nbcount++;
	}
	return safe_nbcount;
}

int slow_nbcount( core::pose::Pose const & pose, float DIST, Vec const & rv ){
	int safe_nbcount = 0;
	for(Size ir = 1; ir <= pose.size(); ++ir) {
		for(Size ia = 1; ia <= pose.residue(ir).natoms(); ++ia){
			if( rv.distance_squared(pose.residue(ir).xyz(ia)) <= DIST*DIST ){
				++safe_nbcount;
			}
		}
	}
	return safe_nbcount;
}

void slow_nbset( core::pose::Pose const & pose, float DIST, Vec const & rv, std::set<int> & nbrs ){
	for(Size ir = 1; ir <= pose.size(); ++ir) {
		if( rv.distance_squared(pose.residue(ir).xyz("CA")) <= DIST*DIST ){
			nbrs.insert(ir);
		}
	}
}

bool slow_clash( core::pose::Pose const & pose, float DIST, Vec const & rv ){
	int safe_nbcount = 0;
	for(Size ir = 1; ir <= pose.size(); ++ir) {
		for(Size ia = 1; ia <= pose.residue(ir).natoms(); ++ia){
			if( rv.distance_squared(pose.residue(ir).xyz(ia)) <= DIST*DIST ){
				return true;
			}
		}
	}
	return false;
}


utility::vector1<Vec>
slow_nbget( protocols::sic_dock::xyzStripeHashPose const & xyzhash, float DIST, Vec const & rv ){
	utility::vector1<Vec> nbrs;
	Vec tmp = rv+xyzhash.translation();
	int safe_nbcount = 0;
	for(int j = 0; j < (int)xyzhash.natom(); ++j) {
		Vec const & c = *((Vec*)(xyzhash.grid_atoms()+j));
		if( c.distance_squared(tmp) <= DIST*DIST ){
			nbrs.push_back( c );
		}
	}
	return nbrs;
}

void dump_points_pdb(utility::vector1<Vec> & p, std::string fn) {
	using namespace ObjexxFCL::format;
	std::ofstream o(fn.c_str());
	for(Size i = 1; i <= p.size(); ++i) {
		std::string rn = "VIZ";
		o<<"HETATM"<<I(5,i)<<' '<<" CA "<<' '<<rn<<' '<<"A"<<I(4,i)<<"    "<<F(8,3,p[i].x())<<F(8,3,p[i].y())<<F(8,3,p[i].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
	}
	o.close();
}

struct CollectNeighbors {
	utility::vector1<Vec> nbrs;
	void visit(Vec const &, Vec const & c, float const &){
		nbrs.push_back(c);
	}
};

void
test_find_neighbors(core::pose::Pose const & p){
	float t,to=0,tg=0,ts=0;
	{
		std::cout << "stripe" << std::endl;
		core::conformation::AtomGraphOP pg = new core::conformation::AtomGraph;
		core::conformation::atom_graph_from_conformation(p.conformation(),pg);
		t = time_highres();
		find_neighbors_stripe(pg,6.0);
		ts += time_highres()-t;
	}
	{
		std::cout << "octree" << std::endl;
		core::conformation::AtomGraphOP pg = new core::conformation::AtomGraph;
		core::conformation::atom_graph_from_conformation(p.conformation(),pg);
		t = time_highres();
		find_neighbors_octree(pg,6.0,core::conformation::AUTOMATIC);
		to += time_highres()-t;
	}
	{
		std::cout << "3dgrid" << std::endl;
		core::conformation::AtomGraphOP pg = new core::conformation::AtomGraph;
		core::conformation::atom_graph_from_conformation(p.conformation(),pg);
		t = time_highres();
		find_neighbors_3dgrid(pg,6.0);
		tg += time_highres()-t;
	}
	std::cout << "test find_neighbors " << to << " " << tg << " " << ts << std::endl;
}

int main(int argc, char *argv[]) {

 try {

	using numeric::geometry::hashing::xyzStripeHash;
	register_options();
	devel::init(argc,argv);

	std::cout << "sizeof(utility::vector1< xyzVector<float> >) "
	          << sizeof(utility::vector1< numeric::xyzVector<float> >) << " "
	          << sizeof(xyzStripeHash::ushort2)
	          << std::endl;

	float const DIST(6);

	core::pose::Pose p;
	core::import_pose::pose_from_file(p,basic::options::option[basic::options::OptionKeys::in::file::s]()[1], core::import_pose::PDB_file);
	if(true) {
		for(Size ir = 1; ir <= p.size(); ++ir) {
			if( p.residue(ir).is_lower_terminus() ) core::pose::remove_lower_terminus_type_from_pose_residue(p,ir);
			if( p.residue(ir).is_upper_terminus() ) core::pose::remove_upper_terminus_type_from_pose_residue(p,ir);
		}
	}

	test_find_neighbors(p);
	// utility_exit_with_message("test find_neighbors");

	protocols::scoring::ImplicitFastClashCheck ifc(p,DIST);

	std::cout << "hash mode CA" << std::endl;
	protocols::sic_dock::xyzStripeHashPose xyzhash(p,protocols::sic_dock::ALL,DIST);
	xyzhash.sanity_check();

	// xyzStripeHash<float>::float4 const * j = xyzhash.grid_atoms_;
	// for(xyzStripeHash<float>::xyz_iterator i = xyzhash.xyz_begin(); i != xyzhash.xyz_end(); ++i){
	// 	numeric::xyzVector<float> const & a = *i;
	// 	xyzStripeHash<float>::float4 const & b = *(j++);
	// 	std::cout << a.x() <<" "<< b.x <<" "<< a.y() <<" "<< b.y <<" "<< a.z() <<" "<< b.z << std::endl;
	// 	if( a.x() != b.x && a.y() != b.y && a.z() != b.z ){
	// 		utility_exit_with_message("BAD!!!!!");
	// 	}
	// }

	// utility_exit_with_message("FOO");

	if(basic::options::option[basic::options::OptionKeys::dump_hash]()) {
		using namespace ObjexxFCL::format;
		std::cout << xyzhash.grid_size() << std::endl;
		std::cout << xyzhash.natom() << std::endl;
		std::cout << xyzhash.xdim() << " " << xyzhash.ydim() << " " << xyzhash.zdim() << std::endl;
		for(int j = 0; j < (int)xyzhash.natom(); ++j) {
			std::cout << F(12,7,xyzhash.grid_atoms()[j].x()) << " ";
			std::cout << F(12,7,xyzhash.grid_atoms()[j].y()) << " ";
			std::cout << F(12,7,xyzhash.grid_atoms()[j].z()) << std::endl;
		}
		for(int j = 0; j < (int)xyzhash.xdim()*xyzhash.ydim()*xyzhash.zdim(); ++j) {
			std::cout << xyzhash.grid_stripe()[j].x << " ";
			std::cout << xyzhash.grid_stripe()[j].y << std::endl;
		}
		return 0;
	}

	// for(Size ir = 1; ir <= p.size(); ++ir) {
	// 	for(Size ia = 1; ia <= p.residue_type(ir).natoms(); ++ia) {
	// 		p.set_xyz(AtomID(ia,ir),p.xyz(AtomID(ia,ir)) + (numeric::xyzVector<core::Real>)xyzhash.translation() );
	// 	}
	// }
	// p.dump_pdb("input_trans.pdb");


	Vec mn(9e9),mx(-9e9);
	int real_natom = 0.0;
	for(Size ir = 1; ir <= p.size(); ++ir) {
		for(Size ia = 1; ia <= p.residue(ir).natoms(); ++ia) {
			real_natom++;
			// mn.x() = numeric::min(mn.x(),p.xyz(AtomID(ia,ir)).x());
			// mn.y() = numeric::min(mn.y(),p.xyz(AtomID(ia,ir)).y());
			// mn.z() = numeric::min(mn.z(),p.xyz(AtomID(ia,ir)).z());
			// mx.x() = numeric::max(mx.x(),p.xyz(AtomID(ia,ir)).x());
			// mx.y() = numeric::max(mx.y(),p.xyz(AtomID(ia,ir)).y());
			// mx.z() = numeric::max(mx.z(),p.xyz(AtomID(ia,ir)).z());
		}
	}
	// mn -= Vec(15,15,15);
	// mx += Vec(15,15,15);
	mn = Vec(-10,-10,-10);
	mx = Vec( 50, 50, 50);

	std::cout << "ncells: " << xyzhash.xdim()*xyzhash.ydim()*xyzhash.zdim() << std::endl;
	std::cout << "Bounds: " << mn << " " << mx << std::endl;
	std::cout << "Natom: " << real_natom << " " << xyzhash.natom() << std::endl;
	{
		utility::vector1<Vec> hashpts;
		for(int i = 0; i < (int)xyzhash.natom(); ++i) {
			hashpts.push_back( Vec( xyzhash.grid_atoms()[i].x(), xyzhash.grid_atoms()[i].y(), xyzhash.grid_atoms()[i].z() ) );
		}
		// dump_points_pdb(hashpts,"input_xyzhash.pdb");
	}

	protocols::sic_dock::xyzStripeHashPose cahash(p,protocols::sic_dock::CA,DIST);

	float tifc=0.0,th=0.0,ts=0.0,t=0.0;
	int tot = 0;
	for(int iter = 0; iter < (int)basic::options::option[basic::options::OptionKeys::out::nstruct](); ++iter) {

			Vec rv( numeric::random::uniform(),numeric::random::uniform(),numeric::random::uniform() );
			rv.x() = (mx.x()-mn.x()) * rv.x() + mn.x();
			rv.y() = (mx.y()-mn.y()) * rv.y() + mn.y();
			rv.z() = (mx.z()-mn.z()) * rv.z() + mn.z();


			utility::vector1<std::pair<int,int> > hpairs;
			std::set<int> resset;
			cahash.fill_pairs(rv,0,hpairs);
			slow_nbset(p,DIST,rv,resset);
			if( hpairs.size() != resset.size() ){
				utility_exit_with_message("biarist");
			}
			for(utility::vector1<std::pair<int,int> >::const_iterator i = hpairs.begin(); i != hpairs.end(); ++i){
				if(resset.count(i->second)!= 1) utility_exit_with_message("BAD nbrs");
			}

			numeric::geometry::hashing::Counter counter;

			int hash_nbcount=0,safe_nbcount=0;
			t = time_highres();
			hash_nbcount  = xyzhash.clash(rv);
			th += time_highres()-t;

			t = time_highres();
			safe_nbcount  = slow_clash(p,DIST,rv);
			ts += time_highres()-t;

			t = time_highres();
			int ifc_nbcount = ifc.clash_count(rv);
			tifc += time_highres()-t;

			tot += hash_nbcount;
			if( safe_nbcount != hash_nbcount /*|| safe_nbcount != ifc_nbcount*/ ) {
				std::cout << rv << std::endl;
				std::cout << safe_nbcount << " " << hash_nbcount << " " << ifc_nbcount << std::endl;


				CollectNeighbors cn; xyzhash.visit(rv,cn);
				utility::vector1<Vec> ref(slow_nbget(xyzhash,DIST,rv)), tst(cn.nbrs);

				for(utility::vector1<Vec>::const_iterator i = tst.begin(); i != tst.end(); ++i){
					bool i_in = false;
					for(utility::vector1<Vec>::const_iterator j = ref.begin(); j != ref.end(); ++j){
						if( i->distance_squared(*j) < 0.0001 ) i_in = true;
					}
					if(!i_in) std::cout << "TST NOT IN REF: " << ObjexxFCL::format::F(20,18,i->distance(rv+xyzhash.translation())) << " " << (*i) << std::endl;
				}
				for(utility::vector1<Vec>::const_iterator i = ref.begin(); i != ref.end(); ++i){
					bool i_in = false;
					for(utility::vector1<Vec>::const_iterator j = tst.begin(); j != tst.end(); ++j){
						if( i->distance_squared(*j) < 0.0001 ) i_in = true;
					}
					if(!i_in) std::cout << "REF NOT IN TST: " << ObjexxFCL::format::F(20,18,i->distance(rv+xyzhash.translation())) << " " << (*i) << std::endl;
				}


				utility_exit_with_message("FAIL after "+ObjexxFCL::format::I(10,iter)+"!!!");
			}
	}
	std::cout << ObjexxFCL::format::I(10,basic::options::option[basic::options::OptionKeys::out::nstruct]())+" nb counts of random points match. Woot! " << tot << std::endl;
	std::cout << "hash speedup over N^2 count: " << ts  /th << std::endl;
	std::cout << "hash speedup over IFC count: " << tifc/th << std::endl;

	return 0;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
