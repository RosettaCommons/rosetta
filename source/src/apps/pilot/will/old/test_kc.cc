// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /src/apps/pilat/will/genmatch.cc
/// @brief ???

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <core/conformation/Residue.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <numeric/random/random.hh>
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>
#include <numeric/kinematic_closure/bridgeObjects.hh>
#include <numeric/kinematic_closure/kinematic_closure_helpers.hh>

#include <sstream>
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>

#include <core/kinematics/Jump.hh>
#include <numeric/conversions.hh>

// #include <devel/init.hh>

// #include <core/scoring/constraints/LocalCoordinateConstraint.hh>

using core::Real;
using core::Size;
using utility::vector1;
typedef numeric::xyzVector<Real> Vec;
typedef numeric::xyzMatrix<Real> Mat;
typedef utility::vector1<Size> sizes;
typedef utility::vector1<int> ints;
typedef utility::vector1<Real> reals;
typedef utility::vector1<Vec>  vecs;

static basic::Tracer TR( "test_kc" );

vector1<reals> vecs2vv(vecs const & v) {
	vector1<reals> vv;
	for ( vecs::const_iterator i = v.begin(); i != v.end(); ++i ) {
		reals r(3);
		r[1] = i->x();
		r[2] = i->y();
		r[3] = i->z();
		vv.push_back(r);
	}
	return vv;
}


void test_bridgeObjects() {
	using namespace numeric::kinematic_closure;
	// input arguments
	int N=8;
	utility::vector1<Real> t_ang0 (N), b_ang0 (N), b_len0 (N), dt (N), db (N), da (N);
	utility::vector1<int> pivots (3), order (3);
	// input argument values
	// DJM [ if inputting torsions, uncomment next 4 lines ]
	Real t_ang0vals[] = {36.0213, 0.0, 253.8843, 52.8441, 36.0366, 359.8985, 253.9602, 52.8272};
	Real b_ang0vals[] = {114.9908, 114.9805, 115.0112, 114.9923, 114.9582, 115.0341, 114.9920, 115.0594};
	Real b_len0vals[] = {1.5200, 1.5202, 1.5197, 1.5205, 1.5198, 1.5201, 1.5197, 1.5196};
	Real dtvals[] = {180.0000, 180.0000, -5.3651, 180.0000, 180.0000, 180.0000, 180.0000, -5.0393};

	utility::vector1<utility::vector1<Real> > atoms (N);
	// output arguments
	int nsol;
	utility::vector1<utility::vector1<Real> > t_ang, b_ang, b_len;
	/*
	* DJM uncomment if inputting torsions
	*/
	int ind=0;
	for ( int i=1; i<=8; i++ ) {
		t_ang0[i]=t_ang0vals[ind++];
	}

	ind=0;
	for ( int i=1; i<=8; i++ ) {
		b_ang0[i]=b_ang0vals[ind++];
	}

	ind=0;
	for ( int i=1; i<=8; i++ ) {
		b_len0[i]=b_len0vals[ind++];
	}

	ind=0;
	for ( int i=1; i<=8; i++ ) {
		dt[i]=dtvals[ind++];
	}

	for ( int i=1; i<=8; i++ ) {
		db[i]=1.52;
	}

	for ( int i=1; i<=8; i++ ) {
		da[i]=115;
	}
	//*/
	pivots[1] = 2;
	pivots[2] = 5;
	pivots[3] = 7;
	order[1] = 1;
	order[2] = 2;
	order[3] = 3;

	//bridgeObjects(N, t_ang0, b_ang0, b_len0, dt, db, da, pivots, order, t_ang, b_ang, b_len, nsol);
	bridgeObjects(atoms, dt, db, da, pivots, order, t_ang, b_ang, b_len, nsol);
	//bridgeObjects(atoms, pivots, order, t_ang, b_ang, b_len, nsol);
	printMatrix(t_ang);
	printMatrix(b_ang);
	printMatrix(b_len);
	std::cout << nsol << std::endl;
}

/* // DJM: this was for benchmarking current version of bridgeObjects against earlier (slower) version.
// v2 is now current version. 6-8-08
void test_bridgeObjects_v2() {
// input arguments
int N=27;
utility::vector1<Real> dt (N), db (N), da (N);
utility::vector1<int> pivots (3), order (3);
// input argument values
// DJM [ if inputting torsions, uncomment next 4 lines ]
Real dt_ang_vals[] = {348.798, 337.654, 181.609, -170, 140, 179.743, -140, 80, 181.769, -130, 10, 180.419, -90, -30, 177.171, -90, -20, 184.154, -90, -40, 182.472, -80, 170, 184.776, 310.274, 150.675, 290.245};
Real db_ang_vals[] = {88.2267, 112.669, 115.978, 122.109, 110.442, 117.607, 121.606, 110.194, 115.415, 120.286, 110.888, 117.314, 125.907, 110.287, 117.773, 121.692, 109.691, 115.437, 122.169, 109.827, 119.212, 120.735, 111.479, 115.289, 121.921, 114.467, 85.1688};
Real db_len_vals[] = {1.41632, 1.48472, 1.32645, 1.46204, 1.54295, 1.31506, 1.45302, 1.51076, 1.28081, 1.46153, 1.51705, 1.33131, 1.46801, 1.5172, 1.34797, 1.45556, 1.49902, 1.31172, 1.41681, 1.49235, 1.36072, 1.4561, 1.52821, 1.31126, 1.44372, 1.48442, 9.90623};
// initialize the atoms vector
utility::vector1<utility::vector1<Real> > atoms (N);
Real atoms_vals[] = {13.92, 13.881, 12.607, 11.925, 10.671, 9.481, 9.657, 8.597, 8.771, 8.01, 8.116, 6.871, 5.803, 4.504, 3.703, 2.616, 1.752, 0.851, 0.206, -0.747, -0.251, 0.965, 1.458, 2.652, 2.912, 4.066, 4.232, 43.318, 44.382, 44.429, 43.293, 43.179, 43.781, 44.046, 44.619, 46.119, 46.785, 48.241, 48.954, 48.205, 48.696, 47.567, 47.894, 46.878, 46.229, 47.052, 46.598, 45.73, 45.989, 45.197, 44.279, 44.041, 43.277, 41.984, 25.637, 24.703, 23.942, 23.88, 23.137, 23.913, 25.189, 26.001, 26.047, 25.261, 25.191, 25.684, 25.95, 26.426, 27.047, 27.774, 28.357, 27.35, 26.558, 25.613, 24.505, 23.952, 22.834, 23.093, 24.356, 24.767, 24.057};
for (int i=1; i<=N; i++) {
atoms[i].resize(3);
}
int ind=0;
for (int j=1; j<=3; j++) {
for (int i=1; i<=N; i++) {
atoms[i][j] = atoms_vals[ind++];
}
}
ind=0;
// initialize the designed torsions, angles, and length
for (int i=1; i<=N; i++) {
dt[i]=dt_ang_vals[ind++];
}
ind=0;
for (int i=1; i<=N; i++) {
da[i]=db_ang_vals[ind++];
}
ind=0;
for (int i=1; i<=N; i++) {
db[i]=db_len_vals[ind++];
}
// these pivots specifically correspond to the designed inputs to produce the expected outputs
pivots[1] = 5;
pivots[2] = 14;
pivots[3] = 23;
order[1] = 1;
order[2] = 2;
order[3] = 3;

// DJM: debug
//std::cout << "dt: " << std::endl;
//printVector(dt);
//std::cout << "db: " << std::endl;
//printVector(db);
//std::cout << "da: " << std::endl;
//printVector(da);
// output arguments
int nits=100000; // number times to call bridgeObjects for performance benchmarking
int nsol;
utility::vector1<utility::vector1<Real> > t_ang, b_ang, b_len;

std::cout << "calling bridge2" << std::endl;
std::time_t const start_time_v2( std::time( NULL ) );
for (int d=1; d<=nits; d++) {
bridgeObjects_v2(atoms, dt, da, db, pivots, order, t_ang, b_ang, b_len, nsol);
}
std::time_t const end_time_v2( std::time( NULL ) );
std::cout << "start time: " << start_time_v2 << std::endl;
std::cout << "end time: " << end_time_v2 << std::endl;
std::cout << "start time minus end time v2: " << end_time_v2 - start_time_v2 << std::endl;
std::cout << "from bridgeObjects v2: " << std::endl;
printMatrix(t_ang);

nsol=0;

std::cout << "calling bridge1" << std::endl;
std::time_t const start_time_v1( std::time( NULL ) );
for (int d=1; d<=nits; d++) {
bridgeObjects(atoms, dt, da, db, pivots, order, t_ang, b_ang, b_len, nsol);
}
std::time_t const end_time_v1( std::time( NULL ) );
std::cout << "start time: " << start_time_v1 << std::endl;
std::cout << "end time: " << end_time_v1 << std::endl;
std::cout << "start time minus end time v1: " << end_time_v1 - start_time_v1 << std::endl;
std::cout << "from bridgeObjects v1: " << std::endl;
printMatrix(t_ang);

//std::cout << "output t_ang: " << std::endl;
//printMatrix(t_ang);
//std::cout << "output b_ang: " << std::endl;
//printMatrix(b_ang);
//std::cout << "output b_len: " << std::endl;
//printMatrix(b_len);
////TR.Debug << nsol << std::endl;
}
*/


void test_kc() {


	//scoring::Ramachandran const & rama( core::scoring::ScoringManager::get_instance()->get_Ramachandran() );
	using numeric::conversions::radians;
	using numeric::conversions::degrees;

	using core::id::AtomID;
	using core::pose::Pose;
	using namespace numeric::kinematic_closure;
	using basic::options::option;
	using namespace basic::options::OptionKeys;

	Pose pose;
	core::import_pose::pose_from_file(pose,option[in::file::s]()[1], core::import_pose::PDB_file);


	Size start_res_  = 4;
	Size middle_res_ = 6;
	Size end_res_    = 8;

	// DJM: debug
	//TR << "from " << start_res_ << " to " << end_res_ << std::endl;

	// inputs to loop closure
	utility::vector1<utility::vector1<Real> > atoms;
	utility::vector1<Size> pivots (3), order (3);
	// outputs from loop closure
	utility::vector1<utility::vector1<Real> > t_ang, b_ang, b_len;
	int nsol=0;
	// for eliminating identical solutions
	utility::vector1<utility::vector1<Real> > Q0 (3);
	utility::vector1<Real> dt_ang, db_len, db_ang, R0 (3);

	Size middle_offset = middle_res_ - start_res_; // is used to set central pivot atom
	Size seg_len = 1 + end_res_ - start_res_; // length of the closure chain
	atoms.resize((seg_len + 2) * 3); // one extra residue on each side to establish the geometric frame


	// Get the current coords of the loop closure segment (if a pivot is terminal we took the coords above so skip)
	Size ind = 1;
	for ( Size i = start_res_-1; i <= end_res_+1; i++ ) {
		core::conformation::Residue res=pose.residue(i);
		for ( Size j=1; j<=3; j++ ) { // DJM: just keeping N, CA, C atoms. We assume these are always the first 3.  BAD -- PROTEIN ONLY ASSUMPTION -- How about metal ions with only 1 atom?
			atoms[ind].resize(3);
			atoms[ind][1] = static_cast<Real> (res.xyz(j).x());//+numeric::random::gaussian());
			atoms[ind][2] = static_cast<Real> (res.xyz(j).y());//+numeric::random::gaussian());
			atoms[ind][3] = static_cast<Real> (res.xyz(j).z());//+numeric::random::gaussian());
			ind++;
		}
	}
	order[1]=1;
	order[2]=2;
	order[3]=3;
	Size pvatom1=5; // second C-alpha
	Size pvatom2=5 + (3 * middle_offset); // middle res C-alpha
	Size pvatom3=(3 * (seg_len+1)) - 1; // second-to-last C-alpha
	pivots[1]=pvatom1;
	pivots[2]=pvatom2;
	pivots[3]=pvatom3;
	chainTORS(atoms.size(), atoms, dt_ang, db_ang, db_len, R0, Q0);

	//  if(i%10000==0) std::cout << i << std::endl;

	// for(Size i = 1; i <= 1000000; ++i) {
	bridgeObjects(atoms, dt_ang, db_ang, db_len, pivots, order, t_ang, b_ang, b_len, nsol);
	//std::cout << "NSOL " << nsol << std::endl;
	//}

	TR << "nsol " << nsol << std::endl;

	for ( int isol = 1; isol <= nsol; isol++ ) {
		utility::vector1<utility::vector1<core::Real> > atm_out;
		numeric::kinematic_closure::chainXYZ(atoms.size(),b_len[isol],b_ang[isol],t_ang[isol],false,R0,Q0,atm_out);
		utility::io::ozstream out("test_"+ObjexxFCL::string_of(isol)+".pdb");
		for ( Size i = 1; i <= atm_out.size(); ++i ) {
			using namespace ObjexxFCL::format;
			out<<"HETATM"<<I(5,i)<<' '<<"VIZ "<<' '<<"VIZ"<<' '<<"A"<<I(4,i)<<"    "<<F(8,3,atm_out[i][1])<<F(8,3,atm_out[i][2])<<F(8,3,atm_out[i][3])<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
		}
		out.close();
	}

}


void test_kc2() {
	using basic::options::option;
	using namespace basic::options::OptionKeys;
	using core::pose::Pose;

	vecs atoms;
	reals dt,da,db,R0(3);
	ints pivots,order;
	vector1<reals> t_ang,b_ang,b_len,Q0(3);
	int nsol = -1;
	order.push_back(1);
	order.push_back(2);
	order.push_back(3);

	// N      CA     CB     SG    CEN     SG     CB      CA      N
	//    1.5    1.5    1.8   3.89   3.89    1.8    1.5     1.5
	//        109.5  109.5  100.0 109.5  100.0   109.5   109.5

	Pose pose;
	core::import_pose::pose_from_file(pose,option[in::file::s]()[1], core::import_pose::PDB_file);


	for ( int i = 1; i <= 8; ++i ) {
		atoms.push_back(pose.residue(i).xyz("N" ));
		atoms.push_back(pose.residue(i).xyz("CA"));
		atoms.push_back(pose.residue(i).xyz("C" ));
	}
	// {
	//  utility::io::ozstream out("start.pdb");
	//  using namespace ObjexxFCL::format;
	//  for(int i = 0; i < 8; ++i) {
	//   out<<"ATOM  "<<I(5,3*i+1)<<' '<<"  N "<<' '<<"ALA"<<' '<<"A"<<I(4,i+1)<<"    "<<F(8,3,atoms[3*i+1].x())<<F(8,3,atoms[3*i+1].y())<<F(8,3,atoms[3*i+1].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
	//   out<<"ATOM  "<<I(5,3*i+2)<<' '<<" CA "<<' '<<"ALA"<<' '<<"A"<<I(4,i+1)<<"    "<<F(8,3,atoms[3*i+2].x())<<F(8,3,atoms[3*i+2].y())<<F(8,3,atoms[3*i+2].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
	//   out<<"ATOM  "<<I(5,3*i+3)<<' '<<"  C "<<' '<<"ALA"<<' '<<"A"<<I(4,i+1)<<"    "<<F(8,3,atoms[3*i+3].x())<<F(8,3,atoms[3*i+3].y())<<F(8,3,atoms[3*i+3].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
	//  }
	//  out.close();
	// }


	pivots.push_back(5);
	pivots.push_back(8);
	pivots.push_back(20);

	numeric::kinematic_closure::chainTORS(atoms.size(), vecs2vv(atoms), dt, da, db, R0, Q0);
	// db[12] = 3.0;
	numeric::kinematic_closure::bridgeObjects( vecs2vv(atoms), dt, da, db, pivots, order, t_ang, b_ang, b_len, nsol );
	std::cout << nsol << std::endl;
	for ( int isol = 1; isol <= nsol; isol++ ) {
		utility::vector1<utility::vector1<core::Real> > atm_out;
		numeric::kinematic_closure::chainXYZ(atoms.size(),b_len[isol],b_ang[isol],t_ang[isol],false,R0,Q0,atm_out);
		utility::io::ozstream out("test2_"+ObjexxFCL::string_of(isol)+".pdb");
		using namespace ObjexxFCL::format;
		for ( int i = 0; i < 8; ++i ) {
			out<<"ATOM  "<<I(5,3*i+1)<<' '<<"  N "<<' '<<"ALA"<<' '<<"A"<<I(4,i+1)<<"    "<<F(8,3,atm_out[3*i+1][1])<<F(8,3,atm_out[3*i+1][2])<<F(8,3,atm_out[3*i+1][3])<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
			out<<"ATOM  "<<I(5,3*i+2)<<' '<<" CA "<<' '<<"ALA"<<' '<<"A"<<I(4,i+1)<<"    "<<F(8,3,atm_out[3*i+2][1])<<F(8,3,atm_out[3*i+2][2])<<F(8,3,atm_out[3*i+2][3])<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
			out<<"ATOM  "<<I(5,3*i+3)<<' '<<"  C "<<' '<<"ALA"<<' '<<"A"<<I(4,i+1)<<"    "<<F(8,3,atm_out[3*i+3][1])<<F(8,3,atm_out[3*i+3][2])<<F(8,3,atm_out[3*i+3][3])<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
		}
		out.close();
	}
}


int main (int argc, char *argv[]) {

	try {

		devel::init(argc,argv);
		// numeric::kinematic_closure::test_bridgeObjects();
		test_kc2();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
