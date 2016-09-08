// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#define CONTACT_D2 20.25
#define CONTACT_TH 0
#define NSS 672 // 1812 8192 17282
#define MIN_HELEX_RES 20
#define MAX_CYS_RES 3
#define MAX_NRES 99999

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/kinematics/Stub.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/packing/compute_holes_score.hh>
#include <core/scoring/rms_util.hh>
#include <numeric/conversions.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <numeric/xyzVector.hh>
#include <protocols/sic_dock/SICFast.hh>
#include <protocols/sic_dock/designability_score.hh>

#include <apps/pilot/will/will_util.ihh>

#include <numeric/xyzTransform.hh>

typedef numeric::xyzVector<core::Real> Vec;
typedef numeric::xyzMatrix<core::Real> Mat;
typedef numeric::xyzTransform<core::Real> Xform;
Vec Ux(1,0,0),Uy(0,1,0),Uz(0,0,1),V0(0,0,0);

using core::id::AtomID;
using basic::options::option;
using core::pose::Pose;
using core::Real;
using core::scoring::ScoreFunctionOP;
using core::Size;
using numeric::max;
using numeric::min;
using numeric::random::gaussian;
using numeric::random::uniform;
using numeric::rotation_matrix_degrees;
using numeric::conversions::radians;
using numeric::conversions::degrees;
using ObjexxFCL::format::F;
using ObjexxFCL::format::I;
using ObjexxFCL::string_of;
using std::cerr;
using std::cout;
using std::string;
using utility::io::izstream;
using utility::io::ozstream;
using utility::vector1;
using std::endl;
using core::import_pose::pose_from_file;
using core::kinematics::Stub;

static THREAD_LOCAL basic::Tracer TR( "CXdock" );
static core::io::silent::SilentFileData sfd;

OPT_1GRP_KEY( FileVector, cxdock, bench2 )
OPT_1GRP_KEY( FileVector, cxdock, bench3 )
OPT_1GRP_KEY( FileVector, cxdock, bench4 )
OPT_1GRP_KEY( FileVector, cxdock, bench5 )
OPT_1GRP_KEY( FileVector, cxdock, bench6 )
OPT_1GRP_KEY( IntegerVector, cxdock, nfold )

void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	NEW_OPT( cxdock::bench2 , "benchmark dimers (must be aligned-z)" , "" );
	NEW_OPT( cxdock::bench3 , "benchmark trimers (must be aligned-z)" , "" );
	NEW_OPT( cxdock::bench4 , "benchmark tetra (must be aligned-z)" , "" );
	NEW_OPT( cxdock::bench5 , "benchmark penta (must be aligned-z)" , "" );
	NEW_OPT( cxdock::bench6 , "benchmark hexa (must be aligned-z)" , "" );
	NEW_OPT( cxdock::nfold  , "C syms to test", 3 );
}

#ifdef USE_OPENMP
#include <omp.h>
#endif
int num_threads() {
	#ifdef USE_OPENMP
		return omp_get_max_threads();
	#else
		return 1;
	#endif
}
int thread_num() {
	#ifdef USE_OPENMP
		return omp_get_thread_num();
	#else
		return 0;
	#endif
}

// info about a "hit"
struct Hit {
	int iss,irt,sym;
	Real cbc,xscore,rmsd;
	Xform s1,s2;
	Hit(int is, int ir, Real cb, int sm) : iss(is),irt(ir),sym(sm),cbc(cb) {}
};
bool cmpcbc (Hit i,Hit j) { return i.cbc    > j.cbc   ; }
bool cmpxsc (Hit i,Hit j) { return i.xscore > j.xscore; }
bool cmprmsd(Hit i,Hit j) { return i.rmsd   < j.rmsd; }


inline Vec get_rot_center(Xform const & xsym, int sym){
	Vec cen(0,0,0),pt(0,0,0);
	for(int i = 1; i <= sym; ++i){ cen += pt; pt = xsym*pt; }
	cen /= Real(sym);
	return cen;
}
inline Vec get_rot_center(Xform const & x1, Xform const & x2, int sym){
	return get_rot_center( x2 * ~x1, sym );
}

inline
Xform
get_cx_xform(Hit const & h){
	// Real halfside = (h.s1.t-h.s2.t).length()/2.0;
	// Real r = halfside / tan(numeric::constants::d::pi/(Real)h.sym); // tan(a/2) * edge_len/2
	// cout << halfside << " " << r << endl;
	// core::kinematics::Stub x(h.s2);
	// x.t = x.t + Vec(0,r,halfside);
	// Mat R = rotation_matrix_degrees(Vec(1,0,1),180.0);
	// return Xform( R*x.R, R*x.t );
	Xform x1(h.s1.R,h.s1.t),x2(h.s2.R,h.s2.t);
	Vec cen = get_rot_center(x1,x2, h.sym);
	// cout << cen << endl;
	return Xform( h.s2.R, -cen );
}


// dock pose against rotated version of itself
// rotate all different ways by looping over rotation axes on unit sphere point
// samples in ssamp, then loopint over rotation magnitudes from 1-180°
// picks cases with most CB contacts in helices
void
dock(
	Pose const & init_pose,
	std::string const & fn
){
	using namespace basic::options;
	using namespace OptionKeys;

	// set up one SIC per thread
	protocols::sic_dock::SICFast sic;
	sic.init(init_pose);

	utility::vector1<Vec> init_ca;
	for(Size ir = 1; ir <= init_pose.size(); ++ir){
		if(init_pose.residue(ir).has("CA")){
			init_ca.push_back(init_pose.residue(ir).xyz("CA"));
		}
	}

	init_pose.dump_pdb(utility::file_basename(fn)+"_start.pdb");
	utility::io::ozstream xforms_out(utility::file_basename(fn)+"_xforms.dat");

	//store initial coords & angle samples 1-180 (other distribution of samps would be better...)
	// vector1<double> asamp; for(Real i = 1.0; i < 180.0; i += 5.0) asamp.push_back(i);

	vector1<int> syms = option[cxdock::nfold]();
	vector1<Mat> Rsym(syms.size());
	for(int ic = 1; ic <= (int)syms.size(); ic++) Rsym[ic] = rotation_matrix_degrees(Ux,360.0/Real(syms[ic]));

	// utility::vector1<core::kinematics::Stub> stubs;
	// utility::vector1<char> ss;
	// for(Size ir = 1; ir <= init_pose.size(); ++ir){
	// 	if( !init_pose.residue(ir).is_protein() ) continue;
	// 	if( !init_pose.residue(ir).has("CB") ) continue;
	// 	if( init_pose.secstruct(ir) == 'L' ) continue;
	// 	stubs.push_back( core::kinematics::Stub(init_pose.residue(ir).xyz("CB"),init_pose.residue(ir).xyz("CA"),init_pose.residue(ir).xyz("N")));
	// 	ss.push_back(init_pose.secstruct(ir));
	// }

	// loop over axes on sphere points & rotaton amounts 1-180
	// #ifdef USE_OPENMP
	// #pragma omp parallel for schedule(dynamic,1)
	// #endif
	// for(int iss = 1; iss <= (int)ssamp.size(); ++iss){
	// 	if(iss%1000==0) TR << iss << " of " << NSS << std::endl;
		// for(int irt = 1; irt <= (int)asamp.size(); ++irt) {
		// 	Mat R = rotation_matrix_degrees(ssamp[iss],(Real)asamp[irt]);


	int S1 = 3;
	int S2 = 2;
	Vec A1 = Ux;
	Vec A2 = Uy;
	Mat swapXZ = rotation_matrix_degrees(Vec(1,0,1),180.0); // swap axis
	Real ANG1 = 360.0/Real(S1);
	Real ANG2 = 360.0/Real(S2);

	Mat Rsym1 = rotation_matrix_degrees(Ux,ANG1);
	Mat Rsym2 = rotation_matrix_degrees(Uy,ANG2);

	//int irt=1, iss=1, ic=1;  // unused ~Labonte
	for(Real az = 0.0; az < 1.0; az += 5.0){
		Mat Rz = rotation_matrix_degrees(Uz,az);
	for(Real ax = 0.0; ax < 1.0; ax += 5.0){
		Mat Rx = rotation_matrix_degrees(Ux,ax);
	for(Real ay = 0.0; ay < 360.0; ay += 5.0){
		Mat Ry = rotation_matrix_degrees(Uy,ay);

		Mat R = Ry*Rx*Rz;
		Xform const x1(       R, V0);
		Xform const x2( Rsym1*R, V0 );

		Real t = sic.slide_into_contact(x1,x2,Uz);

		Hit h(1,1,1,S1); h.s1 = x1; h.s2 = x2; h.s1.t += t*Uz;

		Xform x = get_cx_xform(h);

		Pose tmp0(init_pose),tmp1(init_pose),tmp2(init_pose);
		xform_pose(tmp0,Stub(swapXZ*x.R,swapXZ*x.t));
		xform_pose(tmp1,h.s1);
		xform_pose(tmp2,h.s2);
		tmp0.dump_pdb("tmp0.pdb");
		tmp1.dump_pdb("tmp1.pdb");
		tmp2.dump_pdb("tmp2.pdb");
		utility_exit_with_message("test c2");

		for(Real dx = 0; dx < ANG1; dx += 5.0){
			Real t = 0.0;
			Xform xd = rotation_matrix_degrees(A1,dx) * x;
			Stub xmx1,xmx2;
			Real angmax = -1;
			for(int it = 1; it <= S1; ++it){
				Mat Rd = rotation_matrix_degrees(A1,Real(it-1)*ANG1);
				Stub x3(         xd.R,         xd.t);
				Stub x4(Rd*Rsym2*xd.R,Rd*Rsym2*xd.t);
				Real tmp = sic.slide_into_contact_DEPRICATED(x3,x4,A1);
				if(tmp > 9e8) continue;
				if( tmp > t && tmp < 9e8 ){
					t = tmp;
					x3.v += t*A1;
					// x3.v += t/2.0*A1;
					// x4.v -= t/2.0*A1;
					xmx1 = x3;
					xmx2 = x4;
					angmax = Real(it-1)*ANG1;
				}
				cout << t << " " << tmp << " " << angmax << endl;
			}
			cout << t << " " << angmax << endl;
			Xform x3(xmx1.M,xmx1.v);
			Xform x4(xmx2.M,xmx2.v);
			Vec cen = get_rot_center(x3,x4,S2);
			// Xform xout( x3.R, x3.t-cen-t/2.0*A1 );

			Xform xout( x4.R, x4.t - cen );
			xout = swapXZ * xout;
			xforms_out << ((numeric::xyzMatrix<Real>)(xout.R)) << " " << xout.t << std::endl;

			Pose tmp0(init_pose),tmp1(init_pose),tmp2(init_pose);
			xform_pose(tmp0,Stub(xout.R,xout.t));
			xform_pose(tmp1,Stub(swapXZ*xmx1.M,swapXZ*xmx1.v));
			xform_pose(tmp2,Stub(swapXZ*xmx2.M,swapXZ*xmx2.v));
			tmp0.dump_pdb("tmp0.pdb");
			tmp1.dump_pdb("tmp1.pdb");
			tmp2.dump_pdb("tmp2.pdb");
			cout << swapXZ*cen << endl;
			utility_exit_with_message("test d3");

		}

		utility_exit_with_message("testd3");
	}
	}
	}

	xforms_out.close();


}

utility::vector1<Vec>
align_native_state(
	core::pose::Pose & pose,
	int nfold
){
	// align so that axis is X and COM is along Z
	// for convenience of RMSD calc
	if(nfold<2) return utility::vector1<Vec>();
	Vec cen(0,0,0);
	Real ncen = 0.0;
	for(Size ir = 1; ir <= pose.size(); ++ir){
		if(pose.residue(ir).has("CA")){
			cen += pose.residue(ir).xyz("CA");
			ncen += 1.0;
		}
	}
	cen /= ncen;
	rot_pose( pose, Uz, -dihedral_degrees(Ux,V0,Uz,cen) );
	rot_pose( pose, Vec(1,0,1), 180.0 ); // symZ to symX
	utility::vector1<Vec> CAs;
	for(Size ir = 1; ir <= pose.size(); ++ir){
		if(pose.residue(ir).has("CA")){
			CAs.push_back(pose.residue(ir).xyz("CA"));
		}
	}
	return CAs;
}

void
read_sphere(
	vector1<Vec> & ssamp
){
	izstream is;
	if(!basic::database::open(is,"sampling/spheres/sphere_"+str(NSS)+".dat"))
		utility_exit_with_message("can't open sphere data");
	for(int i = 1; i <= NSS; ++i) {
		double x,y,z;
		is >> x >> y >> z;
		ssamp[i] = Vec(x,y,z);
	}
	is.close();
}

void
get_tasks_from_command_line(
	utility::vector1<std::pair<string,int> > & tasks
){
	using namespace basic::options;
	using namespace OptionKeys;
	if(option[in::file::s].user()){
		for(Size ifn = 1; ifn <= option[in::file::s]().size(); ++ifn){
			tasks.push_back(std::make_pair(option[in::file::s]()[ifn],1));
		}
	}
	if(option[cxdock::bench2].user()){
		for(Size ifn = 1; ifn <= option[cxdock::bench2]().size(); ++ifn){
			tasks.push_back(std::make_pair(option[cxdock::bench2]()[ifn],2));
		}
	}
	if(option[cxdock::bench3].user()){
		for(Size ifn = 1; ifn <= option[cxdock::bench3]().size(); ++ifn){
			tasks.push_back(std::make_pair(option[cxdock::bench3]()[ifn],3));
		}
	}
	if(option[cxdock::bench4].user()){
		for(Size ifn = 1; ifn <= option[cxdock::bench4]().size(); ++ifn){
			tasks.push_back(std::make_pair(option[cxdock::bench4]()[ifn],4));
		}
	}
	if(option[cxdock::bench5].user()){
		for(Size ifn = 1; ifn <= option[cxdock::bench5]().size(); ++ifn){
			tasks.push_back(std::make_pair(option[cxdock::bench5]()[ifn],5));
		}
	}
	if(option[cxdock::bench6].user()){
		for(Size ifn = 1; ifn <= option[cxdock::bench6]().size(); ++ifn){
			tasks.push_back(std::make_pair(option[cxdock::bench6]()[ifn],6));
		}
	}
}

int main(int argc, char *argv[]) {

	try {

	register_options();
	devel::init(argc,argv);
	using namespace basic::options::OptionKeys;

	vector1<Vec> ssamp(NSS); read_sphere(ssamp);

	// protocols::sic_dock::XfoxmScore const xfs("/Users/sheffler/project/designability_stats/results/");

	vector1<std::pair<string,int> > tasks;
	get_tasks_from_command_line(tasks);

	for(vector1<std::pair<string,int> >::const_iterator i = tasks.begin(); i != tasks.end(); ++i){
		string fn = i->first;
		int nfold = i->second;
		TR << "searching " << fn << " " << nfold << std::endl;

		Pose pnat;
		core::import_pose::pose_from_file(pnat,fn, core::import_pose::PDB_file);
		core::scoring::dssp::Dssp dssp(pnat);
		dssp.insert_ss_into_pose(pnat);

		// center on z
		trans_pose(pnat,Vec(0,0,-center_of_geom(pnat,1,pnat.size()).z()));

		// Pose olig; make_native_olig(pnat,olig,nfold);

		// // sym axis now along X, rotated around X so CA com in XZ plane,
		// // this makes computing the RMSD of docked states easier
		// vector1<Vec> native_ca = align_native_state(pnat,nfold); // align CA com
		// pnat.dump_pdb("native_align.pdb");

		Vec cen = center_of_geom(pnat,1,pnat.size());
		trans_pose(pnat,-Vec(0,cen.y(),cen.z()));

		if( pnat.size() > MAX_NRES ) continue;
		//Size cyscnt=0, nhelix=0;  // unused ~Labonte
		// for(Size ir = 2; ir <= pnat.size()-1; ++ir) {
		// 	if(pnat.secstruct(ir) == 'H') nhelix++;
		// 	//if(!pnat.residue(ir).is_protein()) goto cont1;
		// 	if(pnat.residue(ir).is_lower_terminus()) remove_lower_terminus_type_from_pose_residue(pnat,ir);//goto cont1;
		// 	if(pnat.residue(ir).is_upper_terminus()) remove_upper_terminus_type_from_pose_residue(pnat,ir);//goto cont1;
		// 	if(pnat.residue(ir).name3()=="CYS") { if(++cyscnt > MAX_CYS_RES) goto cont1; }
		// } goto done1; cont1: TR << "skipping " << fn << std::endl; continue; done1:
		// if( nhelix < MIN_HELEX_RES ) continue;

		dock(pnat,fn);
	}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}


// core::Real
// compute_xform_score(
// 	protocols::sic_dock::XfoxmScore const & xfs,
// 	core::kinematics::Stub const & x1,
// 	core::kinematics::Stub const & x2,
// 	utility::vector1<core::kinematics::Stub> const & stubs,
// 	utility::vector1<char> const & ss
// ){
// 	assert(ss.size()==stubs.size());
// 	core::Real score = 0.0;
// 	utility::vector1<char>::const_iterator ssi = ss.begin();
// 	for(utility::vector1<core::kinematics::Stub>::const_iterator i = stubs.begin(); i != stubs.end(); ++i,++ssi){
// 		utility::vector1<char>::const_iterator ssj = ss.begin();
// 		for(utility::vector1<core::kinematics::Stub>::const_iterator j = stubs.begin(); j != stubs.end(); ++j,++ssj){
// 			Stub const a( x1.M * i->M, x1.M * i->v + x1.v );
// 			Stub const b( x2.M * j->M, x2.M * j->v + x2.v );
// 			if( (a.v-b.v).length_squared() > 64.0 ) continue;
// 			score += xfs.score(a,b,*ssi,*ssj);
// 		}
// 	}
// 	return score;
// }

// core::Real
// compute_xform_score(
// 	protocols::sic_dock::XfoxmScore const & xfs,
// 	core::pose::Pose const & pose1,
// 	core::pose::Pose const & pose2
// ){
// 	float tot_score = 0.0;
// 	for(core::Size ir = 1; ir <= pose1.size(); ++ir){
// 		if(!pose1.residue(ir).is_protein()) continue;
// 		if(!pose1.residue(ir).has("CB")) continue;
// 		Vec CBi = pose1.residue(ir).xyz("CB");
// 		Vec CAi = pose1.residue(ir).xyz("CA");
// 		Vec  Ni = pose1.residue(ir).xyz( "N");
// 		core::kinematics::Stub sir(CBi,CAi,Ni);
// 		for(core::Size jr = 1; jr <= pose2.size(); ++jr){
// 			if(!pose2.residue(jr).is_protein()) continue;
// 			if(!pose2.residue(jr).has("CB")) continue;
// 			Vec CBj = pose2.residue(jr).xyz("CB");
// 			Vec CAj = pose2.residue(jr).xyz("CA");
// 			Vec  Nj = pose2.residue(jr).xyz( "N");
// 			if( CBi.distance_squared(CBj) > 64.0 ) continue;
// 			core::kinematics::Stub sjr(CBj,CAj,Nj);
// 			core::Real tmp = xfs.score(sir,sjr,pose1.secstruct(ir),pose2.secstruct(jr));
// 			tot_score += tmp;
// 		}
// 	}
// 	return tot_score;
// }

// void
// make_native_olig(
// 	core::pose::Pose const & pose,
// 	core::pose::Pose       & olig,
// 	int nfold
// ){
// 	olig = pose;
// 	Pose tmp(pose);
// 	for(int i = 2; i <= nfold; ++i){
// 		rot_pose(tmp,Uz,360.0/nfold);
// 		for(Size i = 1; i <= tmp.size(); ++i){
// 			if(olig.residue(i).is_lower_terminus()||olig.residue(i).is_ligand()){
// 				olig.append_residue_by_jump(tmp.residue(i),1);
// 			} else {
// 				olig.append_residue_by_bond(tmp.residue(i));
// 			}
// 		}
// 	}
// }

// void
// make_dock_olig(
// 	core::pose::Pose const & pose,
// 	core::pose::Pose       & olig,
// 	Hit h
// ){
// 	Pose tmp1(pose);
// 	Real halfside = (h.s1.v-h.s2.v).length()/2.0;
// 	Real r = halfside / tan(numeric::constants::d::pi/(Real)h.sym); // tan(a/2) * edge_len/2
// 	xform_pose(tmp1,core::kinematics::Stub(h.s1.M,V0));
// 	trans_pose(tmp1,Vec(0,0,halfside));
// 	trans_pose(tmp1,Vec(0,r,0));
// 	rot_pose(tmp1,Vec(1,0,1),180.0);
// 	olig = tmp1;
// 	// tmp1.dump_pdb(outdir+tag);
// 	// make_native_olig(tmp1,olig,h.sym);
// }


// inline
// Real
// get_rmsd(
// 	utility::vector1<Vec> const & native_ca,
// 	utility::vector1<Vec> const & init_ca,
// 	Hit const & h//,
// 	// Pose pose
// ){
// 	Real halfside = (h.s1.v-h.s2.v).length()/2.0;
// 	Real r = halfside / tan(numeric::constants::d::pi/(Real)h.sym); // tan(a/2) * edge_len/2
// 	core::kinematics::Stub x(h.s1);
// 	x.v = x.v + Vec(0,r,halfside);
// 	{
// 		Real ang = dihedral_degrees(Uz,V0,Ux, x.local2global(V0) );
// 		Mat R = numeric::x_rotation_matrix_degrees( -ang );
// 		x.M = R * x.M;
// 		x.v = R * x.v;
// 	}

// 	Mat Rz = numeric::z_rotation_matrix_degrees(180.0);
// 	Real rmsd=0.0, rev_rmsd=0.0;
// 	utility::vector1<Vec>::const_iterator i=native_ca.begin(), j=init_ca.begin();
// 	for(; i != native_ca.end(); ++i,++j){
// 		Vec v = x.local2global(*j);
// 		rmsd     += i->distance_squared(      v );
// 		rev_rmsd += i->distance_squared( Rz * v );
// 	}
// 	if(rev_rmsd < rmsd){
// 		rmsd = rev_rmsd;
// 	}
// 	rmsd = sqrt( rmsd/(Real)native_ca.size() );

// 	return rmsd;
// }

// inline
// Real
// get_rmsd_debug(
// 	utility::vector1<Vec> const & native_ca,
// 	utility::vector1<Vec> const & init_ca,
// 	Hit const & h,
// 	Pose pose
// ){
// 	Real halfside = (h.s1.v-h.s2.v).length()/2.0;
// 	Real r = halfside / tan(numeric::constants::d::pi/(Real)h.sym); // tan(a/2) * edge_len/2
// 	core::kinematics::Stub x(h.s1);
// 	x.v = x.v + Vec(0,r,halfside);
// 	{
// 		Real ang = dihedral_degrees(Uz,V0,Ux, x.local2global(V0) );
// 		Mat R = numeric::x_rotation_matrix_degrees( -ang );
// 		x.M = R * x.M;
// 		x.v = R * x.v;
// 	}

// 	xform_pose(pose,x);

// 	Mat Rz = numeric::z_rotation_matrix_degrees(180.0);
// 	Real rmsd=0.0, rev_rmsd=0.0;
// 	utility::vector1<Vec>::const_iterator i=native_ca.begin(), j=init_ca.begin();
// 	for(; i != native_ca.end(); ++i,++j){
// 		Vec v = x.local2global(*j);
// 		rmsd     += i->distance_squared(      v );
// 		rev_rmsd += i->distance_squared( Rz * v );
// 	}
// 	if(rev_rmsd < rmsd){
// 		rmsd = rev_rmsd;
// 		rot_pose(pose,Uz,180.0);
// 	}
// 	rmsd = sqrt( rmsd/(Real)native_ca.size() );


// 	pose.dump_pdb("align_rmsd.pdb");
// 	// utility_exit_with_message("foo");

// 	return rmsd;
// }

