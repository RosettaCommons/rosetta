// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:

#define CONTACT_D2 20.25
#define CONTACT_TH 0
#define NSS 1812 // 672 8192 17282
#define MIN_HELEX_RES 20
#define MAX_CYS_RES 3
#define MAX_NRES 200

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

typedef numeric::xyzVector<core::Real> Vec;
typedef numeric::xyzMatrix<core::Real> Mat;

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
using ObjexxFCL::fmt::F;
using ObjexxFCL::fmt::I;
using ObjexxFCL::string_of;
using std::cerr;
using std::cout;
using std::string;
using utility::io::izstream;
using utility::io::ozstream;
using utility::vector1;
using std::endl;
using core::import_pose::pose_from_pdb;
using core::kinematics::Stub;

static basic::Tracer TR("CXdock");
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
	NEW_OPT( cxdock::nfold  , "C syms to test", 1 );
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
	core::kinematics::Stub s1,s2;
	Hit(int is, int ir, Real cb, int sm) : iss(is),irt(ir),sym(sm),cbc(cb) {}
};
bool cmpcbc (Hit i,Hit j) { return i.cbc    > j.cbc   ; }
bool cmpxsc (Hit i,Hit j) { return i.xscore > j.xscore; }
bool cmprmsd(Hit i,Hit j) { return i.rmsd   < j.rmsd; }


core::Real
compute_xform_score(
	protocols::sic_dock::XfoxmScore const & xfs,
	core::kinematics::Stub const & x1,
	core::kinematics::Stub const & x2,
	utility::vector1<core::kinematics::Stub> const & stubs,
	utility::vector1<char> const & ss
){
	assert(ss.size()==stubs.size());
	core::Real score = 0.0;
	utility::vector1<char>::const_iterator ssi = ss.begin();
	for(utility::vector1<core::kinematics::Stub>::const_iterator i = stubs.begin(); i != stubs.end(); ++i,++ssi){
		utility::vector1<char>::const_iterator ssj = ss.begin();
		for(utility::vector1<core::kinematics::Stub>::const_iterator j = stubs.begin(); j != stubs.end(); ++j,++ssj){
			Stub const a( x1.M * i->M, x1.M * i->v + x1.v );
			Stub const b( x2.M * j->M, x2.M * j->v + x2.v );
			if( (a.v-b.v).length_squared() > 64.0 ) continue;
			score += xfs.score(a,b,*ssi,*ssj);
		}
	}
	return score;
}

core::Real
compute_xform_score(
	protocols::sic_dock::XfoxmScore const & xfs,
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2	
){
	float tot_score = 0.0;
	for(core::Size ir = 1; ir <= pose1.n_residue(); ++ir){
		if(!pose1.residue(ir).is_protein()) continue;
		if(!pose1.residue(ir).has("CB")) continue;
		Vec CBi = pose1.residue(ir).xyz("CB");
		Vec CAi = pose1.residue(ir).xyz("CA");
		Vec  Ni = pose1.residue(ir).xyz( "N");
		core::kinematics::Stub sir(CBi,CAi,Ni);
		for(core::Size jr = 1; jr <= pose2.n_residue(); ++jr){
			if(!pose2.residue(jr).is_protein()) continue;
			if(!pose2.residue(jr).has("CB")) continue;
			Vec CBj = pose2.residue(jr).xyz("CB");
			Vec CAj = pose2.residue(jr).xyz("CA");
			Vec  Nj = pose2.residue(jr).xyz( "N");
			if( CBi.distance_squared(CBj) > 64.0 ) continue;
			core::kinematics::Stub sjr(CBj,CAj,Nj);
			core::Real tmp = xfs.score(sir,sjr,pose1.secstruct(ir),pose2.secstruct(jr));
			tot_score += tmp;
		}
	}
	return tot_score;
}

void
make_native_olig(
	core::pose::Pose const & pose,
	core::pose::Pose       & olig,
	int nfold
){
	olig = pose;
	Pose tmp(pose);
	for(int i = 2; i <= nfold; ++i){
		rot_pose(tmp,Vec(0,0,1),360.0/nfold);
		for(Size i = 1; i <= tmp.n_residue(); ++i){
			if(olig.residue(i).is_lower_terminus()||olig.residue(i).is_ligand()){
				olig.append_residue_by_jump(tmp.residue(i),1);
			} else {
				olig.append_residue_by_bond(tmp.residue(i));
			}
		}
	}
}

void
make_dock_olig(
	core::pose::Pose const & pose,
	core::pose::Pose       & olig,
	Hit h
){
	Pose tmp1(pose);
	Real halfside = (h.s1.v-h.s2.v).length()/2.0;
	Real r = halfside / tan(numeric::constants::d::pi/(Real)h.sym); // tan(a/2) * edge_len/2
	xform_pose(tmp1,h.s1);
	trans_pose(tmp1,Vec(0,0,halfside));
	trans_pose(tmp1,Vec(0,r,0));
	rot_pose(tmp1,Vec(1,0,1),180.0);
	olig = tmp1;
	// tmp1.dump_pdb(outdir+tag);
	// make_native_olig(tmp1,olig,h.sym);
}

inline
Real
get_rmsd(
	utility::vector1<Vec> const & native_ca,
	utility::vector1<Vec> const & init_ca,
	Hit const & h//,
	// Pose pose
){
	Real halfside = (h.s1.v-h.s2.v).length()/2.0;
	Real r = halfside / tan(numeric::constants::d::pi/(Real)h.sym); // tan(a/2) * edge_len/2
	core::kinematics::Stub x(h.s1);
	x.v = x.v + Vec(0,r,halfside);
	{
		Real ang = dihedral_degrees(Vec(0,0,1),Vec(0,0,0),Vec(1,0,0), x.local2global(Vec(0,0,0)) );
		Mat R = numeric::x_rotation_matrix_degrees( -ang );
		x.M = R * x.M;
		x.v = R * x.v;
	}

	Mat Rz = numeric::z_rotation_matrix_degrees(180.0);
	Real rmsd=0.0, rev_rmsd=0.0;
	utility::vector1<Vec>::const_iterator i=native_ca.begin(), j=init_ca.begin();
	for(; i != native_ca.end(); ++i,++j){
		Vec v = x.local2global(*j);
		rmsd     += i->distance_squared(      v );
		rev_rmsd += i->distance_squared( Rz * v );
	}
	if(rev_rmsd < rmsd){
		rmsd = rev_rmsd;
	}
	rmsd = sqrt( rmsd/(Real)native_ca.size() );

	return rmsd;
}

inline
Real
get_rmsd_debug(
	utility::vector1<Vec> const & native_ca,
	utility::vector1<Vec> const & init_ca,
	Hit const & h,
	Pose pose
){
	Real halfside = (h.s1.v-h.s2.v).length()/2.0;
	Real r = halfside / tan(numeric::constants::d::pi/(Real)h.sym); // tan(a/2) * edge_len/2
	core::kinematics::Stub x(h.s1);
	x.v = x.v + Vec(0,r,halfside);
	{
		Real ang = dihedral_degrees(Vec(0,0,1),Vec(0,0,0),Vec(1,0,0), x.local2global(Vec(0,0,0)) );
		Mat R = numeric::x_rotation_matrix_degrees( -ang );
		x.M = R * x.M;
		x.v = R * x.v;
	}

	xform_pose(pose,x);

	Mat Rz = numeric::z_rotation_matrix_degrees(180.0);
	Real rmsd=0.0, rev_rmsd=0.0;
	utility::vector1<Vec>::const_iterator i=native_ca.begin(), j=init_ca.begin();
	for(; i != native_ca.end(); ++i,++j){
		Vec v = x.local2global(*j);
		rmsd     += i->distance_squared(      v );
		rev_rmsd += i->distance_squared( Rz * v );
	}
	if(rev_rmsd < rmsd){
		rmsd = rev_rmsd;
		rot_pose(pose,Vec(0,0,1),180.0);
	}
	rmsd = sqrt( rmsd/(Real)native_ca.size() );


	pose.dump_pdb("align_rmsd.pdb");
	// utility_exit_with_message("foo");

	return rmsd;
}

// dock pose against rotated version of itself
// rotate all different ways by looping over rotation axes on unit sphere point
// samples in ssamp, then loopint over rotation magnitudes from 1-180°
// picks cases with most CB contacts in helices
void
dock(
	protocols::sic_dock::XfoxmScore const & xfs,
	Pose const & init_pose,
	std::string const & fn,
	vector1<Vec> const & ssamp,

	Pose const & native_olig,
	int native_nfold = 1,
	vector1<Vec> const & native_ca = vector1<Vec>()

){
	using namespace basic::options;
	using namespace OptionKeys;

	// set up one SIC per thread
	std::vector<protocols::sic_dock::SICFast*> sics_;
	sics_.resize(num_threads());
	for(int i = 0; i < num_threads(); ++i) sics_[i] = new protocols::sic_dock::SICFast;
	for(int i = 0; i < num_threads(); ++i) sics_[i]->init(init_pose);

	utility::vector1<Vec> init_ca;
	for(Size ir = 1; ir <= init_pose.n_residue(); ++ir){
		if(init_pose.residue(ir).has("CA")){
			init_ca.push_back(init_pose.residue(ir).xyz("CA"));
		}
	}

	//store initial coords & angle samples 1-180 (other distribution of samps would be better...)
	vector1<double> asamp; for(Real i = 1.0; i < 180.0; i += 5.0) asamp.push_back(i);

	vector1<int> syms;
	if(native_nfold==1) syms = option[cxdock::nfold]();
	else syms.push_back(native_nfold);

	// setup rotation matricies
	vector1<Mat> Rsym(syms.size());
	for(int ic = 1; ic <= (int)syms.size(); ic++) Rsym[ic] = rotation_matrix_degrees(Vec(1,0,0),360.0/Real(syms[ic]));
	vector1<vector1<Hit> > hits(syms.size());

	utility::vector1<core::kinematics::Stub> stubs;
	utility::vector1<char> ss;
	for(Size ir = 1; ir <= init_pose.n_residue(); ++ir){
		if( !init_pose.residue(ir).is_protein() ) continue;
		if( !init_pose.residue(ir).has("CB") ) continue;
		if( init_pose.secstruct(ir) == 'L' ) continue;
		stubs.push_back( core::kinematics::Stub(init_pose.residue(ir).xyz("CB"),init_pose.residue(ir).xyz("CA"),init_pose.residue(ir).xyz("N")));
		ss.push_back(init_pose.secstruct(ir));
	}

	// loop over axes on sphere points & rotaton amounts 1-180
	#ifdef USE_OPENMP
	#pragma omp parallel for schedule(dynamic,1)
	#endif
	for(int iss = 1; iss <= (int)ssamp.size(); ++iss){
		if(iss%1000==0) TR << iss << " of " << NSS << std::endl;
		for(int irt = 1; irt <= (int)asamp.size(); ++irt) {
			Mat R = rotation_matrix_degrees(ssamp[iss],(Real)asamp[irt]);
			Stub const x1( R, Vec(0,0,0));
			// loop over symmitreis to check
			for(int ic = 1; ic <= (int)syms.size(); ic++){
				Stub const x2( Rsym[ic]*R, Vec(0,0,0) );
				Real cbc; 
				
				Real t = sics_[thread_num()]->slide_into_contact(x1,x2,Vec(0,0,1),cbc);
				
				Hit h(iss,irt,cbc,syms[ic]);
				h.s1 = x1;
				h.s2 = x2;
				h.s1.v += t*Vec(0,0,1);
				h.xscore = compute_xform_score(xfs,h.s1,h.s2,stubs,ss);
				h.rmsd = get_rmsd( native_ca, init_ca, h ); //, init_pose );
				if(cbc >= CONTACT_TH){
					#ifdef USE_OPENMP
					#pragma omp critical
					#endif
					{
						hits[ic].push_back(h);

						// Pose olig;
						// make_dock_olig(init_pose,olig,h);
						// Real actualrmsd = core::scoring::CA_rmsd(olig,native_olig);
						
							std::cout << "RMSD " 
							          // << actualrmsd << " " 
							          << h.rmsd << " " 
							          << h.cbc << " "
							          << h.xscore << " " 
							          << std::endl;

						// if( fabs(actualrmsd - h.rmsd) > 5.0 && h.rmsd < 5.0 ){
						// 	get_rmsd_debug( native_ca, init_ca, h, init_pose );
						// 	// core::scoring::calpha_superimpose_pose(olig,native_olig);
						// 	olig.dump_pdb("olig.pdb");
						// 	native_olig.dump_pdb("native.pdb");
						// 	utility_exit_with_message("test rmsd");
						// }


					}
				}
			}
		}
	}

	Size nstruct = basic::options::option[basic::options::OptionKeys::out::nstruct]();
	std::string outdir = basic::options::option[basic::options::OptionKeys::out::file::o]()+"/";

	// test
	for(int ic = 1; ic <= (int)syms.size(); ic++) {
		std::sort(hits[ic].begin(),hits[ic].end(),cmpcbc);
		for(int i = 1; i <= (int)hits[ic].size(); ++i) {
			Hit & h(hits[ic][i]);

			// Pose tmp1(init_pose),tmp2(init_pose);
			// xform_pose(tmp1,h.s1);
			// xform_pose(tmp2,h.s2);
			// tmp1.dump_pdb("test1.pdb");
			// tmp2.dump_pdb("test2.pdb");

			// Stub const a( h.s1.M * stubs[1].M, h.s1.M * stubs[1].v + h.s1.v );
			// Stub const b( h.s2.M * stubs[1].M, h.s2.M * stubs[1].v + h.s2.v );			
			// std::cout << core::kinematics::Stub(tmp1.residue(1).xyz("CB"),tmp1.residue(1).xyz("CA"),tmp1.residue(1).xyz("N")) << std::endl;
			// std::cout << a << std::endl;
			// std::cout << core::kinematics::Stub(tmp2.residue(1).xyz("CB"),tmp2.residue(1).xyz("CA"),tmp2.residue(1).xyz("N")) << std::endl;
			// std::cout << b << std::endl;

			// std::cout << compute_xform_score(xfs,tmp1,tmp2) << " == "
			//           << compute_xform_score(xfs,h.s1,h.s2,stubs,ss) << " " 
			//           << h.cbc << " " << h.xscore << std::endl;;

			// utility_exit_with_message("test");

			// std::cout <<"SCORES " << h.cbc <<" "<< compute_xform_score(xfs,h.s1,h.s2,stubs,ss) << " " << h.xscore << std::endl;;

		}
	}


	// dump top results
	for(int ic = 1; ic <= (int)syms.size(); ic++) {
		std::sort(hits[ic].begin(),hits[ic].end(),cmpcbc);
		for(int i = 1; i <= (int)min(nstruct,hits[ic].size()); ++i) {
			Hit & h(hits[ic][i]);
			std::string tag = utility::file_basename(fn);
			if(tag.substr(tag.size()-3)==".gz" ) tag = tag.substr(0,tag.size()-3);
			if(tag.substr(tag.size()-4)==".pdb") tag = tag.substr(0,tag.size()-4);			
			tag += "_C"+ObjexxFCL::string_of(h.sym)+"_cbc"+ObjexxFCL::string_of(i)+".pdb";
			cout << "RESULT " << h.sym << " " << h.iss << " " << NSS  << " " << h.irt << " " << h.cbc << " " << tag << endl;
			Pose tmp;
			make_dock_olig(init_pose,tmp,h);
			tmp.dump_pdb( option[out::file::o]()+"/"+tag);
			// utility_exit_with_message("test");
		}
	}
	// dump top results
	for(int ic = 1; ic <= (int)syms.size(); ic++) {
		std::sort(hits[ic].begin(),hits[ic].end(),cmpxsc);
		for(int i = 1; i <= (int)min(nstruct,hits[ic].size()); ++i) {
			Hit & h(hits[ic][i]);
			std::string tag = utility::file_basename(fn);
			if(tag.substr(tag.size()-3)==".gz" ) tag = tag.substr(0,tag.size()-3);
			if(tag.substr(tag.size()-4)==".pdb") tag = tag.substr(0,tag.size()-4);			
			tag += "_C"+ObjexxFCL::string_of(h.sym)+"_xsc"+ObjexxFCL::string_of(i)+".pdb";
			cout << "RESULT " << h.sym << " " << h.iss << " " << NSS  << " " << h.irt << " " << h.cbc << " " << tag << endl;
			Pose tmp;
			make_dock_olig(init_pose,tmp,h);
			tmp.dump_pdb( option[out::file::o]()+"/"+tag);
			// utility_exit_with_message("test");
		}
	}
	// dump top results
	for(int ic = 1; ic <= (int)syms.size(); ic++) {
		std::sort(hits[ic].begin(),hits[ic].end(),cmprmsd);
		for(int i = 1; i <= (int)min(nstruct,hits[ic].size()); ++i) {
			Hit & h(hits[ic][i]);
			std::string tag = utility::file_basename(fn);
			if(tag.substr(tag.size()-3)==".gz" ) tag = tag.substr(0,tag.size()-3);
			if(tag.substr(tag.size()-4)==".pdb") tag = tag.substr(0,tag.size()-4);			
			tag += "_C"+ObjexxFCL::string_of(h.sym)+"_rmsd"+ObjexxFCL::string_of(i)+".pdb";
			cout << "RESULT " << h.sym << " " << h.iss << " " << NSS  << " " << h.irt << " " << h.cbc << " " << tag << endl;
			Pose tmp;
			make_dock_olig(init_pose,tmp,h);
			tmp.dump_pdb( option[out::file::o]()+"/"+tag);
			// utility_exit_with_message("test");
		}
	}

	if(native_nfold != 1){ // benchmark
		// Real rms_goodcbc=999.0, rms_goodxsc=999.0;
		std::cout << "RESULT RMSD for best my contact count";
		std::sort(hits[1].begin(),hits[1].end(),cmpcbc);
		for(int i = 1; i <= (int)min(20ul,hits[1].size()); ++i) std::cout << " " << hits[1][i].rmsd;
		std::cout << endl;

		std::cout << "RESULT RMSD for best my xform score  ";
		std::sort(hits[1].begin(),hits[1].end(),cmpxsc);
		for(int i = 1; i <= (int)min(20ul,hits[1].size()); ++i) std::cout << " " << hits[1][i].rmsd;
		std::cout << endl;

		Real rms_goodcbc=999.0, rms_goodxsc=999.0;
		std::sort(hits[1].begin(),hits[1].end(),cmpcbc);
		for(int i = 1; i <= (int)min(20ul,hits[1].size()); ++i) rms_goodcbc = min(rms_goodcbc,hits[1][i].rmsd);
		Real best_cbc = hits[1][1].rmsd;
		std::sort(hits[1].begin(),hits[1].end(),cmpxsc);
		for(int i = 1; i <= (int)min(20ul,hits[1].size()); ++i) rms_goodxsc = min(rms_goodxsc,hits[1][i].rmsd);
		Real best_xsc = hits[1][1].rmsd;
		// rms_goodcbc /= 10.0;
		// rms_goodxsc /= 10.0;
		std::cout << "Best RMSD in top 20 by num. contacts: " << rms_goodcbc << " top " << best_cbc << std::endl;
		std::cout << "Best RMSD in top 20 by xform score:   " << rms_goodxsc << " top " << best_xsc << std::endl;		

	}


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
	for(Size ir = 1; ir <= pose.n_residue(); ++ir){
		if(pose.residue(ir).has("CA")){
			cen += pose.residue(ir).xyz("CA");
			ncen += 1.0;
		}
	}
	cen /= ncen;
	rot_pose( pose, Vec(0,0,1), -dihedral_degrees(Vec(1,0,0),Vec(0,0,0),Vec(0,0,1),cen) );
	rot_pose( pose, Vec(1,0,1), 180.0 ); // symZ to symX
	utility::vector1<Vec> CAs;
	for(Size ir = 1; ir <= pose.n_residue(); ++ir){
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
	register_options();
	devel::init(argc,argv);
	using namespace basic::options::OptionKeys;

	vector1<Vec> ssamp(NSS); read_sphere(ssamp);

	protocols::sic_dock::XfoxmScore const xfs("/Users/sheffler/project/designability_stats/results/");

	vector1<std::pair<string,int> > tasks;
	get_tasks_from_command_line(tasks);

	for(vector1<std::pair<string,int> >::const_iterator i = tasks.begin(); i != tasks.end(); ++i){
		string fn = i->first;
		int nfold = i->second;
		TR << "searching " << fn << " " << nfold << std::endl;

		Pose pnat;
		core::import_pose::pose_from_pdb(pnat,fn);
		core::scoring::dssp::Dssp dssp(pnat);
		dssp.insert_ss_into_pose(pnat);

		// center on z
		trans_pose(pnat,Vec(0,0,-center_of_geom(pnat,1,pnat.n_residue()).z()));

		Pose olig; make_native_olig(pnat,olig,nfold);

		// sym axis now along X, rotated around X so CA com in XZ plane,
		// this makes computing the RMSD of docked states easier
		vector1<Vec> native_ca = align_native_state(pnat,nfold); // align CA com
		pnat.dump_pdb("native_align.pdb");

		Vec cen = center_of_geom(pnat,1,pnat.n_residue());
		trans_pose(pnat,-Vec(0,cen.y(),cen.z()));

		if( pnat.n_residue() > MAX_NRES ) continue;
		Size cyscnt=0, nhelix=0;
		// for(Size ir = 2; ir <= pnat.n_residue()-1; ++ir) {
		// 	if(pnat.secstruct(ir) == 'H') nhelix++;
		// 	//if(!pnat.residue(ir).is_protein()) goto cont1;
		// 	if(pnat.residue(ir).is_lower_terminus()) remove_lower_terminus_type_from_pose_residue(pnat,ir);//goto cont1;
		// 	if(pnat.residue(ir).is_upper_terminus()) remove_upper_terminus_type_from_pose_residue(pnat,ir);//goto cont1;
		// 	if(pnat.residue(ir).name3()=="CYS") { if(++cyscnt > MAX_CYS_RES) goto cont1; }
		// } goto done1; cont1: TR << "skipping " << fn << std::endl; continue; done1:
		// if( nhelix < MIN_HELEX_RES ) continue;

		dock(xfs,pnat,fn,ssamp,olig,nfold,native_ca);
	}
}

