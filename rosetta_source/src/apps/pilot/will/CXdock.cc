// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:

#define sym_Clo 3
#define sym_Chi 12
#define CONTACT_D2 20.25
#define CONTACT_TH 5
#define NSS 8192 // 672 8192 17282
#define MIN_HELEX_RES 20
#define MAX_CYS_RES 3
#define MAX_NRES 200

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
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
#include <protocols/sic_dock/sic_dock.hh>

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
	Real cbc;
	core::kinematics::Stub s1,s2;
	Hit(int is, int ir, Real cb, int sm) : iss(is),irt(ir),sym(sm),cbc(cb) {}
};
bool cmp(Hit i,Hit j) { return i.cbc > j.cbc; }


// dock pose against rotated version of itself
// rotate all different ways by looping over rotation axes on unit sphere point
// samples in ssamp, then loopint over rotation magnitudes from 1-180°
// picks cases with most CB contacts in helices
void dock(Pose const init_pose, std::string const & fn, vector1<Vec> const & ssamp) {

	// set up one SIC per thread
	std::vector<protocols::sic_dock::SICFast*> sics_;
	sics_.resize(num_threads());
	for(int i = 0; i < num_threads(); ++i) sics_[i] = new protocols::sic_dock::SICFast;
	for(int i = 0; i < num_threads(); ++i) sics_[i]->init(init_pose);

	//store initial coords & angle samples 1-180 (other distribution of samps would be better...)
	vector1<double> asamp; for(Real i = 0; i < 180; ++i) asamp.push_back(i);

	// setup rotation matricies
	vector1<Mat> Rsym(sym_Chi);
	for(int ic = sym_Clo; ic <= sym_Chi; ic++) Rsym[ic] = rotation_matrix_degrees(Vec(1,0,0),360.0/Real(ic));
	vector1<vector1<Hit> > hits(sym_Chi);

	// loop over axes on sphere points & rotaton amounts 1-180
	#ifdef USE_OPENMP
	#pragma omp parallel for schedule(dynamic,1)
	#endif
	for(int iss = 1; iss <= (int)ssamp.size(); ++iss) {
		if(iss%1000==0) TR << iss << " of " << NSS << std::endl;
		for(int irt = 1; irt <= (int)asamp.size(); ++irt) {
			Mat R = rotation_matrix_degrees(ssamp[iss],(Real)asamp[irt]);
			Stub const x1( R, Vec(0,0,0));
			// loop over symmitreis to check
			for(int ic = sym_Clo; ic <= sym_Chi; ic++) {
				Stub const x2( Rsym[ic]*R, Vec(0,0,0) );
				Real cbc; Real t = sics_[thread_num()]->slide_into_contact(x1,x2,Vec(0,0,1),cbc);
				if(cbc >= CONTACT_TH) {
					#ifdef USE_OPENMP
					#pragma omp critical
					#endif
					{
						Hit h(iss,irt,cbc,ic);
						h.s1 = x1;
						h.s2 = x2;
						h.s1.v += t*Vec(0,0,1);
						hits[ic].push_back(h);
					}
				}
			}
		}
	}
	// dump top results
	Size nstruct = basic::options::option[basic::options::OptionKeys::out::nstruct]();
	std::string outdir = basic::options::option[basic::options::OptionKeys::out::file::o]()+"/";
	for(int ic = sym_Clo; ic <= sym_Chi; ic++) {
		std::sort(hits[ic].begin(),hits[ic].end(),cmp);
		for(int i = 1; i <= (int)min(nstruct,hits[ic].size()); ++i) {
			Hit & h(hits[ic][i]);
			std::string tag = utility::file_basename(fn);
			if(tag.substr(tag.size()-3)==".gz" ) tag = tag.substr(0,tag.size()-3);
			if(tag.substr(tag.size()-4)==".pdb") tag = tag.substr(0,tag.size()-4);			
			tag += "_C"+ObjexxFCL::string_of(h.sym)+"_"+ObjexxFCL::string_of(i)+".pdb";
			cout << "RESULT " << h.sym << " " << h.iss << " " << NSS  << " " << h.irt << " " << h.cbc << " " << tag << endl;
			Pose tmp1(init_pose),tmp2(init_pose);
			Real halfside = (h.s1.v-h.s2.v).length()/2.0;
			Real r = halfside / tan(numeric::constants::d::pi/(Real)h.sym); // tan(a/2) * edge_len/2
			xform_pose(tmp1,h.s1);
			trans_pose(tmp1,Vec(0,0,halfside));
			trans_pose(tmp1,Vec(0,r,0));
			rot_pose(tmp1,Vec(0,1,0),90.0);
			tmp1.dump_pdb(outdir+tag);
			// utility_exit_with_message("test");
		}
	}

}


int main(int argc, char *argv[]) {
	devel::init(argc,argv);
	using namespace basic::options::OptionKeys;

	// read in sphere poients
	vector1<Vec> ssamp(NSS); {
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
	// loop over input files, do some checks, call dock
	for(Size ifn = 1; ifn <= option[in::file::s]().size(); ++ifn) {
		string fn = option[in::file::s]()[ifn];
		Pose pnat;
		TR << "searching " << fn << std::endl;
		core::import_pose::pose_from_pdb(pnat,fn);
		trans_pose(pnat,-center_of_geom(pnat,1,pnat.n_residue()));
		core::scoring::dssp::Dssp dssp(pnat);
		dssp.insert_ss_into_pose(pnat);
		if( pnat.n_residue() > MAX_NRES ) continue;
		Size cyscnt=0, nhelix=0;
		for(Size ir = 2; ir <= pnat.n_residue()-1; ++ir) {
			if(pnat.secstruct(ir) == 'H') nhelix++;
			//if(!pnat.residue(ir).is_protein()) goto cont1;
			if(pnat.residue(ir).is_lower_terminus()) remove_lower_terminus_type_from_pose_residue(pnat,ir);//goto cont1;
			if(pnat.residue(ir).is_upper_terminus()) remove_upper_terminus_type_from_pose_residue(pnat,ir);//goto cont1;
			if(pnat.residue(ir).name3()=="CYS") { if(++cyscnt > MAX_CYS_RES) goto cont1; }
		} goto done1; cont1: TR << "skipping " << fn << std::endl; continue; done1:
		if( nhelix < MIN_HELEX_RES ) continue;
		Pose pala(pnat);
		dock(pala,fn,ssamp);
	}
}

